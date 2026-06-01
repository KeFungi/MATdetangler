"""Enumerate candidate alleles as walks through each k's SPAdes assembly graph.

Replaces `bubble_cands.py` + `repeat_stitch.py`. The merit of working on a GFA is to
recover alleles SPAdes never produced as a single contig — e.g., an allele whose path
crosses a repeat (MITE) that the kmer can't span. The previous design picked contigs
from `contigs.fasta` (SPAdes' linear output) and stitched pairs of contigs across
repeat-flanking matches, which was indirect: it could only recover alleles whose two
"halves" had already been correctly separated as contigs by SPAdes.

Here we go straight to the graph:

  1. Parse the GFA into segments (with sequence + depth) and directed L-line links.
  2. Label each segment by content via tblastn(variable_proteins) and blastn(flankL/R,
     repeats, known_degHD).
  3. DFS-enumerate simple paths from any flankL-bearing segment to any flankR-bearing
     segment, with constraints (max bp, max nodes, must touch a variable-gene segment).
  4. For each path, reconstruct the allele sequence by concatenating segment sequences
     with proper k-mer overlap (from L-line CIGAR) and reverse-complementing as the
     orientation tags dictate.
  5. Annotate the reconstructed sequence: which variable genes it covers, whether it
     has both flanks, mean depth along the path, whether it includes a repeat/degHD
     segment.

Per-k outputs:
  bubble_alleles_<k>.fasta
  bubble_alleles_<k>.ann.tsv  — name, len, k, segments, variable_genes,
                              flankL, flankR, cov, is_degHD, has_repeat
                              (flankL/flankR split lets downstream see open-bubble
                               arms — both T = complete; one T = partial / open)

Combined outputs (one row per k):
  bubble_alleles.fasta
  bubble_alleles.ann.tsv
"""
from __future__ import annotations
import os, sys, subprocess, argparse, tempfile, collections, re
from . import blast_utils as bu
from .paths import spades_k_paths

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def _rc(s: str) -> str: return s.translate(_COMP)[::-1]

# ---------- GFA parsing ----------

def parse_segments(gfa: str) -> dict[str, tuple[str, float]]:
    """Return {seg_id: (seq, depth)}. Depth from DP:f: tag, falling back to KC:i: / len."""
    out: dict[str, tuple[str, float]] = {}
    for ln in open(gfa):
        if ln[0] != "S": continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 3: continue
        sid, seq = f[1], f[2]
        depth = 0.0
        for tag in f[3:]:
            if tag.startswith("DP:f:"):
                try: depth = float(tag[5:])
                except ValueError: pass
                break
            if tag.startswith("KC:i:"):
                try: depth = int(tag[5:]) / max(1, len(seq))
                except ValueError: pass
                break
        out[sid] = (seq, depth)
    return out

def parse_links(gfa: str) -> list[tuple[str, str, str, str, int]]:
    """Return list of (seg1, orient1, seg2, orient2, overlap_bp) from L-lines."""
    out = []
    for ln in open(gfa):
        if ln[0] != "L": continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 6: continue
        m = re.match(r"(\d+)M", f[5])
        ov = int(m.group(1)) if m else 0
        out.append((f[1], f[2], f[3], f[4], ov))
    return out

def build_directed_adj(links: list[tuple[str, str, str, str, int]]
                       ) -> dict[tuple[str, str], list[tuple[str, str, int]]]:
    """SPAdes L-line semantics:
        L A +  B +  ov  =>  A+ -> B+ (overlap ov)
                            B- -> A- (overlap ov, reverse-traverse)
        L A +  B -  ov  =>  A+ -> B-, B+ -> A-
        L A -  B +  ov  =>  A- -> B+, B- -> A+
        L A -  B -  ov  =>  A- -> B-, B+ -> A+
    Build a directed graph keyed by (segment, orientation) for DFS path enumeration."""
    def opp(o): return "+" if o == "-" else "-"
    adj: dict[tuple[str, str], list[tuple[str, str, int]]] = collections.defaultdict(list)
    for s1, o1, s2, o2, ov in links:
        adj[(s1, o1)].append((s2, o2, ov))
        adj[(s2, opp(o2))].append((s1, opp(o1), ov))
    return adj

# ---------- content labeling ----------

def label_segments(segs: dict[str, tuple[str, float]], queries_dir: str,
                   repeats: str | None, known_degHD: str | None,
                   tblastn_pid: float = 30.0, tblastn_aa: int = 50,
                   blastn_var_pid: float = 80.0, blastn_var_minlen: int = 200,
                   flank_pid: float = 85.0, flank_minlen: int = 100,
                   degHD_pid: float = 95.0, degHD_minlen: int = 1000,
                   repeat_pid: float = 80.0, repeat_minlen: int = 50,
                   min_seg_len_for_label: int = 100,
                   threads: int = 1,
                   hd_labels_from_step22: dict[str, set[str]] | None = None,
                   blast_out_dir: str | None = None,
                   blast_tag: str = "",
                   ) -> tuple[dict[str, str], dict[str, set[str]]]:
    """Return (labels, var_per) where labels are "+"-joined feature tags per segment and
    var_per is {seg_id: {variable_gene_names}}.

    Search strategy:
    - Variable genes (HD-core): tblastn(variable_proteins) ∪ blastn(variable_nt)
      tblastn rescues highly-divergent allele variants whose nucleotide identity is too
      low to register; blastn catches the conserved nucleotide stretches the tblastn
      threshold may underweight. The union gives the more sensitive labeling.
    - Flanks (flankL, flankR): blastn only. Flanks are conserved within a species — the
      nucleotide signal is reliable and there's no protein query to translate.
    - Repeats / degHD: blastn (nucleotide refs).
    """
    proteins   = os.path.join(queries_dir, "variable_proteins.fasta")
    variable_nt = os.path.join(queries_dir, "variable_nt.fasta")
    flankL     = os.path.join(queries_dir, "flankL.fasta")
    flankR     = os.path.join(queries_dir, "flankR.fasta")
    feats: dict[str, list[str]] = collections.defaultdict(list)
    var_per: dict[str, set[str]] = collections.defaultdict(set)
    # Default min_seg_len_for_label=0: label EVERY segment regardless of length
    # (2026-05-30 — user-requested). Previously short connector segs (50-300 bp)
    # were filtered out for perf, but that left them unlabeled in the per-seg
    # labels file used by pick_alleles' pair-topology classifier. Keep the param
    # configurable so a caller can still raise the cut if blast wall-time becomes
    # an issue on huge GFAs.
    blast_segs = {sid: seq for sid, (seq, _) in segs.items()
                  if len(seq) >= min_seg_len_for_label}
    # Build out_tsv path for cached blast calls (or None to skip caching).
    # All cached blast tsvs share the "blast_" prefix so the wrapper's --re-blast
    # flag can wipe them with a single glob.
    def _out(prog: str, query: str) -> str | None:
        if not blast_out_dir or not blast_tag: return None
        return os.path.join(blast_out_dir, f"blast_gfa_{query}_{prog}_{blast_tag}.tsv")
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "n.fa")
        with open(sf, "w") as o:
            for sid, seq in blast_segs.items(): o.write(f">{sid}\n{seq}\n")
        db = bu.fasta_to_db(sf, t, name="db")
        # ----- VARIABLE GENES (HD) -----
        # Fast path: if the caller already ran HD-search against the GFA segments
        # (step 2.2 anchor_search --source segments) and supplied the per-segment
        # HD labels, reuse them — skip the expensive tblastn/blastn-variable-nt
        # re-blast on the BFS neighborhood. Otherwise run the full HD-blast pass.
        if hd_labels_from_step22 is not None:
            # Use precomputed HD labels (subset of segs that step 2.2 found HD genes on);
            # segments not in the dict are assumed to have NO HD genes (i.e. step 2.2
            # already searched everything and these are below threshold).
            for sid, genes in hd_labels_from_step22.items():
                if sid not in blast_segs: continue  # not in our BFS neighborhood; skip
                for gene in genes:
                    feats[sid].append(gene); var_per[sid].add(gene)
        else:
            # variable genes — tblastn primary, blastn fallback only when needed
            # (same pattern as bubble_cands: skip the expensive blastn call when tblastn alone
            # already labels at least one segment per variable gene)
            if os.path.exists(proteins):
                for r in bu.tblastn_hits(proteins, db, min_pid=tblastn_pid,
                                          min_aa=tblastn_aa, threads=threads,
                                          out_tsv=_out("tblastn", "hd")):
                    feats[r[0]].append(r[1]); var_per[r[0]].add(r[1])
            # how many variable genes did tblastn surface across all segments?
            genes_via_tblastn = set().union(*var_per.values()) if var_per else set()
            if os.path.exists(variable_nt):
                # only fall back to nucleotide blastn if tblastn missed at least one variable gene
                nvar_total = sum(1 for ln in open(proteins) if ln.startswith(">")) if os.path.exists(proteins) else 0
                if len(genes_via_tblastn) < nvar_total:
                    for r in bu.blastn_hits(variable_nt, db, min_pid=blastn_var_pid,
                                             min_len=blastn_var_minlen, threads=threads,
                                             out_tsv=_out("blastn", "hd_nt")):
                        sid, gene = r[0], r[1]
                        if gene not in var_per[sid]:
                            feats[sid].append(gene); var_per[sid].add(gene)
        # flanks: blastn only
        if os.path.exists(flankL):
            for sid in bu.hits_ids(bu.blastn_hits(flankL, db, min_pid=flank_pid,
                                                  min_len=flank_minlen, threads=threads,
                                                  out_tsv=_out("blastn", "flankL"))):
                feats[sid].append("flankL")
        if os.path.exists(flankR):
            for sid in bu.hits_ids(bu.blastn_hits(flankR, db, min_pid=flank_pid,
                                                  min_len=flank_minlen, threads=threads,
                                                  out_tsv=_out("blastn", "flankR"))):
                feats[sid].append("flankR")
        if known_degHD and os.path.exists(known_degHD):
            for sid in bu.hits_ids(bu.blastn_hits(known_degHD, db, min_pid=degHD_pid,
                                                  min_len=degHD_minlen, threads=threads,
                                                  out_tsv=_out("blastn", "degHD"))):
                feats[sid].append("degHD")
        if repeats and os.path.exists(repeats):
            for sid in bu.hits_ids(bu.blastn_hits(repeats, db, min_pid=repeat_pid,
                                                  min_len=repeat_minlen, threads=threads,
                                                  out_tsv=_out("blastn", "repeat"))):
                feats[sid].append("repeat")
    labels: dict[str, str] = {}
    for sid in segs:
        fs = list(dict.fromkeys(feats.get(sid, [])))
        labels[sid] = "degHD" if "degHD" in fs else "+".join(fs)
    return labels, var_per

# ---------- path enumeration + sequence reconstruction ----------

def _flip(o: str) -> str: return "-" if o == "+" else "+"


def mirror_path(path: list[tuple[str, str, int]]) -> tuple[tuple[str, str, int], ...]:
    """The same physical walk traversed in the opposite direction.
        mirror[0]   = (P[-1].seg, flip(P[-1].o), 0)
        mirror[i>0] = (P[-1-i].seg, flip(P[-1-i].o), P[-i].ov)
    Overlaps shift one slot because in P the overlap at index i sits BETWEEN
    P[i-1] and P[i]; in the mirror that same junction is at the corresponding
    index from the other end.
    """
    n = len(path)
    if n == 0: return tuple()
    out = [(path[-1][0], _flip(path[-1][1]), 0)]
    for i in range(1, n):
        seg, orient, _ = path[n - 1 - i]
        ov = path[n - i][2]
        out.append((seg, _flip(orient), ov))
    return tuple(out)


def enumerate_paths(adj, segs, starts: set[str], ends: set[str],
                    must_visit_any: set[str],
                    max_bp: int, max_nodes: int, max_paths: int = 50000
                    ) -> list[list[tuple[str, str, int]]]:
    """Iterative DFS. Each path element is (seg, orient, overlap_to_prev).
    The first element's overlap is 0. The path may start at any orientation of each `start`.
    A path is reported when it reaches any `end` segment (in either orientation) AND has
    touched at least one segment in `must_visit_any` (so flank-only walks are rejected).

    RC-mirror dedup (pure graph topology, no sequence): when the SAME physical
    walk could emerge twice with reverse-complemented orientation — possible
    when starts and ends overlap or the graph is symmetric enough that both
    members of a mirror pair get DFS-enumerated — we keep the first emitted
    and drop the second. If a walk's mirror is never enumerated (the typical
    HD case with disjoint flankL/flankR seeds vs ends), nothing is dropped.
    """
    paths: list[list[tuple[str, str, int]]] = []
    emitted: set[tuple[tuple[str, str, int], ...]] = set()
    for s0 in sorted(starts):
        for o0 in ("+", "-"):
            seq0_len = len(segs[s0][0])
            stack = [(s0, o0, [(s0, o0, 0)], {s0}, seq0_len)]
            while stack:
                if len(paths) >= max_paths: return paths
                seg, orient, path, visited, bp = stack.pop()
                if seg in ends and len(path) >= 1 and (visited & must_visit_any):
                    tp = tuple(path)
                    if mirror_path(path) not in emitted:
                        emitted.add(tp); paths.append(list(path))
                    continue
                # `max_bp <= 0` is the no-bp-limit sentinel — skip the budget checks entirely
                if len(path) >= max_nodes: continue
                if max_bp > 0 and bp >= max_bp: continue
                for next_seg, next_orient, ov in adj.get((seg, orient), []):
                    if next_seg in visited: continue
                    nb = bp + len(segs[next_seg][0]) - ov
                    if max_bp > 0 and nb > max_bp: continue
                    stack.append((next_seg, next_orient, path + [(next_seg, next_orient, ov)],
                                  visited | {next_seg}, nb))
    return paths

def reconstruct(path: list[tuple[str, str, int]], segs: dict[str, tuple[str, float]]) -> str:
    out = []
    for i, (sid, orient, ov) in enumerate(path):
        seq = segs[sid][0]
        if orient == "-": seq = _rc(seq)
        out.append(seq if i == 0 else seq[ov:])
    return "".join(out)

def path_depth(path: list[tuple[str, str, int]], segs) -> float:
    """Length-weighted mean depth across the path's segments."""
    num, den = 0.0, 0
    for sid, _, _ in path:
        seq, dp = segs[sid]
        num += dp * len(seq); den += len(seq)
    return num / den if den else 0.0

# ---------- main ----------

def run_one_k(sample: str, k: str, gfa: str, queries_dir: str, outdir: str,
              repeats: str | None, known_degHD: str | None,
              max_locus_len: int, min_allele_len: int,
              max_nodes: int, max_paths: int, threads: int) -> tuple[str, str]:
    segs = parse_segments(gfa)
    links = parse_links(gfa)
    adj = build_directed_adj(links)
    labels, var_per = label_segments(segs, queries_dir, repeats, known_degHD, threads=threads)
    # which segments carry what?
    starts = {sid for sid, l in labels.items() if "flankL" in l.split("+")}
    ends   = {sid for sid, l in labels.items() if "flankR" in l.split("+")}
    var_segs = {sid for sid, vs in var_per.items() if vs}
    fa_path  = os.path.join(outdir, f"bubble_alleles_{k}.fasta")
    tsv_path = os.path.join(outdir, f"bubble_alleles_{k}.ann.tsv")
    if not starts or not ends or not var_segs:
        open(fa_path, "w").close(); open(tsv_path, "w").close()
        print(f"  [{k}] no flankL/flankR/variable-gene segments — no candidates")
        return fa_path, tsv_path
    paths = enumerate_paths(adj, segs, starts, ends, var_segs,
                             max_bp=max_locus_len, max_nodes=max_nodes, max_paths=max_paths)
    # dedup paths whose reconstructed sequence is identical (incl. revcomp). On a tie we keep
    # the path with MORE nodes — it preserves more graph structure for the visualization, and
    # the picker's length tiebreaker already prefers longer within --max-locus-len. Picking the
    # shortest would silently collapse multi-segment walks into single-segment views.
    seen: dict[str, list[tuple[str, str, int]]] = {}
    for p in paths:
        s = reconstruct(p, segs)
        key = min(s, _rc(s))
        if key not in seen or len(p) > len(seen[key]): seen[key] = p
    # write outputs
    n_written = 0
    with open(fa_path, "w") as fa, open(tsv_path, "w") as tsv:
        for p in seen.values():
            seq = reconstruct(p, segs)
            L = len(seq)
            if L < min_allele_len: continue
            seg_str = ",".join(f"{sid}{o}" for sid, o, _ in p)
            var_hits = set().union(*(var_per.get(sid, set()) for sid, _, _ in p))
            has_flankL = any("flankL" in labels[sid].split("+") for sid, _, _ in p)
            has_flankR = any("flankR" in labels[sid].split("+") for sid, _, _ in p)
            is_degHD = any(labels[sid] == "degHD" for sid, _, _ in p)
            has_repeat = any("repeat" in labels[sid].split("+") for sid, _, _ in p)
            dp = path_depth(p, segs)
            nm = f"{sample}__path_{k}_n{len(p)}_L{L}_d{dp:.1f}_p{n_written}"
            fa.write(f">{nm}\n")
            for i in range(0, L, 80): fa.write(seq[i:i + 80] + "\n")
            tsv.write(f"{nm}\t{L}\t{k}\t{seg_str}\t{','.join(sorted(var_hits)) or '-'}"
                      f"\t{('T' if (has_flankL and has_flankR) else 'F')}"
                      f"\t{dp:.2f}\t{'T' if is_degHD else 'F'}\t{'T' if has_repeat else 'F'}\n")
            n_written += 1
    print(f"  [{k}] enumerated {len(paths)} raw paths -> {len(seen)} unique -> {n_written} kept (>= {min_allele_len} bp)")
    return fa_path, tsv_path

def run(sample: str, spades_dir: str, queries_dir: str, ks: list[str], outdir: str,
        repeats: str | None = None, known_degHD: str | None = None,
        max_locus_len: int = 12000, min_allele_len: int = 2000,
        max_nodes: int = 15, max_paths: int = 50000, threads: int = 4) -> str:
    os.makedirs(outdir, exist_ok=True)
    cand_fa = os.path.join(outdir, "bubble_alleles.fasta")
    cand_tsv = os.path.join(outdir, "bubble_alleles.ann.tsv")
    open(cand_fa, "w").close(); open(cand_tsv, "w").close()
    for k in ks:
        _, gfa = spades_k_paths(spades_dir, k)
        if not gfa:
            print(f"  [{k}] no GFA at {spades_dir}/k{k.lstrip('kK')}/; skipping"); continue
        fa, tsv = run_one_k(sample, k, gfa, queries_dir, outdir,
                             repeats, known_degHD, max_locus_len, min_allele_len,
                             max_nodes, max_paths, threads)
        with open(cand_fa, "a") as o:
            for ln in open(fa): o.write(ln)
        with open(cand_tsv, "a") as o:
            for ln in open(tsv): o.write(ln)
    nrec = sum(1 for ln in open(cand_fa) if ln.startswith(">"))
    print(f"[segment_alleles] {sample}: {nrec} candidate paths across {len(ks)} k's -> {cand_fa}")
    return cand_fa

def _cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--sample", required=True)
    p.add_argument("--spades-dir", required=True)
    p.add_argument("--queries-dir", required=True)
    p.add_argument("--ks", default="k21,k33,k55")
    p.add_argument("--outdir", required=True)
    p.add_argument("--repeats", default=None)
    p.add_argument("--known-degHD", default=None)
    p.add_argument("--max-locus-len", type=int, default=12000)
    p.add_argument("--min-allele-len", type=int, default=2000)
    p.add_argument("--max-nodes", type=int, default=15)
    p.add_argument("--max-paths", type=int, default=50000)
    p.add_argument("--threads", type=int, default=4)
    a = p.parse_args(argv)
    run(a.sample, a.spades_dir, a.queries_dir, a.ks.split(","), a.outdir,
        a.repeats, a.known_degHD, a.max_locus_len, a.min_allele_len,
        a.max_nodes, a.max_paths, a.threads)

if __name__ == "__main__":
    _cli()
