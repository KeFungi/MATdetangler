"""Pick the two locus alleles per sample from the combined candidate pool.

Completeness-aware + coverage-aware. Generalizes the Pcub pick_primary.py to any number of variable
genes (BED-driven) and any locus annotation.

Inputs:
  --candidates FASTA           — combined candidate pool (anchor_contig.fasta + bubble_alleles.fasta;
                                 cat'd by the wrapper as cand_combined.fasta)
  queries_dir/                 — flankL.fasta, flankR.fasta, variable_proteins.fasta
  (optional) known_degHD       — degenerated-element reference fasta
  expected_count               — 1 (haploid) or 2 (dikaryon)
  genome_coverage              — used for the coverage-aware tiebreak

Output:
  primary_alleles.fasta    A1 (+ A2 if dikaryon), headers carry origin/k/type tags
  picks.tsv                sample,name,origin,k,type,len,from_contig,cov,both_vars,is_degHD
"""
from __future__ import annotations
import os, sys, re, subprocess, argparse, tempfile, collections
from .input_process import read_fasta
from . import blast_utils as bu

DUP_PID = 0.95
DUP_AF  = 0.80

def _kof(name: str) -> str | None:
    m = re.search(r"_(?:path|bubble)_(k\d+)_", name); return m.group(1) if m else None

def _cov(name: str) -> float:
    # segment_alleles names: "...__path_<k>_n<N>_L<len>_d<depth>_p<idx>"
    # contig-style names:     "...NODE_<i>_length_<L>_cov_<depth>"
    m = re.search(r"_d([\d.]+)_p\d+$", name)
    if m: return float(m.group(1))
    m = re.search(r"cov_([\d.]+)", name); return float(m.group(1)) if m else 0.0

def _vars_in(contig: str, seq_for_db: dict[str, str], variable_proteins_fa: str,
             tblastn_pid: float, tblastn_aa: int) -> dict[str, set[str]]:
    """For each contig name, which variable genes does it carry?"""
    out = collections.defaultdict(set)
    if not seq_for_db: return out
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "s.fa")
        with open(sf, "w") as o:
            for n, s in seq_for_db.items(): o.write(f">{n}\n{s}\n")
        db = bu.fasta_to_db(sf, t, name="db")
        for r in bu.tblastn_hits(variable_proteins_fa, db, min_pid=tblastn_pid, min_aa=tblastn_aa):
            out[r[0]].add(r[1])
    return out

def _vars_aa_coverage(seq_for_db: dict[str, str], variable_proteins_fa: str,
                      tblastn_pid: float, tblastn_aa: int) -> dict[str, int]:
    """For each candidate name, return the total tblastn aa-coverage across all variable genes.

    Higher = the candidate covers more of the variable proteins, which is the right
    discriminator between "complete but truncated" candidates and "fully covers both alleles".
    For two candidates that both pass `has_all_vars` (some hit to every variable gene), the
    one with higher aa coverage actually contains a larger fraction of the protein lengths.
    Per-gene merged intervals would be ideal, but the tblastn hit lengths give a useful
    approximation: a candidate that covers all 735 aa of HD1 has a hit length around 735,
    while a truncated one hits only ~300 aa.
    """
    out: dict[str, int] = collections.defaultdict(int)
    if not seq_for_db: return out
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "s.fa")
        with open(sf, "w") as o:
            for n, s in seq_for_db.items(): o.write(f">{n}\n{s}\n")
        db = bu.fasta_to_db(sf, t, name="db")
        for r in bu.tblastn_hits(variable_proteins_fa, db, min_pid=tblastn_pid, min_aa=tblastn_aa):
            try: out[r[0]] += int(r[3])  # r[3] = aa-length of the tblastn hit
            except (IndexError, ValueError): pass
    return out

def _flank_hits(seq_for_db: dict[str, str], flank_fa: str,
                min_len: int = 100, min_pid: float = 85.0) -> set[str]:
    if not seq_for_db or not os.path.exists(flank_fa): return set()
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "s.fa")
        with open(sf, "w") as o:
            for n, s in seq_for_db.items(): o.write(f">{n}\n{s}\n")
        db = bu.fasta_to_db(sf, t, name="db")
        return bu.hits_ids(bu.blastn_hits(flank_fa, db, min_pid=min_pid, min_len=min_len))

def _flank_coverage_bp(seq_for_db: dict[str, str], flank_fas: list[str],
                       min_pid: float = 85.0, min_len: int = 100) -> dict[str, int]:
    """For each candidate, total nucleotide blastn coverage against the flank queries (sum of
    aligned lengths across flankL.fasta + flankR.fasta). This discriminates truncated alleles
    (covers only part of the flanks) from fully-extended alleles (covers the entire flank).
    Higher = more complete flank reach. Critical because `has_both_flanks` is binary and only
    requires a 100 bp hit — a 200-bp hit and a 5000-bp hit are equal under that test."""
    out: dict[str, int] = collections.defaultdict(int)
    if not seq_for_db: return out
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "s.fa")
        with open(sf, "w") as o:
            for n, s in seq_for_db.items(): o.write(f">{n}\n{s}\n")
        db = bu.fasta_to_db(sf, t, name="db")
        for fa in flank_fas:
            if not fa or not os.path.exists(fa): continue
            for r in bu.blastn_hits(fa, db, min_pid=min_pid, min_len=min_len):
                # blastn_hits returns: sseqid, qseqid, pident, length, ...
                try: out[r[0]] += int(r[3])
                except (IndexError, ValueError): pass
    return out

def _degHD_hits(seq_for_db: dict[str, str], known_degHD: str | None,
                min_pid: float = 95.0, min_len: int = 1000) -> set[str]:
    if not known_degHD or not os.path.exists(known_degHD) or not seq_for_db: return set()
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "s.fa")
        with open(sf, "w") as o:
            for n, s in seq_for_db.items(): o.write(f">{n}\n{s}\n")
        db = bu.fasta_to_db(sf, t, name="db")
        return bu.hits_ids(bu.blastn_hits(known_degHD, db, min_pid=min_pid, min_len=min_len))

def _trim_locus(seq: str, variable_proteins_fa: str, pad: int = 2500) -> str:
    """Trim a long contig to the locus window (variable-gene hits +/- pad bp)."""
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "s.fa")
        open(sf, "w").write(f">s\n{seq}\n")
        db = bu.fasta_to_db(sf, t, name="db")
        rows = bu.tblastn_hits(variable_proteins_fa, db, extra_outfmt="sstart send")
    pts = []
    for r in rows:
        # sstart, send are extra cols
        try: a, b = int(r[4]), int(r[5])
        except (IndexError, ValueError): continue
        pts += [a, b]
    if not pts: return seq
    lo = max(0, min(pts) - pad); hi = min(len(seq), max(pts) + pad)
    return seq[lo:hi]

def _hd_core_pairwise_identities(seq_map: dict[str, str],
                                    queries_dir: str,
                                    locus_ref_fa: str | None
                                    ) -> dict[tuple[str, str], tuple[float, float]]:
    """One-shot HD-core pairwise identity for ALL pairs in seq_map.

    Uses the SAME machinery as step 3's `cluster_alleles`: one MAFFT MSA on
    (locus_ref + candidates trimmed to HD-core ± 1000 bp), HD-core columns
    projected via the reference's ungapped HD-protein-hit span. Identity is
    computed ONLY on those columns. Matches the pipeline's other MAFFT tests
    (step 3 dedup, cluster step) instead of relying on whole-allele MAFFT
    (which is dominated by shared flanks in closed-bubble samples).

    Returns {sorted-tuple(a, b): (alnid, aln_frac)} for every pair.
    Empty dict if locus_ref_fa is missing — caller should fall back to
    whole-allele `mafft_pair`.
    """
    if not locus_ref_fa or not os.path.exists(locus_ref_fa): return {}
    if len(seq_map) < 2: return {}
    from .graph_path_search import _align_cores
    from .pairwise_identity import _score_columns
    proteins = os.path.join(queries_dir, "variable_proteins.fasta")
    named = [(name, s) for name, s in seq_map.items()]
    try:
        aligned, hd_cols = _align_cores(named, locus_ref_fa, proteins)
    except Exception as e:
        print(f"[pick_alleles] _align_cores failed ({e}); HD-core pair identities skipped")
        return {}
    out: dict[tuple[str, str], tuple[float, float]] = {}
    names = sorted(aligned.keys())
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            sa, sb = aligned[a], aligned[b]
            if not hd_cols:
                alnid, aln_frac = 0.0, 0.0
            else:
                # use the shorter HD-core ungapped length as the denominator
                # (same convention as pairwise_identity.mafft_pair_core)
                hd_len_a = sum(1 for ci in hd_cols if ci < len(sa) and sa[ci] != '-')
                hd_len_b = sum(1 for ci in hd_cols if ci < len(sb) and sb[ci] != '-')
                denom = min(hd_len_a, hd_len_b) if (hd_len_a and hd_len_b) else 0
                alnid, aln_frac = _score_columns(sa, sb, hd_cols, denom_len=denom)
            out[(a, b)] = (alnid, aln_frac)
    return out


_SEG_TOKEN_RE = re.compile(r"^(\d+)([+-])?$")

def _parse_segs_string(s: str) -> list[tuple[str, str]]:
    """Parse "22399580-,4622205-,31967890+" into [('22399580','-'), ...]."""
    out: list[tuple[str, str]] = []
    if not s or s == "-": return out
    for tok in s.split(","):
        tok = tok.strip()
        if not tok: continue
        m = _SEG_TOKEN_RE.match(tok)
        if m: out.append((m.group(1), m.group(2) or "+"))
    return out


def _read_seg_labels(outdir: str, k: str) -> dict[str, str]:
    """Load per-segment labels written by step 3 (graph_path_search) as
    seg_labels_<k>.tsv. Returns {seg_id: label_str} where label_str is "+"-joined
    feature tags (e.g. "flankL", "flankR+HD1", "HD2") or "" for unlabeled.
    Returns {} if the file is missing — pair topology classification then degrades
    to rank 0 (other) for everything."""
    out: dict[str, str] = {}
    f = os.path.join(outdir, f"seg_labels_{k}.tsv")
    if not os.path.exists(f): return out
    for ln in open(f):
        fld = ln.rstrip("\n").split("\t")
        if len(fld) < 2: continue
        out[fld[0]] = fld[1]
    return out


def _segments_via_contigs_paths(cand_name: str, spades_dir: str | None,
                                  paths_cache: dict[str, dict]) -> list[tuple[str, str]]:
    """Resolve a bubble-origin candidate's segment walk by looking up the
    underlying SPAdes contig in <spades_dir>/<k>/contigs.paths.

    cand_name format from anchor_search: "<sample>__bubble_<k>_<NODE_*_length_*_cov_*>"
    contigs.paths key:                   "<NODE_*_length_*_cov_*>" (drops trailing _<component>)

    paths_cache is a {k: {contig_name: [(seg, orient), ...]}} cache so we don't
    re-parse contigs.paths for every candidate.
    """
    if not spades_dir: return []
    m = re.match(r"^.+?__bubble_(k\d+)_(NODE_\d+_length_\d+_cov_[\d.]+)", cand_name)
    if not m: return []
    k, contig_key = m.group(1), m.group(2)
    if k not in paths_cache:
        from .graph_path_search import parse_contigs_paths
        paths_file = os.path.join(spades_dir, k, "contigs.paths")
        paths_cache[k] = parse_contigs_paths(paths_file)
    return paths_cache[k].get(contig_key, [])


def _pair_topology_rank(segs_a: list[tuple[str, str]], segs_b: list[tuple[str, str]],
                          labels: dict[str, str]) -> tuple[int, str]:
    """Pair-level topology rank. Compares the two candidates' segment SETS and
    classifies based on which (if any) flank-bearing segments they share.

    Returns (rank, label):
        (3, "closed_bubble") — share ≥1 flankL-bearing seg AND ≥1 flankR-bearing seg
        (2, "open_bubble")   — share ≥1 flank-bearing seg (one side only)
        (1, "detached")      — share segs but none flank-bearing
        (0, "other")         — share no segs (or labels unavailable)
    """
    ids_a = {s for s, _ in segs_a}
    ids_b = {s for s, _ in segs_b}
    shared = ids_a & ids_b
    if not shared: return (0, "other")
    shared_flankL = any("flankL" in labels.get(s, "") for s in shared)
    shared_flankR = any("flankR" in labels.get(s, "") for s in shared)
    if shared_flankL and shared_flankR: return (3, "closed_bubble")
    if shared_flankL or shared_flankR:  return (2, "open_bubble")
    return (1, "detached")


def _flank_completeness_rank(a_in_L: bool, a_in_R: bool,
                               b_in_L: bool, b_in_R: bool) -> int:
    """Joint flank-completeness rank for a pair. Ordering:
        (2,2)            -> 22  ("both both")
        (2,1) or (1,2)   -> 12  ("one both + one one")
        (1,1)            -> 11  ("each ≥1, no both")
        (2,0) or (0,2)   ->  2  ("one both + one none")
        (1,0) or (0,1)   ->  1  ("one one + one none")
        (0,0)            ->  0  ("both none")
    Naturally orders: 22 > 12 > 11 > 2 > 1 > 0. Respects the "balanced over
    lopsided" intuition: (1,1) > (2,0) — having BOTH alleles each with one flank
    beats one allele with both + other with none.
    """
    na = (1 if a_in_L else 0) + (1 if a_in_R else 0)
    nb = (1 if b_in_L else 0) + (1 if b_in_R else 0)
    lo, hi = (na, nb) if na <= nb else (nb, na)
    return lo * 10 + hi


def _var_completeness_rank(nvars_a: int, nvars_b: int) -> int:
    """Joint variable-gene completeness rank for a pair, analogous to
    _flank_completeness_rank but for variable genes. Encoded as lo*10+hi so that
    a balanced pair (e.g., 1+1) outranks a lopsided one (2+0).

    For HD locus with 2 variable genes:
        (2,2) -> 22, (2,1)/(1,2) -> 12, (1,1) -> 11, (2,0) -> 2, (1,0) -> 1, (0,0) -> 0
    Naturally orders: 22 > 12 > 11 > 2 > 1 > 0.
    """
    lo, hi = (nvars_a, nvars_b) if nvars_a <= nvars_b else (nvars_b, nvars_a)
    return lo * 10 + hi


def _score_pair(a: str, b: str,
                  seqs: dict[str, str],
                  seg_of_cand: dict[str, list[tuple[str, str]]],
                  labels: dict[str, str],
                  hL: set[str], hR: set[str],
                  var_per: dict[str, set[str]],
                  var_aa_cov: dict[str, int],
                  flank_cov: dict[str, int],
                  pair_alnid: float) -> tuple:
    """Lexicographic / sequential pair score — first criterion is primary; ties
    broken by second; etc.

    (var_compl_rank, topo_rank, flank_compl_rank, pair_divergence,
     joint_HD_aa_cov, joint_flank_bp_cov, -|len_a - len_b|)

    Variable-gene completeness ranks HIGHEST: a pair where each allele carries
    the variable genes wins over a topologically-cleaner pair where one allele
    is missing them. Topology and flank completeness are next; HD-core
    divergence is only a tiebreaker within the same completeness/topology
    layer. This preserves the biallelic-pair preference for true closed
    bubbles only when both alleles are individually complete.

    pair_divergence = -pair_alnid (lower MAFFT identity = more divergent = higher score).
    """
    topo, _ = _pair_topology_rank(seg_of_cand.get(a, []), seg_of_cand.get(b, []), labels)
    vc = _var_completeness_rank(len(var_per.get(a, set())), len(var_per.get(b, set())))
    fc = _flank_completeness_rank(a in hL, a in hR, b in hL, b in hR)
    aa = var_aa_cov.get(a, 0) + var_aa_cov.get(b, 0)
    fb = flank_cov.get(a, 0) + flank_cov.get(b, 0)
    diff = abs(len(seqs[a]) - len(seqs[b]))
    pair_divergence = -pair_alnid    # negative MAFFT id → more divergent ranks higher
    return (vc, topo, fc, pair_divergence, aa, fb, -diff)


def _is_dup(a: str, b: str,
            queries_dir: str | None = None,
            locus_ref_fa: str | None = None) -> bool:
    """Same biological allele?

    Prefers HD-core MAFFT identity (matching step 3's `cluster_alleles` and the
    new pair-divergence scoring) when queries_dir + locus_ref_fa are supplied.
    Falls back to whole-allele `mafft_pair` otherwise. Both paths apply the
    same DUP_PID / DUP_AF thresholds.

    Note: the modern picker (`run()` → `best_pair_in_k`) bypasses this function
    by precomputing the full HD-core pairwise matrix once per K and reusing it
    for both the dup-filter and the divergence-score. This function exists for
    backward compatibility and standalone use."""
    if not a or not b: return False
    if queries_dir and locus_ref_fa and os.path.exists(locus_ref_fa):
        d = _hd_core_pairwise_identities({"_a": a, "_b": b}, queries_dir, locus_ref_fa)
        if ("_a", "_b") in d:
            alnid, aln_frac = d[("_a", "_b")]
            return alnid >= DUP_PID and aln_frac >= DUP_AF
    try:
        from .pairwise_identity import mafft_pair
        alnid, aln_frac = mafft_pair(a, b)
    except Exception:
        return False
    return alnid >= DUP_PID and aln_frac >= DUP_AF

def _trim_outer_flanks(picks: list[str],
                        seqs: dict[str, str],
                        seg_of_cand_list: dict[str, list[tuple[str, str]]],
                        seg_labels_by_k: dict[str, dict[str, str]],
                        spades_dir: str | None,
                        sample: str) -> None:
    """Drop outer flanks from the picked allele(s):
      * PAIR — if two flank-only fragments share the same GFA segment ID
        across the two arms, drop everything OUTSIDE those anchors. For each
        arm: trim to [leftmost shared flankL : rightmost shared flankR].
      * SINGLETON — drop the very-first segment if it's flank-only AND the
        very-last if it's flank-only.

    Modifies `seqs` and `seg_of_cand_list` in place — both downstream emitters
    (primary_alleles.fasta, picks.tsv) read those structures.

    Skips silently for cross-K pair picks (the two k's would need separate GFA
    loads with no shared segment IDs by definition).
    """
    if not picks or not spades_dir: return
    from .GFA_search import parse_segments, parse_links, _rc

    def is_flank_only(sid: str, labels: dict[str, str]) -> bool:
        """True if the segment's label has ONLY flank tokens (no variable gene).
        Treats flankL and flankR as a single 'flank-only' category — the walk
        position determines whether a given anchor is a left or right boundary,
        not the flank tag itself (so reversed walks still trim correctly).
        """
        lab = labels.get(sid, "")
        toks = [t for t in lab.split("+") if t]
        return bool(toks) and all(t.startswith("flank") for t in toks)

    def has_var(sid: str, labels: dict[str, str]) -> bool:
        lab = labels.get(sid, "")
        toks = [t for t in lab.split("+") if t]
        return any(not t.startswith("flank") for t in toks)

    # Lazy per-k GFA loader; each k's segment set + overlap map is parsed once.
    gfa_cache: dict[str, tuple[dict, dict] | tuple[None, None]] = {}
    def _opp(o): return "+" if o == "-" else "-"
    def load_k(k: str):
        if k in gfa_cache: return gfa_cache[k]
        gfa = os.path.join(spades_dir, k, "assembly_graph_after_simplification.gfa")
        if not os.path.exists(gfa): gfa_cache[k] = (None, None); return gfa_cache[k]
        segs = parse_segments(gfa)
        ov_map: dict[tuple[str, str, str, str], int] = {}
        for sa, oa, sb, ob, ov in parse_links(gfa):
            ov_map[(sa, oa, sb, ob)] = ov
            ov_map[(sb, _opp(ob), sa, _opp(oa))] = ov
        gfa_cache[k] = (segs, ov_map)
        return gfa_cache[k]

    def reconstruct(walk: list[tuple[str, str]], segs: dict, ov_map: dict) -> str:
        out: list[str] = []
        for i, (sid, orient) in enumerate(walk):
            if sid not in segs: return ""
            seq = segs[sid][0]
            if orient == "-": seq = _rc(seq)
            if i == 0:
                out.append(seq)
            else:
                ps, po = walk[i - 1]
                ov = ov_map.get((ps, po, sid, orient), 0)
                out.append(seq[ov:])
        return "".join(out)

    def apply_trim(c: str, new_walk: list[tuple[str, str]]) -> None:
        k = _kof(c)
        if not k: return
        segs, ov_map = load_k(k)
        if not segs: return
        new_seq = reconstruct(new_walk, segs, ov_map)
        if not new_seq: return
        old_len = len(seqs.get(c, ""))
        # Safety: trim must always SHRINK the sequence. Bubble-origin picks
        # are SPAdes contigs already trimmed by _trim_locus to the locus
        # envelope; the underlying GFA segments can be much larger (the
        # contig spans only a portion of them). Reconstructing from a
        # walk subset would inflate beyond the original. Skip in that case.
        if len(new_seq) > old_len:
            print(f"[pick_alleles] {sample}: trim ABORTED on {c}: "
                  f"reconstruct {len(new_seq)} > original {old_len} bp "
                  f"(contig spans a subset of GFA segments)")
            return
        seqs[c] = new_seq
        seg_of_cand_list[c] = new_walk
        print(f"[pick_alleles] {sample}: trimmed outer flanks on {c}: "
              f"{old_len} -> {len(new_seq)} bp ({len(new_walk)} segs)")

    # SINGLETON (only 1 pick): no other arm to anchor a shared-rule trim
    # against, so we use a "drop outermost flank-only fragment" rule. Drop the
    # very-first segment IFF it's flank-only AND another flank-only exists
    # elsewhere in the walk (so we never strip the only flank-anchor). Same
    # logic for the very-last segment.
    if len(picks) == 1:
        c = picks[0]
        walk = seg_of_cand_list.get(c, [])
        if not walk: return
        k = _kof(c)
        if not k: return
        labels = seg_labels_by_k.get(k, {})
        # Aggressive singleton trim: keep only the span bounded by the
        # OUTERMOST flank-tagged segments on each side (any tag containing
        # "flankL" or "flankR" — composite labels included). Drop everything
        # before the leftmost flank-tagged segment AND everything after the
        # rightmost flank-tagged segment, EVEN if those outer segments carry
        # variable-gene tags. Orientation-agnostic.
        flank_positions = [
            i for i, (sid, _) in enumerate(walk)
            if ("flankL" in labels.get(sid, "").split("+")
                or "flankR" in labels.get(sid, "").split("+"))
        ]
        if not flank_positions: return
        left = min(flank_positions)
        right = max(flank_positions)
        new_walk = walk[left:right + 1]
        if new_walk and new_walk != list(walk):
            apply_trim(c, new_walk)
        return

    # PAIR: drop everything OUTSIDE shared flank-only anchors (treating
    # flankL / flankR as one category — position in walk determines the
    # boundary side, so reversed walks work correctly). Trim only within
    # the outer tails (positions before the first variable-gene segment /
    # after the last), so we never cross a var-gene boundary.
    if len(picks) != 2: return
    a, b = picks[0], picks[1]
    wa = seg_of_cand_list.get(a, [])
    wb = seg_of_cand_list.get(b, [])
    ka, kb = _kof(a), _kof(b)
    if ka != kb or not ka: return
    labels = seg_labels_by_k.get(ka, {})
    fa = {sid for sid, _ in wa if is_flank_only(sid, labels)}
    fb = {sid for sid, _ in wb if is_flank_only(sid, labels)}
    shared = fa & fb
    if not shared: return
    for c, walk in ((a, wa), (b, wb)):
        var_idx = [i for i, (sid, _) in enumerate(walk) if has_var(sid, labels)]
        if not var_idx: continue
        first_var, last_var = var_idx[0], var_idx[-1]
        # Pick the INNERMOST shared anchor in each outer tail (closest to the
        # var-gene block). Cutting at the outermost shared anchor would be a
        # no-op when the walk starts/ends ON a shared flank segment — and on
        # complex topologies where the two arms share a long flank chain
        # before diverging (e.g. AG17), we want the trim to drop that whole
        # shared chain, not preserve it.
        left_anchors = [i for i in range(first_var) if walk[i][0] in shared]
        left_cut = left_anchors[-1] if left_anchors else 0
        right_anchors = [i for i in range(last_var + 1, len(walk)) if walk[i][0] in shared]
        right_cut = (right_anchors[0] + 1) if right_anchors else len(walk)
        new_walk = walk[left_cut:right_cut]
        if new_walk and new_walk != list(walk):
            apply_trim(c, new_walk)


def _normalize_orientation(picks: list[str],
                             seqs: dict[str, str],
                             seg_of_cand_list: dict[str, list[tuple[str, str]]],
                             seg_labels_by_k: dict[str, dict[str, str]],
                             sample: str) -> None:
    """For each picked walk, if flankL-only segments sit RIGHTWARD of
    flankR-only segments, the walk is in reverse orientation. Flip the walk
    (reverse order + invert each segment's +/-) and reverse-complement the
    sequence so every emitted allele points flankL → HD → flankR uniformly.
    Modifies seqs and seg_of_cand_list in place.
    """
    from .GFA_search import _rc
    def opp(o): return "-" if o == "+" else "+"
    for c in picks:
        walk = seg_of_cand_list.get(c, [])
        if not walk: continue
        k = _kof(c)
        if not k: continue
        labels = seg_labels_by_k.get(k, {})
        # Use ALL segments carrying the flank tag — including mixed labels
        # like "HD1+flankL" — to compute orientation. A composite segment
        # still carries that flank signal in the walk.
        L_positions, R_positions = [], []
        for i, (sid, _) in enumerate(walk):
            toks = [t for t in labels.get(sid, "").split("+") if t]
            if "flankL" in toks: L_positions.append(i)
            if "flankR" in toks: R_positions.append(i)
        if not L_positions or not R_positions: continue
        if sum(L_positions) / len(L_positions) <= sum(R_positions) / len(R_positions): continue
        new_walk = [(sid, opp(o)) for sid, o in reversed(walk)]
        new_seq = _rc(seqs[c]) if c in seqs else ""
        seg_of_cand_list[c] = new_walk
        if new_seq: seqs[c] = new_seq
        print(f"[pick_alleles] {sample}: flipped {c} "
              f"(flankR→flankL detected; now flankL→flankR)")


def run(sample: str, anchor_contig_fa: str, queries_dir: str, outdir: str,
        known_degHD: str | None,
        expected_count: int, genome_coverage: float,
        tblastn_pid: float = 30.0, tblastn_aa: int = 50,
        flank_pid: float = 85.0, flank_minlen: int = 100,
        degHD_pid: float = 95.0, degHD_minlen: int = 1000,
        max_locus_len: int = 12000, min_allele_len: int = 2000,
        cov_dev_frac: float = 0.5, cov_multicopy_x: float = 1.6,
        cand_ann_tsv: str | None = None,
        spades_dir: str | None = None,
        locus_ref_fa: str | None = None) -> dict:
    """Returns dict with picks list + scoring info.

    `anchor_contig_fa` is the candidate pool (typically cand_combined.fasta — cat of
    anchor_contig.fasta + bubble_alleles.fasta). `cand_ann_tsv` is the matching
    ann.tsv whose `segments` column we propagate into picks so downstream tools
    can use the path directly without re-aligning.
    """
    os.makedirs(outdir, exist_ok=True)
    seqs = read_fasta(anchor_contig_fa)
    # parse the candidate annotation tsv if provided — gives us segments per candidate
    # ann.tsv columns from graph_path_search.py (10 cols):
    #   name, len, k, segments, variable_genes, flankL, flankR, cov, is_degHD, has_repeat
    seg_of_cand: dict[str, str] = {}
    if cand_ann_tsv and os.path.exists(cand_ann_tsv):
        for ln in open(cand_ann_tsv):
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 4: seg_of_cand[f[0]] = f[3]
    # pre-trim oversized contigs to the locus before scoring.
    # `max_locus_len <= 0` is the no-cap sentinel — skip trim entirely.
    proteins = os.path.join(queries_dir, "variable_proteins.fasta")
    flankL  = os.path.join(queries_dir, "flankL.fasta")
    flankR  = os.path.join(queries_dir, "flankR.fasta")
    if max_locus_len > 0:
        for c, s in list(seqs.items()):
            if len(s) > max_locus_len:
                seqs[c] = _trim_locus(s, proteins)
    # NOTE: drop the --min-allele-len hard filter at the picker layer (per user
    # 2026-05-30) — pick_alleles considers all candidates regardless of length;
    # the score tuple's "length similarity" component handles preferences.
    if not seqs:
        print(f"[pick_alleles] {sample}: empty pool"); return dict(picks=[])
    # annotations
    var_per = _vars_in(seqs, seqs, proteins, tblastn_pid, tblastn_aa)
    var_aa_cov = _vars_aa_coverage(seqs, proteins, tblastn_pid, tblastn_aa)
    hL = _flank_hits(seqs, flankL, min_len=flank_minlen, min_pid=flank_pid)
    hR = _flank_hits(seqs, flankR, min_len=flank_minlen, min_pid=flank_pid)
    flank_cov = _flank_coverage_bp(seqs, [flankL, flankR],
                                    min_pid=flank_pid, min_len=flank_minlen)
    deg = _degHD_hits(seqs, known_degHD, min_pid=degHD_pid, min_len=degHD_minlen)
    # Fragment pre-filter has moved upstream: the wrapper now drops
    # flankL-only / flankR-only / pure-junk anchor contigs from
    # `cand_combined.fasta` using `anchor_contig.ann.tsv`. The picker just
    # consumes the filtered pool.
    # how many variable genes total?
    nvar_total = sum(1 for _ in open(proteins) if _.startswith(">"))
    expected_single_cov = (genome_coverage / 2.0) if expected_count == 2 else float(genome_coverage)
    # --- NEW pair-based picker (2026-05-30 rewrite) ---
    #
    # Hard filter (drop pair only if):
    #   * pair is _is_dup (MAFFT HD-core id >= dup_id AND aln_frac >= dup_frac)
    #
    # Sequential pair score (highest tuple wins, Python tuple comparison):
    #   1. var_compl_rank      : variable-gene joint completeness (lo*10+hi); 22 > 12 > 11 > 2 > 1 > 0
    #   2. pair_topology_rank  : closed_bubble (3) > open_bubble (2) > detached (1) > other (0)
    #   3. flank_compl_rank    : (2,2)=22 > (2,1)=12 > (1,1)=11 > (2,0)=2 > (1,0)=1 > (0,0)=0
    #   4. pair_divergence     : -HDcore_alnid (more divergent wins)
    #   5. joint_HD_aa_cov     : sum of HD aa-cov across pair
    #   6. joint_flank_bp_cov  : sum of flank bp-cov across pair
    #   7. -|len_a - len_b|    : prefer similar-length alleles
    #
    # Level 1: per-K, enumerate C(N,2) pairs; drop dups; pick highest-scoring pair.
    # Level 2: across K's, rank each K's winner by the SAME tuple, pick the K.
    # Fallback: if no K yields a valid divergent pair, return 1 allele (best single
    # by topology-implied per-cand surrogate: max(aa_cov + flank_cov) ).

    # candidate name encodes source k as "<sample>__path_<k>_..." (graph_path_search)
    # or "<sample>__bubble_<k>_..." (anchor_search contigs).
    def _k_of(c):
        m = re.search(r"__(?:path|bubble)_(k\d+)_", c); return m.group(1) if m else None

    # Resolve segments per candidate: path-origin from ann.tsv col 3; bubble-origin via contigs.paths.
    paths_cache: dict[str, dict] = {}
    seg_of_cand_list: dict[str, list[tuple[str, str]]] = {}
    for c in seqs:
        listed = _parse_segs_string(seg_of_cand.get(c, ""))
        if listed:
            seg_of_cand_list[c] = listed
        elif "__bubble_" in c:
            seg_of_cand_list[c] = _segments_via_contigs_paths(c, spades_dir, paths_cache)
        else:
            seg_of_cand_list[c] = []

    by_k: dict[str, list[str]] = collections.defaultdict(list)
    for c in seqs:
        k = _k_of(c) or "k?"
        by_k[k].append(c)
    # Per-K segment-label maps (from step 3's seg_labels_<k>.tsv)
    seg_labels_by_k: dict[str, dict[str, str]] = {k: _read_seg_labels(outdir, k) for k in by_k}

    def best_pair_in_k(k: str, cands: list[str]) -> tuple[tuple[str, str], tuple] | None:
        """Enumerate all C(N,2) pairs in cands, drop _is_dup pairs, return the
        (pair, score_tuple) with the highest score. None if no valid pair.

        Completeness-first optimization: var_completeness_rank is just counting
        (cheap — uses cached var_per). Compute it for all pairs first. The score
        tuple has var_compl as its primary key, so a pair can only WIN if no
        other pair has higher var_compl. Run the expensive HD-core MAFFT MSA
        only on the SUBSET of candidates participating in the leading layer.
        Descend through completeness layers if the top layer yields no
        divergent pair (all _is_dup).
        """
        from .pairwise_identity import mafft_pair as _mafft_pair_fallback
        labels = seg_labels_by_k.get(k, {})
        n = len(cands)
        # Phase 1: compute var_completeness rank for every pair (cheap)
        pair_vc: list[tuple[str, str, int]] = []
        for i in range(n):
            for j in range(i + 1, n):
                a, b = cands[i], cands[j]
                vc = _var_completeness_rank(
                    len(var_per.get(a, set())),
                    len(var_per.get(b, set())))
                pair_vc.append((a, b, vc))
        if not pair_vc: return None
        # Phase 2: descend var_compl layers from highest to lowest. For each
        # layer, run ONE MAFFT MSA only on candidates in that layer's pairs,
        # apply _is_dup filter, score remaining pairs, return the best. If a
        # layer has no surviving (non-dup) pair, descend to the next layer.
        layers = sorted({v for _, _, v in pair_vc}, reverse=True)
        for layer in layers:
            cur_pairs = [(a, b) for a, b, v in pair_vc if v == layer]
            subset_names: set[str] = set()
            for a, b in cur_pairs: subset_names.add(a); subset_names.add(b)
            subset_seqs = {c: seqs[c] for c in subset_names}
            hd_pair_id = _hd_core_pairwise_identities(subset_seqs, queries_dir, locus_ref_fa)
            best: tuple[tuple[str, str], tuple] | None = None
            for a, b in cur_pairs:
                key = (a, b) if a < b else (b, a)
                if key in hd_pair_id:
                    alnid, aln_frac = hd_pair_id[key]   # HD-core (preferred)
                else:
                    try:
                        alnid, aln_frac = _mafft_pair_fallback(seqs[a], seqs[b])
                    except Exception:
                        alnid, aln_frac = 0.0, 0.0
                # Hard filter: drop pair if it's a duplicate.
                if alnid >= DUP_PID and aln_frac >= DUP_AF: continue
                sc = _score_pair(a, b, seqs, seg_of_cand_list, labels,
                                   hL, hR, var_per, var_aa_cov, flank_cov,
                                   pair_alnid=alnid)
                if best is None or sc > best[1]: best = ((a, b), sc)
            if best is not None: return best
        return None

    pair_by_k: dict[str, tuple[tuple[str, str], tuple]] = {}
    for k, kc in by_k.items():
        pair = best_pair_in_k(k, kc)
        if pair is not None: pair_by_k[k] = pair

    picks: list[str] = []
    best_k_choice: str | None = None
    if pair_by_k:
        # Level 2: pick the K whose pair has the highest score tuple
        best_k_choice = max(pair_by_k, key=lambda k: pair_by_k[k][1])
        picks = list(pair_by_k[best_k_choice][0])
        # Score tuple: (var_compl, topo, flank_compl, -alnid, joint_aa, joint_fl, -len_diff)
        print(f"[pick_alleles] {sample}: per-k best pair scores: "
              + "; ".join(
                  f"{k}=var{pair_by_k[k][1][0]}/topo{pair_by_k[k][1][1]}/flank{pair_by_k[k][1][2]}/div{-pair_by_k[k][1][3]:.3f}"
                  for k in pair_by_k)
              + f"  -> chose {best_k_choice} (var_compl={pair_by_k[best_k_choice][1][0]}, "
                f"topo={pair_by_k[best_k_choice][1][1]}, "
                f"flank_compl={pair_by_k[best_k_choice][1][2]}, "
                f"div_HDcore={-pair_by_k[best_k_choice][1][3]:.3f})")
    elif expected_count == 1 and seqs:
        # Haploid mode: single allele requested
        picks = [max(seqs.keys(),
                     key=lambda c: var_aa_cov.get(c, 0) + flank_cov.get(c, 0))]
        print(f"[pick_alleles] {sample}: haploid mode -> 1 allele picked")
    elif seqs:
        # Fallback: no valid divergent pair in any K. Return 1 allele —
        # the highest-scoring SINGLE by (HD aa-cov + flank bp-cov) as surrogate.
        picks = [max(seqs.keys(),
                     key=lambda c: var_aa_cov.get(c, 0) + flank_cov.get(c, 0))]
        print(f"[pick_alleles] {sample}: no divergent pair in any K -> "
              f"fallback to single allele")
    else:
        print(f"[pick_alleles] {sample}: empty pool -> no picks")

    # Drop outer flanks (2026-05-31). Pair case: trim each arm to the span
    # between the leftmost / rightmost SHARED flank-only segment IDs.
    # Singleton case: drop the very-first segment if it's flank-only and the
    # very-last if it's flank-only. Operates on the picked records only.
    if picks and spades_dir:
        _trim_outer_flanks(picks, seqs, seg_of_cand_list, seg_labels_by_k, spades_dir, sample)
    # Normalize orientation: if a walk's flankL-only segments sit RIGHTWARD of
    # its flankR-only segments, reverse the walk + RC the sequence so every
    # allele points flankL → HD → flankR uniformly across samples.
    if picks:
        _normalize_orientation(picks, seqs, seg_of_cand_list, seg_labels_by_k, sample)
    # emit
    fa = os.path.join(outdir, "primary_alleles.fasta")
    tsv = os.path.join(outdir, "picks.tsv")
    with open(fa, "w") as o_fa, open(tsv, "w") as o_tsv:
        o_tsv.write("sample\tallele\torigin\tk\ttype\tlen\tfrom_contig\tsegments\tcov\tn_variable_genes\thas_both_flanks\tis_degHD\n")
        for i, c in enumerate(picks, 1):
            s = seqs[c]
            nvars = len(var_per.get(c, set()))
            is_degHD = c in deg
            has_all_vars = (nvars == nvar_total)
            has_both_flanks = (c in hL) and (c in hR)
            # eased: "complete" iff all variable genes present (flanks NOT required); the score tuple
            # uses (has_all_vars AND has_both_flanks) as a higher sub-tier so a flanked complete still wins
            typ = ("element-like" if is_degHD
                   else "complete" if has_all_vars
                   else "partial")
            origin = ("path" if "__path_" in c else "bubble")
            k = _kof(c) or "-"
            # `seg_of_cand_list[c]` is the canonical list-of-tuples walk —
            # populated for path-origin from cand_combined.ann.tsv col-3 and
            # for bubble-origin via contigs.paths. The trim + orientation
            # normalization steps update this in place, so we always re-serialize
            # FROM it to make picks.tsv consistent with the emitted sequence.
            # Fall back to the raw `seg_of_cand` string if the list is empty.
            if seg_of_cand_list.get(c):
                segs = ",".join(f"{sid}{ori}" for sid, ori in seg_of_cand_list[c])
            else:
                segs = seg_of_cand.get(c, "-")
            hdr = (f"Pcub_{sample}_allele{i} sample={sample} allele={i} origin={origin} k={k} "
                   f"type={typ} len={len(s)} from={c} cov={_cov(c):.1f}")
            o_fa.write(f">{hdr}\n")
            for j in range(0, len(s), 80): o_fa.write(s[j:j+80] + "\n")
            o_tsv.write(f"{sample}\tallele{i}\t{origin}\t{k}\t{typ}\t{len(s)}\t{c}\t{segs}\t{_cov(c):.1f}\t{nvars}/{nvar_total}\t{has_both_flanks}\t{is_degHD}\n")
    print(f"[pick_alleles] {sample}: picked {len(picks)} allele(s) -> {fa}")
    return dict(picks=[(seqs[c], c) for c in picks], picks_names=picks, fasta=fa, tsv=tsv,
                var_per=var_per, hL=hL, hR=hR, deg=deg)

def _cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--sample", required=True)
    p.add_argument("--candidates", required=True)
    p.add_argument("--queries-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--known-degHD", default=None)
    p.add_argument("--expected-count", type=int, choices=[1, 2], required=True)
    p.add_argument("--genome-coverage", type=float, default=0)
    p.add_argument("--tblastn-pid", type=float, default=30.0)
    p.add_argument("--tblastn-aa", type=int, default=50)
    p.add_argument("--flank-pid", type=float, default=85.0)
    p.add_argument("--flank-minlen", type=int, default=100)
    p.add_argument("--degHD-id", dest="degHD_pid", type=float, default=95.0)
    p.add_argument("--degHD-minlen", type=int, default=1000)
    p.add_argument("--max-locus-len", type=int, default=12000)
    p.add_argument("--min-allele-len", type=int, default=2000)
    p.add_argument("--cand-ann", default=None,
                   help="bubble_alleles.ann.tsv from segment_alleles.py; gives us the segment "
                        "path per candidate so picks.tsv can record it")
    p.add_argument("--spades-dir", default=None,
                   help="parent dir of per-k SPAdes outputs. Used to resolve "
                        "bubble-origin (whole-contig) candidates' segments via "
                        "<spades-dir>/<k>/contigs.paths so pair-topology "
                        "classification has the segment list for both origins.")
    p.add_argument("--locus-ref", default=None,
                   help="locus reference fasta. When supplied, pair divergence "
                        "is computed on HD-core columns (one MAFFT MSA per K, "
                        "matching step 3's cluster_alleles metric) instead of "
                        "whole-allele MAFFT.")
    a = p.parse_args(argv)
    run(a.sample, a.candidates, a.queries_dir, a.outdir, a.known_degHD,
        a.expected_count, a.genome_coverage, a.tblastn_pid, a.tblastn_aa,
        a.flank_pid, a.flank_minlen, a.degHD_pid, a.degHD_minlen,
        a.max_locus_len, a.min_allele_len, cand_ann_tsv=a.cand_ann,
        spades_dir=a.spades_dir, locus_ref_fa=a.locus_ref)

if __name__ == "__main__":
    _cli()
