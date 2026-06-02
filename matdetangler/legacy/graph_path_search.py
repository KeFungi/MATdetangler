"""Step 3 — graph path search + bubble typing.

Always runs. Consumes anchor records produced by step 2 (`anchor_search`) and walks
the GFA from them to enumerate candidate allele paths.

Anchor inputs (one or both):
    --anchors-contigs  anchor_contig.fasta     from anchor_search --source contigs
    --anchors-segments anchor_segments.fasta   from anchor_search --source segments

The two fastas use distinct record-name templates so we can split them per k and
per kind without ambiguity:
    <sample>__bubble_<k>_<contig_id>   = a contig anchor (must be mapped to GFA segments)
    <sample>__seg_<k>_<segment_id>     = a GFA segment anchor (already an S-line id)

Per k:
  1. Collect this-k anchors from both inputs (contigs need to be mapped to GFA segs
     via SPAdes' contigs.paths primary + blastn-vs-segments fallback; segment anchors
     are already segment IDs).
  2. BFS-expand the GFA neighborhood from the union of anchor segments.
  3. Label every neighborhood segment (tblastn variable_proteins + blastn flankL /
     flankR / repeats / degHD).
  4. Enumerate simple paths from any flankL-bearing segment to any flankR-bearing
     segment inside the neighborhood, bounded by --max-locus-len / --max-nodes /
     --max-paths.
  5. Completeness filter: keep paths that cover every variable gene AND include both
     flanks.
  6. If genome coverage was supplied, a two-pass DFS deprioritizes high-coverage
     segments (likely repeats / paralogs / SD-collapsed) in pass 1; pass 2 fallback
     allows them back in if pass 1 produced zero complete paths.
  7. Widen the BFS by +1 hop and retry up to --max-hops.

Output per call:
    bubble_alleles.fasta + bubble_alleles.ann.tsv

Return value: number of complete candidate paths written across all k's. The orchestrator
uses this to decide whether to fall back to step 2.2 (segment anchors) and re-invoke
step 3.
"""
from __future__ import annotations
import os, sys, argparse, subprocess, tempfile, collections, re
from . import blast_utils as bu
from .paths import spades_k_paths
from .GFA_search import (
    _rc, parse_segments, parse_links, build_directed_adj,
    reconstruct, path_depth, enumerate_paths, label_segments,
)
from .pairwise_identity import mafft_pair, mafft_pair_core, detect_core_span


# ---------- anchor input parsing ----------

_RX_BUBBLE = re.compile(r"^.+?__bubble_([^_]+)_(.+)$")   # → (k, contig_id)
_RX_SEG    = re.compile(r"^.+?__seg_([^_]+)_(.+)$")      # → (k, segment_id)


def _read_fasta_names(fa: str) -> list[str]:
    if not fa or not os.path.exists(fa) or os.path.getsize(fa) == 0:
        return []
    out: list[str] = []
    for ln in open(fa):
        if ln.startswith(">"): out.append(ln[1:].split()[0])
    return out


def _read_ann_hd_carriers(ann_tsv: str | None) -> set[str]:
    """Return the set of anchor record names whose ann.tsv `hd_genes` column is
    non-"-" — i.e. records that carry at least one variable gene. Returns an empty
    set if the ann.tsv is missing or empty.
    Schema (cand_*.ann.tsv): name, len, k, kind, hd_genes, flanks
    """
    out: set[str] = set()
    if not ann_tsv or not os.path.exists(ann_tsv): return out
    for ln in open(ann_tsv):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 5: continue
        if f[4] and f[4] != "-": out.add(f[0])
    return out


def collect_anchors_per_k(anchors_contigs_fa: str | None,
                          anchors_segments_fa: str | None,
                          seeds_from: str = "hd"
                          ) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Group anchor record names by k.

    Returns ({k: contig_anchor_names}, {k: segment_anchor_ids}).
    `contig_anchor_names` use SPAdes contig names (matches the contigs.fasta header
    and contigs.paths keys). `segment_anchor_ids` are GFA S-line ids directly.

    seeds_from = "hd" (default): filter to records whose ann.tsv `hd_genes` column
                                  is non-"-" (i.e. records that carry at least one
                                  variable gene). Flank-only records are dropped —
                                  the BFS reaches flank-bearing segs from HD seeds
                                  via L-line adjacency anyway, and pure-flank big
                                  contigs would otherwise pull in noise.
    seeds_from = "all":           keep every record. Use under --debug.
    """
    contig_per_k: dict[str, set[str]] = collections.defaultdict(set)
    seg_per_k:    dict[str, set[str]] = collections.defaultdict(set)
    # Build the per-input HD-carrier filter if requested.
    hd_contigs: set[str] = set()
    hd_segs:    set[str] = set()
    if seeds_from == "hd":
        if anchors_contigs_fa:
            ann_ctg = anchors_contigs_fa.rsplit(".fasta", 1)[0] + ".ann.tsv"
            hd_contigs = _read_ann_hd_carriers(ann_ctg)
        if anchors_segments_fa:
            ann_seg = anchors_segments_fa.rsplit(".fasta", 1)[0] + ".ann.tsv"
            hd_segs = _read_ann_hd_carriers(ann_seg)
    for nm in _read_fasta_names(anchors_contigs_fa):
        if seeds_from == "hd" and nm not in hd_contigs: continue
        m = _RX_BUBBLE.match(nm)
        if m: contig_per_k[m.group(1)].add(m.group(2))
    for nm in _read_fasta_names(anchors_segments_fa):
        if seeds_from == "hd" and nm not in hd_segs: continue
        m = _RX_SEG.match(nm)
        if m: seg_per_k[m.group(1)].add(m.group(2))
    return dict(contig_per_k), dict(seg_per_k)


# ---------- contigs.paths and contig→segment mapping (kept) ----------

_PATHS_NAME_RX = re.compile(r"^(.*_cov_[\d.]+)(?:_\d+)?'?$")


def parse_contigs_paths(paths_file: str) -> dict[str, list[tuple[str, str]]]:
    """SPAdes contigs.paths formats — accept both:
        Modern (3.15+):
            NODE_<id>_length_<L>_cov_<cov>          (forward walk follows)
            <seg_id><orient>,<seg_id><orient>,...
            NODE_<id>_length_<L>_cov_<cov>'         (RC marker `'`; reverse walk follows)
            <seg_id><orient>,...
        Older:
            NODE_<id>_length_<L>_cov_<cov>_<component>
            <seg_id><orient>,...
        (blank line between records)

    Returns {contig_base_name: [(seg_id, orient), ...]}. RC walks are ignored
    (the forward walk plus orientations is enough). For older multi-component
    contigs, walks from all components concatenate under the base name.
    """
    out: dict[str, list[tuple[str, str]]] = collections.defaultdict(list)
    if not paths_file or not os.path.exists(paths_file): return {}
    cur: str | None = None
    cur_is_rc: bool = False
    for ln in open(paths_file):
        ln = ln.rstrip("\n")
        if not ln: cur = None; continue
        if ln.startswith("NODE_") or ln.startswith("EDGE_"):
            cur_is_rc = ln.endswith("'")
            m = _PATHS_NAME_RX.match(ln)
            cur = m.group(1) if m else (ln[:-1] if cur_is_rc else ln)
        elif cur is not None and not cur_is_rc:
            for tok in ln.replace(";", ",").split(","):
                tok = tok.strip()
                if not tok: continue
                m = re.match(r"^(\d+)([+-])?$", tok)
                if m: out[cur].append((m.group(1), m.group(2) or "+"))
    return dict(out)




def map_anchors_to_segments(anchors: set[str], contigs_paths_file: str | None,
                            contigs_fa: str, gfa_segs: dict[str, tuple[str, float]],
                            threads: int = 1) -> dict[str, list[tuple[str, str]]]:
    """For each anchor contig, return its list of (seg_id, orient) from SPAdes'
    contigs.paths. Contigs not found in contigs.paths are silently dropped
    (a warning is printed); they cannot be mapped to oriented segment walks
    without contigs.paths. SPAdes emits contigs.paths alongside contigs.fasta
    by default, so this is the normal case.

    `contigs_fa`, `gfa_segs`, `threads` are kept for API symmetry but unused
    now that the blastn fallback has been retired.
    """
    out: dict[str, list[tuple[str, str]]] = {}
    if not contigs_paths_file:
        print(f"[map_anchors_to_segments] no contigs.paths file -> 0 of "
              f"{len(anchors)} anchors mapped")
        return out
    paths_map = parse_contigs_paths(contigs_paths_file)
    missing: list[str] = []
    for c in anchors:
        if c in paths_map: out[c] = paths_map[c]; continue
        base = c.rsplit("_", 1)[0]
        if base in paths_map: out[c] = paths_map[base]; continue
        missing.append(c)
    if missing:
        print(f"[map_anchors_to_segments] {len(missing)} of {len(anchors)} "
              f"anchor contigs missing from contigs.paths -> dropped: "
              f"{missing[:5]}{'...' if len(missing) > 5 else ''}")
    return out


# ---------- BFS / path helpers ----------

def undirected_adj_from_links(links) -> dict[str, set[str]]:
    adj: dict[str, set[str]] = collections.defaultdict(set)
    for s1, _, s2, _, _ in links:
        adj[s1].add(s2); adj[s2].add(s1)
    return adj


def bfs_expand(adj_und: dict[str, set[str]], seeds: set[str], hops: int,
               repeat_segs: set[str] | None = None) -> set[str]:
    """Hop-limited BFS over undirected adjacency.

    Asymmetric repeat rule:
        - A NON-repeat segment in the frontier expands normally — all its
          neighbors (repeat OR non-repeat) are added to the visited set.
        - A REPEAT segment in the frontier does NOT expand — none of its
          neighbors are added.

    Effect: paths CAN enter a repeat segment (it shows up in the neighborhood
    via a non-repeat neighbor), but the repeat doesn't drag in its many
    neighbors and blow up the search space at high-copy regions.
    """
    repeat_segs = repeat_segs or set()
    visited = set(seeds); frontier = set(seeds)
    for _ in range(max(0, hops)):
        nxt: set[str] = set()
        for s in frontier:
            if s in repeat_segs: continue   # absorb but don't expand
            nxt |= adj_und.get(s, set())
        nxt -= visited
        if not nxt: break
        visited |= nxt; frontier = nxt
    return visited


def restrict_adj_to_subgraph(adj, keep: set[str]):
    out: dict[tuple[str, str], list[tuple[str, str, int]]] = {}
    for (s, o), edges in adj.items():
        if s not in keep: continue
        kept = [(s2, o2, ov) for (s2, o2, ov) in edges if s2 in keep]
        if kept: out[(s, o)] = kept
    return out


def classify_neighborhood_topology(
        nhood: set[str], adj_und: dict[str, set[str]],
        segs: dict[str, tuple[str, float]], labels: dict[str, str],
        var_per: dict[str, set[str]], min_core_len: int = 1000) -> dict:
    """Bubble topology on the BFS NEIGHBORHOOD (not the whole GFA).

    Same classification scheme as the old bubble_topo.py, but it reuses the
    already-computed BFS neighborhood + segment labels from run_one_k. No extra
    tblastn.

    Returns dict with: type, n_main_seg, n_path, n_shared_anchor,
    n_shared_flank_anchor, main_lens, main_ids.

    Classification:
        n_paths == 0                              -> no_main
        n_paths == 1                              -> single        (haploid / collapsed)
        n_paths == 2  + shared flank anchors >= 2 -> closed_bubble (one flankL + one flankR)
        n_paths == 2  + shared flank anchor == 1  -> open_bubble
        n_paths == 2  + 0 shared flank anchors    -> detached
        n_paths >  2                              -> complexed
    """
    main_set = {sid for sid in nhood
                if var_per.get(sid) and len(segs[sid][0]) >= min_core_len}
    if not main_set:
        return dict(type="no_main", n_main_seg=0, n_path=0,
                    n_shared_anchor=0, n_shared_flank_anchor=0,
                    main_lens=[], main_ids=[])
    flank_segs = {sid for sid in nhood
                  if "flankL" in labels.get(sid, "").split("+")
                  or "flankR" in labels.get(sid, "").split("+")}
    # Connected components within main_set, induced by L-line adjacency
    seen, comps = set(), []
    for m in main_set:
        if m in seen: continue
        comp, stack = set(), [m]
        while stack:
            x = stack.pop()
            if x in seen: continue
            seen.add(x); comp.add(x)
            for y in adj_und.get(x, set()) & main_set:
                if y not in seen: stack.append(y)
        comps.append(comp)
    n = len(comps); nm = len(main_set)
    main_ids  = sorted(main_set)
    main_lens = sorted([len(segs[m][0]) for m in main_set], reverse=True)
    if n == 1:
        return dict(type="single", n_main_seg=nm, n_path=1, n_shared_anchor=0,
                    n_shared_flank_anchor=0, main_lens=main_lens, main_ids=main_ids)
    if n > 2:
        return dict(type="complexed", n_main_seg=nm, n_path=n,
                    n_shared_anchor=0, n_shared_flank_anchor=0,
                    main_lens=main_lens, main_ids=main_ids)
    # n == 2: count shared external anchors and shared flank anchors.
    def ext(p): return set().union(*(adj_und.get(x, set()) for x in p)) - main_set
    shared = ext(comps[0]) & ext(comps[1])
    ns = len(shared)
    shared_flank = shared & flank_segs
    nsf = len(shared_flank)
    t = ("closed_bubble" if nsf >= 2 else
         "open_bubble"   if nsf == 1 else
         "detached")
    return dict(type=t, n_main_seg=nm, n_path=2, n_shared_anchor=ns,
                n_shared_flank_anchor=nsf, main_lens=main_lens, main_ids=main_ids)


def is_complete_path(path, labels: dict[str, str],
                     var_per: dict[str, set[str]], nvar_total: int) -> bool:
    seen_genes: set[str] = set()
    sL = sR = False
    for sid, _, _ in path:
        seen_genes |= var_per.get(sid, set())
        toks = labels.get(sid, "").split("+")
        if "flankL" in toks: sL = True
        if "flankR" in toks: sR = True
    return (len(seen_genes) == nvar_total) and sL and sR


# ---------- per-k driver ----------

def run_one_k(sample: str, k: str, gfa: str, contigs_fa: str, queries_dir: str, outdir: str,
              anchor_contigs: set[str], anchor_segments_direct: set[str],
              repeats: str | None, known_degHD: str | None,
              max_locus_len: int, min_allele_len: int,
              max_nodes: int, max_paths: int,
              max_walk_bp: int = 0,
              init_hops: int = 5, max_hops: int = 10,
              threads: int = 1,
              genome_cov: float = 0.0, cov_repeat_factor: float = 2.0,
              min_core_len: int = 1000,
              asymmetric_bfs: bool = False,
              expected_count: int = 2,
              dup_id: float = 0.95,
              dup_frac: float = 0.80,
              locus_ref_fa: str | None = None,
              blast_out_dir: str | None = None) -> tuple[str, str, int, dict]:
    """Return (fa_path, tsv_path, n_complete_paths_written, topology_dict).
    n_complete_paths_written = 0 signals the caller can fall back to step 2.2.

    BFS widens hop-by-hop. Each hop runs enumerate_paths (which already
    deduplicates RC mirrors at emit time via path-topology canonicalization),
    filters to is_complete_path, then clusters the complete paths via MAFFT
    on the HD-core columns at id >= dup_id AND aln-fraction >= dup_frac.
    The loop breaks when MAFFT clusters >= expected_count truly distinct
    alleles (decoration-variant paths collapse into one cluster, so we never
    quit on false-distinct paths), or when max_hops is hit.
    """
    fa_path  = os.path.join(outdir, f"bubble_alleles_{k}.fasta")
    tsv_path = os.path.join(outdir, f"bubble_alleles_{k}.ann.tsv")
    open(fa_path, "w").close(); open(tsv_path, "w").close()
    if not anchor_contigs and not anchor_segments_direct:
        print(f"  [{k}] no anchors supplied for this k -> skipping")
        return fa_path, tsv_path, 0, dict(type="no_main", n_main_seg=0, n_path=0, n_shared_anchor=0, n_shared_flank_anchor=0, main_lens=[], main_ids=[])
    segs = parse_segments(gfa)
    links = parse_links(gfa)
    adj_dir = build_directed_adj(links)
    adj_und = undirected_adj_from_links(links)
    # Map contig anchors to GFA segments (contigs.paths or blastn fallback).
    spades_k_dir = os.path.dirname(gfa)
    paths_file = None
    for cand in ("contigs.paths", "final_contigs.paths", "scaffolds.paths"):
        p = os.path.join(spades_k_dir, cand)
        if os.path.exists(p): paths_file = p; break
    contig_to_segs: dict[str, list[tuple[str, str]]] = {}
    if anchor_contigs:
        contig_to_segs = map_anchors_to_segments(
            anchor_contigs, paths_file, contigs_fa, segs, threads=threads)
    anchor_segs_from_contigs = {sid for lst in contig_to_segs.values() for sid, _ in lst}
    # Direct segment anchors (already GFA seg ids).
    anchor_segs_direct = {sid for sid in anchor_segments_direct if sid in segs}
    anchor_segs = anchor_segs_from_contigs | anchor_segs_direct
    print(f"  [{k}] anchors: {len(anchor_contigs)} contig + {len(anchor_segments_direct)} segment "
          f"-> {len(anchor_segs)} GFA segments "
          f"(from contigs: {len(anchor_segs_from_contigs)} via "
          f"{'contigs.paths' if paths_file else 'blastn fallback'}; direct: {len(anchor_segs_direct)})")
    if not anchor_segs:
        print(f"  [{k}] no GFA segments resolved from the supplied anchors -> skipping")
        return fa_path, tsv_path, 0, dict(type="no_main", n_main_seg=0, n_path=0, n_shared_anchor=0, n_shared_flank_anchor=0, main_lens=[], main_ids=[])
    # Steps 3-6 iteration
    proteins = os.path.join(queries_dir, "variable_proteins.fasta")
    nvar_total = sum(1 for ln in open(proteins) if ln.startswith(">"))
    # Repeat-aware BFS: coverage-based detection only.
    # We used to also blastn `--repeats` against ALL segments (`all_seg_db`),
    # but that required building a 100k-300k segment DB per k for marginal
    # benefit — coverage > cov_repeat_factor × genome_cov already catches the
    # high-copy hubs that drive the BFS blow-up. Removed for efficiency.
    high_cov_segs: set[str] = set()
    if genome_cov > 0:
        cutoff = cov_repeat_factor * genome_cov
        high_cov_segs = {sid for sid, (_, dp) in segs.items() if dp > cutoff}
    repeat_segs = high_cov_segs
    rule_msg = ("ON: BFS will absorb but not expand from repeats"
                if asymmetric_bfs else
                "OFF (default): symmetric BFS — repeats expand normally")
    print(f"  [{k}] repeats detected: {len(repeat_segs)} segments "
          f"(coverage-based; sequence-based blastn vs --repeats was retired)  "
          f"asymmetric_bfs={rule_msg}")
    hops = init_hops
    complete_paths: list = []
    paths: list = []
    labels: dict[str, str] = {}
    var_per: dict[str, set[str]] = collections.defaultdict(set)
    bfs_repeats_arg = repeat_segs if asymmetric_bfs else None
    # Auto-detect step 2.2 output (anchor_segments_<k>.ann.tsv) and reuse its HD
    # labels in label_segments — step 2.2 already tblastn'd HD-proteins against
    # every GFA segment for this k, so we should never re-blast in step 3.
    hd_labels_step22: dict[str, set[str]] = {}
    seg_hd_ann = os.path.join(outdir, f"anchor_segments_{k}.ann.tsv")
    if os.path.exists(seg_hd_ann):
        for ln in open(seg_hd_ann):
            f = ln.rstrip("\n").split("\t")
            # anchor_segments ann schema: name, len, k, kind, hd_genes, flanks
            if len(f) < 5: continue
            nm = f[0]
            # name looks like "<sample>__seg_<k>_<segment_id>" — recover segment_id
            m = re.match(r".*__seg_k\d+_(.+)$", nm)
            sid = m.group(1) if m else nm
            hd_genes_str = f[4]
            if hd_genes_str and hd_genes_str != "-":
                hd_labels_step22[sid] = set(g for g in hd_genes_str.split(",") if g)
        print(f"  [{k}] using step 2.2 HD-labels: {len(hd_labels_step22)} segments "
              f"have precomputed HD-gene tags (skipping HD-tblastn re-blast in label_segments)")
    hd_arg = hd_labels_step22 if hd_labels_step22 else None
    while hops <= max_hops:
        nhood = bfs_expand(adj_und, anchor_segs, hops, repeat_segs=bfs_repeats_arg)
        nhood_segs = {sid: segs[sid] for sid in nhood if sid in segs}
        labels, var_per = label_segments(nhood_segs, queries_dir,
                                          repeats=repeats, known_degHD=known_degHD,
                                          threads=threads,
                                          hd_labels_from_step22=hd_arg,
                                          blast_out_dir=blast_out_dir,
                                          blast_tag=f"nhood_{k}_hops{hops}")
        starts = {sid for sid, l in labels.items() if "flankL" in l.split("+")}
        ends   = {sid for sid, l in labels.items() if "flankR" in l.split("+")}
        var_segs = {sid for sid in nhood if var_per.get(sid)}
        if not (starts and ends and var_segs):
            print(f"  [{k}] hops={hops}: |nhood|={len(nhood)}, missing "
                  f"{'flankL ' if not starts else ''}{'flankR ' if not ends else ''}"
                  f"{'var ' if not var_segs else ''}-> widening")
            hops += 1; continue
        # Single-pass path enumeration over the full neighborhood. Paths CAN
        # traverse repeat-labeled segments (rule 2). Repeats can't blow up the
        # search because they didn't pull in their own neighbors during BFS
        # (rule 3).
        adj_sub = restrict_adj_to_subgraph(adj_dir, nhood)
        paths = enumerate_paths(adj_sub, nhood_segs, starts, ends, var_segs,
                                  max_bp=max_walk_bp, max_nodes=max_nodes,   # --max-walk-bp is the DFS cap; --max-locus-len is for pick_alleles trim only
                                  max_paths=max_paths)
        complete_paths = [p for p in paths if is_complete_path(p, labels, var_per, nvar_total)]
        # MAFFT-core distinctness is the sole stopping criterion. RC mirrors are
        # already deduped at emission time inside enumerate_paths.
        n_repeat_in_nhood = len(repeat_segs & nhood)
        expand_note = ("none expanded" if asymmetric_bfs else "expanded normally")
        if not complete_paths:
            print(f"  [{k}] hops={hops}: |nhood|={len(nhood)} "
                  f"(includes {n_repeat_in_nhood} repeat segs, {expand_note}) -> "
                  f"{len(paths)} paths, 0 complete (need {expected_count})")
            hops += 1; continue
        recs = [(f"p{i}", reconstruct(p, segs)) for i, p in enumerate(complete_paths)]
        clusters = cluster_alleles(recs, queries_dir=queries_dir,
                                     locus_ref_fa=locus_ref_fa,
                                     id_thresh=dup_id, frac_thresh=dup_frac)
        n_distinct = len(clusters)
        # Classify topology on the current neighborhood — also break early on
        # closed_bubble even if n_distinct < expected_count (the bubble has
        # converged; widening further won't change biology).
        topology_now = (classify_neighborhood_topology(
                          set(nhood), adj_und, segs, labels, var_per,
                          min_core_len=min_core_len)
                        if labels else dict(type="no_main"))
        print(f"  [{k}] hops={hops}: |nhood|={len(nhood)} "
              f"(includes {n_repeat_in_nhood} repeat segs, {expand_note}) -> "
              f"{len(paths)} paths, {len(complete_paths)} complete, "
              f"{n_distinct} MAFFT-core distinct "
              f"@ id>={dup_id:.2f},frac>={dup_frac:.2f} "
              f"(need {expected_count})  topology={topology_now['type']}")
        # Break when `n_distinct >= expected_count` AND topology has converged
        # to either `closed_bubble` (the ideal) or `complexed` (the bubble has
        # settled into a shape that won't change with more BFS — main segs
        # don't share flank anchors, and widening just adds noise). For other
        # topologies (no_main, etc.) keep widening until max_hops.
        if n_distinct >= expected_count and topology_now['type'] in ("closed_bubble", "complexed"):
            break
        hops += 1
    # Final widen-loop stats — visible at end-of-step-3 even when buffered logs
    # haven't flushed the per-hop lines yet.
    final_nhood = len(nhood) if 'nhood' in dir() else 0
    print(f"  [{k}] widen-loop stats: last_hops={hops}, last_nhood={final_nhood}, "
          f"max_hops={max_hops}")
    # Persist per-segment labels from the final widen-loop iteration so the
    # downstream pick_alleles step can classify pair-level topology without
    # re-blasting. Format: seg_id<TAB>label_str (label_str is "+"-joined feature
    # tags like "flankL", "flankR+HD1", "HD2"; empty string = unlabeled).
    # label_segments now uses min_seg_len_for_label=0 by default so every seg in
    # the BFS nhood gets a label entry — no anchor-path patch needed.
    if labels:
        seg_labels_tsv = os.path.join(outdir, f"seg_labels_{k}.tsv")
        with open(seg_labels_tsv, "w") as o:
            for sid in sorted(labels):
                o.write(f"{sid}\t{labels[sid]}\n")
        print(f"  [{k}] wrote per-segment labels -> {seg_labels_tsv} "
              f"({len(labels)} segments)")
    if os.environ.get("HDD_DEBUG_BFS"):
        nhood_members = sorted(nhood) if 'nhood' in dir() else []
        print(f"  [{k}] DEBUG nhood members ({len(nhood_members)}): "
              f"{','.join(nhood_members)}")
    n_complete_written = len(complete_paths)
    if not complete_paths:
        complete_paths = paths  # best-effort fallback
    # Stage 1 dedup: canonical fwd/rc sequence — collapses byte-equal walks
    # (and walk + its reverse). Keeps the LONGER path per canonical seq.
    canon: dict[str, list] = {}
    for p in complete_paths:
        s = reconstruct(p, segs); key = min(s, _rc(s))
        if key not in canon or len(p) > len(canon[key]): canon[key] = p
    # Stage 2 dedup: MAFFT-core clustering across the canonical-distinct survivors.
    # Two paths whose HD-core MAFFT identity is >= dup_id AND aligned fraction
    # is >= dup_frac are the SAME allele (decoration variants — different walks
    # through repeat segments or different overlap choices, same biological allele).
    # Emit ONE representative per cluster: the LONGEST canonical sequence in the
    # cluster (most graph structure captured).
    canon_recs = [(f"c{i}", reconstruct(p, segs), p) for i, p in enumerate(canon.values())]
    if len(canon_recs) > 1:
        named_seqs = [(n, s) for n, s, _ in canon_recs]
        clusters = cluster_alleles(named_seqs, queries_dir=queries_dir,
                                     locus_ref_fa=locus_ref_fa,
                                     id_thresh=dup_id, frac_thresh=dup_frac)
        path_of = {n: p for n, _, p in canon_recs}
        seq_of  = {n: s for n, s, _ in canon_recs}
        seen: dict[str, list] = {}
        for cluster in clusters:
            rep_name = max(cluster, key=lambda n: len(seq_of[n]))
            rep_path = path_of[rep_name]
            rep_seq  = seq_of[rep_name]
            seen[min(rep_seq, _rc(rep_seq))] = rep_path
        print(f"  [{k}] emission dedup: {len(complete_paths)} paths -> "
              f"{len(canon)} canonical -> {len(clusters)} MAFFT-clusters; "
              f"writing {len(seen)} representative(s)")
    else:
        seen = canon
        if canon:
            print(f"  [{k}] emission dedup: {len(complete_paths)} paths -> "
                  f"{len(canon)} canonical (single — MAFFT clustering skipped)")
    n_written = 0
    with open(fa_path, "w") as fa, open(tsv_path, "w") as tsv:
        for p in seen.values():
            seq = reconstruct(p, segs); L = len(seq)
            if L < min_allele_len: continue
            seg_str = ",".join(f"{sid}{o}" for sid, o, _ in p)
            var_hits = set().union(*(var_per.get(sid, set()) for sid, _, _ in p))
            has_flankL = any("flankL" in labels.get(sid, "").split("+") for sid, _, _ in p)
            has_flankR = any("flankR" in labels.get(sid, "").split("+") for sid, _, _ in p)
            is_degHD   = any(labels.get(sid) == "degHD" for sid, _, _ in p)
            has_repeat = any("repeat" in labels.get(sid, "").split("+") for sid, _, _ in p)
            dp = path_depth(p, segs)
            nm = f"{sample}__path_{k}_n{len(p)}_L{L}_d{dp:.1f}_p{n_written}"
            fa.write(f">{nm}\n")
            for i in range(0, L, 80): fa.write(seq[i:i + 80] + "\n")
            tsv.write(f"{nm}\t{L}\t{k}\t{seg_str}\t{','.join(sorted(var_hits)) or '-'}"
                      f"\t{'T' if has_flankL else 'F'}\t{'T' if has_flankR else 'F'}"
                      f"\t{dp:.2f}\t{'T' if is_degHD else 'F'}\t{'T' if has_repeat else 'F'}\n")
            n_written += 1
    # Bubble topology on the BFS neighborhood (uses labels + var_per from the final loop iteration;
    # falls back to "no_main" if BFS never produced a usable neighborhood).
    topology = (classify_neighborhood_topology(
                  set(nhood_segs), adj_und, segs, labels, var_per,
                  min_core_len=min_core_len)
                if labels else
                dict(type="no_main", n_main_seg=0, n_path=0,
                     n_shared_anchor=0, n_shared_flank_anchor=0,
                     main_lens=[], main_ids=[]))
    print(f"  [{k}] topology: {topology['type']} "
          f"(main_segs={topology['n_main_seg']}, paths={topology['n_path']}, "
          f"shared={topology['n_shared_anchor']}, shared_flank={topology['n_shared_flank_anchor']})")
    print(f"  [{k}] -> {n_written} kept candidates (>= {min_allele_len} bp); complete-path coverage signal = {n_complete_written}")
    return fa_path, tsv_path, n_complete_written, topology


# ---------- entry point ----------

def _read_fasta_records(fa: str) -> list[tuple[str, str]]:
    """Return [(name, seq), ...]."""
    out, name, body = [], None, []
    for ln in open(fa):
        ln = ln.rstrip()
        if not ln: continue
        if ln.startswith(">"):
            if name is not None: out.append((name, "".join(body)))
            name = ln[1:].split()[0]; body = []
        else:
            body.append(ln)
    if name is not None: out.append((name, "".join(body)))
    return out


def _trim_to_locus_envelope(seq: str, locus_ref_fa: str, pad: int = 100) -> str:
    """Blastn locus_ref → seq; trim seq to the candidate's [min(sstart)-pad,
    max(send)+pad] envelope. No qualifying hits → returns seq unchanged.

    Uses subprocess directly (no blast_cache) because the subject DB content
    is per-sequence but its basename is fixed — the cache would falsely hit.
    """
    if not seq or not os.path.exists(locus_ref_fa):
        return seq
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "s.fa")
        with open(sf, "w") as o: o.write(f">x\n{seq}\n")
        db = os.path.join(t, "s_db")
        subprocess.run(["makeblastdb", "-in", sf, "-dbtype", "nucl", "-out", db],
                        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out = subprocess.run(
            ["blastn", "-query", locus_ref_fa, "-db", db, "-dust", "no",
             "-outfmt", "6 sseqid pident length sstart send"],
            check=True, capture_output=True, text=True).stdout
    pts: list[int] = []
    for ln in out.splitlines():
        f = ln.rstrip("\n").split("\t")
        if len(f) < 5: continue
        try:
            if float(f[1]) < 80.0 or int(f[2]) < 100: continue   # min_pid=80, min_len=100
            pts += [int(f[3]), int(f[4])]
        except ValueError: continue
    if not pts: return seq
    lo = max(0, min(pts) - pad)
    hi = min(len(seq), max(pts) + pad)
    return seq[lo:hi]


def _read_first_fasta(fa: str) -> tuple[str, str]:
    name = ""; buf: list[str] = []
    for ln in open(fa):
        if ln.startswith(">"):
            if name: break
            name = ln[1:].strip().split()[0]
        else:
            buf.append(ln.strip())
    return name, "".join(buf).upper()


_REF_KEY = "__REF__"   # reserved name inside the MSA; must not collide with candidates


def _hd_cols_from_aligned_ref(aligned_ref: str,
                               ref_hd_lo: int, ref_hd_hi: int,
                               flipped: bool) -> set[int]:
    """Project an HD-core span (1-based inclusive, in ORIGINAL forward-ref
    ungapped coords) onto column indices of a gapped, MAFFT-aligned reference.

    When MAFFT --adjustdirection reverse-complements the reference (signalled by
    a "_R_" prefix on its output record), position 1 of the aligned ref maps to
    the LAST nucleotide of the original ungapped ref, so the span has to be
    remapped to the RC frame before column-walking:

        scan_lo = ref_ungap_len - ref_hd_hi + 1
        scan_hi = ref_ungap_len - ref_hd_lo + 1

    Returns the set of 0-based column indices in `aligned_ref` covered by the
    span (after any RC remap).
    """
    ref_ungap_len = sum(1 for c in aligned_ref if c != '-')
    if flipped:
        scan_lo = ref_ungap_len - ref_hd_hi + 1
        scan_hi = ref_ungap_len - ref_hd_lo + 1
    else:
        scan_lo, scan_hi = ref_hd_lo, ref_hd_hi
    hd_cols: set[int] = set()
    ungap = 0
    for ci, ch in enumerate(aligned_ref):
        if ch != '-':
            ungap += 1
            if scan_lo <= ungap <= scan_hi:
                hd_cols.add(ci)
    return hd_cols


def _trim_by_hd_core(seq: str, hd_proteins_fa: str, pad: int = 1000) -> tuple[str, int, int]:
    """tblastn `hd_proteins_fa` → seq. Trim seq to [min(sstart) - pad,
    max(send) + pad] (clamped to seq bounds). Returns (trimmed_seq, hd_lo, hd_hi)
    where (hd_lo, hd_hi) are the HD-core span in the UN-GAPPED TRIMMED coords
    (1-based inclusive). If no qualifying tblastn hits, returns (seq, 1, len(seq))
    — i.e. no trim, treat whole sequence as core.
    """
    if not seq or not (hd_proteins_fa and os.path.exists(hd_proteins_fa)):
        return seq, 1, len(seq) if seq else 0
    lo_orig, hi_orig = detect_core_span(seq, hd_proteins_fa)
    if lo_orig == 1 and hi_orig == len(seq):
        return seq, 1, len(seq)
    L = len(seq)
    cut_lo = max(0, lo_orig - 1 - pad)            # 0-based slice start
    cut_hi = min(L, hi_orig + pad)                # 0-based slice end (exclusive)
    trimmed = seq[cut_lo:cut_hi]
    new_lo = lo_orig - cut_lo                      # HD-core start in trimmed (1-based)
    new_hi = hi_orig - cut_lo                      # HD-core end   in trimmed (1-based)
    return trimmed, new_lo, new_hi


def _align_cores(seqs: list[tuple[str, str]],
                 locus_ref_fa: str | None,
                 hd_proteins_fa: str | None,
                 pad: int = 1000,
                 fast: bool = False,
                 ) -> tuple[dict[str, str], set[int]]:
    """1. TRIM each candidate to HD-core span ± `pad` bp (tblastn HD proteins
          on the candidate). Generous pad keeps enough flank for MSA anchoring.
       2. MSA — run ONE MAFFT with the REFERENCE LOCUS (similarly trimmed) in
          the input alongside the trimmed candidates. The reference is the
          coordinate ruler: everything is aligned to the same frame.
       3. HD-CORE COLUMNS — the trimmed reference has a known HD-core span
          (ungapped). Walk the reference's aligned string; the MSA columns
          where the reference's un-gapped position lies in that span are the
          HD-core columns. ONE FIXED SET, shared across all candidate pairs.

    Returns (aligned, hd_cols):
      aligned[name] = MSA-aligned uppercase string for each CANDIDATE
                      (reference dropped before returning)
      hd_cols       = set of MSA column indices for HD-core (reference frame)

    Identity is computed by the caller over hd_cols, counting only positions
    where BOTH candidates are non-gap. Flanks anchor the MSA but don't
    contribute to identity.

    Tests mock this function to skip mafft+blast binary dependencies.
    """
    # 1. tblastn-trim each candidate to HD-core ± pad bp.
    trimmed: list[tuple[str, str]] = []
    for name, seq in seqs:
        if hd_proteins_fa and os.path.exists(hd_proteins_fa):
            t, _lo, _hi = _trim_by_hd_core(seq, hd_proteins_fa, pad=pad)
            trimmed.append((name, t or seq))
        else:
            trimmed.append((name, seq))
    # 2. Trim the reference the same way, so it sits in the same coord frame.
    ref_seq = ""
    ref_hd_lo = ref_hd_hi = 0
    if locus_ref_fa and os.path.exists(locus_ref_fa):
        _ref_name, full_ref = _read_first_fasta(locus_ref_fa)
        if full_ref and hd_proteins_fa and os.path.exists(hd_proteins_fa):
            ref_seq, ref_hd_lo, ref_hd_hi = _trim_by_hd_core(full_ref, hd_proteins_fa, pad=pad)
        else:
            ref_seq, ref_hd_lo, ref_hd_hi = full_ref, 1, len(full_ref)
    # 3. MSA: trimmed reference + trimmed candidates together.
    msa_input = ([(_REF_KEY, ref_seq)] if ref_seq else []) + trimmed
    aligned: dict[str, str] = {}
    with tempfile.TemporaryDirectory() as t:
        in_fa  = os.path.join(t, "in.fa")
        out_fa = os.path.join(t, "msa.fa")
        with open(in_fa, "w") as o:
            for n, s in msa_input: o.write(f">{n}\n{s}\n")
        with open(out_fa, "w") as o:
            # --adjustdirection: MAFFT tries both strands per record and picks
            # the one that aligns best. Without it, candidates on opposite
            # strands of the same allele (a routine outcome of graph traversal)
            # align poorly and get scored as falsely distinct.
            # fast=True selects FFT-NS-1 (single progressive alignment, no
            # iterative refinement) — ~10× faster than --auto on N>50 inputs.
            # Good enough for dedup at id≥0.95; not used for final identity.
            mafft_args = (["mafft", "--adjustdirection", "--retree", "1",
                           "--maxiterate", "0", in_fa]
                           if fast else
                           ["mafft", "--adjustdirection", "--auto", in_fa])
            subprocess.run(mafft_args, check=True, stdout=o, stderr=subprocess.DEVNULL)
        cur_n = None; buf = []
        flipped: set[str] = set()    # records MAFFT reverse-complemented
        for ln in open(out_fa):
            if ln.startswith(">"):
                if cur_n is not None: aligned[cur_n] = "".join(buf).upper()
                nm = ln[1:].strip().split()[0]
                # --adjustdirection prefixes reverse-complemented records with
                # "_R_". Strip the prefix but REMEMBER which records got flipped.
                if nm.startswith("_R_"):
                    nm = nm[3:]; flipped.add(nm)
                cur_n = nm; buf = []
            else:
                buf.append(ln.strip())
        if cur_n is not None: aligned[cur_n] = "".join(buf).upper()
    # 4. Find HD-core MSA columns via the (trimmed) reference's aligned positions.
    hd_cols: set[int] = set()
    if _REF_KEY in aligned and ref_hd_lo > 0:
        hd_cols = _hd_cols_from_aligned_ref(
            aligned[_REF_KEY], ref_hd_lo, ref_hd_hi,
            flipped=(_REF_KEY in flipped))
        del aligned[_REF_KEY]
    else:
        # No reference available — degraded mode: union of all non-gap cols.
        for s in aligned.values():
            for ci, ch in enumerate(s):
                if ch != '-': hd_cols.add(ci)
    return aligned, hd_cols


def _cluster_alleles_one_msa(seqs: list[tuple[str, str]],
                              queries_dir: str | None = None,
                              locus_ref_fa: str | None = None,
                              id_thresh: float = 0.95,
                              frac_thresh: float = 0.80,
                              fast: bool = False) -> list[list[str]]:
    """Single-shot cluster_alleles. Runs ONE MAFFT MSA on all inputs and
    greedy-clusters by HD-core identity. See the wrapper `cluster_alleles`
    (below) for the production entry point — it batches large input pools
    through this function to keep individual MSAs manageable.

    fast=True uses lossy MAFFT (--retree 1 --maxiterate 0); enough for dedup
    at id≥0.95 over HD-core columns, ~10× faster than --auto on N>50 inputs.
    """
    if not seqs: return []
    if len(seqs) == 1: return [[seqs[0][0]]]
    hd_proteins_fa = (os.path.join(queries_dir, "variable_proteins.fasta")
                       if queries_dir else None)
    aligned, hd_cols = _align_cores(seqs, locus_ref_fa, hd_proteins_fa, fast=fast)
    n_hd = len(hd_cols) or 1
    def pair(name_a: str, name_b: str) -> tuple[float, float]:
        sa, sb = aligned[name_a], aligned[name_b]
        match = aln = 0
        for ci in hd_cols:
            ca, cb = sa[ci], sb[ci]
            if ca == '-' or cb == '-': continue
            aln += 1
            if ca == cb: match += 1
        pid = (match / aln) if aln else 0.0
        frac = aln / n_hd   # fraction of REFERENCE HD-core covered by both
        return pid, frac
    clusters: list[list[str]] = []
    reps: list[str] = []
    for name, _seq in seqs:
        if name not in aligned:
            clusters.append([name]); reps.append(name); continue
        matched = False
        for i, rep in enumerate(reps):
            pid, frac = pair(name, rep)
            if pid >= id_thresh and frac >= frac_thresh:
                clusters[i].append(name); matched = True; break
        if not matched:
            clusters.append([name]); reps.append(name)
    return clusters


def cluster_alleles(seqs: list[tuple[str, str]],
                    queries_dir: str | None = None,
                    locus_ref_fa: str | None = None,
                    id_thresh: float = 0.95,
                    frac_thresh: float = 0.80,
                    batch_size: int = 20) -> list[list[str]]:
    """Cluster candidates by HD-core MAFFT identity.

    Pipeline:
      0. CANONICAL DEDUP — collapse byte-equal fwd/rc sequences. Often the
         widen loop emits many walks that share a sequence (different orders
         of identical segments); these are merged before any MAFFT runs.
      1. If canonical-distinct count ≤ batch_size: ONE accurate MAFFT MSA
         (--auto) on all canonical reps. Result clusters expand back to
         original names.
      2. Otherwise BATCHED: sort canonical reps by length desc (longer
         anchors more flank context for MAFFT). Round 0 runs lossy MAFFT
         (--retree 1 --maxiterate 0) on the first `batch_size` reps and
         keeps one rep per cluster (the longest). Round i runs MAFFT on
         (carryover rep survivors + next `batch_size` new reps), recomputes
         clusters, keeps one rep per cluster. Repeats until input is
         exhausted. Lossy MAFFT is fine for dedup at id≥0.95; the final
         pairwise_identity step (8a/8b) uses --auto on the picked pair.

    Returns a list of clusters, each = list of ORIGINAL record names.
    """
    if not seqs: return []
    if len(seqs) == 1: return [[seqs[0][0]]]

    # 0. Canonical dedup
    seq_of = dict(seqs)
    canon: dict[str, list[str]] = {}      # canonical_key -> all original names
    rep_of: dict[str, str] = {}           # canonical_key -> longest name in group
    for n, s in seqs:
        key = min(s, _rc(s))
        canon.setdefault(key, []).append(n)
        if key not in rep_of or len(s) > len(seq_of[rep_of[key]]):
            rep_of[key] = n
    canon_reps = [(rep_of[k], seq_of[rep_of[k]]) for k in canon]
    rep_to_key = {rep_of[k]: k for k in canon}

    if len(canon_reps) == 1:
        return [list(canon.values())[0]]

    def _expand(rep_clusters: list[list[str]]) -> list[list[str]]:
        """Expand a clustering over canonical reps back to original names."""
        out = []
        for rc in rep_clusters:
            members: list[str] = []
            for rep_name in rc:
                members.extend(canon[rep_to_key[rep_name]])
            out.append(members)
        return out

    # 1. Small pool: single-shot accurate MAFFT
    if len(canon_reps) <= batch_size:
        rep_clusters = _cluster_alleles_one_msa(canon_reps, queries_dir, locus_ref_fa,
                                                  id_thresh, frac_thresh, fast=False)
        return _expand(rep_clusters)

    # 2. Large pool: batched lossy MAFFT, longest-first
    canon_reps.sort(key=lambda x: -len(x[1]))
    survivors: list[str] = []                      # current cluster reps (canonical-rep names)
    survivor_to_members: dict[str, list[str]] = {} # rep -> all canonical-rep names in its cluster
    n_rounds = 0
    for i in range(0, len(canon_reps), batch_size):
        batch = canon_reps[i:i + batch_size]
        round_input = [(n, seq_of[n]) for n in survivors] + batch
        carryover = len(survivors)
        rc_clusters = _cluster_alleles_one_msa(round_input, queries_dir, locus_ref_fa,
                                                 id_thresh, frac_thresh, fast=True)
        new_survivors: list[str] = []
        new_to_members: dict[str, list[str]] = {}
        for rc in rc_clusters:
            rep = max(rc, key=lambda n: len(seq_of[n]))
            new_survivors.append(rep)
            members: list[str] = []
            for n in rc:
                if n in survivor_to_members:
                    members.extend(survivor_to_members[n])
                else:
                    members.append(n)
            new_to_members[rep] = members
        survivors, survivor_to_members = new_survivors, new_to_members
        n_rounds += 1
        print(f"  [cluster_alleles] batched dedup round {n_rounds}: input "
              f"{len(round_input)} ({carryover} carryover + {len(batch)} new) -> "
              f"{len(survivors)} survivors")
    return _expand(list(survivor_to_members.values()))


def run(sample: str, spades_dir: str, queries_dir: str, ks: list[str], outdir: str,
        anchors_contigs_fa: str | None = None,
        anchors_segments_fa: str | None = None,
        repeats: str | None = None, known_degHD: str | None = None,
        max_locus_len: int = 0, min_allele_len: int = 2000,
        max_nodes: int = 15, max_paths: int = 50000,
        max_walk_bp: int = 0,
        init_hops: int = 5, max_hops: int = 10, threads: int = 4,
        genome_cov: float = 0.0, cov_repeat_factor: float = 2.0,
        expected_count: int = 2,
        dup_id: float = 0.95, dup_frac: float = 0.80,
        asymmetric_bfs: bool = False,
        seeds_from: str = "hd",
        locus_ref_fa: str | None = None,
        blast_out_dir: str | None = None) -> int:
    """Return TOTAL number of DISTINCT complete alleles found across all k's.

    "Distinct" = greedy single-link MAFFT clustering at id >= `dup_id` AND aligned
    fraction (over the shorter) >= `dup_frac`. Two complete paths at id >= 0.95
    over >= 80% length are the same allele; they count as ONE.

    Early termination: as soon as the number of distinct allele clusters reaches
    `expected_count` (default 2 = dikaryon), the remaining k's are skipped — the
    first k that delivers two divergent complete alleles is enough.
    """
    os.makedirs(outdir, exist_ok=True)
    cand_fa  = os.path.join(outdir, "bubble_alleles.fasta")
    cand_tsv = os.path.join(outdir, "bubble_alleles.ann.tsv")
    open(cand_fa, "w").close(); open(cand_tsv, "w").close()
    contig_per_k, seg_per_k = collect_anchors_per_k(anchors_contigs_fa, anchors_segments_fa,
                                                      seeds_from=seeds_from)
    complete_pool: list[tuple[str, str]] = []   # accumulated (name, seq) of complete paths
    n_complete_total = 0
    n_distinct = 0
    topology_per_k: dict[str, dict] = {}
    for k in ks:
        ctg, gfa = spades_k_paths(spades_dir, k)
        if not (ctg and gfa):
            print(f"  [{k}] no contigs+gfa under {spades_dir}; skipping")
            topology_per_k[k] = dict(type="NA_no_gfa", n_main_seg=0, n_path=0,
                                       n_shared_anchor=0, n_shared_flank_anchor=0,
                                       main_lens=[], main_ids=[])
            continue
        fa, tsv, n_complete, topology = run_one_k(
            sample, k, gfa, ctg, queries_dir, outdir,
            anchor_contigs=contig_per_k.get(k, set()),
            anchor_segments_direct=seg_per_k.get(k, set()),
            repeats=repeats, known_degHD=known_degHD,
            max_locus_len=max_locus_len, min_allele_len=min_allele_len,
            max_walk_bp=max_walk_bp,
            max_nodes=max_nodes, max_paths=max_paths,
            init_hops=init_hops, max_hops=max_hops, threads=threads,
            genome_cov=genome_cov, cov_repeat_factor=cov_repeat_factor,
            asymmetric_bfs=asymmetric_bfs,
            expected_count=expected_count,
            dup_id=dup_id, dup_frac=dup_frac,
            locus_ref_fa=locus_ref_fa,
            blast_out_dir=blast_out_dir)
        topology_per_k[k] = topology
        n_complete_total += n_complete
        with open(cand_fa, "a") as o:
            for ln in open(fa): o.write(ln)
        with open(cand_tsv, "a") as o:
            for ln in open(tsv): o.write(ln)
        # Pull this k's complete paths (the FIRST n_complete records in fa) into the
        # accumulated pool and re-cluster by MAFFT divergence.
        recs_this_k = _read_fasta_records(fa)[:n_complete]
        complete_pool.extend(recs_this_k)
        clusters = cluster_alleles(complete_pool, queries_dir=queries_dir,
                                     id_thresh=dup_id, frac_thresh=dup_frac)
        n_distinct = len(clusters)
        print(f"  [{k}] running totals: {len(complete_pool)} complete path(s), "
              f"{n_distinct} distinct allele(s) at id>={dup_id:.2f}, aln>={dup_frac:.2f} "
              f"(MAFFT aligns full allele but identity is computed on HD-core columns only)")
        # Early termination: enough divergent complete alleles -> skip remaining k's.
        if n_distinct >= expected_count:
            i = ks.index(k)
            skipped = ks[i + 1:]
            if skipped:
                print(f"[graph_path_search] {n_distinct} distinct >= expected {expected_count} after k={k} "
                      f"-> skipping remaining k's: {','.join(skipped)}")
            break
    n_records = sum(1 for ln in open(cand_fa) if ln.startswith(">"))
    # Aggregate per-k bubble topology into one TSV (replaces bubble_topo.py output).
    topology_tsv = os.path.join(outdir, "bubble_topology.tsv")
    with open(topology_tsv, "w") as o:
        o.write("sample\tk\ttype\tn_main_seg\tn_path\tn_shared_anchor"
                "\tn_shared_flank_anchor\thd_seg_lens\n")
        for k, t in topology_per_k.items():
            lens = ",".join(str(x) for x in t["main_lens"]) or "-"
            o.write(f"{sample}\t{k}\t{t['type']}\t{t['n_main_seg']}\t{t['n_path']}"
                    f"\t{t['n_shared_anchor']}\t{t['n_shared_flank_anchor']}\t{lens}\n")
    print(f"[graph_path_search] {sample}: {n_records} candidate path(s), "
          f"{n_complete_total} complete, {n_distinct} distinct allele(s) -> {cand_fa}")
    print(f"[graph_path_search] bubble topology per k -> {topology_tsv}")
    # Emit the picker-input pool: filtered anchor contigs ∪ all bubble paths.
    _emit_picker_candidates(sample, outdir, ks)
    return n_distinct


def _emit_picker_candidates(sample: str, outdir: str, ks: list[str]) -> None:
    """Build picker_candidates.fasta + picker_candidates.ann.tsv = filtered
    anchor contigs ∪ all bubble paths. The filter drops anchor contigs whose
    `anchor_contig_<k>.ann.tsv` row shows 0 variable genes AND lacks at least
    one of (flankL, flankR) — flankL-only / flankR-only / pure-junk fragments
    that step 3 needs for BFS seeding but that the picker can never use as a
    real allele. Bubble paths (end-to-end walks by construction) pass through
    unconditionally. The picker consumes this pool with no further filter
    logic — see METHODS §4 / "Pool assembly for the picker".
    """
    pick_fa  = os.path.join(outdir, "picker_candidates.fasta")
    pick_tsv = os.path.join(outdir, "picker_candidates.ann.tsv")
    keep_names: set[str] = set()
    n_total = n_kept = n_dropped = 0
    open(pick_fa, "w").close(); open(pick_tsv, "w").close()
    for k in ks:
        ann_p = os.path.join(outdir, f"anchor_contig_{k}.ann.tsv")
        fa_p  = os.path.join(outdir, f"anchor_contig_{k}.fasta")
        if not (os.path.exists(ann_p) and os.path.exists(fa_p)): continue
        for ln in open(ann_p):
            fields = ln.rstrip("\n").split("\t")
            if len(fields) < 6: continue
            name, vars_col, flanks_col = fields[0], fields[4], fields[5]
            n_total += 1
            has_var = vars_col not in ("", "-")
            flank_set = set() if flanks_col in ("", "-") else set(flanks_col.split(","))
            has_both_flanks = ("flankL" in flank_set) and ("flankR" in flank_set)
            if has_var or has_both_flanks:
                keep_names.add(name); n_kept += 1
                with open(pick_tsv, "a") as o: o.write(ln if ln.endswith("\n") else ln + "\n")
            else:
                n_dropped += 1
        # stream sequences from fa_p, keeping only names in keep_names
        emit = False
        with open(fa_p) as fh, open(pick_fa, "a") as o:
            for line in fh:
                if line.startswith(">"):
                    nm = line[1:].rstrip().split()[0]
                    emit = nm in keep_names
                if emit: o.write(line)
    # append all bubble path records unconditionally
    bb_fa  = os.path.join(outdir, "bubble_alleles.fasta")
    bb_tsv = os.path.join(outdir, "bubble_alleles.ann.tsv")
    n_paths = 0
    if os.path.exists(bb_fa):
        with open(bb_fa) as fh, open(pick_fa, "a") as o:
            for line in fh:
                if line.startswith(">"): n_paths += 1
                o.write(line)
    if os.path.exists(bb_tsv):
        with open(bb_tsv) as fh, open(pick_tsv, "a") as o:
            for line in fh: o.write(line)
    print(f"[graph_path_search] {sample}: picker_candidates -> "
          f"{n_kept}/{n_total} anchor contig(s) kept ({n_dropped} fragments dropped) "
          f"+ {n_paths} bubble path(s)  -> {pick_fa}")


def _cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--sample", required=True)
    p.add_argument("--spades-dir", required=True)
    p.add_argument("--queries-dir", required=True)
    p.add_argument("--ks", default="k21,k33,k55")
    p.add_argument("--outdir", required=True)
    p.add_argument("--anchors-contigs",  default=None,
                   help="anchor_contig.fasta produced by step 2.1 (anchor_search --source contigs)")
    p.add_argument("--anchors-segments", default=None,
                   help="anchor_segments.fasta produced by step 2.2 (anchor_search --source segments)")
    p.add_argument("--repeats", default=None)
    p.add_argument("--known-degHD", default=None)
    p.add_argument("--locus-ref", default=None,
                   help="locus reference fasta (single record). Used by cluster_alleles "
                        "to anchor the MSA and project HD-core columns from the reference's "
                        "ungapped span. Without it, clustering falls back to whole-allele "
                        "identity, which misses HD divergence in samples with long flanks "
                        "or wandering walks through repeats.")
    p.add_argument("--max-locus-len", type=int, default=0,
                   help="bp budget for path enumeration. Default 0 = NO LIMIT — let "
                        "downstream MAFFT clustering dedup the long paths instead of "
                        "pre-pruning them here. (Long flank-bearing GFA segments are "
                        "common, and a tight budget cuts off real allele walks.)")
    p.add_argument("--min-allele-len", type=int, default=2000)
    p.add_argument("--max-nodes", type=int, default=15)
    p.add_argument("--max-paths", type=int, default=50000)
    p.add_argument("--max-walk-bp", type=int, default=0,
                   help="bp budget for DFS path enumeration (independent of "
                        "--max-locus-len, which only triggers pick_alleles' "
                        "pre-MAFFT trim). 0 = no cap. The wrapper sets a "
                        "sensible default = 10x manifest's derived_max_locus_len.")
    p.add_argument("--init-hops", type=int, default=5)
    p.add_argument("--max-hops", type=int, default=10,
                   help="Widen BFS up to this many hops. Loop also breaks early "
                        "on closed_bubble topology OR n_distinct >= --expected-count.")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--blast-out-dir", default=None,
                   help="If set, the heavy GFA-search blasts inside label_segments (per-hop "
                        "neighborhood flankL/R, degHD, repeat) write their hit tables to this "
                        "directory as named TSVs (e.g. gfa_flankL_blastn_nhood_k45_hops5.tsv). "
                        "On re-run, an existing file is reused as-is (file-existence-only "
                        "check). The wrapper exposes --re-blast to wipe stale files when "
                        "blast parameters change.")
    p.add_argument("--genome-coverage", type=float, default=0.0,
                   help="mean genome read depth (estimated upstream from flank-region depth). "
                        "Drives the two-pass DFS: segments with depth >= cov-repeat-factor x "
                        "this value are deprioritized in pass 1. <=0 disables.")
    p.add_argument("--cov-repeat-factor", type=float, default=2.0,
                   help="multiplier on --genome-coverage above which a segment is treated as a "
                        "coverage-by-repeat (MITE/degHD/SD-collapsed). Default 2.0 — "
                        "segments with depth > 2x single-copy genome cov are treated as "
                        "collapsed two-copy or higher.")
    p.add_argument("--expected-count", type=int, default=2,
                   help="exit 0 only if AT LEAST this many DISTINCT complete alleles "
                        "(complete = flank + HD genes + flank; distinct = below the dup-id / "
                        "dup-frac MAFFT threshold from any other) were produced across all k's. "
                        "Default 2 (dikaryon). Drives BOTH the exit code AND the early-termination "
                        "loop inside run() — the loop stops at the first k that delivers two "
                        "divergent complete alleles.")
    p.add_argument("--dup-id",   type=float, default=0.95,
                   help="MAFFT identity threshold above which two complete paths are considered "
                        "the SAME allele (default 0.95).")
    p.add_argument("--dup-frac", type=float, default=0.80,
                   help="MAFFT aligned-fraction threshold above which two complete paths are "
                        "considered the SAME allele (default 0.80, over the shorter sequence).")
    p.add_argument("--seeds-from", choices=("hd", "all"), default="hd",
                   help="Which anchor records seed the BFS. Default 'hd' = only records "
                        "whose ann.tsv hd_genes column is non-'-' (i.e. carry a variable "
                        "gene). 'all' = use every anchor including flank-only. The BFS "
                        "reaches flank-bearing segs from HD seeds via L-line adjacency "
                        "anyway, so flank-only seeds typically just pull in noise.")
    p.add_argument("--asymmetric-bfs", action="store_true",
                   help="OPTIONAL: BFS may ABSORB repeat segments into the neighborhood but "
                        "does NOT expand neighbors FROM them. Repeats are detected via "
                        "depth > --cov-repeat-factor x --genome-coverage. Default OFF: "
                        "symmetric BFS (repeats expand normally). Turn on to dampen "
                        "combinatorial blow-up at high-copy regions if you see the BFS "
                        "neighborhood explode.")
    a = p.parse_args(argv)
    if not a.anchors_contigs and not a.anchors_segments:
        p.error("at least one of --anchors-contigs / --anchors-segments must be supplied")
    n_distinct = run(a.sample, a.spades_dir, a.queries_dir, a.ks.split(","), a.outdir,
                     anchors_contigs_fa=a.anchors_contigs,
                     anchors_segments_fa=a.anchors_segments,
                     repeats=a.repeats, known_degHD=a.known_degHD,
                     max_locus_len=a.max_locus_len, min_allele_len=a.min_allele_len,
                     max_nodes=a.max_nodes, max_paths=a.max_paths,
                     max_walk_bp=a.max_walk_bp,
                     init_hops=a.init_hops, max_hops=a.max_hops, threads=a.threads,
                     genome_cov=a.genome_coverage, cov_repeat_factor=a.cov_repeat_factor,
                     expected_count=a.expected_count,
                     dup_id=a.dup_id, dup_frac=a.dup_frac,
                     asymmetric_bfs=a.asymmetric_bfs,
                     seeds_from=a.seeds_from,
                     blast_out_dir=a.blast_out_dir,
                     locus_ref_fa=a.locus_ref)
    if n_distinct >= a.expected_count:
        print(f"[graph_path_search] distinct alleles {n_distinct} >= expected {a.expected_count} -> OK")
        sys.exit(0)
    print(f"[graph_path_search] distinct alleles {n_distinct} < expected {a.expected_count} "
          f"-> exit 2 (orchestrator should fall back to step 2.2 if not already done)",
          file=sys.stderr)
    sys.exit(2)

if __name__ == "__main__":
    _cli()
