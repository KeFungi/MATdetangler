"""Per-k allele caller — one GFA + one label-hits TSV in, one allele set out.

Contract: a single k (Suilu k33 / k45 / k53; Pcub k55 / k77). Cross-k
integration is a downstream selection step that picks the best per-k
verdict + allele set.

Algorithm
=========

Phase 1 (main BFS — var-seeded):

    while nhop <= max_nhop:
        1. neighborhood = N-hop BFS from var segments in the GFA.
        2. Apply coverage filter [lo*genome_cov, hi*genome_cov] on segments.
        3. Classify the resulting bubble.
           If closed_bubble AND arms are divergent: break.
        4. Fallback: remove the coverage filter, classify again.
           If closed_bubble AND arms are divergent: break.
        5. Else: nhop += 1, loop.

Phase 2 (flank fallback — var + flank seeded):

    Same loop, but seed BFS from var nodes AND flank nodes.

Emission rule (uniform dedup across all topologies)
===================================================

Every topology produces a list of candidate sequences (one per arm /
component / simple-path). The candidate list is then collapsed by
pairwise divergence — sequences that are NOT divergent from an
already-kept representative drop out, and the longest representative
wins. The verdict **always preserves the classifier's origin**; only
the emitted record naming reflects the post-dedup count `n`:

    n == 1  → allele1                                    (verdict = origin)
    n == 2  → allele1 + allele2                          (verdict = origin)
    n >= 3  → chimera1, chimera2, … (one per dedup seq)  (verdict = origin)

Examples:
    single        (1 candidate)                          → allele1            verdict=single
    closed_bubble (2 arms, dedup→1)                      → allele1            verdict=closed_bubble
    closed_bubble (2 arms, dedup→2)                      → allele1+allele2    verdict=closed_bubble
    open_bubble   (2 arms, dedup→2)                      → allele1+allele2    verdict=open_bubble
    separate      (3 components, dedup→2)                → allele1+allele2    verdict=separate
    separate      (3 components, dedup→3)                → chimera1..3        verdict=separate
    complexed     (5 paths, dedup→1)                     → allele1            verdict=complexed
    complexed     (5 paths, dedup→2)                     → allele1+allele2    verdict=complexed
    complexed     (5 paths, dedup→4)                     → chimera1..4        verdict=complexed

Divergence
==========

Two reconstructed arm sequences are "divergent" if pairwise identity is
below `1 - divergence_threshold`. The default threshold is 0.05 (i.e.
arms must differ by ≥5% to count as distinct alleles).
"""
from __future__ import annotations
from collections import deque
from .seg_processor import directional_split
from .bubble_classifier import classify
from .bubble_bfs import build_adj
from .arm_sequences import parse_gfa_sequences, reconstruct_arms, reverse_complement
from .labeler import read_seg_label_hits, read_segment_lengths_from_hits


# -----------------------------------------------------------------------------
# GFA-level helpers
# -----------------------------------------------------------------------------

def parse_gfa(gfa_path: str
              ) -> tuple[dict[str, str], dict[str, float],
                          set[frozenset], dict[tuple, tuple[str, str]]]:
    """Parse a GFA. Returns (sequences, depths, edges, edge_endpoints)."""
    sequences: dict[str, str] = {}
    depths: dict[str, float] = {}
    edges: set[frozenset] = set()
    endpoints: dict[tuple, tuple[str, str]] = {}
    with open(gfa_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f: continue
            if f[0] == "S":
                if len(f) >= 3:
                    sequences[f[1]] = f[2]
                    for tag in f[3:]:
                        if tag.startswith("DP:f:"):
                            try: depths[f[1]] = float(tag[5:])
                            except ValueError: pass
                            break
            elif f[0] == "L" and len(f) >= 5:
                a, oa, b, ob = f[1], f[2], f[3], f[4]
                end_a = "R" if oa == "+" else "L"
                end_b = "L" if ob == "+" else "R"
                edges.add(frozenset((a, b)))
                key = tuple(sorted([a, b]))
                if a == key[0]:
                    endpoints[key] = (end_a, end_b)
                else:
                    endpoints[key] = (end_b, end_a)
    return sequences, depths, edges, endpoints


def parse_contigs_paths(path: str) -> dict[str, list[str]]:
    """Parse SPAdes contigs.paths. Returns {contig_name: [seg_id, ...]}.

    Format (modern SPAdes):
        NODE_1_length_5000_cov_45.6
        12+,3-,45+;
        78+,9-;
        NODE_1_length_5000_cov_45.6'                   ← RC variant
        45-,3+,12-
        ...

    Each contig has 1+ subpaths (semicolon-separated) of `seg_id±` tokens.
    We strip orientation marks and return the union of all visited seg IDs
    (orientation-insensitive). RC entries (trailing `'`) are merged into
    their parent contig.
    """
    walks: dict[str, list[str]] = {}
    cur: str | None = None
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line: continue
            if cur is None:
                cur = line.rstrip("'")                          # strip RC marker
                walks.setdefault(cur, [])
            else:
                for chunk in line.replace(";", ",").split(","):
                    seg = chunk.rstrip("+-").strip()
                    if seg:
                        walks[cur].append(seg)
                cur = None
    return walks


def contig_seeds_for_locus(
        contig_blast_tsvs: list[str],
        contigs_paths_file: str,
        min_pid: float = 80.0,
        min_alnlen: int = 100,
) -> set[str]:
    """Identify GFA segments anchored by any contig that BLAST-hits the locus.

    Reads outfmt-6 contig-BLAST TSVs (qseqid sseqid pident length ...),
    keeps subjects (sseqid = contig names) with pident ≥ min_pid AND
    length ≥ min_alnlen. Maps those contigs → segments via contigs.paths.
    Returns the union of segments across all locus-bearing contigs.
    """
    import os
    if not os.path.exists(contigs_paths_file):
        return set()
    walks = parse_contigs_paths(contigs_paths_file)
    locus_contigs: set[str] = set()
    for path in contig_blast_tsvs:
        if not os.path.exists(path): continue
        with open(path) as fh:
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) < 4: continue
                try:
                    pid = float(f[2]); ln_ = int(f[3])
                except ValueError: continue
                if pid >= min_pid and ln_ >= min_alnlen:
                    locus_contigs.add(f[1])
    segs: set[str] = set()
    for c in locus_contigs:
        # SPAdes contigs.paths may name contigs without prefix; also
        # strip the trailing _cov_X.X if header has it. Try both forms.
        if c in walks:
            segs.update(walks[c])
            continue
        # Heuristic: tolerate suffix differences
        cn = c.rstrip("'")
        for k, v in walks.items():
            if k == cn or k.startswith(cn + "_") or cn.startswith(k + "_"):
                segs.update(v)
                break
    return segs


def bfs_expand_segments(seeds: set[str], adj_und: dict[str, set[str]],
                         n_hops: int) -> set[str]:
    """N-hop BFS in the GFA undirected adjacency, returns the reachable set."""
    if n_hops <= 0: return set(seeds)
    seen = set(seeds)
    frontier = set(seeds)
    for _ in range(n_hops):
        nxt = set()
        for n in frontier:
            for m in adj_und.get(n, ()):
                if m not in seen:
                    nxt.add(m)
        if not nxt: break
        seen |= nxt
        frontier = nxt
    return seen


def filter_segs_by_depth(segs: set[str], depths: dict[str, float],
                          genome_cov: float,
                          lo_mult: float = 0.25,
                          hi_mult: float = 2.0) -> set[str]:
    """Keep only segments whose depth is within [lo_mult * genome_cov,
    hi_mult * genome_cov]. Segments without a depth value are kept (we don't
    want to throw them away just because we couldn't look up depth)."""
    lo, hi = lo_mult * genome_cov, hi_mult * genome_cov
    return {s for s in segs if s not in depths or lo <= depths[s] <= hi}


def restrict_subgraph(nodes: set[str], edges: set[frozenset],
                       endpoints: dict[tuple, tuple[str, str]]
                       ) -> tuple[set[frozenset], dict[tuple, tuple[str, str]]]:
    """Drop edges whose endpoints aren't both in `nodes`."""
    kept = {e for e in edges if all(x in nodes for x in tuple(e))}
    kept_ep = {k: v for k, v in endpoints.items()
                if k[0] in nodes and k[1] in nodes}
    return kept, kept_ep


# -----------------------------------------------------------------------------
# Divergence
# -----------------------------------------------------------------------------

def _identity(a: str, b: str) -> float:
    """Pairwise nucleotide identity via edlib infix (HW) edit distance,
    strand-agnostic. Query = shorter seq, target = longer (end gaps free).
    We try the target both forward AND reverse-complemented and take the
    closer alignment — graph walks through the same bubble can be emitted
    in either orientation depending on which anchor BFS started from, so
    two alleles that are reverse complements of each other must collapse
    in dedup. Returns 1.0 for identical (or fully contained), lower for
    divergent. Handles indels and end-fragmentation properly."""
    import edlib
    if not a or not b: return 0.0
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    r_fwd = edlib.align(q, t, task="distance", mode="HW")
    r_rc  = edlib.align(q, reverse_complement(t), task="distance", mode="HW")
    best_ed = min(r_fwd["editDistance"], r_rc["editDistance"])
    return 1.0 - best_ed / max(len(q), 1)


def is_divergent(seq_a: str, seq_b: str, threshold: float = 0.05) -> bool:
    """True if the two sequences differ by at least `threshold` (edit-distance
    identity below 1 − threshold). Default 0.05 → pairs with ≥ 5% edit distance
    count as divergent alleles, < 5% collapse as duplicates."""
    if not seq_a or not seq_b: return False
    return (1.0 - _identity(seq_a, seq_b)) >= threshold


# -----------------------------------------------------------------------------
# Orchestrator
# -----------------------------------------------------------------------------

def _try_one_pass(seeds: set[str], all_edges: set[frozenset],
                   endpoints: dict[tuple, tuple[str, str]],
                   adj_und: dict[str, set[str]], n_hops: int,
                   seg_labels: dict[str, list], seg_length: dict[str, int],
                   depths: dict[str, float], gfa_sequences: dict[str, str],
                   genome_cov: float | None,
                   divergence_threshold: float,
                   apply_cov_filter: bool,
                   lo_mult: float, hi_mult: float,
                   max_paths: int = 50,
                   max_path_length: int = 15) -> dict:
    """One pass: BFS-expand seeds N hops, optionally cov-filter, classify,
    return verdict dict with arm sequences attached if it classified."""
    nhood = bfs_expand_segments(seeds, adj_und, n_hops)
    if apply_cov_filter and genome_cov:
        nhood = filter_segs_by_depth(nhood, depths, genome_cov, lo_mult, hi_mult)
    pass_edges, pass_endpoints = restrict_subgraph(nhood, all_edges, endpoints)
    pass_seg_length = {s: seg_length.get(s, len(gfa_sequences.get(s, "")))
                       for s in nhood}
    pass_seg_labels = {s: seg_labels.get(s, []) for s in nhood}
    nodes, edges_pp, labels, var_per, provenance = directional_split(
        pass_seg_labels, pass_seg_length, pass_edges, pass_endpoints,
    )
    res = classify(nodes, edges_pp, labels, var_per,
                    max_paths=max_paths, max_path_length=max_path_length)
    res["_bfs_limits"] = res.get("bfs_limits", {})
    res["_provenance"] = provenance
    res["_nhood"] = nhood
    res["_var_per_node"] = var_per
    res["_label_per_node"] = labels
    # Adjacency over post-P1 nodes — needed by _emit_result to look up the
    # exterior neighbors of each candidate path's endpoints (for the
    # innermost-flank record).
    _adj_pp: dict[str, set[str]] = {}
    for e in edges_pp:
        t = tuple(e)
        if len(t) == 1: continue
        a, b = t
        _adj_pp.setdefault(a, set()).add(b)
        _adj_pp.setdefault(b, set()).add(a)
    res["_adj"] = _adj_pp
    return res


def find_alleles(
        seg_label_hits_tsv: str,
        gfa_path: str,
        genome_cov: float | None = None,
        init_nhop: int = 3,
        max_nhop: int = 10,
        divergence_threshold: float = 0.01,
        lo_mult: float = 0.2,
        hi_mult: float = 2.0,
        k: int | str | None = None,
        var_proteins_ref: str | None = None,
        expected_var_tags: set[str] | None = None,
        locus_padding: int = 4000,
        contig_seeds: set[str] | None = None,
        queries_dir: str | None = None,
        out_candidate_fa: str | None = None,
        seeds_mode: str = "both",            # "flank" | "var" | "both"
        cov_filter: bool = True,             # cov filter ON by default
        max_paths: int = 50,
        max_path_length: int = 15,
        min_allele_bp: int = 3000,           # hard floor on per-allele length
) -> dict:
    """Run the full orchestrator. Returns a dict with the final classification
    and any emitted allele/chimera sequences:

      {
        "verdict": "closed_bubble" | "open_bubble" | "single" | "separate" | "complexed",
        "divergent": bool,
        "n_hops_used": int,
        "phase": "main" | "flank_fallback",
        "alleles": list[(name, sequence)],
        "info": <last classify info>,
      }
    """
    gfa_seqs, depths, all_edges, endpoints = parse_gfa(gfa_path)
    seg_labels = read_seg_label_hits(seg_label_hits_tsv)
    seg_length = read_segment_lengths_from_hits(seg_label_hits_tsv)
    for s, seq in gfa_seqs.items():
        seg_length.setdefault(s, len(seq))

    adj_und: dict[str, set[str]] = {}
    for e in all_edges:
        t = tuple(e)
        if len(t) == 1: continue                            # self-loop in GFA
        a, b = t
        adj_und.setdefault(a, set()).add(b)
        adj_und.setdefault(b, set()).add(a)

    # Seeds: derived from seeds_mode argument.
    #   "var"   → seed BFS from var-labeled GFA segments only
    #   "flank" → seed from flank-labeled segments only
    #   "both"  → union (default — broadest entry into the locus)
    # contig_seeds is added regardless (when populated by upstream).
    contig_seeds = set(contig_seeds or ())
    var_seg_set = {s for s, hits in seg_labels.items()
                    if any(h.kind == "var" for h in hits)}
    flank_seg_set = {s for s, hits in seg_labels.items()
                      if any(h.kind == "flank" for h in hits)}
    if seeds_mode == "var":
        seeds = var_seg_set | contig_seeds
    elif seeds_mode == "flank":
        seeds = flank_seg_set | contig_seeds
    else:  # "both"
        seeds = var_seg_set | flank_seg_set | contig_seeds

    def _log(nhop: int, res: dict) -> None:
        n_arms = res.get("n_arms", 0)
        n_var  = res.get("n_var", 0)
        nhood  = len(res.get("_nhood", ()))
        limits = res.get("_bfs_limits", {})
        lim_str = ""
        if limits.get("max_paths_hit") or limits.get("max_path_length_hit"):
            lim_str = (f"  ⚠ limits: max_paths_hit={limits.get('max_paths_hit',0)}"
                       f" max_path_length_hit={limits.get('max_path_length_hit',0)}")
        print(f"  [nhop={nhop} seeds={seeds_mode} cov={'on' if cov_filter else 'off'}] "
              f"|nhood|={nhood:<6} var={n_var:<3} cls={res['class']:<14} arms={n_arms}{lim_str}",
              flush=True)

    # New design (loop simplified per user spec):
    #   - drop phase 2 (flank_fallback was a no-op in 144/144 picks)
    #   - drop cov-off pass (cov-on won 144/144 picks; cov-off only added noise)
    #   - seeds + cov_filter are configuration, not loop axes
    # The loop is just `for nhop in init_nhop..max_nhop`, with a single
    # variant per nhop. Hard short-circuit on the first complete closed_bubble.
    # All emitted candidates ranked at the end by the unified 4-tier key.

    candidates: list[dict] = []
    iter_metadata: dict[str, dict] = {}
    short_circuit = False

    for nhop in range(init_nhop, max_nhop + 1):
        if short_circuit: break
        iter_id = f"h{nhop}"
        res = _try_one_pass(
            seeds, all_edges, endpoints, adj_und, nhop,
            seg_labels, seg_length, depths, gfa_seqs,
            genome_cov, divergence_threshold,
            apply_cov_filter=cov_filter,
            lo_mult=lo_mult, hi_mult=hi_mult,
            max_paths=max_paths, max_path_length=max_path_length,
        )
        res["_phase"]      = seeds_mode
        res["_n_hops"]     = nhop
        res["_cov_filter"] = cov_filter
        _log(nhop, res)

        if res.get("n_var", 0) == 0 or res["class"] == "no_var":
            continue

        pools_out = _emit_result(
            res, gfa_seqs, depths=depths,
            divergence_threshold=divergence_threshold,
            var_proteins_ref=var_proteins_ref,
            locus_padding=locus_padding,
            expected_var_tags=expected_var_tags,
            genome_cov=genome_cov,
            lo_mult=lo_mult, hi_mult=hi_mult,
            queries_dir=queries_dir,
            return_pools=True,
            min_allele_bp=min_allele_bp,
        )
        pools = pools_out["pools"]
        iter_metadata[iter_id] = {
            "res": res, "pools_out": pools_out,
            "phase": seeds_mode, "nhop": nhop, "cov": cov_filter,
        }

        for net_i, pool in enumerate(pools, start=1):
            if not pool["alleles"]: continue
            net_in_iter = net_i if len(pools) > 1 else 0
            cand = {
                "iter_id":   iter_id,
                "net_in_iter": net_in_iter,
                "phase":     seeds_mode,
                "nhop":      nhop,
                "cov":       cov_filter,
                "verdict":   pool.get("verdict") or res["class"],
                "alleles":   pool["alleles"],
                "allele_cov": pool["allele_cov"],
                "extend_bounds": pool["extend_bounds"],
                "allele_segments": pool["allele_segments"],
                "n_dedup":   pool["n_dedup"],
                "n_raw":     pool["n_raw"],
                # Tri-state at the NETWORK level: 0=none, 1=some, 2=all
                "complete_var":   int(pool.get("complete_var") or 0),
                "complete_locus": int(pool.get("complete_locus") or 0),
                "basepair":  pool["basepair"],
                "diploid_dist": pool["diploid_dist"],
                "found_tags": pool["found_tags_surviving"],
            }
            candidates.append(cand)

            # Acceptance: closed_bubble n≥2 with FULLY complete locus + var
            # (both tri-states at level 2).
            if (cand["verdict"] == "closed_bubble"
                    and cand["n_dedup"] >= 2
                    and cand["complete_locus"] >= 2
                    and cand["complete_var"] >= 2):
                print(f"  [accept] {iter_id} net={net_in_iter}: "
                      f"closed_bubble n={cand['n_dedup']} complete — short-circuiting",
                      flush=True)
                short_circuit = True
                break

    # Write candidate_allele.fasta if requested — every emitted sequence
    # across every (iter, network), tagged with its provenance in the header.
    if out_candidate_fa and candidates:
        _write_candidate_fasta(candidates, out_candidate_fa)

    if not candidates:
        # No iteration produced any emissible content — return an empty result
        # so the caller can write a "no result" row.
        return _empty_find_alleles_result(k, genome_cov)

    # Unified ranking — 4 tiers (smaller value = better):
    #   0. bubble_priority      (K-picker order)
    #   1. complete_locus       DESC   tri-state 2 > 1 > 0
    #   2. complete_var         DESC   tri-state 2 > 1 > 0
    #   3. diploid_dist         ASC    |mean(allele_cov)/D_k − 0.5|
    def _rank_key(c: dict) -> tuple:
        return (
            _bubble_priority(c["verdict"], c["n_dedup"]),    # 0
            -int(c["complete_locus"]),                        # 1 (negate so 2 sorts first)
            -int(c["complete_var"]),                          # 2
            c["diploid_dist"],                                # 3
        )

    best = min(candidates, key=_rank_key)
    print(f"  [pick] best={best['iter_id']} net={best['net_in_iter']} "
          f"verdict={best['verdict']} n={best['n_dedup']} "
          f"complete_locus={best['complete_locus']}", flush=True)

    # Gather every same-iteration sibling (different networks of the same
    # iteration). Cross-iteration mixing is disallowed.
    siblings = [c for c in candidates if c["iter_id"] == best["iter_id"]]

    return _finalize_candidates(
        siblings, iter_metadata[best["iter_id"]],
        gfa_seqs=gfa_seqs, depths=depths,
        divergence_threshold=divergence_threshold,
        expected_var_tags=expected_var_tags,
        genome_cov=genome_cov,
        queries_dir=queries_dir,
        seg_labels=seg_labels,
        k=k,
        min_allele_bp=min_allele_bp,
    )


_BUBBLE_PRIORITY = {
    # Order matches the cross-K picker (pick_k.py): smaller = better.
    # Diploid signatures (closed/open with n=2) preferred over chimeras
    # (complexed with n!=2) and singletons.
    ("closed_bubble", "any"): 1,
    ("open_bubble",   "any"): 2,
    ("separate",      "div2"): 3,
    ("single",        "any"): 4,
    ("complexed",     "div2"): 5,
    ("separate",      "other"): 6,
    ("complexed",     "other"): 7,
    ("no_var",        "any"): 8,
}


def _bubble_priority(verdict: str, n_dedup: int) -> int:
    """Priority tier for cross-iteration ranking — smaller = better."""
    if verdict in ("closed_bubble", "open_bubble", "single", "no_var"):
        return _BUBBLE_PRIORITY[(verdict, "any")]
    if verdict in ("separate", "complexed"):
        cls = "div2" if n_dedup == 2 else "other"
        return _BUBBLE_PRIORITY[(verdict, cls)]
    return 99


def build_longest_alleles_fasta(candidate_fa: str, out_fa: str,
                                  divergence_threshold: float = 0.01) -> int:
    """Read every emission from `candidate_fa` (= candidate_allele.fasta),
    run length-first RC-aware dedup at `divergence_threshold` (5% default),
    and write the surviving sequences — the LONGEST representative of each
    edit-distance equivalence class — to `out_fa`. Returns # records written.

    Different from the primary picker:
      - Operates on the FULL candidate pool (all iterations, all networks)
      - Length-first rank (keep longest per class), not completeness-first
      - Output is the "longest unique walks ever seen" — useful for
        downstream analyses that want a wide net of variants.
    """
    import os
    if not os.path.exists(candidate_fa): return 0
    seqs: list[tuple[str, str]] = []
    cur = None; buf: list[str] = []
    with open(candidate_fa) as fh:
        for ln in fh:
            ln = ln.rstrip()
            if ln.startswith(">"):
                if cur is not None: seqs.append((cur, "".join(buf)))
                cur = ln[1:].split()[0]; buf = []
            else:
                buf.append(ln)
        if cur is not None: seqs.append((cur, "".join(buf)))
    seqs.sort(key=lambda x: -len(x[1]))
    kept: list[tuple[str, str]] = []
    for n, s in seqs:
        if not s: continue
        if any(not is_divergent(s, ks, threshold=divergence_threshold)
                for _, ks in kept):
            continue
        kept.append((n, s))
    with open(out_fa, "w") as fh:
        for n, s in kept:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")
    return len(kept)


def _write_candidate_fasta(candidates: list[dict], out_path: str) -> None:
    """Write every candidate's alleles to a single FASTA. Header carries
    iter_id + network + verdict + cov filter state so the file can be
    forensically diffed against the picker's choice."""
    with open(out_path, "w") as fh:
        for ci, c in enumerate(candidates, start=1):
            cov_tag = "con" if c["cov"] else "coff"
            for name, seq in c["alleles"]:
                hdr = (f">cand{ci:04d}_{c['iter_id']}_n{c['net_in_iter']}"
                       f"_{cov_tag}_{c['verdict']}_{name}")
                fh.write(hdr + "\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")


def _empty_find_alleles_result(k, genome_cov) -> dict:
    """Returned when no iteration produced any emissible candidate at all."""
    return {
        "verdict": "no_var", "topology": "no_var",
        "complete_var": False, "complete_locus": False,
        "locus_coverage": 0.0, "found_var_tags": [],
        "n_candidates": 0, "n_after_dedup": 0,
        "divergent": False,
        "n_hops_used": None, "phase": None, "cov_filter_used": None,
        "alleles": [], "allele_cov": [],
        "extend_bounds": [], "allele_segments": [],
        "subnode_seqs": {},
        "info": {},
        "k": k, "genome_cov": genome_cov,
        "bubble_type": "no_var",
        "component_list": [], "segments": [], "segments_labeled": [],
        "basepair": 0,
    }


def _finalize_candidates(siblings: list[dict], iter_meta: dict,
                          gfa_seqs: dict[str, str],
                          depths: dict[str, float] | None,
                          divergence_threshold: float,
                          expected_var_tags: set[str] | None,
                          genome_cov: float | None,
                          queries_dir: str | None,
                          seg_labels: dict,
                          k,
                          min_allele_bp: int = 3000) -> dict:
    """Take all same-iteration sibling candidates, run cross-network dedup
    (RC-aware, completeness-first ranking), and build the final output dict
    in the same shape that the legacy _emit_result+_finalize path produced.

    Sample-level verdict:
      * surviving sequences from >= 2 distinct networks → "separate"
      * all from 1 network → that network's sub-verdict
    """
    # Flatten: each sibling contributes its (allele, cov, bound, segs)
    # tuples — keep network-of-origin so we can detect cross-network survival.
    # Apply the same hard min_allele_bp floor as in-pool dedup — drops any
    # sub-min-bp fragment that came in via a sibling (e.g. KYH069 n1's 308 bp
    # pair sneaking through because they pass at the sibling-dedup level).
    items = []
    for s in siblings:
        for (name, seq), cov, bnd, segs in zip(
                s["alleles"], s["allele_cov"], s["extend_bounds"], s["allele_segments"]):
            if min_allele_bp > 0 and len(seq) < min_allele_bp: continue
            items.append({
                "name": name, "seq": seq,
                "net": s["net_in_iter"], "cov": cov,
                "bnd": bnd, "segs": segs,
                "sib_complete_var":   s["complete_var"],
                "sib_complete_locus": s["complete_locus"],
                "diploid_dist": s["diploid_dist"],
            })

    # Sort by completeness-first key (same as the in-pool dedup), then walk.
    # The cross-network dedup may collapse two networks' alleles into one if
    # they're within 5% edit distance (RC-aware) — which would mean they're
    # actually duplicate calls of the same allele in two different graph
    # components.
    def _sort_key(it):
        # Unified 4-tier (no bubble_priority — siblings share the iter's
        # verdict so tier 0 is constant; tiers 1–3 differentiate).
        # Tri-state DESC → negate.
        return (-int(it["sib_complete_locus"]),
                -int(it["sib_complete_var"]),
                it["diploid_dist"])
    items.sort(key=_sort_key)

    final = []
    for it in items:
        if any(not is_divergent(it["seq"], kept["seq"], threshold=divergence_threshold)
               for kept in final):
            continue
        final.append(it)

    surviving_nets = {it["net"] for it in final}
    if len(surviving_nets) >= 2:
        sample_verdict = "separate"
    elif final:
        # Find the sibling whose net matches the surviving net to get its verdict
        net = next(iter(surviving_nets))
        matching = next((s for s in siblings if s["net_in_iter"] == net), siblings[0])
        sample_verdict = matching["verdict"]
    else:
        sample_verdict = siblings[0]["verdict"] if siblings else "no_var"

    # Re-emit names: when sample_verdict == "separate" we keep network
    # prefixes; otherwise we strip them (cleaner output) and re-number as
    # allele1/allele2 (or chimera if n>=3).
    final_alleles = []
    final_cov = []; final_bnd = []; final_segs = []
    n_final = len(final)
    for i, it in enumerate(final, start=1):
        if sample_verdict == "separate":
            name = it["name"]                    # keep network-prefixed
        else:
            if n_final == 1: nm = "allele1"
            elif n_final == 2: nm = f"allele{i}"
            else: nm = f"chimera{i}"
            name = nm
        final_alleles.append((name, it["seq"]))
        final_cov.append(it["cov"])
        final_bnd.append(it["bnd"])
        final_segs.append(it["segs"])

    # Pull provenance + sub-seqs from the chosen iteration
    res = iter_meta["res"]
    prov = res.get("_provenance", {})

    surviving_tags = set()
    for s in siblings:
        surviving_tags |= s["found_tags"]

    # Sample-level tri-states: report the MAX across surviving sibling
    # networks (the best per-network completeness this sample achieved).
    if siblings:
        complete_var   = max(int(s.get("complete_var",   0)) for s in siblings)
        complete_locus = max(int(s.get("complete_locus", 0)) for s in siblings)
    else:
        complete_var = 0; complete_locus = 0
    if expected_var_tags:
        locus_coverage = (len(expected_var_tags & surviving_tags)
                          / max(1, len(expected_var_tags)))
    else:
        locus_coverage = None

    # Emitted sub-node IDs → materialized sub-seqs (for run_per_k's
    # subnode_seqs.fasta side file)
    sub_seqs: dict[str, str] = {}
    referenced_ids: set[str] = set()
    for seg_list in final_segs:
        referenced_ids.update(seg_list)
    for sid in referenced_ids:
        if sid not in prov: continue
        parent, sst, eend, strand = prov[sid]
        if parent not in gfa_seqs: continue
        sub = gfa_seqs[parent][sst:eend]
        if not sub: continue
        if strand == "-": sub = reverse_complement(sub)
        sub_seqs[sid] = sub

    # Segments + segment-label list, from the surviving candidates' walks
    emitted_segs: set[str] = set()
    for sl in final_segs:
        for sid in sl:
            if sid in prov: emitted_segs.add(prov[sid][0])
    labels_by_seg: dict[str, str] = {}
    for s in emitted_segs:
        tags = sorted({h.tag for h in seg_labels.get(s, [])})
        if tags: labels_by_seg[s] = "+".join(tags)

    return {
        "verdict": sample_verdict,
        "topology": res.get("class", sample_verdict),
        "complete_var": complete_var,       # tri-state 0/1/2 (max across networks)
        "complete_locus": complete_locus,   # tri-state 0/1/2 (max across networks)
        "locus_coverage": locus_coverage,
        "found_var_tags": sorted(surviving_tags),
        "n_candidates": sum(s["n_raw"] for s in siblings),
        "n_after_dedup": n_final,
        "divergent": n_final >= 2,
        "n_hops_used": iter_meta["nhop"],
        "phase": iter_meta["phase"],
        "cov_filter_used": iter_meta["cov"],
        "alleles": final_alleles,
        "allele_cov": final_cov,
        "extend_bounds": final_bnd,
        "allele_segments": final_segs,
        "subnode_seqs": sub_seqs,
        "info": {},
        "k": k,
        "genome_cov": genome_cov,
        "bubble_type": sample_verdict,
        "component_list": [n for n, _ in final_alleles],
        "segments": sorted(emitted_segs),
        "segments_labeled": [(s, labels_by_seg.get(s, "")) for s in sorted(emitted_segs)],
        "basepair": sum(len(s) for _, s in final_alleles),
    }


def _arm_sequence_for(arm_path: list[str], provenance: dict,
                       gfa_seqs: dict[str, str]) -> str:
    """Build a sequence from one arm path, without joint-skipping.
    Used for the orchestrator's _accept divergence check (full arm content)."""
    parts = []
    for n in arm_path:
        if n not in provenance: continue
        seg, start, end, strand = provenance[n]
        sub = gfa_seqs.get(seg, "")[start:end]
        if strand == "-":
            sub = reverse_complement(sub)
        parts.append(sub)
    return "".join(parts)


def _arm_sequence_for_skip(arm_path: list[str], provenance: dict,
                            gfa_seqs: dict[str, str],
                            skip_nodes: set[str]) -> str:
    """Build a sequence from one arm path, omitting nodes in `skip_nodes`
    (joint/fork nodes shared with other candidate paths). Used by emission
    so that each candidate's sequence reflects only its distinguishing
    (non-shared) content."""
    parts = []
    for n in arm_path:
        if n in skip_nodes: continue
        if n not in provenance: continue
        seg, start, end, strand = provenance[n]
        sub = gfa_seqs.get(seg, "")[start:end]
        if strand == "-":
            sub = reverse_complement(sub)
        parts.append(sub)
    return "".join(parts)


def _dedup_sequences(seqs: list[str],
                      divergence_threshold: float = 0.01) -> list[str]:
    """Collapse non-divergent sequences. Walk in length-desc order so the
    longest representative of each equivalence class is kept."""
    nonempty = [s for s in seqs if s]
    nonempty.sort(key=len, reverse=True)
    kept: list[str] = []
    for s in nonempty:
        if any(not is_divergent(s, k, threshold=divergence_threshold) for k in kept):
            continue
        kept.append(s)
    return kept


def _tblastn_trim_each(
        named_seqs: list[tuple[str, str]],
        var_proteins_ref: str,
        padding: int,
        min_pid: float = 30.0,
        min_aa: int = 50,
) -> tuple[list[tuple[str, str]], dict[str, set[str]], dict[str, tuple[int, int]]]:
    """Run a single tblastn(var proteins → all candidates concat-as-multi-fasta).
    For each candidate, trim to [min hit start − padding, max hit end + padding].
    Drop candidates with no hits.

    Returns (kept_named_seqs, found_tags_by_name, hd_span_in_padded).

    `hd_span_in_padded[name] = (lo, hi)` — coordinates of the HD-only region
    INSIDE the padded `kept_named_seqs` slice (so callers can do
    `seq[lo:hi]` to get the un-padded HD span for var-region-only
    divergence in dedup).
    """
    import subprocess, tempfile, os
    if not named_seqs: return [], {}, {}
    with tempfile.TemporaryDirectory() as td:
        sfa = os.path.join(td, "all_cands.fa")
        with open(sfa, "w") as fh:
            for name, seq in named_seqs:
                fh.write(f">{name}\n{seq}\n")
        out = subprocess.run(
            ["tblastn", "-query", var_proteins_ref, "-subject", sfa,
             "-evalue", "1e-5", "-outfmt",
             "6 qseqid sseqid pident length sstart send"],
            capture_output=True, text=True).stdout
        spans_by: dict[str, list[tuple[int, int]]] = {}
        tags_by:  dict[str, set[str]]               = {}
        for ln in out.splitlines():
            f = ln.split("\t")
            if len(f) >= 6 and float(f[2]) >= min_pid and int(f[3]) >= min_aa:
                qid, sid = f[0], f[1]
                ss, se = int(f[4]), int(f[5])
                spans_by.setdefault(sid, []).append((min(ss, se), max(ss, se)))
                tags_by .setdefault(sid, set()).add(qid)
        kept: list[tuple[str, str]] = []
        hd_span: dict[str, tuple[int, int]] = {}
        for name, seq in named_seqs:
            spans = spans_by.get(name)
            if not spans: continue
            hd_lo = min(s for s, _ in spans)
            hd_hi = max(e for _, e in spans)
            lo = max(0, hd_lo - padding)
            hi = min(len(seq), hd_hi + padding)
            padded = seq[lo:hi]
            kept.append((name, padded))
            # HD coords in the padded slice: (hd_lo - lo, hd_hi - lo)
            hd_span[name] = (hd_lo - lo, hd_hi - lo)
    return kept, tags_by, hd_span


def _dedup_named(named_seqs: list[tuple[str, str]],
                  divergence_threshold: float = 0.01
                  ) -> list[tuple[str, str]]:
    """Length-desc walk, drop seqs whose k-mer Jaccard distance < threshold from
    any already-kept seq. Like _dedup_sequences but preserves the (name, seq) pairing."""
    nonempty = [(n, s) for n, s in named_seqs if s]
    nonempty.sort(key=lambda ns: len(ns[1]), reverse=True)
    kept: list[tuple[str, str]] = []
    for n, s in nonempty:
        if any(not is_divergent(s, ks, threshold=divergence_threshold)
               for _, ks in kept):
            continue
        kept.append((n, s))
    return kept


def _dedup_named_ranked(named_seqs: list[tuple[str, str]],
                         divergence_threshold: float,
                         key_fn,
                         compare_seqs: dict[str, str] | None = None) -> list[tuple[str, str]]:
    """Same as _dedup_named, but the walk order is determined by `key_fn`
    (smaller value = higher priority) instead of length-desc. This lets us
    rank by completeness FIRST so the more-complete representative of each
    equivalence class survives, not just the longest one.

    `compare_seqs`: optional {name: seq_for_comparison}. When given, the
    divergence check uses these sequences instead of the candidate's full
    emitted sequence — e.g. the HD-only span sliced from the padded
    candidate, so dedup focuses on the var region and ignores conserved
    flank context that would otherwise dilute the divergence signal.
    The kept tuple is still (name, emit_seq) — only the comparison changes.
    """
    nonempty = [(n, s) for n, s in named_seqs if s]
    nonempty.sort(key=lambda ns: key_fn(ns[0], ns[1]))
    def _cmp(name: str, seq: str) -> str:
        if compare_seqs is None: return seq
        return compare_seqs.get(name, seq)
    kept: list[tuple[str, str, str]] = []   # (name, emit_seq, cmp_seq)
    for n, s in nonempty:
        cmp_s = _cmp(n, s)
        if any(not is_divergent(cmp_s, k_cmp, threshold=divergence_threshold)
               for _, _, k_cmp in kept):
            continue
        kept.append((n, s, cmp_s))
    return [(n, s) for n, s, _ in kept]


def _blastn_flank_presence(named_seqs: list[tuple[str, str]],
                            queries_dir: str | None
                            ) -> dict[str, tuple[bool, bool]]:
    """Run a single blastn pass with queries/flankL.fasta + flankR.fasta
    against the candidate sequences (one DB built per call). Returns
    {name: (has_flankL, has_flankR)} where each bool is True iff at least
    one hit at pid >= 85% and aln length >= 100 bp survives.

    Returns {name: (False, False)} for all candidates when queries_dir is
    missing, the flank FASTAs are absent, or BLAST exits non-zero — so
    completeness ranking degrades gracefully (every candidate looks
    "incomplete at flanks" and the next tiers do the ranking)."""
    import subprocess, tempfile, os
    presence = {n: (False, False) for n, _ in named_seqs}
    if not queries_dir or not named_seqs: return presence
    flankL_q = os.path.join(queries_dir, "flankL.fasta")
    flankR_q = os.path.join(queries_dir, "flankR.fasta")
    if not (os.path.exists(flankL_q) and os.path.exists(flankR_q)):
        return presence
    with tempfile.TemporaryDirectory() as td:
        sfa = os.path.join(td, "cands.fa")
        with open(sfa, "w") as fh:
            for n, s in named_seqs:
                fh.write(f">{n}\n{s}\n")
        db = os.path.join(td, "candb")
        r = subprocess.run(
            ["makeblastdb", "-in", sfa, "-dbtype", "nucl", "-out", db],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0: return presence
        for q_path, idx in [(flankL_q, 0), (flankR_q, 1)]:
            r = subprocess.run(
                ["blastn", "-query", q_path, "-db", db, "-outfmt", "6",
                 "-evalue", "1e-10", "-dust", "no"],
                capture_output=True, text=True)
            if r.returncode != 0: continue
            for ln in r.stdout.splitlines():
                f = ln.split("\t")
                if len(f) < 12: continue
                sname = f[1]
                try:
                    pid = float(f[2]); aln = int(f[3])
                except ValueError: continue
                if pid >= 85.0 and aln >= 100 and sname in presence:
                    cur = list(presence[sname]); cur[idx] = True
                    presence[sname] = tuple(cur)
    return presence


def _path_mean_cov(path: list[str], provenance: dict,
                    depths: dict[str, float]) -> float:
    """Length-weighted mean DP:f: depth across a path's nodes."""
    total_bp = 0
    total_w  = 0.0
    for n in path:
        if n not in provenance: continue
        seg, s, e, _ = provenance[n]
        L = e - s
        d = depths.get(seg)
        if d is None or L <= 0: continue
        total_bp += L
        total_w  += L * d
    return (total_w / total_bp) if total_bp else 0.0


def _emit_result(res: dict, gfa_seqs: dict[str, str],
                  depths: dict[str, float] | None = None,
                  divergence_threshold: float = 0.01,
                  var_proteins_ref: str | None = None,
                  locus_padding: int = 1500,
                  expected_var_tags: set[str] | None = None,
                  genome_cov: float | None = None,
                  lo_mult: float = 0.25,
                  hi_mult: float = 2.0,
                  queries_dir: str | None = None,
                  return_pools: bool = False,
                  min_allele_bp: int = 3000,
                  force: bool = False) -> dict:
    """Pipe every var-bearing candidate through trim → dedup → emit.

    1. Build raw candidate sequences (full path content, no node-level trim).
    2. If `var_proteins_ref`: tblastn(var_proteins → each candidate), trim each
       to [min_hit − padding, max_hit + padding]. Candidates with no var hits
       are dropped here. This is the locus-trim step, inline with dedup.
    3. Dedup the trimmed sequences (k-mer Jaccard at `divergence_threshold`).
    4. Emit by post-dedup count `n`:
         n == 1 → allele1
         n == 2 → allele1 + allele2
         n >= 3 → chimera1, chimera2, …

    Candidate sources (the classifier populates whichever apply):
      - closed_arms + dangling_arms   (single, open/closed_bubble, complexed)
      - var_components                (separate — disjoint var subgraphs)
    """
    cls = res["class"]
    prov = res.get("_provenance", {})
    var_per_node = res.get("_var_per_node", {})
    var_nodes = {n for n, v in var_per_node.items() if v}
    # Label tokens per post-P1 node (set by directional_split). Carries
    # the per-node tag string like "HD1+HD2" or "flankL"; empty for Uvars.
    label_per_node = res.get("_label_per_node", {})
    adj_pp = res.get("_adj", {})

    def _is_flankL(n: str) -> bool:
        return "flankL" in (label_per_node.get(n, "")).split("+")

    def _is_flankR(n: str) -> bool:
        return "flankR" in (label_per_node.get(n, "")).split("+")

    def _innermost_flank_bounds(p: list[str]) -> tuple[str | None, str | None]:
        """Find innermost flank-labeled nodes bracketing the candidate path.
        Path nodes are bubble-internal (var + unlabeled); flank-labeled
        nodes sit OUTSIDE the bubble, adjacent to its boundary. So we look
        at each endpoint's exterior neighbors:
            L_bound = a flankL-labeled neighbor of path[0]   (or None)
            R_bound = a flankR-labeled neighbor of path[-1]  (or None)
        Emission stays var-trimmed; this is a metadata record only."""
        if not p: return (None, None)
        L_bound = R_bound = None
        for nb in adj_pp.get(p[0], ()):
            if _is_flankL(nb): L_bound = nb; break
        for nb in adj_pp.get(p[-1], ()):
            if _is_flankR(nb): R_bound = nb; break
        # Symmetric check (in case the path's "first" is actually flankR-side
        # due to enumeration order from a different anchor): try swapping.
        if L_bound is None and R_bound is None:
            for nb in adj_pp.get(p[-1], ()):
                if _is_flankL(nb): L_bound = nb; break
            for nb in adj_pp.get(p[0], ()):
                if _is_flankR(nb): R_bound = nb; break
        return (L_bound, R_bound)

    def _seq_full(p): return _arm_sequence_for(p, prov, gfa_seqs)

    def _trim_path_to_var(p: list[str]) -> list[str]:
        """Node-level trim: ordered path → sub-path [first var, last var]."""
        if not p: return p
        var_idx = [i for i, n in enumerate(p) if n in var_nodes]
        if not var_idx: return p
        return p[var_idx[0]:var_idx[-1] + 1]

    def _trim_component_to_var(c: list[str]) -> list[str]:
        """Node-level trim: unordered component → only var nodes."""
        return [n for n in c if n in var_nodes]

    var_per_node = res.get("_var_per_node", {})

    # The component cov-band check is GATED by the same `cov_filter` state
    # the BFS iteration used (res["_cov_filter"]). When the BFS ran with
    # cov-filter on, we apply it to var_components too; when the BFS ran
    # with cov-filter off (the lenient fallback), we don't.
    cov_filter_active = bool(res.get("_cov_filter", False))

    def _component_passes_filter(comp_nodes: list[str]) -> bool:
        """Filter for raw `var_components` emission (rarely used now that
        we recurse classify per network).
        (a) ≥ 2 DISTINCT var tags — always applies (the orphan-fragment guard).
        (b) mean depth in cov-filter band — applied only when the BFS this
            iteration also ran with the cov filter on (consistent state)."""
        tags: set[str] = set()
        for n in comp_nodes:
            t = var_per_node.get(n)
            if t: tags |= set(t)
        if len(tags) < 2:
            return False
        if cov_filter_active and genome_cov and depths:
            lo, hi = lo_mult * genome_cov, hi_mult * genome_cov
            total_bp = 0; weighted = 0.0
            for n in comp_nodes:
                if n not in prov: continue
                seg, s, e, _ = prov[n]
                L = e - s
                d = depths.get(seg)
                if d is None or L <= 0: continue
                total_bp += L; weighted += L * d
            mean_dp = (weighted / total_bp) if total_bp else 0.0
            if mean_dp and not (lo <= mean_dp <= hi):
                return False
        return True

    # ---- Per-network pool builder -------------------------------------------
    #
    # For a `separate` verdict the classifier recurses into each disjoint
    # network and returns its sub-results in res["sub_results"]. We then run
    # the full trim → dedup → emit pipeline INDEPENDENTLY per network and
    # tag each network's emitted alleles with an "n{i}_" prefix, so the
    # output for an N-network separate sample is the concatenation of N
    # per-network allele sets. Dedup does NOT cross networks (alleles from
    # disjoint graph regions stay as distinct records even if their
    # sequences happen to be similar — they represent different loci).
    #
    # For non-separate verdicts there's a single "pool" with prefix = "".

    def _has_flank_in_or_near(path: list[str], side: str) -> bool:
        """Check whether `path` carries a flank label of `side` ('flankL'
        or 'flankR') — either ON one of the path's own nodes, or on one
        of the immediate (post-P1) neighbors of the path's endpoints.

        The latter case is the common one for closed_bubble: after var-
        trim, the path's first/last nodes are var-bearing; the flank-
        labeled bubble anchor sits OUTSIDE the path as a neighbor."""
        for n in path:
            if side in (label_per_node.get(n, "")).split("+"):
                return True
        if path:
            for ep in (path[0], path[-1]):
                for m in adj_pp.get(ep, set()):
                    if side in (label_per_node.get(m, "")).split("+"):
                        return True
        return False

    def _process_pool(closed_arms, dangling_arms, var_components, prefix,
                       pool_verdict: str = None):
        """Build raw candidates → locus-trim → dedup → emit names.
        Returns dict with keys: alleles, allele_cov, allele_segments,
        extend_bounds, n_raw, n_dedup, found_tags_surviving,
        verdict (per-pool sub-verdict), complete_var, complete_locus,
        basepair, diploid_dist (for cross-iteration ranking)."""
        local_raw: list[tuple[str, str]] = []
        local_path: dict[str, list[str]] = {}
        # Keep pre-trim paths (for graph-level flank-presence checks). After
        # the var-trim, the path's nodes are all var-bearing — flank labels
        # only show up on the original arm (and via post-P1 neighbors).
        local_path_pre: dict[str, list[str]] = {}
        for j, p in enumerate(closed_arms):
            tp = _trim_path_to_var(list(p))
            nm = f"{prefix}cl{j}"
            local_raw.append((nm, _seq_full(tp))); local_path[nm] = tp
            local_path_pre[nm] = list(p)
        for j, p in enumerate(dangling_arms):
            tp = _trim_path_to_var(list(p))
            nm = f"{prefix}da{j}"
            local_raw.append((nm, _seq_full(tp))); local_path[nm] = tp
            local_path_pre[nm] = list(p)
        for j, c in enumerate(var_components or []):
            tc = _trim_component_to_var(list(c))
            if not _component_passes_filter(tc): continue
            nm = f"{prefix}vc{j}"
            local_raw.append((nm, _arm_sequence_for(tc, prov, gfa_seqs)))
            local_path[nm] = tc
            local_path_pre[nm] = list(c)

        n_raw_pool = len(local_raw)

        # Hard minimum-length floor — applied as the FIRST dedup filter,
        # BEFORE tblastn locus-trim. Drops fragment candidates (e.g. the
        # KYH069 n1 308bp segments) before paying the tblastn cost.
        if min_allele_bp > 0:
            local_raw = [(n, s) for n, s in local_raw if len(s) >= min_allele_bp]

        # Locus trim — drops no-HD candidates and trims survivors to
        # [min_hit − padding, max_hit + padding]. Also returns the
        # HD-only span coords WITHIN each padded slice so dedup can
        # compare on the var-region-only sequence (flanks excluded).
        local_found: dict[str, set[str]] = {}
        local_hd_span: dict[str, tuple[int, int]] = {}
        if var_proteins_ref:
            local_raw, local_found, local_hd_span = _tblastn_trim_each(
                local_raw, var_proteins_ref, locus_padding)

        # Second min_allele_bp pass — after locus-trim, in case the trim
        # window narrowed a survivor below the floor.
        if min_allele_bp > 0:
            local_raw = [(n, s) for n, s in local_raw if len(s) >= min_allele_bp]

        depths_d = depths or {}
        _cov  = lambda rn: _path_mean_cov(local_path.get(rn, []), prov, depths_d)
        _bnds = lambda rn: _innermost_flank_bounds(local_path.get(rn, []))

        # Completeness ranking for dedup. Per candidate, compute:
        #   complete_var   = expected_var_tags ⊆ found_var_tags(this candidate)
        #   complete_locus = complete_var AND has_flankL AND has_flankR
        # Flank presence comes from the path's GRAPH-LEVEL labels (the
        # pre-trim arm carries flank-labeled nodes as endpoints/anchors),
        # NOT from a blastn over the locus-trimmed sequence — the trim
        # window strips off the flanks by design, so blastn would always
        # report False here.
        flank_presence = {
            rn: (_has_flank_in_or_near(local_path_pre.get(rn, []), "flankL"),
                 _has_flank_in_or_near(local_path_pre.get(rn, []), "flankR"))
            for rn, _ in local_raw
        }

        def _rank_key(name: str, seq: str) -> tuple:
            """Unified 4-tier rank, tri-state form. bubble_priority lives at
            the candidate (pool) level — within a pool all entries share the
            verdict so tier 0 is constant; tiers 1–3 differentiate.

            Per-CANDIDATE tri-states (computed here):
              cv_level  = 0/1/2 of expected_var_tags found in this candidate
              cl_level  = 0/1/2 of {flankL, flankR} present in this candidate
            """
            fv = local_found.get(name, set())
            if expected_var_tags:
                hit = fv & expected_var_tags
                if not hit:                          cv_level = 0
                elif hit == set(expected_var_tags):  cv_level = 2
                else:                                cv_level = 1
            else:
                cv_level = 0
            hL, hR = flank_presence.get(name, (False, False))
            cl_level = int(hL) + int(hR)
            cov = _cov(name)
            diploid_dist = abs(cov / genome_cov - 0.5) if genome_cov else 1.0
            return (-cl_level, -cv_level, diploid_dist)

        # Dedup within this pool — divergence is computed on the HD-ONLY
        # slice (var region without padding) so flank conservation can't
        # dilute the var-region divergence signal. The emitted sequences
        # remain the full padded versions.
        hdonly_by = {
            name: padded_seq[lo:hi]
            for (name, padded_seq), (lo, hi) in (
                (ns, local_hd_span.get(ns[0], (0, len(ns[1]))))
                for ns in local_raw
            )
        }
        local_dedup = _dedup_named_ranked(
            local_raw, divergence_threshold, _rank_key,
            compare_seqs=hdonly_by,
        )
        nd = len(local_dedup)

        # Joint detection via multi-source BFS in the post-P1 graph.
        # "Joints" are the cycle convergence points where ≥ 2 surviving
        # arms meet when extended outward through the FULL graph — they
        # are the topological boundary of the bubble's cycle, defined by
        # graph structure alone, NOT by labels (so an unlabeled hub-node
        # joint is detected just as well as a flank-labeled one).
        #
        # Each surviving arm's full walk = [joint_L, ..., arm_internal,
        # ..., joint_R]. The emission FASTA uses the first-non-joint to
        # last-non-joint slice (joints stripped at the ends only), then
        # tblastn-trim. Dedup ranking + divergence compare still use
        # `local_raw` (var-trimmed → HD-only slice) — emission is decoupled
        # from dedup.
        survivor_pre: dict[str, list[str]] = {
            rn: list(local_path_pre.get(rn, []) or [])
            for rn, _ in local_dedup
        }
        arm_internals: dict[str, set[str]] = {
            rn: set(p) for rn, p in survivor_pre.items()
        }

        def _bfs_outward(src: str, exclude: set[str]) -> dict[str, tuple[int, str | None]]:
            """BFS in adj_pp from src, never entering nodes in `exclude`.
            Returns {node: (dist, prev)} for every node reached (incl. src)."""
            visited: dict[str, tuple[int, str | None]] = {src: (0, None)}
            front: list[str] = [src]
            while front:
                nxt: list[str] = []
                for u in front:
                    for v in adj_pp.get(u, ()):
                        if v in visited or v in exclude: continue
                        visited[v] = (visited[u][0] + 1, u)
                        nxt.append(v)
                front = nxt
            return visited

        # Per-arm BFS from each endpoint. For a single-node arm, both
        # endpoints are the same node — one BFS, shared as L and R.
        arm_bfs: dict[str, tuple[dict, dict]] = {}
        for rn, p in survivor_pre.items():
            if not p:
                arm_bfs[rn] = ({}, {}); continue
            excl_l = arm_internals[rn] - {p[0]}
            vl = _bfs_outward(p[0], excl_l)
            if len(p) == 1:
                arm_bfs[rn] = (vl, vl)
            else:
                excl_r = arm_internals[rn] - {p[-1]}
                vr = _bfs_outward(p[-1], excl_r)
                arm_bfs[rn] = (vl, vr)

        # Candidate joints: nodes reached by BFSes of ≥ 2 distinct arms,
        # excluding any arm's internal nodes.
        all_internal = set().union(*arm_internals.values()) if arm_internals else set()
        joint_visits: dict[str, dict[str, int]] = {}
        for rn, (vl, vr) in arm_bfs.items():
            for n_id in (set(vl.keys()) | set(vr.keys())) - all_internal:
                d_l = vl.get(n_id, (10**9, None))[0]
                d_r = vr.get(n_id, (10**9, None))[0]
                joint_visits.setdefault(n_id, {})[rn] = min(d_l, d_r)
        joint_cands: dict[str, dict[str, int]] = {
            n_id: dists for n_id, dists in joint_visits.items()
            if len(dists) >= 2
        }

        def _trace(visits: dict[str, tuple[int, str | None]],
                    target: str) -> list[str] | None:
            """Reconstruct path from BFS source to `target` via `visits`."""
            if target not in visits: return None
            out: list[str] = []; cur: str | None = target
            while cur is not None:
                out.append(cur); cur = visits[cur][1]
            return list(reversed(out))   # source at out[0], target at out[-1]

        # Pick 2 joints by minimizing TOTAL walk length across all arms over
        # the joint-pair. Picking simply by min(max_dist) tie-breaks badly
        # in true closed-bubble graphs where ≥ 3 nodes hit the same max_dist
        # (e.g., vietnam_BD1417 k45: {flankL_anchor, flankR_anchor,
        # arm2_flankL_subseg} all at max=2). The flank anchors sit on
        # opposite sides of the cycle (short walks); the false candidate is
        # one arm's own flank sub-seg (only reachable from the OTHER arm by
        # walking all the way around → long walks). Pair-total-length
        # selection naturally prefers the short-walk pair.
        #
        # Restrict candidate set to the K smallest-max_dist nodes to bound
        # the O(k²) pair enumeration; K = 20 is comfortably larger than any
        # real-graph candidate count we've seen.
        K_MAX = 20
        sorted_j = sorted(joint_cands.keys(),
                          key=lambda n: max(joint_cands[n].values()))
        cand_set = sorted_j[:K_MAX]

        def _pair_total(j_x: str, j_y: str) -> tuple[float, dict[str, tuple]]:
            """Total walk length for (j_x, j_y) across all surviving arms.
            Returns (total, per_arm_chosen_options). Each option is
            (j_left_for_arm, j_right_for_arm, pl, pr) for multi-node arms,
            or (j_a, j_b, pa, pb) for single-node arms. inf if any arm
            can't reach both joints."""
            total = 0; per_arm: dict[str, tuple] = {}
            for rn, p in survivor_pre.items():
                if not p:
                    per_arm[rn] = (None, None, None, None); continue
                vl, vr = arm_bfs[rn]
                if len(p) == 1:
                    pa = _trace(vl, j_x); pb = _trace(vl, j_y)
                    if not pa or not pb: return float("inf"), {}
                    per_arm[rn] = (j_x, j_y, pa, pb)
                    total += len(pa) + len(pb) - 1
                    continue
                opts = []
                for jl, jr in ((j_x, j_y), (j_y, j_x)):
                    pl = _trace(vl, jl); pr = _trace(vr, jr)
                    if pl and pr:
                        opts.append((len(pl) + len(pr), jl, jr, pl, pr))
                if not opts: return float("inf"), {}
                opts.sort()
                w, jl, jr, pl, pr = opts[0]
                per_arm[rn] = (jl, jr, pl, pr)
                total += w
            return total, per_arm

        chosen_joints: set[str] = set()
        j_a: str | None = None
        j_b: str | None = None
        best_per_arm: dict[str, tuple] = {}
        if len(survivor_pre) >= 2 and len(cand_set) >= 2:
            best_total = float("inf")
            for i in range(len(cand_set)):
                for j in range(i + 1, len(cand_set)):
                    j_x, j_y = cand_set[i], cand_set[j]
                    total, per_arm = _pair_total(j_x, j_y)
                    if total < best_total:
                        best_total = total; j_a, j_b = j_x, j_y
                        best_per_arm = per_arm
            if j_a is not None:
                chosen_joints.add(j_a); chosen_joints.add(j_b)
        elif len(survivor_pre) >= 2 and len(cand_set) == 1:
            j_a = cand_set[0]; chosen_joints.add(j_a)

        def _flank_extend(p: list[str]) -> tuple[list[str], set[str]]:
            """Legacy fallback when no joint pair exists (single-arm pool
            or disconnected): prepend the L-side flank neighbor and append
            the R-side flank neighbor of the arm endpoints. Returns the
            extended walk AND the set of added neighbors (treated as
            implicit joints for the non-joint emission slice)."""
            walk = list(p); added: set[str] = set()
            if not p: return walk, added
            Lb_s = next((nb for nb in adj_pp.get(p[0],  ()) if _is_flankL(nb)), None)
            Rb_e = next((nb for nb in adj_pp.get(p[-1], ()) if _is_flankR(nb)), None)
            Lb_e = next((nb for nb in adj_pp.get(p[-1], ()) if _is_flankL(nb)), None)
            Rb_s = next((nb for nb in adj_pp.get(p[0],  ()) if _is_flankR(nb)), None)
            if Lb_s or Rb_e:
                if Lb_s: walk.insert(0, Lb_s); added.add(Lb_s)
                if Rb_e: walk.append(Rb_e);   added.add(Rb_e)
            elif Lb_e or Rb_s:
                if Rb_s: walk.insert(0, Rb_s); added.add(Rb_s)
                if Lb_e: walk.append(Lb_e);   added.add(Lb_e)
            return walk, added

        # Build full_walk for every survivor.
        full_walk: dict[str, list[str]] = {}
        for rn, p in survivor_pre.items():
            if not p:
                full_walk[rn] = []; continue
            if j_a is None:
                walk, added = _flank_extend(p)
                full_walk[rn] = walk
                chosen_joints |= added
                continue
            vl, vr = arm_bfs[rn]
            if j_b is None:
                # Only one joint candidate — extend whichever side reaches it.
                t_l = _trace(vl, j_a)
                t_r = _trace(vr, j_a) if vr is not vl else None
                if t_l and (not t_r or len(t_l) <= len(t_r)):
                    full_walk[rn] = list(reversed(t_l))[:-1] + list(p)
                elif t_r:
                    full_walk[rn] = list(p) + t_r[1:]
                else:
                    walk, added = _flank_extend(p)
                    full_walk[rn] = walk
                    chosen_joints |= added
                continue
            # Two joints: reuse the per-arm pairing chosen by _pair_total
            # so the global "min total walk length" decision drives every
            # arm's individual walk too.
            chosen = best_per_arm.get(rn)
            if chosen is None:
                walk, added = _flank_extend(p)
                full_walk[rn] = walk
                chosen_joints |= added
                continue
            jl, jr, pl, pr = chosen
            if len(p) == 1:
                # pl = [ep, ..., j_a]; pr = [ep, ..., j_b]
                full_walk[rn] = list(reversed(pl)) + pr[1:]
            else:
                # pl = [p[0], ..., j_l]; pr = [p[-1], ..., j_r]
                left_ext  = list(reversed(pl))[:-1]
                right_ext = pr[1:]
                full_walk[rn] = left_ext + list(p) + right_ext

        # seg_processor.directional_split now emits clean position-indexed
        # sub-node IDs ("{parent}#1", "{parent}#2", "{parent}#N"). The
        # provenance dict carries the coords for sequence materialization;
        # the IDs themselves stay coord-free for clean display in bubble.*.
        # `_segs` returns the FULL WALK (joint-to-joint) so result.tsv /
        # bubble.* show the complete bubble cycle. Emission FASTA uses the
        # non-joint slice (joints stripped at the ends) — see below.
        def _segs(rn):
            out, seen = [], set()
            for sub in full_walk.get(rn, []):
                if sub not in prov: continue
                if sub in seen: continue
                seen.add(sub); out.append(sub)
            return out

        def _nonjoint_slice(rn) -> list[str]:
            walk = full_walk.get(rn, [])
            nj = [i for i, n in enumerate(walk) if n not in chosen_joints]
            if not nj: return list(walk)
            return walk[nj[0]: nj[-1] + 1]

        def _emit_seq(rn) -> str:
            return _arm_sequence_for(_nonjoint_slice(rn), prov, gfa_seqs)

        # Pre-build raw emission seqs (no tblastn trim yet) for survivors,
        # then apply tblastn trim once over the survivor set. Falls back to
        # the raw emission seq if tblastn drops the candidate.
        emit_raw: list[tuple[str, str]] = [(rn, _emit_seq(rn))
                                            for rn, _ in local_dedup]
        if var_proteins_ref and emit_raw:
            emit_trimmed, _ef, _es = _tblastn_trim_each(
                emit_raw, var_proteins_ref, locus_padding)
            emit_by = {n: s for n, s in emit_trimmed}
            # Fallback: keep raw emission for any survivor tblastn dropped.
            for rn, rs in emit_raw:
                emit_by.setdefault(rn, rs)
        else:
            emit_by = dict(emit_raw)

        # Per-pool label scheme: n=1→allele1, n=2→allele1/2, n>=3→chimeraN
        if nd <= 1:
            if local_dedup:
                rn, _ = local_dedup[0]
                alleles_p     = [(f"{prefix}allele1", emit_by.get(rn, ""))]
                cov_p         = [_cov(rn)]
                bounds_p      = [_bnds(rn)]
                segs_p        = [_segs(rn)]
            else:
                alleles_p = []; cov_p = []; bounds_p = []; segs_p = []
        elif nd == 2:
            alleles_p = [(f"{prefix}allele1", emit_by.get(local_dedup[0][0], "")),
                          (f"{prefix}allele2", emit_by.get(local_dedup[1][0], ""))]
            cov_p     = [_cov(local_dedup[0][0]), _cov(local_dedup[1][0])]
            bounds_p  = [_bnds(local_dedup[0][0]), _bnds(local_dedup[1][0])]
            segs_p    = [_segs(local_dedup[0][0]), _segs(local_dedup[1][0])]
        else:
            alleles_p = [(f"{prefix}chimera{k+1}", emit_by.get(rn, ""))
                         for k, (rn, _) in enumerate(local_dedup)]
            cov_p    = [_cov(rn) for rn, _ in local_dedup]
            bounds_p = [_bnds(rn) for rn, _ in local_dedup]
            segs_p   = [_segs(rn) for rn, _ in local_dedup]

        surviving = set()
        for rn, _ in local_dedup:
            surviving |= local_found.get(rn, set())

        # Per-ALLELE tri-state, then MIN-aggregated to pool level.
        # 0/1/2 = none/some/all for the SET in question:
        #   complete_var   : expected_var_tags ⊆ THIS allele's found tags
        #   complete_locus : {flankL, flankR}   ⊆ THIS allele's blastn hits
        # MIN across alleles = pool scores "all" only when EVERY emitted allele
        # individually scores "all". Avoids the union-pathology where two
        # single-HD fragments union to look fully complete (the KYH069
        # n1 308 bp pair case).
        def _allele_cv(rn):
            if not expected_var_tags: return 0
            fv = local_found.get(rn, set())
            hit = fv & expected_var_tags
            if not hit: return 0
            if hit == set(expected_var_tags): return 2
            return 1
        def _allele_cl(rn):
            hL, hR = flank_presence.get(rn, (False, False))
            return int(hL) + int(hR)
        if local_dedup:
            per_allele_cv = [_allele_cv(rn) for rn, _ in local_dedup]
            per_allele_cl = [_allele_cl(rn) for rn, _ in local_dedup]
            complete_var_p   = min(per_allele_cv)
            complete_locus_p = min(per_allele_cl)
        else:
            complete_var_p = 0; complete_locus_p = 0
        basepair_p = sum(len(s) for _, s in alleles_p)
        # Diploid signature: closer to ½ × genome_cov is better. Use the
        # mean of all emitted alleles' covs vs genome_cov.
        if genome_cov and cov_p:
            diploid_dist_p = abs((sum(cov_p) / len(cov_p)) / genome_cov - 0.5)
        else:
            diploid_dist_p = 1.0

        return {
            "alleles": alleles_p, "allele_cov": cov_p,
            "extend_bounds": bounds_p, "allele_segments": segs_p,
            "n_raw": n_raw_pool, "n_dedup": nd,
            "found_tags_surviving": surviving,
            "verdict": pool_verdict,                  # per-pool sub-verdict
            "complete_var": complete_var_p,
            "complete_locus": complete_locus_p,
            "basepair": basepair_p,
            "diploid_dist": diploid_dist_p,
            # raw (name, seq) of the dedup-surviving candidates — used to
            # reconstruct flank-presence/coverage info during cross-network
            # dedup at the orchestrator level.
            "raw_dedup": list(local_dedup),
        }

    # Build pools — one per network for `separate`, one global pool otherwise.
    # Each pool carries its sub-verdict (never "separate" at the pool level —
    # that label only emerges at the cross-pool aggregation step).
    sub_results = res.get("sub_results")
    if sub_results:
        pools = []
        for i, sub in enumerate(sub_results):
            pools.append(_process_pool(
                sub.get("closed_arms", []) or [],
                sub.get("dangling_arms", []) or [],
                sub.get("var_components", []) or [],
                prefix=f"n{i+1}_",
                pool_verdict=sub.get("class"),
            ))
    else:
        pools = [_process_pool(
            res.get("closed_arms", []) or [],
            res.get("dangling_arms", []) or [],
            res.get("var_components", []) or [],
            prefix="",
            pool_verdict=cls,
        )]

    # Find-alleles wants to inspect per-pool candidates BEFORE merging — used
    # for cross-iteration ranking. The `return_pools` shortcut returns the
    # raw pool list (already trimmed/deduped) and lets the orchestrator
    # build its own emission.
    if return_pools:
        return {"pools": pools, "topology": cls, "sub_results": sub_results,
                "_provenance": prov, "_var_per_node": var_per_node,
                "_label_per_node": label_per_node}

    alleles         = [a for p in pools for a in p["alleles"]]
    allele_cov      = [c for p in pools for c in p["allele_cov"]]
    extend_bounds   = [b for p in pools for b in p["extend_bounds"]]
    allele_segments = [s for p in pools for s in p["allele_segments"]]
    n_raw           = sum(p["n_raw"]   for p in pools)
    n               = sum(p["n_dedup"] for p in pools)
    surviving_tags: set[str] = set()
    for p in pools: surviving_tags |= p["found_tags_surviving"]

    depths = depths or {}

    if expected_var_tags is not None:
        complete_var = expected_var_tags <= surviving_tags
        locus_coverage = len(expected_var_tags & surviving_tags) / max(1, len(expected_var_tags))
    else:
        complete_var = None
        locus_coverage = None

    # Decide the final verdict from the actual network sourcing of emissions:
    #   * alleles span ≥ 2 distinct networks    → "separate"  (HARD OVERRIDE,
    #                                              wins over any sub-verdict like
    #                                              open_bubble / closed_bubble /
    #                                              single — the sample really is
    #                                              split across disconnected
    #                                              graph regions)
    #   * alleles all from ONE network          → promote that network's
    #                                              sub-verdict (so a separate
    #                                              sample where only one network
    #                                              had real alleles is reported
    #                                              as the network's own shape:
    #                                              closed_bubble / open_bubble /
    #                                              complexed / single)
    #   * no networks (no sub_results, normal)  → use the classifier's top-level
    #                                              verdict as-is
    final_verdict = cls
    if sub_results:
        emit_nets = set()
        for name, _ in alleles:
            if name.startswith("n") and "_" in name:
                emit_nets.add(name.split("_", 1)[0])
        if len(emit_nets) >= 2:
            final_verdict = "separate"
        elif len(emit_nets) == 1:
            net_idx = int(next(iter(emit_nets))[1:]) - 1
            if 0 <= net_idx < len(sub_results):
                final_verdict = sub_results[net_idx].get("class", cls)
        # else: no network emitted anything — keep classifier's "separate"

    # Materialize each unique sub-node ID referenced by any emitted allele.
    # The IDs themselves are coord-free ("{parent}#N"); coords live in `prov`.
    # Returned as a small {sub_id: sub_sequence} dict so run_per_k can
    # write a per-k FASTA that graph_paths reads when assembling
    # bubble.txt / bubble.gfa / bubble.png.
    sub_seqs: dict[str, str] = {}
    referenced_ids: set[str] = set()
    for seg_list in allele_segments:
        referenced_ids.update(seg_list)
    for sid in referenced_ids:
        if sid not in prov: continue
        parent, s, e, strand = prov[sid]
        if parent not in gfa_seqs: continue
        sub = gfa_seqs[parent][s:e]
        if not sub: continue
        if strand == "-": sub = reverse_complement(sub)
        sub_seqs[sid] = sub

    return {
        "verdict": final_verdict,
        "topology": cls,                                                # raw classifier top-level (always "separate" when sub_results present)
        "complete_var": complete_var,
        "complete_locus": complete_var,                                # same notion: did we recover all var tags
        "locus_coverage": locus_coverage,
        "found_var_tags": sorted(surviving_tags),
        "n_candidates": n_raw,                                         # before dedup (pre-trim count)
        "n_after_dedup": n,
        "divergent": n >= 2,                                           # at least 2 distinct allele sequences
        "n_hops_used": res.get("_n_hops"),
        "phase": res.get("_phase"),
        "cov_filter_used": res.get("_cov_filter"),
        "alleles": alleles,
        "allele_cov": allele_cov,                                       # length-weighted DP:f: mean per allele
        "extend_bounds": extend_bounds,                                 # per allele (innermost_flankL, innermost_flankR)
        "allele_segments": allele_segments,                             # per allele: list of GFA seg IDs walked (path order, deduped)
        "subnode_seqs": sub_seqs,                                       # {sub_id: materialized sub-region sequence} for split parents
        "info": {k: v for k, v in res.items() if not k.startswith("_")
                 and k not in ("closed_arms", "dangling_arms", "var_components")},
    }


def write_fasta(alleles: list[tuple[str, str]], out_path: str,
                 sample: str = "") -> int:
    """Write the emitted alleles/chimeras to a FASTA file."""
    with open(out_path, "w") as fh:
        for name, seq in alleles:
            fh.write(f">{sample + '_' if sample else ''}{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    return len(alleles)
