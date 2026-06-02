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
                   lo_mult: float, hi_mult: float) -> dict:
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
    res = classify(nodes, edges_pp, labels, var_per)
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
        init_nhop: int = 5,
        max_nhop: int = 10,
        divergence_threshold: float = 0.05,
        lo_mult: float = 0.25,
        hi_mult: float = 2.0,
        k: int | str | None = None,
        var_proteins_ref: str | None = None,
        expected_var_tags: set[str] | None = None,
        locus_padding: int = 1500,
        contig_seeds: set[str] | None = None,
        queries_dir: str | None = None,
        out_candidate_fa: str | None = None,
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

    # Seeds for the two phases. Optional contig_seeds = GFA segments walked
    # by SPAdes contigs that BLAST-hit the locus (HD/flank); used to anchor
    # the BFS in samples where the labeler's segment-level hits split the
    # locus across disconnected GFA components.
    contig_seeds = set(contig_seeds or ())
    var_seeds = {s for s, hits in seg_labels.items()
                 if any(h.kind == "var" for h in hits)} | contig_seeds
    flank_seeds = {s for s, hits in seg_labels.items()
                   if any(h.kind == "flank" for h in hits)} | contig_seeds

    def _log(phase: str, nhop: int, use_cov: bool, res: dict) -> None:
        n_arms = res.get("n_arms", 0)
        n_var  = res.get("n_var", 0)
        nhood  = len(res.get("_nhood", ()))
        print(f"  [nhop={nhop} {phase} cov={'on' if use_cov else 'off':>3}] "
              f"|nhood|={nhood:<6} var={n_var:<3} cls={res['class']:<14} arms={n_arms}",
              flush=True)

    # New design (item 1 in BFS proposal):
    #   outer  = nhop                            (5..10)
    #   middle = phase                           (main, flank_fallback)
    #   inner  = use_cov                         (True, False)
    # Each (nhop, phase, cov) iteration produces per-network candidates.
    # Acceptance check (α): closed_bubble + 2 divergent arms + complete_locus.
    # If ANY candidate accepts at ANY iteration, all remaining iterations
    # are skipped — fast path for clean diploid samples.
    #
    # When loop exhausts without acceptance, pick best candidate by the
    # 6-tier ranking (bubble_priority > complete_locus > complete_var >
    # cov_on > basepair > diploid_dist), gather all SAME-iteration
    # siblings, then run a cross-network dedup to produce final alleles.
    # Sample verdict = "separate" iff ≥ 2 surviving networks; else the
    # single surviving network's sub-verdict.

    candidates: list[dict] = []   # one per (iter, network)
    iter_metadata: dict[str, dict] = {}
    short_circuit = False

    for nhop in range(init_nhop, max_nhop + 1):
        if short_circuit: break
        for phase, seeds in [("main",            var_seeds),
                              ("flank_fallback",  var_seeds | flank_seeds)]:
            if short_circuit: break
            for use_cov in (True, False):
                iter_id = f"h{nhop}_{phase[0]}_c{'on' if use_cov else 'off'}"
                res = _try_one_pass(
                    seeds, all_edges, endpoints, adj_und, nhop,
                    seg_labels, seg_length, depths, gfa_seqs,
                    genome_cov, divergence_threshold,
                    apply_cov_filter=use_cov,
                    lo_mult=lo_mult, hi_mult=hi_mult,
                )
                res["_phase"]      = phase
                res["_n_hops"]     = nhop
                res["_cov_filter"] = use_cov
                _log(phase, nhop, use_cov, res)

                # Skip if classifier produced no var-bearing content
                if res.get("n_var", 0) == 0 or res["class"] == "no_var":
                    continue

                # Run trim + dedup + pool-build (per-network) — return only the
                # pool list without merging.
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
                )
                pools = pools_out["pools"]
                iter_metadata[iter_id] = {
                    "res": res, "pools_out": pools_out,
                    "phase": phase, "nhop": nhop, "cov": use_cov,
                }

                # Build a candidate per non-empty pool (per network)
                for net_i, pool in enumerate(pools, start=1):
                    if not pool["alleles"]: continue
                    # When there's only one pool (non-separate case), the
                    # network index is 0 by convention (no per-network split).
                    net_in_iter = net_i if len(pools) > 1 else 0
                    cand = {
                        "iter_id":   iter_id,
                        "net_in_iter": net_in_iter,
                        "phase":     phase,
                        "nhop":      nhop,
                        "cov":       use_cov,
                        "verdict":   pool.get("verdict") or res["class"],
                        "alleles":   pool["alleles"],
                        "allele_cov": pool["allele_cov"],
                        "extend_bounds": pool["extend_bounds"],
                        "allele_segments": pool["allele_segments"],
                        "n_dedup":   pool["n_dedup"],
                        "n_raw":     pool["n_raw"],
                        "complete_var":   bool(pool.get("complete_var")),
                        "complete_locus": bool(pool.get("complete_locus")),
                        "basepair":  pool["basepair"],
                        "diploid_dist": pool["diploid_dist"],
                        "found_tags": pool["found_tags_surviving"],
                    }
                    candidates.append(cand)

                    # New acceptance check (α): hard short-circuit on gold-
                    # standard candidate. Pool-level "verdict" is the per-
                    # network sub-class (never "separate").
                    if (cand["verdict"] == "closed_bubble"
                            and cand["n_dedup"] >= 2
                            and cand["complete_locus"]):
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

    # Pick best candidate by the 6-tier ranking.
    def _rank_key(c: dict) -> tuple:
        return (
            _bubble_priority(c["verdict"], c["n_dedup"]),    # 0
            not c["complete_locus"],                          # 1
            not c["complete_var"],                            # 2
            not c["cov"],                                     # 3
            -c["basepair"],                                   # 4
            c["diploid_dist"],                                # 5
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
                          k) -> dict:
    """Take all same-iteration sibling candidates, run cross-network dedup
    (RC-aware, completeness-first ranking), and build the final output dict
    in the same shape that the legacy _emit_result+_finalize path produced.

    Sample-level verdict:
      * surviving sequences from >= 2 distinct networks → "separate"
      * all from 1 network → that network's sub-verdict
    """
    # Flatten: each sibling contributes its (allele, cov, bound, segs)
    # tuples — keep network-of-origin so we can detect cross-network survival.
    items = []
    for s in siblings:
        for (name, seq), cov, bnd, segs in zip(
                s["alleles"], s["allele_cov"], s["extend_bounds"], s["allele_segments"]):
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
        return (not it["sib_complete_locus"],
                not it["sib_complete_var"],
                -len(it["seq"]),
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

    if expected_var_tags is not None:
        complete_var = expected_var_tags <= surviving_tags
        locus_coverage = (len(expected_var_tags & surviving_tags)
                          / max(1, len(expected_var_tags)))
    else:
        complete_var = None
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
        "complete_var": complete_var,
        "complete_locus": complete_var,
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
                      divergence_threshold: float = 0.05) -> list[str]:
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
) -> tuple[list[tuple[str, str]], dict[str, set[str]]]:
    """Run a single tblastn(var proteins → all candidates concat-as-multi-fasta).
    For each candidate, trim to [min hit start − padding, max hit end + padding].
    Drop candidates with no hits.

    Returns (kept_named_seqs, found_tags_by_name).
    """
    import subprocess, tempfile, os
    if not named_seqs: return [], {}
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
        for name, seq in named_seqs:
            spans = spans_by.get(name)
            if not spans: continue
            lo = max(0, min(s for s, _ in spans) - padding)
            hi = min(len(seq), max(e for _, e in spans) + padding)
            kept.append((name, seq[lo:hi]))
    return kept, tags_by


def _dedup_named(named_seqs: list[tuple[str, str]],
                  divergence_threshold: float = 0.05
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
                         key_fn) -> list[tuple[str, str]]:
    """Same as _dedup_named, but the walk order is determined by `key_fn`
    (smaller value = higher priority) instead of length-desc. This lets us
    rank by completeness FIRST so the more-complete representative of each
    equivalence class survives, not just the longest one. The divergence
    check itself (5% edit-distance, RC-aware) is unchanged."""
    nonempty = [(n, s) for n, s in named_seqs if s]
    nonempty.sort(key=lambda ns: key_fn(ns[0], ns[1]))
    kept: list[tuple[str, str]] = []
    for n, s in nonempty:
        if any(not is_divergent(s, ks, threshold=divergence_threshold)
               for _, ks in kept):
            continue
        kept.append((n, s))
    return kept


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
                  divergence_threshold: float = 0.05,
                  var_proteins_ref: str | None = None,
                  locus_padding: int = 1500,
                  expected_var_tags: set[str] | None = None,
                  genome_cov: float | None = None,
                  lo_mult: float = 0.25,
                  hi_mult: float = 2.0,
                  queries_dir: str | None = None,
                  return_pools: bool = False,
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

        # Locus trim — drops no-HD candidates
        local_found: dict[str, set[str]] = {}
        if var_proteins_ref:
            local_raw, local_found = _tblastn_trim_each(
                local_raw, var_proteins_ref, locus_padding)

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
            fv = local_found.get(name, set())
            cv = bool(expected_var_tags) and (expected_var_tags <= fv)
            hL, hR = flank_presence.get(name, (False, False))
            cl = cv and hL and hR
            bp = len(seq)
            cov = _cov(name)
            diploid_dist = abs(cov / genome_cov - 0.5) if genome_cov else 1.0
            # sort ascending: smaller is better → negate for DESC tiers
            return (not cl, not cv, -bp, diploid_dist)

        # Dedup within this pool only (no cross-network collapse), now with
        # completeness-first ranking. is_divergent threshold unchanged.
        local_dedup = _dedup_named_ranked(local_raw, divergence_threshold, _rank_key)
        nd = len(local_dedup)

        # seg_processor.directional_split now emits clean position-indexed
        # sub-node IDs ("{parent}#1", "{parent}#2", "{parent}#N"). The
        # provenance dict carries the coords for sequence materialization;
        # the IDs themselves stay coord-free for clean display in bubble.*.
        def _segs(rn):
            out, seen = [], set()
            for sub in local_path.get(rn, ()):
                if sub not in prov: continue
                if sub in seen: continue
                seen.add(sub); out.append(sub)
            return out

        # Per-pool label scheme: n=1→allele1, n=2→allele1/2, n>=3→chimeraN
        if nd <= 1:
            if local_dedup:
                rn, rs = local_dedup[0]
                alleles_p     = [(f"{prefix}allele1", rs)]
                cov_p         = [_cov(rn)]
                bounds_p      = [_bnds(rn)]
                segs_p        = [_segs(rn)]
            else:
                alleles_p = []; cov_p = []; bounds_p = []; segs_p = []
        elif nd == 2:
            alleles_p = [(f"{prefix}allele1", local_dedup[0][1]),
                          (f"{prefix}allele2", local_dedup[1][1])]
            cov_p     = [_cov(local_dedup[0][0]), _cov(local_dedup[1][0])]
            bounds_p  = [_bnds(local_dedup[0][0]), _bnds(local_dedup[1][0])]
            segs_p    = [_segs(local_dedup[0][0]), _segs(local_dedup[1][0])]
        else:
            alleles_p = [(f"{prefix}chimera{k+1}", s)
                         for k, (_, s) in enumerate(local_dedup)]
            cov_p    = [_cov(rn) for rn, _ in local_dedup]
            bounds_p = [_bnds(rn) for rn, _ in local_dedup]
            segs_p   = [_segs(rn) for rn, _ in local_dedup]

        surviving = set()
        for rn, _ in local_dedup:
            surviving |= local_found.get(rn, set())

        # Pool-level metrics (used by find_alleles for cross-iteration ranking)
        complete_var_p = (expected_var_tags is not None
                           and bool(expected_var_tags)
                           and expected_var_tags <= surviving)
        # Flank presence — taken from the per-candidate blastn we already ran
        # above. A pool is flank-complete iff at least one surviving allele has
        # BOTH flanks (the canonical L→R locus walk). For pools with multiple
        # alleles, "both flanks at the pool level" means each emitted allele
        # individually has both flanks — needed for closed_bubble acceptance.
        all_have_flanks = bool(alleles_p) and all(
            flank_presence.get(rn, (False, False))[0]
            and flank_presence.get(rn, (False, False))[1]
            for rn, _ in local_dedup
        )
        complete_locus_p = complete_var_p and all_have_flanks
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
