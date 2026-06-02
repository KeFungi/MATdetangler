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
    """Pairwise nucleotide identity via edlib infix (HW) edit distance.
    Query is the shorter seq, target the longer — end gaps in the target are
    free, so terminal length differences don't penalize identity. Returns
    1.0 for identical (or one fully contained in the other), lower for
    divergent. Handles indels and end-fragmentation properly."""
    import edlib
    if not a or not b: return 0.0
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    r = edlib.align(q, t, task="distance", mode="HW")
    return 1.0 - r["editDistance"] / max(len(q), 1)


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
    return res


def find_alleles(
        seg_label_hits_tsv: str,
        gfa_path: str,
        genome_cov: float | None = None,
        init_nhop: int = 3,
        max_nhop: int = 10,
        divergence_threshold: float = 0.05,
        lo_mult: float = 0.25,
        hi_mult: float = 2.0,
        k: int | str | None = None,
        var_proteins_ref: str | None = None,
        expected_var_tags: set[str] | None = None,
        locus_padding: int = 1500,
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

    # Seeds for the two phases
    var_seeds = {s for s, hits in seg_labels.items()
                 if any(h.kind == "var" for h in hits)}
    flank_seeds = {s for s, hits in seg_labels.items()
                   if any(h.kind == "flank" for h in hits)}

    def _accept(res: dict) -> bool:
        if res["class"] != "closed_bubble": return False
        arms = res.get("closed_arms", [])
        if len(arms) < 2: return False
        seqs = [_arm_sequence_for(a, res["_provenance"], gfa_seqs) for a in arms[:2]]
        if is_divergent(seqs[0], seqs[1], divergence_threshold):
            res["_emit_arms"] = arms[:2]
            res["_emit_seqs"] = seqs
            res["divergent"] = True
            return True
        return False

    def _log(phase: str, nhop: int, use_cov: bool, res: dict) -> None:
        n_arms = res.get("n_arms", 0)
        n_var  = res.get("n_var", 0)
        nhood  = len(res.get("_nhood", ()))
        print(f"  [{phase} nhop={nhop} cov={'on' if use_cov else 'off':>3}] "
              f"|nhood|={nhood:<6} var={n_var:<3} cls={res['class']:<14} arms={n_arms}",
              flush=True)

    last_res = None
    fallback_cov_on = None                     # cov-on @ max_nhop, for the no-acceptance path

    def _phase_loop(phase: str, seeds: set[str]) -> dict | None:
        nonlocal last_res, fallback_cov_on
        for nhop in range(init_nhop, max_nhop + 1):
            for use_cov in (True, False):      # cov-on first → if no acceptance, we replay cov-on at max below
                res = _try_one_pass(
                    seeds, all_edges, endpoints, adj_und, nhop,
                    seg_labels, seg_length, depths, gfa_seqs,
                    genome_cov, divergence_threshold,
                    apply_cov_filter=use_cov,
                    lo_mult=lo_mult, hi_mult=hi_mult,
                )
                res["_phase"]       = phase
                res["_n_hops"]      = nhop
                res["_cov_filter"]  = use_cov
                _log(phase, nhop, use_cov, res)
                last_res = res
                if use_cov and nhop == max_nhop:
                    fallback_cov_on = res      # remember the cov-filtered @ max for emission fallback
                if _accept(res):
                    return res
        return None

    def _finalize(res: dict, divergence_threshold: float) -> dict:
        # _emit_result now does the locus-trim + completeness check inline.
        out = _emit_result(
            res, gfa_seqs, depths=depths,
            divergence_threshold=divergence_threshold,
            var_proteins_ref=var_proteins_ref,
            locus_padding=locus_padding,
            expected_var_tags=expected_var_tags,
        )
        out["k"] = k
        out["genome_cov"] = genome_cov
        out["bubble_type"] = out["verdict"]                            # alias for TSV clarity
        # Segments + segment-label list, from the result's emitted candidates.
        prov = res.get("_provenance", {})
        emitted_subnodes: set[str] = set()
        for p in res.get("closed_arms", []):    emitted_subnodes.update(p)
        for p in res.get("dangling_arms", []):  emitted_subnodes.update(p)
        for c in res.get("var_components") or []: emitted_subnodes.update(c)
        emitted_segs: set[str] = {prov[n][0] for n in emitted_subnodes if n in prov}
        labels_by_seg: dict[str, str] = {}
        for s in emitted_segs:
            tags = sorted({h.tag for h in seg_labels.get(s, [])})
            if tags: labels_by_seg[s] = "+".join(tags)
        out["component_list"]   = [n for n, _ in out["alleles"]]
        out["segments"]         = sorted(emitted_segs)
        out["segments_labeled"] = [(s, labels_by_seg.get(s, "")) for s in sorted(emitted_segs)]
        out["basepair"]         = sum(len(s) for _, s in out["alleles"])
        return out

    accepted = _phase_loop("main", var_seeds)
    if accepted is not None:
        return _finalize(accepted, divergence_threshold)

    accepted = _phase_loop("flank_fallback", var_seeds | flank_seeds)
    if accepted is not None:
        return _finalize(accepted, divergence_threshold)

    # Loop exhausted — fallback chain:
    #   1) cov-on @ max_nhop (preferred — selective)
    #   2) if cov-on yielded 0 candidates, fall back to cov-off @ max_nhop (last_res)
    def _has_candidates(r: dict) -> bool:
        return (bool(r.get("closed_arms")) or bool(r.get("dangling_arms"))
                or bool(r.get("var_components")))

    emit_from = fallback_cov_on or last_res
    if emit_from is fallback_cov_on and not _has_candidates(fallback_cov_on):
        emit_from = last_res
        print(f"  [exhausted] cov-on @ max_nhop emitted 0 candidates; falling back to cov-off",
              flush=True)
    print(f"  [exhausted] emitting from phase={emit_from.get('_phase')} "
          f"nhop={emit_from.get('_n_hops')} cov={emit_from.get('_cov_filter')}",
          flush=True)
    return _finalize(emit_from, divergence_threshold)


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

    # 1. Build candidate sequences. Apply node-level path-trim first (graph
    #    structural), then sequence-level tblastn-trim below (content-based).
    raw: list[tuple[str, str]] = []
    path_by_name: dict[str, list[str]] = {}
    for i, p in enumerate(res.get("closed_arms", [])):
        trimmed_p = _trim_path_to_var(list(p))
        nm = f"cl{i}"
        raw.append((nm, _seq_full(trimmed_p)))
        path_by_name[nm] = trimmed_p
    for i, p in enumerate(res.get("dangling_arms", [])):
        trimmed_p = _trim_path_to_var(list(p))
        nm = f"da{i}"
        raw.append((nm, _seq_full(trimmed_p)))
        path_by_name[nm] = trimmed_p
    for i, c in enumerate(res.get("var_components") or []):
        trimmed_c = _trim_component_to_var(list(c))
        nm = f"vc{i}"
        raw.append((nm, _arm_sequence_for(trimmed_c, prov, gfa_seqs)))
        path_by_name[nm] = trimmed_c

    n_raw = len(raw)

    # 2. Locus trim via tblastn(var proteins → candidate). Drops no-hit candidates.
    found_tags: dict[str, set[str]] = {}
    if var_proteins_ref:
        raw, found_tags = _tblastn_trim_each(raw, var_proteins_ref, locus_padding)

    # 3. Dedup the trimmed sequences (k-mer Jaccard at threshold).
    deduped_named = _dedup_named(raw, divergence_threshold)
    n = len(deduped_named)

    # 4. Emit by post-dedup count. Compute length-weighted mean depth per allele.
    depths = depths or {}
    cov_for = lambda raw_name: _path_mean_cov(path_by_name.get(raw_name, []), prov, depths)
    if n <= 1:
        if deduped_named:
            rn, rs = deduped_named[0]
            alleles = [("allele1", rs)]
            allele_cov = [cov_for(rn)]
        else:
            alleles = []
            allele_cov = []
    elif n == 2:
        alleles = [("allele1", deduped_named[0][1]), ("allele2", deduped_named[1][1])]
        allele_cov = [cov_for(deduped_named[0][0]), cov_for(deduped_named[1][0])]
    else:
        alleles = [(f"chimera{i+1}", s) for i, (_, s) in enumerate(deduped_named)]
        allele_cov = [cov_for(rn) for rn, _ in deduped_named]

    # Aggregate var tags hit by the surviving (post-dedup) candidates.
    surviving_tags: set[str] = set()
    for name, _ in deduped_named:
        surviving_tags |= found_tags.get(name, set())

    if expected_var_tags is not None:
        complete_var = expected_var_tags <= surviving_tags
        locus_coverage = len(expected_var_tags & surviving_tags) / max(1, len(expected_var_tags))
    else:
        complete_var = None
        locus_coverage = None

    return {
        "verdict": cls,                                                # always the classifier origin
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
