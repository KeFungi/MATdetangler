"""New walking-based topology classifier.

MAX_LINKER_PADDING is the single coupled parameter: the maximum number of
unlabeled intermediate nodes allowed between two labeled elements (var or flank)
in well-formed input. Used both by the test-data generator (cap on
`add_linkers` chain length) and by the classifier (flank validity = reachable
to some var node through ≤ MAX_LINKER_PADDING unlabeled intermediates).

Algorithm (operates on (nodes, edges, labels, var_per)):

  1. var-only series = connected components of (nodes \\ flank_nodes), restricted to
     nodes that carry a var gene OR are unlabeled connectors. A series with no var
     gene at all is discarded (pure noise / linker chain).

  2. Check that all var genes belong to ONE connected component of the FULL graph
     (full-graph component check); if not -> 'separate'.

  3. For each series, find:
       - Ls : set of flankL nodes adjacent to any node in the series
       - Rs : set of flankR nodes adjacent to any node in the series
       - clean: bool; the series is a "clean chain" iff
                (a) max degree inside the series is <= 2 (no internal branching)
                (b) any flank attachment goes to a chain endpoint
                    (a deg<=1 node within the series)

  4. Classification:
       - n_series != 2 OR any series not clean    -> 'complexed'
       - both series share the same single flankL  AND
         both series share the same single flankR  -> 'closed_bubble'
       - both series share the same single flank   on exactly one side AND
         the other side is empty or attached to only one of the series at one node
                                                   -> 'open_bubble'
       - else                                      -> 'complexed'

Returns a result dict with the verdict, the series breakdown, and diagnostics
explaining the call.
"""
from __future__ import annotations
from collections import defaultdict


MAX_LINKER_PADDING = 3


def _has_tag(label: str, tag: str) -> bool:
    return tag in label.split("+") if label else False


def _connected_components(nodes: set[str], adj: dict[str, set[str]]) -> list[set[str]]:
    seen: set[str] = set()
    comps: list[set[str]] = []
    for s in nodes:
        if s in seen: continue
        comp: set[str] = set()
        stack = [s]
        while stack:
            x = stack.pop()
            if x in seen: continue
            seen.add(x); comp.add(x)
            for y in adj.get(x, ()):
                if y in nodes and y not in seen:
                    stack.append(y)
        comps.append(comp)
    return comps


def _prune_dead_ends(nodes: set[str], adj: dict[str, set[str]],
                     keep: set[str]) -> set[str]:
    """Iteratively remove leaf nodes (deg <= 1 in the current pruned subgraph)
    that are NOT in `keep`. `keep` should contain all var and flank nodes.
    This strips noise tails — chains of unlabeled nodes hanging off a real node —
    without touching any chain that has var or flank termini."""
    alive = set(nodes)
    while True:
        to_remove = {n for n in alive
                     if n not in keep and len(adj.get(n, set()) & alive) <= 1}
        if not to_remove: break
        alive -= to_remove
    return alive


def classify(nodes: set[str], edges: set[frozenset], labels: dict[str, str],
             var_per: dict[str, set[str]]) -> dict:
    # Build helper sets
    var_nodes = {n for n in nodes if var_per.get(n)}
    flankL_nodes = {n for n in nodes if _has_tag(labels.get(n, ""), "flankL")}
    flankR_nodes = {n for n in nodes if _has_tag(labels.get(n, ""), "flankR")}
    flank_nodes = flankL_nodes | flankR_nodes
    # Pure-flank nodes carry a flank tag and no var gene. Composite nodes
    # (label "flankL+HD1", vars={"HD1"}) are NOT pure-flank — they stay on
    # the var side when we excise flanks to find var-series.
    pure_flank_nodes = flank_nodes - var_nodes

    # Adjacency over the FULL graph
    full_adj: dict[str, set[str]] = defaultdict(set)
    for e in edges:
        a, b = tuple(e)
        full_adj[a].add(b); full_adj[b].add(a)

    # Pre-process: prune noise tails (unlabeled leaves) before classifying.
    # Anchoring set = var nodes + flank nodes (the real graph structure).
    kept = _prune_dead_ends(set(nodes), full_adj, var_nodes | flank_nodes)
    nodes = kept
    adj: dict[str, set[str]] = {n: full_adj[n] & kept for n in kept}
    pure_flank_nodes &= kept
    flankL_nodes &= kept
    flankR_nodes &= kept
    flank_nodes &= kept

    info = {
        "n_nodes": len(nodes), "n_edges": len(edges),
        "n_var": len(var_nodes), "n_flankL": len(flankL_nodes), "n_flankR": len(flankR_nodes),
    }

    if not var_nodes:
        return {"class": "no_var", **info}

    # 1. Separate? — all var genes in ONE full-graph component
    full_comps = _connected_components(set(nodes), dict(adj))
    var_full_comps = [c for c in full_comps if c & var_nodes]
    info["n_var_full_components"] = len(var_full_comps)
    if len(var_full_comps) > 1:
        return {"class": "separate", **info,
                "explain": f"var genes in {len(var_full_comps)} disjoint subgraphs"}

    # 2. Var-only series = components of (nodes - flank_nodes), keeping only those
    #    that contain at least one var node. ALL flank-tagged nodes are treated as
    #    boundary — including composites like "flankL+HD1" — because the gene tag
    #    on the boundary doesn't put that segment INSIDE the var region; it puts
    #    it on the flank/var seam.
    non_flank = set(nodes) - flank_nodes
    non_flank_adj = {n: adj[n] & non_flank for n in non_flank}
    raw_series = _connected_components(non_flank, non_flank_adj)
    var_series = [s for s in raw_series if s & var_nodes]
    info["n_var_series"] = len(var_series)

    # 2b. Flank-region components — generalized: a set of same-side flank nodes
    #    that are mutually reachable through any path that doesn't cross an
    #    INTERIOR var node (var-tagged AND not flank-tagged). This treats
    #    linker subdivisions, fragment-split flanks, and composite boundaries
    #    as transparent to flank-region connectivity.
    interior_var = var_nodes - flank_nodes
    allowed = set(nodes) - interior_var
    allowed_adj = {n: adj[n] & allowed for n in allowed}
    region_comps = _connected_components(allowed, allowed_adj)
    # A flank region is only a VALID ANCHOR if at least one flank node in it
    # can reach some var node through a path with ≤ MAX_LINKER_PADDING unlabeled
    # intermediates. This couples the validity check to the test-data linker
    # cap: any flank that's genuinely part of the locus (with reasonable
    # padding) qualifies; stray `extra_flank` leaves attached via long chains
    # or living in their own disconnected region don't.
    def _reaches_var(start: str, N: int) -> bool:
        """Min-intermediates BFS from start; True iff any var reachable within N."""
        best: dict[str, int] = {start: 0}
        stack: list[tuple[str, int]] = [(start, 0)]
        while stack:
            n, k = stack.pop()
            if k > best[n]:
                continue
            for m in adj.get(n, ()):
                if m in var_nodes:
                    return True
                inc = 0 if m in flank_nodes else 1
                new_k = k + inc
                if new_k > N:
                    continue
                if new_k < best.get(m, 1 << 30):
                    best[m] = new_k
                    stack.append((m, new_k))
        return False

    def _region_has_valid_anchor(region: set[str]) -> bool:
        return any(_reaches_var(f, MAX_LINKER_PADDING) for f in region)

    flankL_comps = [c & flankL_nodes for c in region_comps
                    if c & flankL_nodes and _region_has_valid_anchor(c & flankL_nodes)]
    flankR_comps = [c & flankR_nodes for c in region_comps
                    if c & flankR_nodes and _region_has_valid_anchor(c & flankR_nodes)]
    # Aggregate valid-anchor flank node sets for the internal_flank check below
    valid_flankL_set = set().union(*flankL_comps) if flankL_comps else set()
    valid_flankR_set = set().union(*flankR_comps) if flankR_comps else set()
    valid_flank_set = valid_flankL_set | valid_flankR_set

    def _anchor_idx(series: set[str], comps: list[set[str]],
                    flank_set: set[str]) -> frozenset[int]:
        """Which VALID flank-region components does this series touch?
        A series 'touches' a flank-region if any node in the series is in that
        flank-region (composite case) OR is adjacent to a node in that region."""
        touched: set[int] = set()
        for n in series:
            if n in flank_set:
                for i, c in enumerate(comps):
                    if n in c: touched.add(i)
            for nb in adj[n]:
                if nb in flank_set and nb not in series:
                    for i, c in enumerate(comps):
                        if nb in c: touched.add(i)
        return frozenset(touched)
    # Use the validated flankL/R sets above; stray extras are now invisible.

    # 3. Per-series flank attachments + clean-chain check
    series_info = []
    for idx, s in enumerate(var_series):
        L_anchors = _anchor_idx(s, flankL_comps, valid_flankL_set)
        R_anchors = _anchor_idx(s, flankR_comps, valid_flankR_set)
        # Clean-chain check on the series itself: max degree <= 2, and any
        # PURE flank attachment goes to an endpoint of the series.
        sub_adj = {n: adj[n] & s for n in s}
        deg = {n: len(sub_adj[n]) for n in s}
        endpoints = {n for n, d in deg.items() if d <= 1}
        max_deg = max(deg.values()) if deg else 0
        internal_flank = any((adj[n] & valid_flank_set) and n not in endpoints
                              for n in s)
        clean = (max_deg <= 2) and not internal_flank
        series_info.append({
            "idx": idx, "size": len(s),
            "L_anchors": L_anchors, "R_anchors": R_anchors,
            # Also keep the literal node-set view for verbose debug
            "Ls_nodes": {f for n in s for f in (adj[n] & flankL_nodes)} | (s & flankL_nodes),
            "Rs_nodes": {f for n in s for f in (adj[n] & flankR_nodes)} | (s & flankR_nodes),
            "max_deg": max_deg, "internal_flank": internal_flank, "clean": clean,
            "n_var": len(s & var_nodes),
        })
    info["series"] = series_info

    n_series = len(var_series)
    if n_series == 1:
        # Single var-series — could be a clean haploid-like chain ("single") or
        # a branched/tangled lone series ("complexed").
        if series_info[0]["clean"]:
            return {"class": "single", **info,
                    "explain": "only 1 var-series; clean chain — haploid-like / collapsed dikaryon"}
        return {"class": "complexed", **info,
                "explain": "single var-series but branched or flank attached to interior"}
    if n_series != 2:
        return {"class": "complexed", **info,
                "explain": f"expect 2 var-series, found {n_series}"}
    if not all(si["clean"] for si in series_info):
        return {"class": "complexed", **info,
                "explain": "at least one var-series is not a clean chain"
                           " (branched or pure flank attached to interior)"}

    s1, s2 = series_info
    Ls1, Rs1 = s1["L_anchors"], s1["R_anchors"]
    Ls2, Rs2 = s2["L_anchors"], s2["R_anchors"]

    same_single = lambda A, B: bool(A) and A == B and len(A) == 1

    # Closed: both series share the same single flankL and same single flankR
    if same_single(Ls1, Ls2) and same_single(Rs1, Rs2):
        return {"class": "closed_bubble", **info,
                "explain": f"both series anchor at one flankL ({next(iter(Ls1))}) and one flankR ({next(iter(Rs1))})"}

    # Open: both series share the same single flank on EXACTLY one side, and the
    # other side is either empty, or only one series touches it (and at a single node)
    def is_open_on(LsA, LsB, RsA, RsB):
        """L-anchored open: both series share the L flank; R side may be partial."""
        if not same_single(LsA, LsB): return False
        if not RsA and not RsB: return True                # both R-dangling
        all_Rs = RsA | RsB
        if len(all_Rs) <= 1 and (bool(RsA) ^ bool(RsB)):
            return True                                    # one arm reaches R, other doesn't
        return False

    if is_open_on(Ls1, Ls2, Rs1, Rs2) or is_open_on(Rs1, Rs2, Ls1, Ls2):
        return {"class": "open_bubble", **info,
                "explain": "two series share one flank; other side is dangling or partial"}

    return {"class": "complexed", **info,
            "explain": "flank attachments don't match closed_bubble or open_bubble patterns"}


def classify_network(net) -> dict:
    """Convenience wrapper for a TestNetwork."""
    return classify(set(net.nodes), set(net.edges), dict(net.labels),
                    {k: set(v) for k, v in net.var_per.items()})
