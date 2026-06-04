"""Bubble classifier — R1 (separate), R2 (anchors), R3 (arms), R4 (verdict).

Operates on the POST-P1 graph. Assumes every node has at most one role
(pure flank, pure var, or unlabeled connector). Composites are removed
upstream by the seg_processor.

Inputs
======
nodes:          set[node_id]
edges:          set[frozenset({a, b})]
label_per_node: dict[node_id, str]   (flankL / flankR / gene-name / empty)
var_per_node:   dict[node_id, set[str]]   (which var genes the node carries)

Output
======
{class: ..., n_arms: ..., n_closed: ..., n_dangling: ..., explain: ...}
class is one of: closed_bubble, open_bubble, single, complexed, separate, no_var.
"""
from __future__ import annotations
from .bubble_bfs import bubble_bfs, build_adj, connected_components


def _enum_paths(adj, start, ends, allowed, max_paths=1000, max_path_length=50,
                 limit_counter=None,
                 node_bp=None, var_nodes=None, max_bp_since_var=5000):
    """Simple paths from `start` to any node in `ends`, using only `allowed`.

    SHORTEST-FIRST enumeration (BFS by path length). Three hard limits:
      - max_paths          : abort once this many simple paths are emitted
      - max_path_length    : drop any path whose node count would exceed this
                              (the partial path is also dropped from extension)
      - max_bp_since_var   : BP-AWARE cap. When > 0, track per-path the bp
                              accumulated since the last VAR node in the path
                              (or since `start` if no var node yet). Drop
                              extensions where the post-extension bp count
                              would exceed this. Reaching a VAR node OR an
                              END node always extends (we want to reach
                              productive content even at high bp). Set 0 to
                              disable. Saves enumeration time on chains of
                              long unlabeled connectors that would never
                              reach productive content.

    `limit_counter` (optional dict) tallies each cap hit so callers can log:
        {"max_paths_hit": int, "max_path_length_hit": int, "max_bp_hit": int}
    """
    bp_enabled = (max_bp_since_var > 0 and node_bp is not None
                  and var_nodes is not None)
    def _bp(n): return node_bp.get(n, 0) if node_bp else 0
    paths = []
    # State: (path, vis, bp_since_var). bp_since_var starts at start's own
    # bp if start is not var, else 0.
    start_bp = 0 if (bp_enabled and start in var_nodes) else _bp(start)
    frontier = [([start], frozenset({start}), start_bp)]
    while frontier and len(paths) < max_paths:
        next_frontier = []
        for path, vis, bp_acc in frontier:
            if len(paths) >= max_paths: break
            curr = path[-1]
            for nxt in adj.get(curr, ()):
                if nxt in vis: continue
                if nxt not in allowed and nxt not in ends: continue
                if len(path) + 1 > max_path_length:
                    if limit_counter is not None:
                        limit_counter["max_path_length_hit"] = \
                            limit_counter.get("max_path_length_hit", 0) + 1
                    continue
                # BP-aware cap, PRE-extension check:
                #   * Var nodes always extend and reset the tail to 0.
                #   * Non-var nodes extend ONLY when the CURRENT accumulated
                #     bp (bp_acc, NOT bp_acc + size of nxt) is within the
                #     cap. This means a path that hasn't blown the budget
                #     yet can still pick up one MORE node of any size —
                #     e.g., a final 10 kb anchor — and the cap only kicks
                #     in for FURTHER extensions afterward. Stops chain
                #     extension once the path has clearly wandered too far
                #     from var content, but doesn't lose paths that reach
                #     a productive target on the next hop.
                if bp_enabled:
                    if nxt in var_nodes:
                        new_bp = 0
                    else:
                        if bp_acc > max_bp_since_var:
                            if limit_counter is not None:
                                limit_counter["max_bp_hit"] = \
                                    limit_counter.get("max_bp_hit", 0) + 1
                            continue
                        new_bp = bp_acc + _bp(nxt)
                else:
                    new_bp = 0
                new_path = path + [nxt]
                if nxt in ends and nxt != start:
                    paths.append(new_path)
                    if len(paths) >= max_paths:
                        if limit_counter is not None:
                            limit_counter["max_paths_hit"] = \
                                limit_counter.get("max_paths_hit", 0) + 1
                        break
                else:
                    next_frontier.append((new_path, vis | {nxt}, new_bp))
        frontier = next_frontier
    return paths


def _enum_dangling(adj, start, allowed, var_nodes, exclude_var_subset,
                    max_paths=1000, max_path_length=50, limit_counter=None,
                    node_bp=None, max_bp_since_var=5000):
    """Var-bearing simple paths from `start` that dead-end inside `allowed`.
    Drop a path if its var content is wholly inside `exclude_var_subset`.
    Same three caps as `_enum_paths` (incl. BP-since-var)."""
    bp_enabled = max_bp_since_var > 0 and node_bp is not None
    def _bp(n): return node_bp.get(n, 0) if node_bp else 0
    paths, seen = [], set()
    start_bp = 0 if (bp_enabled and start in var_nodes) else _bp(start)
    frontier = [([start], frozenset({start}), start_bp)]
    while frontier and len(paths) < max_paths:
        next_frontier = []
        for path, vis, bp_acc in frontier:
            if len(paths) >= max_paths: break
            curr = path[-1]
            extended = False
            for nxt in adj.get(curr, ()):
                if nxt in vis: continue
                if nxt not in allowed: continue
                if len(path) + 1 > max_path_length:
                    if limit_counter is not None:
                        limit_counter["max_path_length_hit"] = \
                            limit_counter.get("max_path_length_hit", 0) + 1
                    continue
                if bp_enabled:
                    if nxt in var_nodes:
                        new_bp = 0
                    else:
                        if bp_acc > max_bp_since_var:
                            if limit_counter is not None:
                                limit_counter["max_bp_hit"] = \
                                    limit_counter.get("max_bp_hit", 0) + 1
                            continue
                        new_bp = bp_acc + _bp(nxt)
                else:
                    new_bp = 0
                extended = True
                next_frontier.append((path + [nxt], vis | {nxt}, new_bp))
            if not extended and len(path) > 1:
                vs = set(path) & var_nodes
                if vs and not (vs <= exclude_var_subset):
                    k = tuple(path) if path[0] < path[-1] else tuple(reversed(path))
                    if k not in seen:
                        seen.add(k); paths.append(list(path))
                        if len(paths) >= max_paths:
                            if limit_counter is not None:
                                limit_counter["max_paths_hit"] = \
                                    limit_counter.get("max_paths_hit", 0) + 1
                            break
        frontier = next_frontier
    return paths


def _canonical(p):
    return tuple(p) if (len(p) > 1 and p[0] < p[-1]) else tuple(reversed(p))


def classify(nodes: set[str], edges: set[frozenset],
             label_per_node: dict[str, str],
             var_per_node: dict[str, set[str]],
             max_paths: int = 1000,
             max_path_length: int = 50,
             node_bp: dict[str, int] | None = None,
             max_bp_since_var: int = 5000) -> dict:
    """Appendix-B classifier. Returns verdict dict.

    BP-aware path enumeration: when `node_bp` (per-node bp lengths from
    provenance) AND `max_bp_since_var > 0` are supplied, simple-path
    enumeration drops any extension whose accumulated bp since the last
    var node would exceed the threshold. Saves search effort on chains
    of long unlabeled connectors that have already moved far from
    productive content."""
    var_nodes = {n for n in nodes if var_per_node.get(n)}
    flankL = {n for n in nodes if "flankL" in (label_per_node.get(n, "")).split("+")}
    flankR = {n for n in nodes if "flankR" in (label_per_node.get(n, "")).split("+")}
    unlabeled = {n for n in nodes
                 if not label_per_node.get(n) and not var_per_node.get(n)}

    adj = build_adj(nodes, edges)

    info = {"n_nodes": len(nodes), "n_var": len(var_nodes),
            "n_flankL": len(flankL), "n_flankR": len(flankR)}
    if not var_nodes:
        return {"class": "no_var", **info}

    # R1. SEPARATE — var nodes in disjoint full-graph components.
    # When there are multiple disjoint networks, expand and classify each
    # network independently (recursive call on the network's induced
    # subgraph), then return the per-network sub-results. Downstream
    # `_emit_result` runs the trim → dedup → emit chain PER NETWORK and
    # tags each network's alleles with its network index, so an N-network
    # `separate` sample produces N independent allele sets.
    full_comps = connected_components(set(nodes), adj)
    var_full = [c for c in full_comps if c & var_nodes]
    if len(var_full) > 1:
        sub_results = []
        for comp in var_full:
            comp_nodes = set(comp)
            comp_edges = {e for e in edges if all(x in comp_nodes for x in tuple(e))}
            comp_labels  = {n: label_per_node.get(n, "") for n in comp_nodes}
            comp_var_per = {n: var_per_node.get(n, set()) for n in comp_nodes if var_per_node.get(n)}
            sub = classify(comp_nodes, comp_edges, comp_labels, comp_var_per,
                            max_paths=max_paths, max_path_length=max_path_length,
                            node_bp=node_bp,
                            max_bp_since_var=max_bp_since_var)
            sub_results.append(sub)
        return {"class": "separate", **info,
                "n_var_components": len(var_full),
                "var_components": [list(c) for c in var_full],
                "sub_results": sub_results,
                "explain": f"var genes in {len(var_full)} disjoint subgraphs; "
                           f"each network classified independently"}

    # P2. BUBBLE via BFS from var through unlabeled
    bubble = bubble_bfs(adj, var_nodes, unlabeled)
    info["bubble_size"] = len(bubble)

    # R2. ANCHORS — topological boundary of the bubble.
    #
    # Anchor = bubble node that meets the outside world. Two rules:
    #   (a) has ≥ 1 neighbor OUTSIDE the bubble (the common case — typically
    #       a flank-labeled node).
    #   (b) has bubble-degree ≤ 1 (dead-end leaf of the bubble subgraph) —
    #       FALLBACK ONLY, applied when (a) yields fewer than 2 anchors.
    #
    # (b) recovers anchors for fragmented assemblies / chromosome-boundary
    # samples where rule (a) gives no flank-labeled neighbor (lost-label or
    # truncated-graph cases). Otherwise it can promote tiny graph artifacts
    # (e.g. a 90 bp leaf hanging off a Y-fork) to anchor status, adding a
    # spurious dangling arm and downgrading the verdict from closed_bubble
    # to complexed. The fallback gate keeps the fragmented-recovery case
    # while leaving clean closed-bubble samples alone.
    flank_L_adj = {n for n in bubble if any(m in flankL for m in adj.get(n, ()))}
    flank_R_adj = {n for n in bubble if any(m in flankR for m in adj.get(n, ()))}
    anchors: set[str] = set()
    for n in bubble:
        deg_out = sum(1 for m in adj.get(n, ()) if m not in bubble)
        if deg_out > 0:
            anchors.add(n)
    if len(anchors) < 2:
        # Fallback: admit bubble-leaf nodes as anchors so a fragmented
        # assembly can still produce arm candidates.
        for n in bubble:
            if n in anchors: continue
            deg_in = sum(1 for m in adj.get(n, ()) if m in bubble)
            if deg_in <= 1:
                anchors.add(n)
    # L/R labels for orientation diagnostics; not used for arm membership.
    L_anchors = (anchors & flank_L_adj) or anchors
    R_anchors = (anchors & flank_R_adj) or anchors
    info["n_anchors"]   = len(anchors)
    info["n_L_anchors"] = len(L_anchors & flank_L_adj)
    info["n_R_anchors"] = len(R_anchors & flank_R_adj)

    # R3. ARMS — simple paths between anchors, classified by endpoint flank status.
    #
    # Closed = path with one flankL-adjacent endpoint AND one flankR-adjacent
    #          endpoint (a proper L→R locus walk).
    # Dangling = path with one flank-adjacent endpoint and one bare-leaf
    #            endpoint (one side anchored, other side dangling).
    # Ignored: same-side paths (L→L, R→R) and leaf→leaf paths — not
    #          biologically meaningful as arms.
    def _is_closed_endpoints(p: list[str]) -> bool:
        eps = {p[0], p[-1]}
        return bool(eps & flank_L_adj) and bool(eps & flank_R_adj)

    def _is_dangling_endpoints(p: list[str]) -> bool:
        eps = {p[0], p[-1]}
        any_flank = eps & (flank_L_adj | flank_R_adj)
        any_leaf  = eps - (flank_L_adj | flank_R_adj)
        return bool(any_flank) and bool(any_leaf)

    limits = {"max_paths_hit": 0, "max_path_length_hit": 0}
    closed: dict[tuple, list[str]] = {}
    dangling: dict[tuple, list[str]] = {}
    for s in anchors:
        for p in _enum_paths(adj, s, anchors - {s}, bubble,
                              max_paths=max_paths, max_path_length=max_path_length,
                              limit_counter=limits,
                              node_bp=node_bp, var_nodes=var_nodes,
                              max_bp_since_var=max_bp_since_var):
            if not any(n in var_nodes for n in p): continue
            if _is_closed_endpoints(p):
                closed.setdefault(_canonical(p), p)
            elif _is_dangling_endpoints(p):
                dangling.setdefault(_canonical(p), p)
            # else: same-side or leaf-only — drop

    # Single-node bridge: ONE var node that has BOTH flankL and flankR
    # adjacency (a composite segment spanning the whole locus).
    for n in anchors & flank_L_adj & flank_R_adj & var_nodes:
        closed.setdefault(_canonical([n]), [n])

    closed_var_union: set[str] = set()
    for p in closed.values():
        closed_var_union |= set(p) & var_nodes

    # Also collect _enum_dangling paths from flank anchors (those that
    # dead-end inside the bubble — paths that never reach another anchor).
    flank_anchored = anchors & (flank_L_adj | flank_R_adj)
    for s in flank_anchored:
        for p in _enum_dangling(adj, s, bubble, var_nodes, closed_var_union,
                                  max_paths=max_paths, max_path_length=max_path_length,
                                  limit_counter=limits,
                                  node_bp=node_bp,
                                  max_bp_since_var=max_bp_since_var):
            if p[-1] in flank_anchored and p[-1] != s: continue
            dangling.setdefault(_canonical(p), p)

    info["bfs_limits"] = limits
    n_closed = len(closed); n_dangling = len(dangling)
    n_arms = n_closed + n_dangling
    info["n_closed"] = n_closed; info["n_dangling"] = n_dangling
    info["n_arms"] = n_arms
    info["closed_arms"] = [list(p) for p in closed.values()]
    info["dangling_arms"] = [list(p) for p in dangling.values()]

    # R4. VERDICT
    if n_arms == 0:
        return {"class": "complexed", **info, "explain": "no var-bearing arm"}
    if n_arms == 1:
        return {"class": "single", **info, "explain": "one var-bearing arm"}
    if n_arms == 2 and n_closed == 2:
        return {"class": "closed_bubble", **info, "explain": "2 closed arms"}
    if n_arms == 2:
        return {"class": "open_bubble", **info,
                "explain": f"2 arms; {n_closed} closed, {n_dangling} dangling"}
    return {"class": "complexed", **info, "explain": f"{n_arms} arms"}
