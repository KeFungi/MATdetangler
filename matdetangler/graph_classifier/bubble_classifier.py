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


def _enum_paths(adj, start, ends, allowed, max_paths=200):
    """Simple paths from `start` to any node in `ends`, using only `allowed`."""
    paths = []
    def dfs(curr, path, vis):
        if len(paths) >= max_paths: return
        if curr in ends and curr != start:
            paths.append(list(path)); return
        for nxt in adj.get(curr, ()):
            if nxt in vis: continue
            if nxt not in allowed and nxt not in ends: continue
            vis.add(nxt); path.append(nxt)
            dfs(nxt, path, vis)
            path.pop(); vis.discard(nxt)
    dfs(start, [start], {start})
    return paths


def _enum_dangling(adj, start, allowed, var_nodes, exclude_var_subset, max_paths=200):
    """Var-bearing simple paths from `start` that dead-end in `allowed`.
    Drop a path if its var content is wholly inside `exclude_var_subset`."""
    paths, seen = [], set()
    def dfs(curr, path, vis):
        if len(paths) >= max_paths: return
        extended = False
        for nxt in adj.get(curr, ()):
            if nxt in vis: continue
            if nxt not in allowed: continue
            vis.add(nxt); path.append(nxt)
            dfs(nxt, path, vis)
            path.pop(); vis.discard(nxt)
            extended = True
        if not extended and len(path) > 1:
            vs = set(path) & var_nodes
            if vs and not (vs <= exclude_var_subset):
                k = tuple(path) if path[0] < path[-1] else tuple(reversed(path))
                if k not in seen:
                    seen.add(k); paths.append(list(path))
    dfs(start, [start], {start})
    return paths


def _canonical(p):
    return tuple(p) if (len(p) > 1 and p[0] < p[-1]) else tuple(reversed(p))


def classify(nodes: set[str], edges: set[frozenset],
             label_per_node: dict[str, str],
             var_per_node: dict[str, set[str]]) -> dict:
    """Appendix-B classifier. Returns verdict dict."""
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

    # R1. SEPARATE — var nodes in disjoint full-graph components
    full_comps = connected_components(set(nodes), adj)
    var_full = [c for c in full_comps if c & var_nodes]
    if len(var_full) > 1:
        return {"class": "separate", **info,
                "n_var_components": len(var_full),
                "var_components": [list(c) for c in var_full],
                "explain": f"var genes in {len(var_full)} disjoint subgraphs"}

    # P2. BUBBLE via BFS from var through unlabeled
    bubble = bubble_bfs(adj, var_nodes, unlabeled)
    info["bubble_size"] = len(bubble)

    # R2. ANCHORS — topological boundary of the bubble.
    #
    # Anchor = bubble node that meets the outside world, defined as either:
    #   (a) has ≥1 neighbor OUTSIDE the bubble (covers flank-adjacent nodes
    #       and any node whose neighbor isn't var/unlabeled-reachable), OR
    #   (b) has bubble-degree ≤ 1 (dead-end leaf of the bubble subgraph)
    #
    # This is the universal-leaf rule: it generalizes the previous
    # flank-adjacency rule (those nodes still qualify under (a)) and adds
    # robustness when the BFS-expanded subgraph contains var content but
    # no flank labels reached (chromosome ends, fragmented assemblies,
    # truncated BFS). L/R distinction is decorative — derived from flank
    # labels when present — but arm enumeration uses the universal set.
    flank_L_adj = {n for n in bubble if any(m in flankL for m in adj.get(n, ()))}
    flank_R_adj = {n for n in bubble if any(m in flankR for m in adj.get(n, ()))}
    anchors: set[str] = set()
    for n in bubble:
        deg_in  = sum(1 for m in adj.get(n, ()) if m in bubble)
        deg_out = sum(1 for m in adj.get(n, ()) if m not in bubble)
        if deg_out > 0 or deg_in <= 1:
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

    closed: dict[tuple, list[str]] = {}
    dangling: dict[tuple, list[str]] = {}
    for s in anchors:
        for p in _enum_paths(adj, s, anchors - {s}, bubble):
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
        for p in _enum_dangling(adj, s, bubble, var_nodes, closed_var_union):
            if p[-1] in flank_anchored and p[-1] != s: continue
            dangling.setdefault(_canonical(p), p)

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
