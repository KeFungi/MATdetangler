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
                "explain": f"var genes in {len(var_full)} disjoint subgraphs"}

    # P2. BUBBLE via BFS from var through unlabeled
    bubble = bubble_bfs(adj, var_nodes, unlabeled)
    info["bubble_size"] = len(bubble)

    # R2. ANCHORS — bubble nodes adjacent to flank-region members
    L_anchors = {n for n in bubble if any(m in flankL for m in adj.get(n, ()))}
    R_anchors = {n for n in bubble if any(m in flankR for m in adj.get(n, ()))}
    info["n_L_anchors"] = len(L_anchors)
    info["n_R_anchors"] = len(R_anchors)

    # R3. ARMS — closed (L→R) and dangling (one-sided with novel var content)
    closed: dict[tuple, list[str]] = {}
    for s in L_anchors:
        for p in _enum_paths(adj, s, R_anchors, bubble):
            if any(n in var_nodes for n in p):
                closed.setdefault(_canonical(p), p)
    for n in L_anchors & R_anchors & var_nodes:       # single-node bridge
        closed.setdefault(_canonical([n]), [n])

    closed_var_union: set[str] = set()
    for p in closed.values():
        closed_var_union |= set(p) & var_nodes

    dangling: dict[tuple, list[str]] = {}
    for s in L_anchors | R_anchors:
        if s in L_anchors and s in R_anchors: continue
        opp = R_anchors if s in L_anchors else L_anchors
        for p in _enum_dangling(adj, s, bubble, var_nodes, closed_var_union):
            if p[-1] in opp: continue
            dangling.setdefault(_canonical(p), p)

    n_closed = len(closed); n_dangling = len(dangling)
    n_arms = n_closed + n_dangling
    info["n_closed"] = n_closed; info["n_dangling"] = n_dangling
    info["n_arms"] = n_arms

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
