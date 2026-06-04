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


def _enum_paths(adj_dir, start, ends, allowed, max_paths=1000, max_path_length=50,
                 limit_counter=None,
                 node_bp=None, var_nodes=None, max_bp_since_var=5000,
                 depths=None, genome_cov=None, provenance=None):
    """Directed simple-ish paths from `start` to any node in `ends`.
    `adj_dir[n] = {(side, m, m_side), ...}`
    """
    bp_enabled = (max_bp_since_var > 0 and node_bp is not None
                  and var_nodes is not None)
    def _bp(n): return node_bp.get(n, 0) if node_bp else 0

    def _get_max_visits(n):
        if not depths or not genome_cov or not provenance: return 1
        parent = provenance.get(n, (n,))[0]
        d = depths.get(parent, genome_cov)
        return max(1, round(d / genome_cov))

    paths = []
    # State: (path, vis_counts, bp_acc, last_exit_side)
    # We try starting from both sides of the start node.
    frontier = []
    for side in ("L", "R"):
        frontier.append(([start], {start: 1}, start_bp, side))

    while frontier and len(paths) < max_paths:
        next_frontier = []
        for path, vis_counts, bp_acc, last_side in frontier:
            if len(paths) >= max_paths: break
            curr = path[-1]
            # In a directed graph, if we exited through 'R', we must enter
            # the neighbor's connected side. adj_dir[curr] already stores 
            # edges as (exit_side, neighbor, entry_side).
            for exit_side, nxt, entry_side in adj_dir.get(curr, ()):
                if exit_side != last_side: continue # Must exit from the side we are at
                
                count = vis_counts.get(nxt, 0)
                if count >= _get_max_visits(nxt): continue
                if nxt not in allowed and nxt not in ends: continue
                
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
                
                new_path = path + [nxt]
                new_vis = vis_counts.copy()
                new_vis[nxt] = count + 1
                
                # Exit side for the NEXT step is the opposite of the entry side
                next_exit_side = "R" if entry_side == "L" else "L"
                
                if nxt in ends and nxt != start:
                    paths.append(new_path)
                    if len(paths) >= max_paths:
                        if limit_counter is not None:
                            limit_counter["max_paths_hit"] = \
                                limit_counter.get("max_paths_hit", 0) + 1
                        break
                else:
                    next_frontier.append((new_path, new_vis, new_bp, next_exit_side))
        frontier = next_frontier
    return paths


def _enum_dangling(adj_dir, start, allowed, var_nodes, exclude_var_subset,
                    max_paths=1000, max_path_length=50, limit_counter=None,
                    node_bp=None, max_bp_since_var=5000,
                    depths=None, genome_cov=None, provenance=None):
    """Directed var-bearing paths that dead-end."""
    bp_enabled = max_bp_since_var > 0 and node_bp is not None
    def _bp(n): return node_bp.get(n, 0) if node_bp else 0

    def _get_max_visits(n):
        if not depths or not genome_cov or not provenance: return 1
        parent = provenance.get(n, (n,))[0]
        d = depths.get(parent, genome_cov)
        return max(1, round(d / genome_cov))

    paths, seen = [], set()
    start_bp = 0 if (bp_enabled and start in var_nodes) else _bp(start)
    frontier = []
    for side in ("L", "R"):
        frontier.append(([start], {start: 1}, start_bp, side))

    while frontier and len(paths) < max_paths:
        next_frontier = []
        for path, vis_counts, bp_acc, last_side in frontier:
            if len(paths) >= max_paths: break
            curr = path[-1]
            extended = False
            for exit_side, nxt, entry_side in adj_dir.get(curr, ()):
                if exit_side != last_side: continue
                
                count = vis_counts.get(nxt, 0)
                if count >= _get_max_visits(nxt): continue
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
                new_vis = vis_counts.copy()
                new_vis[nxt] = count + 1
                next_exit_side = "R" if entry_side == "L" else "L"
                next_frontier.append((path + [nxt], new_vis, new_bp, next_exit_side))
            
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


def classify(nodes: set[str], adj_dir: dict[str, set[tuple[str, str, str]]],
             label_per_node: dict[str, str],
             var_per_node: dict[str, set[str]],
             max_paths: int = 1000,
             max_path_length: int = 50,
             node_bp: dict[str, int] | None = None,
             max_bp_since_var: int = 5000,
             depths: dict[str, float] | None = None,
             genome_cov: float | None = None,
             provenance: dict[str, tuple] | None = None) -> dict:
    """Appendix-B classifier (Directed)."""
    var_nodes = {n for n in nodes if var_per_node.get(n)}
    flankL = {n for n in nodes if "flankL" in (label_per_node.get(n, "")).split("+")}
    flankR = {n for n in nodes if "flankR" in (label_per_node.get(n, "")).split("+")}
    unlabeled = {n for n in nodes
                 if not label_per_node.get(n) and not var_per_node.get(n)}

    # Build undirected adjacency only for component and bubble-BFS (label-blind)
    adj_und = defaultdict(set)
    for u, edges in adj_dir.items():
        for _, v, _ in edges:
            adj_und[u].add(v); adj_und[v].add(u)

    info = {"n_nodes": len(nodes), "n_var": len(var_nodes),
            "n_flankL": len(flankL), "n_flankR": len(flankR)}
    if not var_nodes:
        return {"class": "no_var", **info}

    full_comps = connected_components(set(nodes), dict(adj_und))
    var_full = [c for c in full_comps if c & var_nodes]
    if len(var_full) > 1:
        sub_results = []
        for comp in var_full:
            comp_nodes = set(comp)
            comp_adj_dir = {n: {e for e in adj_dir.get(n, ()) if e[1] in comp_nodes} 
                            for n in comp_nodes}
            comp_labels  = {n: label_per_node.get(n, "") for n in comp_nodes}
            comp_var_per = {n: var_per_node.get(n, set()) for n in comp_nodes if var_per_node.get(n)}
            sub = classify(comp_nodes, comp_adj_dir, comp_labels, comp_var_per,
                            max_paths=max_paths, max_path_length=max_path_length,
                            node_bp=node_bp,
                            max_bp_since_var=max_bp_since_var,
                            depths=depths, genome_cov=genome_cov,
                            provenance=provenance)
            sub_results.append(sub)
        return {"class": "separate", **info,
                "n_var_components": len(var_full),
                "var_components": [list(c) for c in var_full],
                "sub_results": sub_results,
                "explain": f"var genes in {len(var_full)} disjoint subgraphs; "
                           f"each network classified independently"}

    bubble = bubble_bfs(dict(adj_und), var_nodes, unlabeled)
    info["bubble_size"] = len(bubble)

    flank_L_adj = {n for n in bubble if any(m in flankL for m in adj_und.get(n, ()))}
    flank_R_adj = {n for n in bubble if any(m in flankR for m in adj_und.get(n, ()))}
    anchors: set[str] = set()
    for n in bubble:
        deg_out = sum(1 for m in adj_und.get(n, ()) if m not in bubble)
        if deg_out > 0:
            anchors.add(n)
    if len(anchors) < 2:
        for n in bubble:
            if n in anchors: continue
            deg_in = sum(1 for m in adj_und.get(n, ()) if m in bubble)
            if deg_in <= 1:
                anchors.add(n)
    
    limits = {"max_paths_hit": 0, "max_path_length_hit": 0}
    closed: dict[tuple, list[str]] = {}
    dangling: dict[tuple, list[str]] = {}
    for s in anchors:
        for p in _enum_paths(adj_dir, s, anchors - {s}, bubble,
                              max_paths=max_paths, max_path_length=max_path_length,
                              limit_counter=limits,
                              node_bp=node_bp, var_nodes=var_nodes,
                              max_bp_since_var=max_bp_since_var,
                              depths=depths, genome_cov=genome_cov,
                              provenance=provenance):
            if not any(n in var_nodes for n in p): continue
            eps = {p[0], p[-1]}
            if bool(eps & flank_L_adj) and bool(eps & flank_R_adj):
                closed.setdefault(_canonical(p), p)
            else:
                any_flank = eps & (flank_L_adj | flank_R_adj)
                any_leaf  = eps - (flank_L_adj | flank_R_adj)
                if bool(any_flank) and bool(any_leaf):
                    dangling.setdefault(_canonical(p), p)

    for n in anchors & flank_L_adj & flank_R_adj & var_nodes:
        closed.setdefault(_canonical([n]), [n])

    closed_var_union: set[str] = set()
    for p in closed.values():
        closed_var_union |= set(p) & var_nodes

    flank_anchored = anchors & (flank_L_adj | flank_R_adj)
    for s in flank_anchored:
        for p in _enum_dangling(adj_dir, s, bubble, var_nodes, closed_var_union,
                                  max_paths=max_paths, max_path_length=max_path_length,
                                  limit_counter=limits,
                                  node_bp=node_bp,
                                  max_bp_since_var=max_bp_since_var,
                                  depths=depths, genome_cov=genome_cov,
                                  provenance=provenance):
            if p[-1] in flank_anchored and p[-1] != s: continue
            dangling.setdefault(_canonical(p), p)

    info["bfs_limits"] = limits
    n_closed = len(closed); n_dangling = len(dangling)
    n_arms = n_closed + n_dangling
    info["n_closed"] = n_closed; info["n_dangling"] = n_dangling
    info["n_arms"] = n_arms
    info["closed_arms"] = [list(p) for p in closed.values()]
    info["dangling_arms"] = [list(p) for p in dangling.values()]

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
