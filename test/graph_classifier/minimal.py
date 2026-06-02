"""Minimal-rules topology classifier.

Three rules, applied in order:

  R1. SEPARATE
        If var-bearing nodes split across 2+ connected components of the
        full graph -> 'separate'.

  R2. BUBBLE SUBGRAPH
        bubble = (nodes - PURE_FLANK_NODES).  Pure-flank = labeled flankL/R
        AND NO var gene. Composite (flank+var) nodes stay in the bubble.

  R3. COUNT ARMS = distinct var-bearing simple paths starting at an L-anchor
      OR R-anchor. Deduplicate by minimal var-node set (drop frankenstein
      supersets). An arm is "closed" if its endpoints touch BOTH flank
      anchors; "dangling" if only one. L-anchor = bubble node carrying
      flankL tag (composite) OR adjacent to pure flankL. R-anchor = same
      for flankR.

  Classification by arm count:
        0 arms but var nodes exist -> 'complexed' (no var-bearing arm)
        1 arm                       -> 'single' (haploid or collapsed)
        2 arms                      -> 'closed_bubble' (both arms closed)
                                    OR 'open_bubble' (at least one dangling)
        3+ arms                     -> 'complexed'

That's it. No flank-region BFS, no validity horizons, no arm-membership
consolidation, no edge-disjoint paths. Three rules.
"""
from __future__ import annotations
from collections import defaultdict


def _has_tag(label: str, tag: str) -> bool:
    return tag in label.split("+") if label else False


def _components(nodes: set[str], adj: dict[str, set[str]]) -> list[set[str]]:
    seen: set[str] = set(); out: list[set[str]] = []
    for s in nodes:
        if s in seen: continue
        comp, stack = set(), [s]
        while stack:
            x = stack.pop()
            if x in seen: continue
            seen.add(x); comp.add(x)
            stack.extend(adj.get(x, set()) - seen)
        out.append(comp)
    return out


def _enum_paths(adj, start, ends, allowed, max_paths=200):
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


def classify(nodes, edges, labels, var_per) -> dict:
    var_nodes = {n for n in nodes if var_per.get(n)}
    flankL = {n for n in nodes if _has_tag(labels.get(n, ""), "flankL")}
    flankR = {n for n in nodes if _has_tag(labels.get(n, ""), "flankR")}
    pure_flankL = flankL - var_nodes
    pure_flankR = flankR - var_nodes
    pure_flank  = pure_flankL | pure_flankR

    adj = defaultdict(set)
    for e in edges:
        a, b = tuple(e); adj[a].add(b); adj[b].add(a)
    adj = dict(adj)

    info = {"n_nodes": len(nodes), "n_var": len(var_nodes),
            "n_flankL": len(flankL), "n_flankR": len(flankR),
            "n_pure_flank": len(pure_flank)}
    if not var_nodes:
        return {"class": "no_var", **info}

    # R1. SEPARATE
    full_comps = _components(set(nodes), adj)
    var_full = [c for c in full_comps if c & var_nodes]
    if len(var_full) > 1:
        return {"class": "separate", **info,
                "explain": f"var genes in {len(var_full)} disjoint subgraphs"}

    # R2. BUBBLE SUBGRAPH
    bubble = set(nodes) - pure_flank

    # R3. ARMS = simple paths from L-anchor to R-anchor through bubble.
    # L-anchor: bubble node carrying flankL tag (composite) OR adjacent to pure flankL.
    # R-anchor: symmetric.
    L_anchors = set(bubble & flankL)
    for f in pure_flankL:
        L_anchors |= adj.get(f, set()) & bubble
    R_anchors = set(bubble & flankR)
    for f in pure_flankR:
        R_anchors |= adj.get(f, set()) & bubble

    info["n_L_anchors"] = len(L_anchors)
    info["n_R_anchors"] = len(R_anchors)

    # Step A: closed paths (L→R or R→L). Dedupe by minimal var-set.
    closed_paths = []
    for s in L_anchors:
        closed_paths.extend(_enum_paths(adj, s, R_anchors, bubble))
    for n in L_anchors & R_anchors & var_nodes:
        closed_paths.append([n])
    closed_by_vs = {}
    for p in closed_paths:
        vs = frozenset(set(p) & var_nodes)
        if not vs: continue
        if vs not in closed_by_vs or len(p) < len(closed_by_vs[vs]):
            closed_by_vs[vs] = p
    sorted_sets = sorted(closed_by_vs.keys(), key=len)
    minimal_closed = []
    for vs in sorted_sets:
        if any(other < vs for other in minimal_closed): continue
        minimal_closed.append(vs)
    closed_arms = [closed_by_vs[vs] for vs in minimal_closed]
    closed_var_union = set().union(*minimal_closed) if minimal_closed else set()

    # Step B: dangling arms (var-bearing paths starting at an anchor, ending
    # at a non-anchor dead-end). Keep ONLY if the dangling var content is not
    # already covered by closed arms (filters "frankenstein" detours).
    def _enum_dangling(adj, start, allowed):
        out, seen_keys = [], set()
        def dfs(curr, path, vis):
            extended = False
            for nxt in adj.get(curr, ()):
                if nxt in vis: continue
                if nxt not in allowed: continue
                vis.add(nxt); path.append(nxt)
                dfs(nxt, path, vis)
                path.pop(); vis.discard(nxt)
                extended = True
            if not extended and any(n in var_nodes for n in path) and len(path) > 1:
                k = tuple(path) if path[0] < path[-1] else tuple(reversed(path))
                if k not in seen_keys:
                    seen_keys.add(k); out.append(list(path))
        dfs(start, [start], {start})
        return out
    dangling_paths = []
    for s in L_anchors | R_anchors:
        if s in L_anchors and s in R_anchors: continue
        for p in _enum_dangling(adj, s, bubble):
            other_side = R_anchors if s in L_anchors else L_anchors
            if p[-1] in other_side: continue                              # actually a closed arm
            dangling_paths.append(p)
    dang_by_vs = {}
    for p in dangling_paths:
        vs = frozenset(set(p) & var_nodes)
        if not vs: continue
        # Only count if there's at least one var node NOT in closed arms
        if vs <= closed_var_union: continue
        if vs not in dang_by_vs or len(p) < len(dang_by_vs[vs]):
            dang_by_vs[vs] = p
    minimal_dangling = []
    for vs in sorted(dang_by_vs.keys(), key=len):
        if any(other < vs for other in minimal_dangling): continue
        minimal_dangling.append(vs)
    dangling_arms = [dang_by_vs[vs] for vs in minimal_dangling]

    n_closed = len(closed_arms)
    n_dangling = len(dangling_arms)
    n_arms = n_closed + n_dangling
    info["n_arms"] = n_arms
    info["n_closed"] = n_closed
    info["n_dangling"] = n_dangling

    if n_arms == 0:
        return {"class": "complexed", **info, "explain": "no var-bearing arm"}
    if n_arms == 1:
        return {"class": "single", **info, "explain": "one var-bearing arm"}
    if n_arms == 2:
        if n_closed == 2:
            return {"class": "closed_bubble", **info, "explain": "2 arms, both closed"}
        return {"class": "open_bubble", **info, "explain": f"2 arms, {n_closed} closed"}
    return {"class": "complexed", **info, "explain": f"{n_arms} arms"}


def classify_network(net) -> dict:
    return classify(set(net.nodes), set(net.edges), dict(net.labels),
                    {k: set(v) for k, v in net.var_per.items()})
