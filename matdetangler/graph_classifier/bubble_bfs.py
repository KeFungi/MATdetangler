"""Bubble BFS — P2.

Bubble = (var nodes) ∪ (unlabeled connectors transitively reachable from
any var node through unlabeled-only neighbors).

The BFS expands from each var seed and walks through unlabeled-tagged
adjacent nodes. When it enters a new unlabeled node, that node is pushed
onto the frontier and its own unlabeled neighbors are explored on
subsequent steps. The walk stops at any flank-tagged node (because flank
nodes are not in the `unlabeled` set passed in).

Result: every node in the returned set is either a var-tagged node or an
unlabeled connector that lies in the var-reachable region. Everything
else is flank-region.
"""
from __future__ import annotations
from collections import defaultdict


def bubble_bfs(adj: dict[str, set[str]],
               var_nodes: set[str],
               unlabeled: set[str]) -> set[str]:
    """Return the bubble — var nodes plus transitively-reachable unlabeled."""
    bubble: set[str] = set(var_nodes)
    stack: list[str] = list(var_nodes)
    while stack:
        n = stack.pop()
        for m in adj.get(n, ()):
            if m in unlabeled and m not in bubble:
                bubble.add(m)
                stack.append(m)                  # transitive: explore m's neighbors next
    return bubble


def build_adj(nodes: set[str], edges: set[frozenset]) -> dict[str, set[str]]:
    """Build undirected adjacency from a set of node IDs and frozenset edges."""
    adj: dict[str, set[str]] = defaultdict(set)
    for e in edges:
        t = tuple(e)
        if len(t) == 1: continue                       # self-loop in GFA
        a, b = t
        adj[a].add(b); adj[b].add(a)
    return dict(adj)


def connected_components(nodes: set[str], adj: dict[str, set[str]]) -> list[set[str]]:
    """Standard connected-components on undirected adjacency."""
    seen: set[str] = set()
    out: list[set[str]] = []
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
        out.append(comp)
    return out
