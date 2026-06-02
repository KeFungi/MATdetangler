"""Segment processing — P1 directional split.

Each input segment may carry multiple BLAST/tblastn hits at different
positions. A segment with hits at distinct, non-overlapping positions is
split into a chain of single-labeled sub-segments along its stored strand.

Input
=====
seg_labels:    dict[seg_id, list[Hit]]
seg_length:    dict[seg_id, int]
edges:         set[frozenset({seg_a, seg_b})]
edge_endpoints (optional): dict[edge, (side_a, side_b)] where each side
                            is "L" or "R" — which end of each segment the
                            edge attaches to (from GFA L-line orientation).

Output
======
(nodes, edges, label_per_node, var_per_node) — graph ready for bubble-BFS.
After P1, every node is one of:
  - pure flank  (label_per_node[n] is "flankL" or "flankR")
  - pure var    (var_per_node[n] non-empty, label is the gene tag)
  - unlabeled   (neither — a connector)
No composites remain.

Direction-axis bookkeeping
==========================
Each segment has a stored strand (its sequence as written in the GFA).
A hit's `start`/`end` are coordinates on this stored strand with
0 <= start < end <= seg_length. The hit's own strand (+/-) is recorded
in `Hit.strand` but is NOT used by P1 — only positions matter for
splitting. Edge endpoints record which end of each segment (L = left of
stored strand, R = right) participates in the L-line junction, derived
from GFA orient flags.
"""
from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class Hit:
    """One BLAST/tblastn hit on a segment.

    tag      : the probe name (flankL / flankR / gene name)
    kind     : "flank" or "var"
    start    : 0-based inclusive position on segment's stored strand
    end      : 0-based exclusive position
    strand   : "+" if hit aligns to stored strand, "-" if reverse-complement
    """
    tag: str
    kind: str
    start: int
    end: int
    strand: str = "+"


def directional_split(
        seg_labels: dict[str, list[Hit]],
        seg_length: dict[str, int],
        edges: set[frozenset],
        edge_endpoints: dict[frozenset, tuple[str, str]] | None = None,
) -> tuple[set[str], set[frozenset], dict[str, str], dict[str, set[str]]]:
    """Apply P1 to the input. Returns the post-split graph data."""
    new_nodes: set[str] = set()
    new_edges: set[frozenset] = set()
    new_labels: dict[str, str] = {}
    new_vars: dict[str, set[str]] = {}
    side_to_sub: dict[tuple[str, str], str] = {}
    seg_subs: dict[str, list[str]] = {}

    # Every segment id referenced anywhere
    all_segs = set(seg_labels.keys()) | {n for e in edges for n in tuple(e)}

    for seg in all_segs:
        hits = sorted(seg_labels.get(seg, []), key=lambda h: h.start)
        if not hits:
            new_nodes.add(seg)
            seg_subs[seg] = [seg]
            side_to_sub[(seg, "L")] = seg
            side_to_sub[(seg, "R")] = seg
            continue
        if len(hits) == 1:
            h = hits[0]
            new_nodes.add(seg)
            new_labels[seg] = h.tag
            if h.kind == "var":
                new_vars.setdefault(seg, set()).add(h.tag)
            seg_subs[seg] = [seg]
            side_to_sub[(seg, "L")] = seg
            side_to_sub[(seg, "R")] = seg
            continue
        # Multi-hit: split into chain of pure sub-segments
        subs: list[str] = []
        for i, h in enumerate(hits):
            sub_id = f"{seg}__sub{i}_{h.tag}"
            new_nodes.add(sub_id)
            new_labels[sub_id] = h.tag
            if h.kind == "var":
                new_vars.setdefault(sub_id, set()).add(h.tag)
            subs.append(sub_id)
            if i > 0:
                new_edges.add(frozenset((subs[i - 1], sub_id)))
        seg_subs[seg] = subs
        side_to_sub[(seg, "L")] = subs[0]
        side_to_sub[(seg, "R")] = subs[-1]

    # Re-attach original edges to the correct sub-segments.
    # Edge endpoints are looked up by the canonical (sorted) tuple, so the
    # mapping is deterministic regardless of frozenset iteration order.
    for e in edges:
        a, b = sorted(tuple(e))                       # canonical order
        if edge_endpoints and (a, b) in edge_endpoints:
            sa, sb = edge_endpoints[(a, b)]
        elif edge_endpoints and e in edge_endpoints:
            sa, sb = edge_endpoints[e]                # legacy frozenset key (may flip)
        else:
            sa, sb = "L", "L"                         # safe default when no orientation
        new_a = side_to_sub.get((a, sa), a)
        new_b = side_to_sub.get((b, sb), b)
        if new_a != new_b:
            new_edges.add(frozenset((new_a, new_b)))

    return new_nodes, new_edges, new_labels, new_vars
