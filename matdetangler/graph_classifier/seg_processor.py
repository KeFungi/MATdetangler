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
) -> tuple[set[str], dict[str, set[tuple[str, str, str]]], dict[str, str], dict[str, set[str]],
           dict[str, tuple[str, int, int, str]]]:
    """Apply P1 to the input. Returns (nodes, directed_adj, labels, var_per, provenance).

    `directed_adj[node_id] = {(side, neighbor_id, neighbor_side), ...}`
    where side is "L" or "R".
    `provenance[sub_id] = (original_seg_id, start, end, strand)`
    """
    new_nodes: set[str] = set()
    new_adj: dict[str, set[tuple[str, str, str]]] = {}
    new_labels: dict[str, str] = {}
    new_vars: dict[str, set[str]] = {}
    side_to_sub: dict[tuple[str, str], str] = {}
    seg_subs: dict[str, list[str]] = {}
    provenance: dict[str, tuple[str, int, int, str]] = {}

    # Every segment id referenced anywhere
    all_segs = set(seg_labels.keys()) | {n for e in edges for n in tuple(e)}

    for seg in all_segs:
        hits = sorted(seg_labels.get(seg, []), key=lambda h: h.start)
        slen = seg_length.get(seg, 0)
        if not hits:
            new_nodes.add(seg)
            seg_subs[seg] = [seg]
            side_to_sub[(seg, "L")] = seg
            side_to_sub[(seg, "R")] = seg
            provenance[seg] = (seg, 0, slen, "+")
            continue

        runs: list[Hit] = []
        for h in hits:
            if runs and runs[-1].tag == h.tag:
                prev = runs[-1]
                runs[-1] = Hit(tag=prev.tag, kind=prev.kind,
                                start=prev.start,
                                end=max(prev.end, h.end),
                                strand=prev.strand)
            else:
                runs.append(h)

        if len(runs) == 1:
            r = runs[0]
            new_nodes.add(seg)
            new_labels[seg] = r.tag
            if r.kind == "var":
                new_vars.setdefault(seg, set()).update(r.tag.split("+"))
            seg_subs[seg] = [seg]
            side_to_sub[(seg, "L")] = seg
            side_to_sub[(seg, "R")] = seg
            provenance[seg] = (seg, 0, slen, r.strand)
            continue

        boundaries = [0]
        for i in range(len(runs) - 1):
            boundaries.append((runs[i].end + runs[i + 1].start) // 2)
        boundaries.append(slen)
        subs: list[str] = []
        for i, r in enumerate(runs):
            sub_id = f"{seg}#{i + 1}"
            new_nodes.add(sub_id)
            new_labels[sub_id] = r.tag
            if r.kind == "var":
                new_vars.setdefault(sub_id, set()).update(r.tag.split("+"))
            subs.append(sub_id)
            if i > 0:
                # Internal link between subsegments of same parent
                new_adj.setdefault(subs[i-1], set()).add(("R", sub_id, "L"))
                new_adj.setdefault(sub_id, set()).add(("L", subs[i-1], "R"))
            provenance[sub_id] = (seg, boundaries[i], boundaries[i + 1], r.strand)
        seg_subs[seg] = subs
        side_to_sub[(seg, "L")] = subs[0]
        side_to_sub[(seg, "R")] = subs[-1]

    for e in edges:
        t = sorted(tuple(e))
        a = t[0]
        b = t[1] if len(t) > 1 else t[0]
        if edge_endpoints and (a, b) in edge_endpoints:
            sa, sb = edge_endpoints[(a, b)]
        elif edge_endpoints and e in edge_endpoints:
            sa, sb = edge_endpoints[e]
        else:
            sa, sb = "L", "L"
        new_a = side_to_sub.get((a, sa), a)
        new_b = side_to_sub.get((b, sb), b)
        # Add directed edge in both directions (GFA L-line is undirected link)
        new_adj.setdefault(new_a, set()).add((sa, new_b, sb))
        new_adj.setdefault(new_b, set()).add((sb, new_a, sa))

    return new_nodes, new_adj, new_labels, new_vars, provenance
