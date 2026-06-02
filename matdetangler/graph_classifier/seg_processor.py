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
) -> tuple[set[str], set[frozenset], dict[str, str], dict[str, set[str]],
           dict[str, tuple[str, int, int, str]]]:
    """Apply P1 to the input. Returns (nodes, edges, labels, var_per, provenance).

    `provenance[sub_id] = (original_seg_id, start, end, strand)` — for each
    post-P1 node, the original GFA segment ID and the sub-region of that
    segment on its stored strand. Used for arm-sequence reconstruction.
    """
    new_nodes: set[str] = set()
    new_edges: set[frozenset] = set()
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

        # Distinct-label rule (matches user spec):
        #   n_label == 1 → 1 piece  (no split; segment as whole, labeled with the tag)
        #   n_label == 2 → up to 3 pieces (left-label, middle, right-label).
        #                  Middle is omitted when the two label spans touch.
        #   n_label >= 3 → up to 3 pieces (leftmost-tag at L end, rightmost-tag
        #                  at R end, all other label spans absorbed into the
        #                  middle, which is left unlabeled).
        # Multiple hits of the same tag collapse to that tag's outer span
        # (min start, max end). Strand for a collapsed span is taken from the
        # tag's longest single hit.
        tag_first: list[str] = []
        tag_span: dict[str, list] = {}   # tag -> [min_start, max_end, strand, longest_aln, kind]
        for h in hits:
            if h.tag not in tag_span:
                tag_first.append(h.tag)
                tag_span[h.tag] = [h.start, h.end, h.strand, h.end - h.start, h.kind]
            else:
                s = tag_span[h.tag]
                s[0] = min(s[0], h.start); s[1] = max(s[1], h.end)
                aln = h.end - h.start
                if aln > s[3]: s[2] = h.strand; s[3] = aln

        n_label = len(tag_span)

        if n_label == 1:
            tag = tag_first[0]; sp = tag_span[tag]
            new_nodes.add(seg)
            new_labels[seg] = tag
            if sp[4] == "var":
                new_vars.setdefault(seg, set()).update(tag.split("+"))
            seg_subs[seg] = [seg]
            side_to_sub[(seg, "L")] = seg
            side_to_sub[(seg, "R")] = seg
            provenance[seg] = (seg, 0, slen, sp[2])
            continue

        # n_label >= 2 — order tags by their span's start coordinate; the
        # leftmost-start tag holds the L-end piece, the rightmost-end tag
        # holds the R-end piece. Any other tags (n_label > 2) get absorbed
        # into an unlabeled middle.
        tags_by_left = sorted(tag_span.items(), key=lambda kv: kv[1][0])
        L_tag, L_sp = tags_by_left[0]
        R_tag, R_sp = max(tag_span.items(), key=lambda kv: kv[1][1])
        if L_tag == R_tag:
            # Two distinct tags shouldn't collapse to the same span endpoint —
            # but if a single tag's span dominates both ends (e.g. one wide hit
            # encompasses another tag entirely), fall back to first/last by start.
            R_tag, R_sp = tags_by_left[-1]

        # Cut points: midpoint between L-end and the next-leftmost label's
        # start (gives a clean L-piece boundary) and between the last
        # non-R label's end and R-end's start.
        # Simplest concrete recipe:
        #   left_cut  = (L_sp.end + R_sp.start) // 2  when only 2 tags
        #               and they touch (L_sp.end >= R_sp.start) → 0-bp middle.
        #   For >=3 tags, left_cut goes between L_sp.end and the next
        #   non-L label start; right_cut between previous non-R label end
        #   and R_sp.start.
        if n_label == 2:
            left_cut  = max(L_sp[1], (L_sp[1] + R_sp[0]) // 2)
            right_cut = min(R_sp[0], (L_sp[1] + R_sp[0]) // 2)
        else:
            other_starts = [sp[0] for tag, sp in tag_span.items() if tag not in (L_tag, R_tag)]
            other_ends   = [sp[1] for tag, sp in tag_span.items() if tag not in (L_tag, R_tag)]
            left_cut  = (L_sp[1] + min(other_starts)) // 2 if other_starts else L_sp[1]
            right_cut = (max(other_ends) + R_sp[0]) // 2 if other_ends else R_sp[0]
        left_cut  = max(L_sp[1], min(left_cut, R_sp[0]))
        right_cut = max(L_sp[1], min(right_cut, R_sp[0]))
        if right_cut < left_cut: right_cut = left_cut  # touching spans → empty middle

        # Sub-node IDs: clean position-indexed form "{seg}#{N}" (1-based).
        # No coordinates encoded in the name — coords live in provenance only.
        # Order: #1 = L-piece, #2 = middle (if present), trailing # = R-piece.
        subs: list[str] = []
        # L piece — always emitted
        L_id = f"{seg}#1"
        new_nodes.add(L_id); new_labels[L_id] = L_tag
        if L_sp[4] == "var":
            new_vars.setdefault(L_id, set()).update(L_tag.split("+"))
        provenance[L_id] = (seg, 0, left_cut, L_sp[2])
        subs.append(L_id)

        # Middle piece — only if non-empty
        if right_cut > left_cut:
            M_id = f"{seg}#2"
            new_nodes.add(M_id)
            # Middle keeps any var labels whose span lies between left_cut
            # and right_cut — losing them would erase real HD content that
            # sits between flank-labeled ends (e.g. a single GFA segment
            # spanning flankL + HD1 + HD2 + flankR would otherwise emit a
            # var-less node and the classifier would report no_var).
            # Flank labels in the middle ARE dropped (the canonical layout
            # puts flanks at the ends of the locus, not in its middle).
            middle_var_tags: set[str] = set()
            middle_label_toks: list[str] = []
            middle_strand = "+"
            best_aln = 0
            for tag, sp in tag_span.items():
                if tag in (L_tag, R_tag): continue
                # treat a label as "in middle" if its span overlaps [left_cut, right_cut)
                if sp[1] <= left_cut or sp[0] >= right_cut: continue
                if sp[4] == "var":
                    middle_var_tags.update(tag.split("+"))
                    if sp[3] > best_aln: best_aln = sp[3]; middle_strand = sp[2]
                middle_label_toks.append(tag)
            if middle_label_toks:
                new_labels[M_id] = "+".join(middle_label_toks)
            if middle_var_tags:
                new_vars[M_id] = middle_var_tags
            provenance[M_id] = (seg, left_cut, right_cut, middle_strand)
            subs.append(M_id)
            new_edges.add(frozenset((subs[-2], M_id)))

        # R piece
        R_id = f"{seg}#{len(subs) + 1}"
        new_nodes.add(R_id); new_labels[R_id] = R_tag
        if R_sp[4] == "var":
            new_vars.setdefault(R_id, set()).update(R_tag.split("+"))
        provenance[R_id] = (seg, right_cut, slen, R_sp[2])
        new_edges.add(frozenset((subs[-1], R_id)))
        subs.append(R_id)

        seg_subs[seg] = subs
        side_to_sub[(seg, "L")] = subs[0]
        side_to_sub[(seg, "R")] = subs[-1]

    # Re-attach original edges to the correct sub-segments.
    # Edge endpoints are looked up by the canonical (sorted) tuple, so the
    # mapping is deterministic regardless of frozenset iteration order.
    for e in edges:
        t = sorted(tuple(e))
        if len(t) == 1: continue                      # self-loop in GFA — skip
        a, b = t                                       # canonical order
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

    return new_nodes, new_edges, new_labels, new_vars, provenance
