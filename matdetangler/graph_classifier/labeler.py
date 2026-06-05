"""Directional labeler — aggregate BLAST cache TSVs + GFA into `seg_label_hits.tsv`.

This is the upstream piece that produces the classifier's preferred input.
It does NOT run new BLAST — it consumes the cached `blast_*.tsv` files
that `anchor_search.py` and `run_per_k.py` already produce.

Inputs
======
- gfa_path:     path to assembly_graph_after_simplification.gfa (for seg lengths)
- blast_tsvs:   list of (path, kind) pairs, where kind ∈ {"flank", "var"}.
                The tag (e.g. "flankL" or "HD1") is read from each row's qseqid.
- out_tsv:      output path

Output format (TSV, one hit per row)
====================================
seg_id    seg_length  tag       kind   start   end    strand

start/end are 0-based, half-open, on the segment's stored strand
(i.e., start < end). strand records whether the hit aligns to the stored
strand (`+`) or the reverse complement (`-`).

Hits with the same (seg_id, tag) that overlap (i.e., their intervals
intersect on the segment) are merged into a single coordinate span (the
union). Distinct non-overlapping hits of the same tag stay as separate
rows — they represent genuinely distinct positions of the same probe on
the segment.
"""
from __future__ import annotations
import os
import re
from dataclasses import dataclass


@dataclass
class _RawHit:
    qseqid: str
    sseqid: str
    sstart: int
    send: int


def parse_segment_lengths(gfa_path: str) -> dict[str, int]:
    """Return {seg_id: length} from S-lines."""
    out: dict[str, int] = {}
    with open(gfa_path) as fh:
        for line in fh:
            if not line.startswith("S\t"): continue
            parts = line.split("\t", 3)
            if len(parts) < 3: continue
            out[parts[1]] = len(parts[2])
    return out


def _iter_blast_tsv(path: str):
    """Iterate rows of a BLAST outfmt-6 TSV. Yields `_RawHit`."""
    if not os.path.isfile(path): return
    with open(path) as fh:
        header_seen = False
        for line in fh:
            if not line.strip() or line.startswith("#"): continue
            cells = line.rstrip("\n").split("\t")
            if not cells: continue
            # Some pipeline TSVs have a header; detect by non-numeric sstart
            if not header_seen:
                try:
                    int(cells[8] if len(cells) >= 12 else cells[2])
                    header_seen = True
                except (ValueError, IndexError):
                    header_seen = True
                    continue
            # outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
            if len(cells) < 12: continue
            try:
                yield _RawHit(qseqid=cells[0], sseqid=cells[1],
                              sstart=int(cells[8]), send=int(cells[9]))
            except ValueError:
                continue


def _merge_overlapping_spans(spans: list[tuple[int, int, str]]) -> list[tuple[int, int, str]]:
    """Stage A — merge same-tag overlapping or touching spans. Strand info is
    kept as the first strand of the merged group."""
    if not spans: return []
    spans = sorted(spans)
    out = [spans[0]]
    for s, e, strand in spans[1:]:
        prev_s, prev_e, prev_strand = out[-1]
        if s <= prev_e:                                   # overlap or touch
            out[-1] = (prev_s, max(prev_e, e), prev_strand)
        else:
            out.append((s, e, strand))
    return out


def _sweepline_merge_kind(spans_by_tag: dict[str, list[tuple[int, int, str]]]
                            ) -> list[tuple[int, int, frozenset, str]]:
    """Stage B — sweepline across multiple tags of the same kind. Returns
    contiguous regions [(start, end, tag_set, strand)] where the active
    tag-set is constant within each region. Adjacent regions with the same
    tag-set are merged."""
    events: list[tuple[int, int, str]] = []
    strand_of: dict[str, str] = {}
    for tag, spans in spans_by_tag.items():
        for s, e, strand in spans:
            if e <= s: continue
            events.append((s, 0, tag))                    # enter (sort first at same pos)
            events.append((e, 1, tag))                    # exit
            strand_of.setdefault(tag, strand)
    if not events: return []
    # At same position: exits before enters so adjacent ranges don't merge.
    events.sort(key=lambda x: (x[0], -x[1]))
    raw: list[tuple[int, int, frozenset]] = []
    active: set[str] = set()
    prev_pos: int | None = None
    for pos, evt, tag in events:
        if active and prev_pos is not None and pos > prev_pos:
            raw.append((prev_pos, pos, frozenset(active)))
        if evt == 0: active.add(tag)
        else: active.discard(tag)
        prev_pos = pos
    # Merge adjacent regions with same tag-set
    out: list[tuple[int, int, frozenset, str]] = []
    for s, e, tags in raw:
        strand = next((strand_of[t] for t in tags if t in strand_of), "+")
        if out:
            ps, pe, ptags, pstrand = out[-1]
            if pe == s and ptags == tags:
                out[-1] = (ps, e, tags, pstrand)
                continue
        out.append((s, e, tags, strand))
    return out


def _clip_flank_against_var(flank_regions: list[tuple[int, int, frozenset, str]],
                             var_regions: list[tuple[int, int, frozenset, str]]
                             ) -> list[tuple[int, int, frozenset, str]]:
    """Stage C — drop or clip flank regions that overlap var regions. Var
    wins on every overlapping base; flank gets clipped to what's NOT
    covered by var."""
    var_intervals = sorted([(s, e) for s, e, _, _ in var_regions])
    out = []
    for fs, fe, tags, strand in flank_regions:
        remaining = [(fs, fe)]
        for vs, ve in var_intervals:
            new_remaining = []
            for s, e in remaining:
                if ve <= s or vs >= e:                    # no overlap
                    new_remaining.append((s, e))
                    continue
                if vs > s:                                # keep left portion
                    new_remaining.append((s, vs))
                if ve < e:                                # keep right portion
                    new_remaining.append((ve, e))
            remaining = new_remaining
        for s, e in remaining:
            if e > s:
                out.append((s, e, tags, strand))
    return out


def emit_seg_label_hits(
        gfa_path: str,
        blast_tsvs: list[tuple[str, str]],
        out_tsv: str,
        min_alnlen: int = 100,
) -> int:
    """Aggregate BLAST cache TSVs into `seg_label_hits.tsv`. Returns row count.

    blast_tsvs is a list of (path, kind, tag_override?). Each path is a
    BLAST outfmt-6 TSV. The tag comes from each row's qseqid by default;
    set `kind` to "flank" or "var" to label each hit's category.
    """
    seg_length = parse_segment_lengths(gfa_path)
    # (seg_id, tag, kind) -> list of (start, end, strand)
    per_key: dict[tuple[str, str, str], list[tuple[int, int, str]]] = {}
    for path, kind in blast_tsvs:
        for hit in _iter_blast_tsv(path):
            s, e = hit.sstart, hit.send
            strand = "+" if s <= e else "-"
            start_0 = min(s, e) - 1
            end_0 = max(s, e)
            if end_0 - start_0 < min_alnlen: continue
            key = (hit.sseqid, hit.qseqid, kind)
            per_key.setdefault(key, []).append((start_0, end_0, strand))

    # Per segment: collect spans by kind, do Stage A (intra-tag merge), then
    # Stage B (sweepline across tags of same kind), then Stage C (var wins
    # over flank where they overlap on the same bases).
    by_seg_kind: dict[tuple[str, str], dict[str, list[tuple[int, int, str]]]] = {}
    for (seg, tag, kind), spans in per_key.items():
        merged = _merge_overlapping_spans(spans)               # Stage A
        by_seg_kind.setdefault((seg, kind), {})[tag] = merged

    # All segs that have any hit
    all_segs = {s for (s, _) in by_seg_kind}
    rows = []
    for seg in sorted(all_segs):
        var_regions = _sweepline_merge_kind(by_seg_kind.get((seg, "var"), {}))
        flank_regions = _sweepline_merge_kind(by_seg_kind.get((seg, "flank"), {}))
        flank_regions = _clip_flank_against_var(flank_regions, var_regions)
        slen = seg_length.get(seg, 0)
        for s, e, tags, strand in var_regions:
            rows.append((seg, slen or e, "+".join(sorted(tags)), "var", s, e, strand))
        for s, e, tags, strand in flank_regions:
            rows.append((seg, slen or e, "+".join(sorted(tags)), "flank", s, e, strand))

    # Unpacked: one row per individual tag. A region carrying N tags emits N
    # rows with the same (start, end) but different `tag` values. Downstream
    # readers aggregate by (seg_id, start, end) to recover the tag-set.
    unpacked = []
    for seg, slen, tag_str, kind, s, e, strand in rows:
        for t in sorted(tag_str.split("+")):
            unpacked.append((seg, slen, t, kind, s, e, strand))
    unpacked.sort(key=lambda r: (r[0], r[4], r[2]))      # seg, start, tag

    with open(out_tsv, "w") as fh:
        fh.write("seg_id\tseg_length\ttag\tkind\tstart\tend\tstrand\n")
        for r in unpacked:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return len(unpacked)


def read_seg_label_hits(path: str) -> dict[str, list]:
    """Inverse: read `seg_label_hits.tsv` and group hits by seg_id.

    The TSV uses long/unpacked format (one row per (tag, position)). Rows
    sharing (seg_id, start, end) represent one positional REGION carrying
    multiple tags. This reader aggregates such rows into ONE Hit object
    per region, with `tag = "+".join(sorted(tags))` and kind = "var" if
    any tag is var else "flank".

    Returns {seg_id: list of Hit}."""
    from .seg_processor import Hit
    # (seg, start, end) -> {"tags": set, "kinds": set, "strand": str}
    grouped: dict[tuple, dict] = {}
    if not os.path.isfile(path): return {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            cells = line.rstrip("\n").split("\t")
            if len(cells) < len(header): continue
            seg = cells[idx["seg_id"]]
            try:
                start = int(cells[idx["start"]]); end = int(cells[idx["end"]])
            except ValueError:
                continue
            key = (seg, start, end)
            g = grouped.setdefault(key, {"tags": set(), "kinds": set(),
                                          "strand": cells[idx["strand"]]
                                          if "strand" in idx else "+"})
            g["tags"].add(cells[idx["tag"]])
            g["kinds"].add(cells[idx["kind"]])
    out: dict[str, list[Hit]] = {}
    for (seg, start, end), g in grouped.items():
        # var wins if mixed (shouldn't happen post-Stage-C but be safe)
        kind = "var" if "var" in g["kinds"] else next(iter(g["kinds"]))
        out.setdefault(seg, []).append(
            Hit(tag="+".join(sorted(g["tags"])), kind=kind,
                start=start, end=end, strand=g["strand"])
        )
    return out


def read_segment_lengths_from_hits(path: str) -> dict[str, int]:
    """Helper: extract seg_id -> seg_length from the hits TSV."""
    out: dict[str, int] = {}
    if not os.path.isfile(path): return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            cells = line.rstrip("\n").split("\t")
            if len(cells) < len(header): continue
            seg = cells[idx["seg_id"]]
            try:
                out[seg] = int(cells[idx["seg_length"]])
            except ValueError: continue
    return out
