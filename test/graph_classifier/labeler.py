"""Directional labeler — aggregate BLAST cache TSVs + GFA into `seg_label_hits.tsv`.

This is the upstream piece that produces the classifier's preferred input.
It does NOT run new BLAST — it consumes the cached `blast_*.tsv` files
that `anchor_search.py` and `graph_path_search.py` already produce.

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
    """Merge overlapping or touching spans (sorted by start). Strand info is
    kept as the majority/first strand of the merged group."""
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

    rows = []
    for (seg, tag, kind), spans in per_key.items():
        for start, end, strand in _merge_overlapping_spans(spans):
            rows.append((seg, seg_length.get(seg, end), tag, kind, start, end, strand))
    rows.sort(key=lambda r: (r[0], r[4]))                # by seg then start
    with open(out_tsv, "w") as fh:
        fh.write("seg_id\tseg_length\ttag\tkind\tstart\tend\tstrand\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return len(rows)


def read_seg_label_hits(path: str) -> dict[str, list]:
    """Inverse: read `seg_label_hits.tsv` and group hits by seg_id.
    Returns {seg_id: list of Hit}. Imports Hit from seg_processor."""
    from .seg_processor import Hit
    out: dict[str, list[Hit]] = {}
    if not os.path.isfile(path): return out
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
            out.setdefault(seg, []).append(
                Hit(tag=cells[idx["tag"]], kind=cells[idx["kind"]],
                    start=start, end=end,
                    strand=cells[idx["strand"]] if "strand" in idx else "+")
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
