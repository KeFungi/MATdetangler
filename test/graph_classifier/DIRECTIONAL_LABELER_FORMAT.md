# Directional segment labeler — algorithm + output format

The Appendix-B classifier needs per-end label coordinates on every GFA
segment so it can split composite segments at their gene/flank boundary.
This document specifies the labeler's job, output format, and the minimal
upstream changes needed.

## What it does

For each segment in the BFS-discovered neighborhood (or, narrower, every
segment that appears in the picker's `bubble.tsv` walks), emit one TSV
row per BLAST/tblastn hit on that segment, with the hit's coordinate
range on the segment's stored strand.

The classifier reads this and applies P1 (directional split): segments
with multiple hits at distinct positions get rewritten into a chain of
single-labeled sub-segments.

## Algorithm

Per segment `S` of length `L`:

1. **tblastn** every variable-gene protein vs. `S`'s nucleotide sequence.
   Collect HSPs (Hit Sequence Pairs) at module-default cutoffs.
2. **blastn** each flank reference (the flankL / flankR queries already
   in use by `anchor_search.py`) vs. `S`'s sequence. Same cutoffs.
3. Merge overlapping HSPs of the same `(query, kind)` into a single
   coordinate span (intra-query collapse — same logic as
   `input_process.py:_hsps_to_spans`).
4. Drop HSPs below the segment's minimum-feature-length threshold
   (default 100 bp; same as `--blastn-minlen`).
5. Emit one row per surviving span.

This is the same scoring/cutoffs that the rest of the pipeline already
uses. The only new artefact is **recording start/end coordinates** instead
of just collapsing every hit into "this segment has label X."

## Output format

A TSV named `seg_label_hits.tsv` co-located with `seg_labels_<k>.tsv`:

```
seg_id    seg_length  tag       kind   start  end    strand
12582487  8843        flankL    flank  0      450    +
12582487  8843        HD1       var    520    1850   +
12582487  8843        HD2       var    1900   2650   +
12582487  8843        flankR    flank  8200   8843   +
4053731   5200        HD1       var    0      1100   +
```

Columns:

| column | type | meaning |
|---|---|---|
| `seg_id` | str  | GFA S-line ID |
| `seg_length` | int  | segment sequence length |
| `tag` | str  | the label name, e.g. `flankL`, `flankR`, `HD1`, `HD2`, or any user-defined gene name |
| `kind` | str  | `flank` or `var` (derived from which probe family the hit came from) |
| `start` | int  | 0-based inclusive start of the hit on the segment |
| `end` | int  | 0-based exclusive end of the hit on the segment |
| `strand` | str  | `+` if hit matches segment's stored strand, `-` if reverse-complement |

Segments with **no hits** are NOT listed in the file. The classifier
infers them from the GFA — any segment in the neighborhood that's
missing from `seg_label_hits.tsv` is an unlabeled connector.

## Edge-end information

The classifier also benefits from knowing **which end of each segment**
a given edge attaches to. GFA L-lines already carry this:

```
L  4053731  +  12582487  -  45M
   ^seg_a   ^orient_a  ^seg_b ^orient_b ^overlap
```

An L-line orient `+` means "the segment's right end participates in this
overlap"; `-` means the left end (after RC). The pipeline already parses
this in `GFA_search.py:parse_links`.

The classifier needs the L-line info propagated alongside the edge —
either as a parallel TSV (`bubble_edge_endpoints.tsv`) or as additional
columns on the existing `bubble.tsv`:

```
node_a   node_b   arm   label_a   label_b   end_a  end_b
12582487 15544054 arm1  flankL    flankR    R      L
...
```

Without this, `directional_split` falls back to a positional heuristic
(neighbor attaches to the sub-segment whose label coordinates are closest
to the segment end on that side). The heuristic works for clean cases
but can misroute neighbors when a composite has multiple labels near the
same end.

## Where this plugs into the pipeline

The data is already computed:

- `matdetangler/input_process.py` runs tblastn for var proteins and
  records `start/end` per HSP — that's `query_hsps_<gene>.tsv`.
- `matdetangler/anchor_search.py` runs blastn for flanks and records
  `start/end` per hit — that's `blast_flank{L,R}_blastn_segments_k<K>.tsv`.

The new labeler step just needs to **aggregate these two cached TSVs
into `seg_label_hits.tsv`** for each k. No new BLAST runs; the
information already exists in the per-step cache.

Pseudocode for the aggregator (`labeler.py`, new module):

```python
def emit_seg_label_hits(blast_cache_dir, gfa_path, out_tsv):
    seg_length = parse_segment_lengths(gfa_path)
    rows = []
    for tsv in find_blast_tsvs(blast_cache_dir):
        kind = "var" if "tblastn" in tsv else "flank"
        for hit in read_tsv(tsv):
            rows.append((hit.sseqid, seg_length[hit.sseqid],
                         hit.qseqid, kind,
                         min(hit.sstart, hit.send) - 1,
                         max(hit.sstart, hit.send),
                         "+" if hit.sstart < hit.send else "-"))
    rows = merge_same_tag_overlaps(rows)
    write_tsv(out_tsv, rows,
              header=["seg_id", "seg_length", "tag", "kind",
                      "start", "end", "strand"])
```

## Minimal upstream change

1. **New module** `matdetangler/labeler.py` (`emit_seg_label_hits`)
2. **Wire it into `MATdetangler run`** after step 2 (anchor search) and
   step 3 (GFA path search) have populated the BLAST cache.
3. **Augment `bubble.tsv`** with `end_a`, `end_b` columns reading from
   the L-line parser in `GFA_search.py:parse_links`.

That's it. The classifier consumes `seg_label_hits.tsv` + the augmented
`bubble.tsv` and the spec from Appendix B runs.

## Backward compatibility fallback

If a sample lacks `seg_label_hits.tsv`, the classifier falls back to:

1. Build per-end label info by parsing the legacy `bubble.tsv` `label_*`
   columns (the `+`-joined composite labels).
2. Heuristically place flank tags at segment ends (start/end) and var
   tags in the middle.
3. Run P1 with the heuristic positions.

The verdict is less precise (composite endpoint routing may be off in
edge cases) but the classifier still produces a sane answer for most
shapes. This lets the spec ship before the upstream plumbing lands.
