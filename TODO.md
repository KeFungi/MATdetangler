# MATdetangler — TODO

## Open

### Per-iteration timeout in `_try_one_pass` to bound classifier explosion

Heavy-combinatorial bubbles can hang `classify()` (specifically
`_enum_paths` DFS) for hours when nhop expands the bubble subgraph to
800+ nodes with many fork points. Example: BPL1195 at k=53 stuck on
nhop=10 cov-off for 20+ minutes burning 99.7% CPU — `max_paths=200`
caps OUTPUT count but DFS keeps exploring fork combinations until it
stumbles into 200 L→R paths.

Fix: wrap each `_try_one_pass` `classify()` call with a per-iteration
deadline (e.g., 60 s). On timeout, treat as `complexed` / abandon, log,
proceed to next nhop. Implementation options:
  - Signal-based: `signal.alarm()` + SIGALRM handler (works on Linux,
    not thread-safe).
  - DFS budget: pass a `node_visits_budget` counter into `_enum_paths`
    that decrements on each `dfs()` call; abort when ≤ 0. Cleaner than
    signal-based, no thread-safety issues, and bounds work deterministically.

Also: bound `_enum_dangling` the same way (it has the same DFS shape).

Quick mitigation while this is open: reduce `max_nhop` to 8 (or even 7)
for the production array, since the explosion almost always happens at
nhop ≥ 9 where the BFS subgraph balloons. Costs us nothing for clean
samples (they accept at nhop=3-4) and prevents the hangs.

### Depth-weighted path decomposition for `n > 2` candidates

When `per_k_caller`'s post-dedup count exceeds 2, fall back to depth-based
ploidy estimation instead of emitting many chimeras.

Use the GFA `DP:f:` field on each segment to estimate ploidy:
  ploidy = round(sum_{u in L_anchors} DP(u) / genome_cov)

Then keep only `ploidy` candidate paths — the ones with the strongest
depth support (sum of `DP(u)` along the path, normalized by length).
Originally rejected when AG5-shape called for "respect the topological
complexity," but the AG17 var-trim run showed the assembly graph
genuinely multi-paths the SAME biological allele (~5 kb cores, all
distinct in k-mer composition due to k=53 graph artifacts). Depth is
the principled way to collapse these to true allele count.

### Change label/text gene order

Make the gene order in allele labels / text output configurable.
Currently `seg_label_hits.tsv` emits multi-tag regions as
`"+".join(sorted(tags))` (alphabetical, e.g. `HD1+HD2`). Decide on a
convention — alphabetical, locus-order, or user-supplied — and apply
consistently across the labeler TSV, FASTA headers, and any downstream
text-based labels.

### Uvar-leaf fallback anchors when no flank is reachable

Currently the bubble classifier requires at least one flank-adjacent
node to act as an L-anchor or R-anchor. If the BFS-expanded
neighborhood contains var content but no flank context (e.g., long
unlabeled chain bounded by Uvars, or a sub-bubble exposed before the
flanks come into range), `L_anchors` and `R_anchors` are both empty,
the classifier returns 0 arms, and the verdict falls through to
`complexed, "no var-bearing arm"` with empty emission.

Fallback rule (~5 LOC in `bubble_classifier.classify`, right after the
existing R2 anchor computation):

```python
# Fallback: promote bubble-leaf Uvars to anchors when no flank reachable
if not L_anchors and not R_anchors:
    bubble_deg = {n: sum(1 for m in adj.get(n, ()) if m in bubble)
                   for n in bubble}
    leaves = {n for n in bubble & unlabeled if bubble_deg.get(n, 0) <= 1}
    L_anchors = R_anchors = leaves
    info["anchor_fallback"] = "uvar_leaves"
```

Path enumeration (`_enum_paths(adj, s, R_anchors, bubble)`) already
handles `L_anchors == R_anchors` correctly because of its
`curr != start` guard — closed paths between any two leaves are found,
duplicates are canonicalized by `_canonical`.

Add a synthetic test case `case_uvar_bounded_bubble`:
`Uvar — var — var — Uvar` with a parallel `Uvar — var — var — Uvar`
and no flank nodes; expect `closed_bubble` with 2 arms.

### graph classifier
develop robust graph classifier assumed flanks and var genes are determined; can be further generalized determine var genes and flank on the flight if var and/or flank is not known later;
make an isolated graph classifier to test algorithm

A. generate data
    1. make a fake network of
        a. close bubble: 
            flankL - HD1a - HD2a - flankR
                   \ HD1b - HD2b/

        b. open bubble:
       case 1
            flankL - HD1a - HD2a - flankR
                   \ HD1b - HD2b
       case 2
            flankL - HD1a - HD2a
                   \ HD1b - HD2b

        c. complexed:
        case 1
                   / HD1b - HD2b \
            flankL - HD1a - HD2a - flankR
                   \ HD1c - HD2c /
        case 2
                   / HD1c - HD2c
            flankL - HD1a - HD2a - flankR
                   \ HD1b - HD2b
        case 3
                   / HD1c - HD2c \
            flankL - HD1a - HD2a - flankR
                   \ HD1b - HD2b
        case 4
                   / HD1c - HD2c 
            flankL - HD1a - HD2a - flankR
                     HD1b - HD2b /
        case 5
                   / HD1c - HD2c 
            flankL - HD1a / HD2a - flankR
                     HD1b - HD2b /

         d. separate:
         flankLa - HD1a - HD2a - flankRa
         flankLb - HD1b - HD2b - flankRb

    2. record the network type; call it clean_network
    3. from clean_network, add random noisy node (random ID), but not disconnecting current network or link it to current network, only outward from clean_network, add up to five distance away; call it noisy_network
    4. from noisy_network, add random linkers, insert 1-5 random linker (random ID) between some connected nodes of clean_network; call it linker_network
    5. from linker_network, combine any random pieces of shared edge in clean_network, separate any random pieces of clean_network (and add edge between the cut pieces); call it reality_network
    6. from linker_network drop one (and only one copy and one var gene); call it degvar_network

    use the data to test the graph algorithm algorism; make the generated data the same format as used in BFS and path enumerator

B. new walking algorism
    1. expand node from HDs as before, start with nhop 1
    2. assess network
       a. not all var gene in one network -> separate
       b. all var gene in one network, inside the network (var series: longest path of var genes not interrupted by any flank) 
           i. exact 2 paths from one flankL to var series to one flankR -> close bubble
           ii. exact 2 path from one flank to var series -> open bubble
           iii. other -> complexed



## OPEN — deferred items
### Picker rewrite (T0-T7 tiered, 9-tuple)

Strategy is **completeness-first, divergence-as-hard-constraint**. Partially
realized 2026-05-31 (var-gene completeness promoted to primary score key). The
full 8-tier T0-T7 ordering and the 9-tuple are still deferred.

Full plan: T0–T7 with T4 (flanked_single_hd) > T2 (flanked_partial) — both-flanks
anchoring locks position on the assembly graph (biological event > technical
truncation). Pair selection via lexicographic 9-tuple:
`(delivered_count, joint_tier_sum, min_tier, joint_aa_cov, joint_fl_cov,
pair_divergence, joint_cov_ok, -|cov(a)-cov(b)|, -length_diff)`.

Priority order when ready:
  1. Pin tier order as T0-T7 constants with ranking comment block.
  2. Rework `best_pair_in_k` for the full lexicographic pair tuple (O(N²)).
  3. Add `--degHD-policy {anchored|allow|forbid}` (default `anchored` = T4+).
  4. Emit `extra_alleles.fasta` + WARN line when N≥3 distinct survive dedup.
  5. Add `picks_debug.tsv` (one row per candidate, not per pick).

Deferred until after (1)–(5): `--prefer-completeness`, `--min-pair-divergence`,
`--cross-k {best-pair|best-per-k-then-vote}`.

### graph_paths whole-GFA labeling fallback (LOW)
Only triggers when `graph_paths` is called standalone without `--cand-ann` and
the picks.tsv has no segments column. Under `--skip-pick` the wrapper always
passes `--cand-ann`, so this fallback never fires.

### User-supplied BED for gene/flank coordinates

- Allow `--locus-bed` to set flank + gene coordinates directly instead of
  inferring from tblastn-on-locus.
- BED rows: `chrom start end name` where `flankL` and `flankR` are reserved
  names and any other name = a variable gene.
- **IMPORTANT — don't translate non-flank features to protein.** Because
  genes have introns, the spliced protein doesn't match the genomic NT.
  Use **tblastx** (6-frame translated NT-vs-NT) when searching gene content
  across all stages where tblastn is currently used. tblastx tolerates
  introns naturally.
- Search-mode standalone flag NOT needed; when `--locus-bed` is supplied,
  tblastx replaces tblastn everywhere automatically.
- Backward compatible: if no `--locus-bed`, keep current tblastn-of-translated
  flow.

## DONE (recent)

- ✅ **2026-05-31 trim_outer_flanks** (pair-pick + singleton-pick):
  Trim each picked allele contig at outer flank-only anchors before writing
  `picks.fasta`. Both modes share `is_flank_only(sid)` = label has tokens and
  every token starts with `flank` (no var gene). Apply through `apply_trim`
  which uses GFA segment reconstruction; if the reconstructed slice would
  exceed the original contig length (SPAdes contig spans only a subset of
  large GFA segments), the trim is **aborted** and the contig is kept intact.
  * **Pair-picks**: per arm, locate `first_var` / `last_var` segment positions
    in the walk. Trim only within the outer tails (prefix before `first_var`,
    suffix after `last_var`). Pick the innermost flank-only segment in each
    tail whose ID is **shared** with the other arm (definite homology). Never
    crosses a variable-gene boundary, so an allele that legitimately has
    flanks-only-no-HD survives untouched.
  * **Singleton-picks** (C1-1, C1-4, C3-1): trim to the span
    `[leftmost flank-tagged segment : rightmost flank-tagged segment]`,
    aggressive — drops everything beyond the outermost flank anchors even if
    var genes are present in the tails. Also emits `bubble_<allele>.png` for
    each singleton.
- ✅ **2026-05-31 _normalize_orientation** flips walks where `flankL` sits
  rightward of `flankR`, so all picked alleles point `flankL → HD → flankR`
  consistently. Counts composite labels too (e.g. `HD1+flankL`). Applies
  before `trim_outer_flanks`. Critical for fragments with mixed var+flank
  content where the orientation is ambiguous from sequence alone.
- ✅ **2026-05-31 max_nodes** default 25 → 15. Tighter DFS bound keeps
  bubble-path enumeration tractable on tangled k33 graphs without losing the
  closed-bubble cases at the original 25 ceiling.
- ✅ **2026-05-31 picks.tsv segments column** always derives from
  `seg_of_cand_list[c]` when available, including for path-origin picks
  (previously only bubble-origin had it filled).
- ✅ **2026-05-31 bubble.png rendering rules**:
  * Coloring: yellow = any variable gene; blue = flank-only (no var gene); white = pure-number/unlabeled.
  * Cross-arm dashed lines: one per shared GFA-segment ID (definite homology) plus a fallback connecting the outermost flank-only node when no shared flank-only-ID exists for that flank type.
- ✅ **2026-05-31 parse_contigs_paths supports modern SPAdes**: was stripping
  `_<last-token>` from every header, which broke modern SPAdes (3.15+) names
  ending in `cov_<float>`. Now keeps the full name (strips trailing `'` RC
  marker) and still tolerates older `_<component>` suffix. Without this fix
  bubble-origin picks couldn't be resolved to GFA walks and graph_paths
  emitted empty bubble.png for them.
- ✅ **2026-05-31 picks.tsv `segments` column** populated for bubble-origin
  picks via `_segments_via_contigs_paths` (was just the kind tag "bubble").
  Lets graph_paths render the PNG without re-deriving the walk.
- ✅ **2026-05-31 wrapper omits empty `--genome-coverage`** from summary_table
  invocation (was passing empty string → invalid-float error).
- ✅ **2026-05-31 batched cluster_alleles** re-implemented with proper
  convergence: canonical pre-dedup → single-shot --auto if pool ≤ batch_size
  → otherwise lossy MAFFT (--retree 1 --maxiterate 0) in longest-first
  batches with carryover survivors. Final pairwise_identity (8a/8b) still
  uses --auto.
- ✅ **2026-05-31 widen-loop topology widen (re-added)**: break condition is now
  `n_distinct >= expected_count AND topology in {closed_bubble, complexed}`.
  Closed-bubble is the ideal; complexed means main segments have settled into
  a stable not-sharing-flanks shape and widening further won't change it.
  Other topologies (`no_main`, ...) keep widening until `max_hops`.
- ✅ **2026-05-31 widen-loop break** initially reverted to `n_distinct >= expected_count`
  alone (was: `AND topology == closed_bubble`), then re-added topology with
  closed_bubble OR complexed as accepted stable shapes. Avoids the SA91 k45
  slowdown the closed-bubble-only requirement caused.
- ✅ **2026-05-31 cluster_alleles** reverted to single-shot MAFFT MSA (was:
  batched 10-at-a-time with carryover survivors, which never converged on
  tangled k33 bubbles).
- ✅ **2026-05-31 max_hops** default 12 → 10.
- ✅ **2026-05-31 persistent BLAST DBs**: `anchor_search` now caches contig
  and segment BLAST DBs (and the GFA-extracted segments fasta) under
  `<outdir>/blast_db/` instead of tempfile.
- ✅ **2026-05-31 picker fragment pre-filter** lives in step 3's
  `_emit_picker_candidates` (drops anchor contigs with 0 variable genes AND
  <2 flanks; emits `picker_candidates.{fasta,ann.tsv}` for step 4).
- ✅ **2026-05-31 picker 7-tuple, completeness-first**:
  `(var_compl_rank, topo_rank, flank_compl_rank, -alnid_HDcore, joint_aa,
  joint_flank_bp, -|len_diff|)`. Each rank is `lo*10+hi` so balanced pairs
  beat lopsided ones. Layered MAFFT by `var_compl` (cheap, just counting).
- ✅ **2026-05-31 bubble.png homology lines**: ONE dashed line per flank type,
  preferring shared-segment-ID anchors, falling back to outermost flank node.
