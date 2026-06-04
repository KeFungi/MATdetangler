# MATdetangler — methods (code-grounded)

This file describes the pipeline as the code at `master` actually implements
it. Where it disagrees with `METHODS.md`, this file is correct — `METHODS.md`
has drifted across recent refactors (GFA-only genome cov, max_nhop=8,
PYTHONHASHSEED lock, MATdetangler-cli rename, step 9 summary.json, etc.).

The top-level entrypoint is `MATdetangler-cli` (bash wrapper, ~770 lines).
It runs 9 steps per sample; steps 1–5 are reads-free, step 6+ need reads
(`--make-consensus`, default OFF). All Python invocations inherit
`PYTHONHASHSEED=0` (locked at `MATdetangler-cli:21`) so set/dict iteration
order is deterministic across runs.

---

## Step 1 — `input_process` → queries/

**Code**: `matdetangler/input_process.py:tblastn_locate()`

Builds the per-sample query set from the user-supplied locus reference and
HD-protein FASTAs. Reads-free.

1. `tblastn` proteins → locus_ref (evalue 1e-5).
2. **Confidence-ordered greedy HSP acceptance**: rank HSPs by `pident × aln_aa` descending; accept each HSP unless it overlaps a previously-accepted HSP **from a different protein** by > 0.5 × shorter HSP length. This prevents weak partial HSPs of HD1 from leaking into HD2's reported span when the two proteins are paralogous.
3. Per-protein span = `[min(start), max(end)]` over that protein's accepted HSPs.
4. **HD envelope** = `[min(start), max(end))` across all proteins. Padded by 500 bp each side internally.
5. **Flanks** = locus sequence OUTSIDE the envelope, trimmed to `--max-flank-len` (default 2000) per side.
6. Sanitize protein IDs (alphanumeric + `_.-` only; auto-suffix on collision).

**Inputs**: `--locus-ref FA`, `--proteins FA`, `--max-flank-len 2000`.

**Outputs** under `<outdir>/<sample>/queries/`:
- `variable_proteins.fasta` — sanitized protein copy
- `variable_nt.fasta` — nt sequence of each tblastn-accepted span (strand-corrected)
- `flankL.fasta`, `flankR.fasta` — derived flank windows
- `manifest.json` — metadata (envelope coords, per-gene spans, `derived_max_locus_len`)

**Defaults** (hardcoded in `tblastn_locate`): `tblastn_pid=30`, `tblastn_aa=50`, `flank_pid=85`, `flank_minlen=100`.

---

## Step 2 — Per-k genome coverage

**Code**: `matdetangler/input_process.py:estimate_genome_cov_from_gfa()` (lines 67–97). Wired up by `MATdetangler-cli:~436`.

Reads-free, **GFA-only** (no contigs.fasta dependency):

1. Parse S-lines of `assembly_graph_after_simplification.gfa[.gz]`.
2. For each segment with sequence length ≥ `--contig-depth-size-cut` (default 5000 bp):
   - Extract depth from `DP:f:` tag if present.
   - Else `KC:i: / seg_length` or `RC:i: / seg_length`.
3. **`D_k` = median(depths)** across qualifying segments. Median (not mean) is robust to the heavy right tail of collapsed-repeat segments.
4. `D_k = 0.0` if no segment qualifies → disables coverage-based filtering downstream (graceful degrade).

`D_k` is bp-equivalent units — same scale as the per-segment `DP:f:` values the per-K caller uses for its depth filter band `[lo_mult × D_k, hi_mult × D_k]`. Cached to `<outdir>/<sample>/genome_cov_spades_k<K>.txt`.

User can override globally with `--genome-coverage FLOAT` (skips the estimator for all k).

---

## Step 3 — `run_per_k` → per-K caller

**Code**: `matdetangler/run_per_k.py:main()` → `matdetangler/graph_classifier/per_k_caller.py:find_alleles()`. The classifier itself lives in `matdetangler/graph_classifier/`.

This step has four sub-stages.

### 3a — Per-segment BLAST + label aggregation

**Code**: `run_per_k.py` lines 121–144 (BLAST) + `matdetangler/graph_classifier/labeler.py:emit_seg_label_hits()`.

1. Extract S-line sequences from the GFA → temp FASTA → `makeblastdb`.
2. `blastn flankL/flankR.fasta` and `tblastn variable_proteins.fasta` against the segment DB.
3. `labeler.py` aggregates the three BLAST TSVs into `seg_label_hits.tsv` via a 3-stage merge:
   - **Stage A**: merge overlapping same-tag intervals per (seg_id, tag).
   - **Stage B**: sweepline across tags of the same kind (flank or var), producing contiguous regions with constant active tag-set.
   - **Stage C**: var wins on overlapping bases (clip flank regions against var regions).
4. Output `seg_label_hits.tsv`: `seg_id, seg_length, tag, kind, start, end, strand`. Multi-tag regions are emitted as multiple rows (one per tag) at the same (start,end).

### 3b — Directional split (P1)

**Code**: `matdetangler/graph_classifier/seg_processor.py:directional_split()`.

A SPAdes segment can carry multiple BLAST labels (e.g. flankL on the left end, HD1 in the middle, HD2 on the right). The classifier needs each post-P1 node to play **at most one role** (pure flank, pure var, or pure unlabeled connector).

Algorithm — for each GFA segment:

1. **No hits** → stays as one node with empty label. `provenance[seg] = (seg, 0, slen, "+")`.
2. **Hits** → sort by `start`, then merge consecutive **same-tag** hits into "runs" (multiple HD1 fragments collapse into one HD1 run). After merging, `runs` is the list of distinct-label spans on the segment.
3. **One run** → stays as one node with that tag. `label_per_node[seg] = tag`. If `kind == "var"`, `var_per_node[seg] = {tag pieces}`.
4. **Multi-run** (composite) → split into N sub-nodes `<seg>#1, #2, …, #N` indexed L→R on the segment's stored strand. Boundaries are at MIDPOINTS between adjacent runs' end/start. Each sub-node gets one tag. New within-segment edges connect consecutive sub-nodes (`#i — #(i+1)`).
5. **Re-attach original GFA L-edges**: each L-line has `edge_endpoints[(a,b)] = (side_a, side_b)` recording which END of seg a/b the edge touched (L or R). Map to `subs[0]` (L side of the parent) or `subs[-1]` (R side).

Provenance: `provenance[sub_id] = (parent_seg_id, start, end, strand)` — for each post-P1 node, the original GFA segment ID and the sub-region on its stored strand. Used by `_arm_sequence_for` downstream to materialize allele sequences.

Returns `(nodes, edges, label_per_node, var_per_node, provenance)`.

### 3c — Two-pass BFS + classify loop

**Code**: `per_k_caller.py:find_alleles()` (lines 306–571).

The orchestrator iterates over **(cov_pass × nhop)** combinations and collects every candidate emission, then ranks across the pool. Two key design choices control the search:

**Cov-pass dimension** (`cov_filter` arg, default `"on"`):
- Pass 1 — **cov-OFF main pass**: BFS neighborhood is NOT depth-filtered (more permissive — recovers fragmented-flank samples where SPAdes broke a flank into too-short segments).
- Pass 2 — **cov-ON fallback**, only if pass 1 didn't short-circuit: BFS neighborhood is filtered to `[lo_mult × D_k, hi_mult × D_k]` (defaults 0.2× to 2.0×) — tightens when cov-off neighborhoods are too noisy.
- `--cov-filter off` → only cov-OFF pass; no fallback.

**Nhop dimension** (`init_nhop=3`, `max_nhop=8`):
- For each cov pass, iterate `nhop` from 3 to 8.
- `nhop` is the BFS hop count radius around seed segments.

**Inner loop body** (per `(cov_pass, nhop)`):
1. **Seed** the BFS from the set determined by `seeds_mode` (default `"both"` = var-labeled + flank-labeled segments).
2. **Expand** N hops → neighborhood node set.
3. If `cov_pass=True`: apply depth filter.
4. **Restrict subgraph** to neighborhood + reachable edges.
5. **Apply P1** (Step 3b) → post-P1 nodes/edges/labels/var/provenance.
6. **Classify** (Step 3d).
7. If verdict has `class != "no_var"` and emission produces a `closed_bubble n>=2 complete_locus=2` candidate → **short-circuit** the entire loop.

### 3d — Bubble classification (R1–R4)

**Code**: `matdetangler/graph_classifier/bubble_classifier.py:classify()`. Helpers: `_enum_paths()`, `_enum_dangling()`, `bubble_bfs()` (in `bubble_bfs.py`).

Inputs (post-P1):
- `nodes` (set), `edges` (frozensets), `label_per_node` (str per node — `"flankL"`, `"flankR"`, gene names, or empty), `var_per_node` (per-node set of var-gene tags)
- From these the classifier derives:
  - `var_nodes = {n : var_per_node[n] is non-empty}`
  - `flankL = {n : "flankL" in label_per_node[n].split("+")}` ; `flankR` similarly
  - `unlabeled = {n : label_per_node[n] empty AND not in var_nodes}`
  - `adj` — undirected adjacency from edges

**R1 — separate.** If `var_nodes` live in ≥2 disjoint full-graph components (`connected_components(nodes, adj)`), recurse `classify()` on each network's induced subgraph and return `{class: "separate", sub_results: [...], var_components: [...]}`.

**P2 — bubble BFS.** `bubble_bfs(adj, var_nodes, unlabeled)`:
- Start with `bubble = set(var_nodes)`.
- DFS through neighbors that are in `unlabeled` (not in bubble yet); add them, continue from there.
- Stops at any flank-labeled node (flanks are NOT in `unlabeled`, so the walk can't cross them).
- Result: `bubble = var_nodes ∪ (unlabeled connectors transitively reachable from var via unlabeled-only steps)`. Everything else is "outside" (flank or unreachable).

**R2 — anchors.** Anchors are bubble nodes that "meet the outside world":
- `flank_L_adj = {n ∈ bubble : any neighbor m ∈ flankL}` — bubble nodes whose outside-the-bubble neighbor is a flankL node. Symmetric for `flank_R_adj`.
- **Rule (a)** anchor: `n ∈ bubble` with `deg_out(n) > 0`, where `deg_out = #{m ∈ adj[n] : m ∉ bubble}`. The common case — a bubble node with at least one neighbor outside the bubble (typically the flank-labeled neighbor itself).
- **Rule (b)** anchor (FALLBACK, only if rule (a) gives `< 2` anchors): admit bubble-leaves — nodes with `deg_in(n) ≤ 1` where `deg_in = #{m ∈ adj[n] : m ∈ bubble}`. Dead-end leaves of the bubble subgraph. Necessary to recover anchors when the assembly broke the flank into short segments that didn't BLAST (lost-label case).

  Rule (b) is gated to avoid spurious anchors on clean closed bubbles. A 90 bp leaf hanging off a Y-fork in a sample with two healthy flank anchors should NOT become a third anchor — that would add a spurious dangling arm.

**R3 — arms.** For each ordered pair of distinct anchors `(s, e)`, enumerate simple paths from `s` to any `e ∈ anchors - {s}` that stay within `bubble`. The path enumerator is `_enum_paths(adj, start=s, ends, allowed=bubble, …)` — see "Path enumeration algorithm" below.

After enumeration, filter each path:

- **Drop** paths that don't carry any var node (`any(n in var_nodes for n in p)` must be true).
- **Closed arm** — both endpoints are flank-adjacent on OPPOSITE sides: `{p[0], p[-1]} ∩ flank_L_adj ≠ ∅ AND ∩ flank_R_adj ≠ ∅`. A proper L↔R walk through the bubble.
- **Dangling arm** — one endpoint is flank-adjacent (L or R), the OTHER is a bare leaf (not flank-adjacent). One side anchored, the other dangling into the bubble interior.
- **Dropped** otherwise — same-side (both endpoints L-adj or both R-adj) and pure leaf-leaf paths.

Canonicalized via `_canonical(p)` (orient so first node-id < last node-id) to dedupe reversed-equivalent walks.

Also handled:

- **Single-node bridge**: a var node that's an anchor with BOTH flankL AND flankR adjacency (a single composite segment spanning the whole locus) → counts as a closed arm with the 1-element path `[n]`.
- **`_enum_dangling` extension**: from each flank-anchored start, enumerate paths that **dead-end inside the bubble** (no neighbor in `bubble` left to extend to) carrying var content NOT entirely covered by the closed arms' var union. Adds to `dangling`. Dedupe via canonical ordering. This catches arms whose other endpoint isn't an anchor at all — assembly truly ends inside the bubble.

**R4 — verdict.**
- `n_arms = n_closed + n_dangling`
- `n_arms == 0` → `complexed` (no var-bearing arm)
- `n_arms == 1` → `single`
- `n_arms == 2 AND n_closed == 2` → `closed_bubble`
- `n_arms == 2` (else) → `open_bubble` (one or both dangling)
- `n_arms ≥ 3` → `complexed`

Note: `n_closed` and `n_dangling` are kept as diagnostic fields in the verdict dict (`closed_arms`, `dangling_arms`); both feed into the picker's "complete_locus" tri-state via flank-presence checks.

### 3d.1 — Path enumeration algorithm (`_enum_paths`)

**Shortest-first BFS by path length.** Frontier is a list of `(path, vis_set, bp_since_var)` triples; iterate until frontier empty or `max_paths` reached.

Initial frontier:
- `start_bp = 0` if `start ∈ var_nodes`, else `node_bp[start]`
- frontier = `[([start], {start}, start_bp)]`

Step (each iteration expands one frontier into next_frontier):

For each `(path, vis, bp_acc)` in frontier, for each `nxt ∈ adj[path[-1]]`:

1. **Simple-path guard**: skip if `nxt ∈ vis` (no revisits).
2. **Allowed set**: skip if `nxt ∉ allowed AND nxt ∉ ends` (only walk through the bubble or land on an anchor).
3. **Path-length cap**: if `len(path) + 1 > max_path_length` → drop the extension (and increment `max_path_length_hit`). The partial path is also pruned from `next_frontier`.
4. **BP-aware cap** (when `node_bp` + `var_nodes` are supplied AND `max_bp_since_var > 0`):
   - If `nxt ∈ var_nodes` → ALWAYS extend; reset `new_bp = 0` (reaching var resets the budget).
   - Else if `bp_acc > max_bp_since_var` → drop (and increment `max_bp_hit`). **Pre-extension check on the CURRENT accumulator, not `bp_acc + bp(nxt)`.** A path that hasn't yet blown the budget can still grab one MORE node of any size (e.g. a final 10 kb anchor); the cap only kicks in for further extensions afterward.
   - Else `new_bp = bp_acc + node_bp[nxt]`.
5. **Termination vs continuation**:
   - If `nxt ∈ ends AND nxt ≠ start` → emit `path + [nxt]` into `paths` list; if `len(paths) ≥ max_paths` → set `max_paths_hit` and abort the whole enumeration.
   - Else → push `(path + [nxt], vis ∪ {nxt}, new_bp)` onto `next_frontier`.

Returns `paths` — a list of node-ID lists, each a simple bubble-internal walk from `start` to some `e ∈ ends`.

### 3d.2 — Dead-end enumeration (`_enum_dangling`)

Similar BFS-by-length, but:

- **No `ends` set** — instead, enumerate paths that DEAD-END (no extending neighbor in `allowed`).
- Per step, track `extended` (bool): set True if any `nxt` was accepted onto `next_frontier`. If a path body had NO valid extension AND `len(path) > 1`:
  - Check var content: `vs = set(path) ∩ var_nodes` must be non-empty AND NOT a subset of `exclude_var_subset` (the closed arms' var union — passed in by the caller). Drops dangling arms whose var content is wholly redundant with already-found closed arms.
  - Canonicalize `tuple(path)` (orient by endpoints), insert into `seen` set, append to `paths`.
- Same three caps: `max_paths`, `max_path_length`, `max_bp_since_var`.

Used in R3 to catch arms whose "outer" endpoint isn't actually another anchor — it's just where the bubble runs out.

### 3e — Emission + dedup + ranking

**Code**: `per_k_caller.py:_emit_result()` + `_finalize_candidates()`.

1. Each `(cov_pass, nhop)` iteration emits a list of "pools" (one per sub_result of a `separate` verdict, otherwise one global pool).
2. Per pool:
   - **Trim each arm to its var span** (first var → last var), then extend through joint-detection BFS to recover the full allele walk (joints = nodes visited by ≥ 2 arms' outward BFS from arm endpoints — see §3.8 of legacy METHODS.md for the multi-source detail).
   - **`samtools consensus`-style sequence reconstruction** from provenance: for each post-P1 sub-node, slice the original GFA segment.
   - **tblastn-trim** with `locus_padding=4000` to clean up flanks beyond the HD envelope.
   - **Dedup** RC-aware via edlib HW edit distance — collapse pairs with identity ≥ (1 − `divergence_threshold`), default threshold 0.01 (1% edit distance).
   - **Min-bp floor**: drop alleles shorter than `min_allele_bp` (default 0 = disabled).
   - Compute per-allele `complete_var` (HD-tag coverage tri-state 0/1/2), `complete_locus` (flank-presence tri-state).
3. **Cross-iteration ranking** (line 506) — 4-tier `_rank_key`, smaller-is-better:
   - 0. `-complete_locus` — both flanks present in EVERY emitted allele (tri-state 2) wins
   - 1. `-complete_var` — all expected HD tags in EVERY allele (tri-state 2) wins
   - 2. `bubble_priority` — `closed_bubble:1 > open_bubble:2 > separate_div2:3 > single:4 > complexed_div2:5 > separate_other:6 > complexed_other:7 > no_var:8`
   - 3. `diploid_dist` — `|mean(allele_cov) / D_k − 0.5|` (closer to 0.5 = better diploid balance)
4. **Short-circuit acceptance** (line 516): the moment any candidate hits `closed_bubble + n_dedup ≥ 2 + complete_locus ≥ 2`, break the (cov_pass, nhop) loop. `complete_var` is NOT required for short-circuit — `complete_locus` is the gate (both flanks reached) and `complete_var` is a secondary tier in ranking.
   - Rationale: a diploid sample with one allele truncated mid-HD (cv=1) is the correct answer, NOT a homozygote-collapsed `single n=1 cv=2`. Forcing cv=2 collapses the diploid signal. See `diploid-preserve-over-cv` memory.

**Outputs** under `<outdir>/<sample>/<k>/`:
- `result.tsv` (22 cols) — sample, k, bubble_type, n_dedup, complete_var, complete_locus, basepair, allele_cov, etc.
- `alleles.fasta` — emitted sequences (`<sample>_<k>_allele1`, `_allele2`, or `_chimera1+` when n > 2)
- `seg_label_hits.tsv`, `flankL_blastn.tsv`, `flankR_blastn.tsv`, `HD_tblastn.tsv` — labeling caches
- `subnode_seqs.fasta` — materialized sub-segment sequences when P1 split anything

---

## Step 4 — `pick_k` cross-K consolidation

**Code**: `matdetangler/pick_k.py:pick_for_sample()` (lines 127–141), scoring at `_score()` (lines 104–124).

Runs only when `--no-skip-pick` is set (default is `--skip-pick`, which leaves the candidate pool from step 3 untouched for QC inspection).

1. Load all per-k `result.tsv` rows + `alleles.fasta` records.
2. Apply the same 4-tier rank key as `per_k_caller._rank_key` (see 3e).
3. Pick the winning k. All alleles from that k get emitted; alleles from losing k's are discarded.

**Outputs**:
- `picks.tsv` — one row per picked allele (legacy 12-col schema: sample, allele, origin, k, type, len, from_contig, segments, cov, n_variable_genes, has_both_flanks, is_degHD)
- `picks_summary.tsv` — sample-level (k_chosen, bubble_type, n_dedup, complete_var, complete_locus, locus_coverage, basepair, genome_cov, allele_cov, n_cand, extend_bounds, components, all_k_tried)
- `primary_alleles.fasta` — concatenated picks (1–N records)
- `longest_alleles.fasta` — length-first dedup of the candidate pool (broader net for downstream cluster work)

---

## Step 5 — `graph_paths` annotated bubble views

**Code**: `matdetangler/graph_paths.py:run()` (entry at the `_cli` block, lines 1000+).

Renders each picked allele's walk on the GFA. Reads-free.

1. For each picked allele:
   - Either use the recorded `segments` from `picks.tsv:segments` verbatim (preferred), or fall back to BLAST-trace alignment.
   - Build the ordered node list through the GFA.
2. **Joint detection**: multi-source BFS from each arm's endpoints (`adj_pp` post-P1) — a node visited by ≥ 2 arms' BFS frontiers is a candidate joint. Pair-search picks the globally-best (j_a, j_b) pair shared between the two alleles.
3. Render five views of the same bubble:
   - `bubble.txt` — ASCII walk (one row per allele)
   - `bubble.gfa` — sub-GFA (loads in Bandage; only the segments + edges that touch the walks)
   - `bubble.dot` — Graphviz DOT
   - `bubble.tsv` — edge list (node_a, node_b, arm, label_a, label_b)
   - `bubble.png` — matplotlib figure
4. **PNG coverage annotation**: below each segment box, a normalized depth tag is drawn — `×<DP:f: / genome_cov>` (e.g. `×0.48`, `×1.02`, `×2.07`). `genome_cov` is read from `<outdir>/genome_cov_spades_k<k>.txt` (or `--genome-cov` override). For a low-het diploid the median GFA depth ≈ 2× haploid coverage, so `×0.5` ≈ a unique bubble-allele segment (only one haploid's reads), `×1.0` ≈ shared-between-alleles or homozygous background, `×2.0` ≈ collapsed two-copy repeat. Helpful for spotting low-coverage tail noise vs real content.

Row layout in the PNG: cross-row solid black lines mark sub-nodes literally shared between rows (the cycle joints). Node face colors: yellow = var-bearing, blue = flank-only, white = unlabeled. Up to 4 rows displayed (extras dropped; title notes truncation).

---

## Step 6 — `map_consensus` (`--make-consensus` only)

**Code**: `matdetangler/map_consensus.sh` (50 lines).

Skipped by default (`SKIP_CONSENSUS=1` in the wrapper). Enable with `--make-consensus` (also requires `--reads-r1`/`--reads-r2`).

1. `bowtie2-build` index ← `primary_alleles.fasta` (the picked references).
2. `bowtie2 --end-to-end --very-sensitive -k 1 --no-unal -q -p <threads>` → `reads.sam`. Competitive mapping: one best hit per read pair, no unmapped records written.
3. `samtools sort reads.sam → reads.sorted.bam` + `samtools index`. The BAM is an **internal intermediate** — used by steps 6.b–6.c and deleted at end of step 6.
4. **Per-allele coverage** via `matdetangler/coverage_core.py` if `queries_dir` has `variable_proteins.fasta` (computes both whole-allele and HD-core-only depths, with optional `--repeats` masking); else fall back to `samtools coverage` (whole-allele only).
5. **Per-allele consensus**: for each record in `primary_alleles.fasta`, run `samtools consensus -r <allele> -f fasta`. Record IDs are **preserved verbatim** from the input reference (`<sample>_<k>_allele1`, etc.) — NOT renamed to `allele_1/2/3`.
6. Cleanup: delete `reads.sorted.bam[.bai]` and the bowtie2 index files. Keep `reads.sam` (human-readable artifact).

**Persisted outputs**:
- `reads.sam` — alignments (re-derive BAM with `samtools sort + index` if needed)
- `consensus_alleles.fasta` — concatenated per-allele consensus
- `coverage.tsv` — `#allele, len, mapped_reads, whole_breadth_pct, whole_meandepth, core_bp, core_meandepth, repeat_bp`

---

## Step 7 — Identity + consensus QC (`--make-consensus` only)

Three parallel sub-steps, all reads-derived (require step 6's outputs).

### 7a — Pairwise identity on the picks

**Code**: `matdetangler/pairwise_identity.py`.

MAFFT-aligns the picked allele pair (`primary_alleles.fasta`):
- `mafft --auto --adjustdirection` — RC-aware alignment.
- `id_pct` = matched aligned columns / total aligned columns (excluding pure-gap columns).
- `aln_frac` = aligned columns / shorter input sequence length.
- `distinct` (bool) = `NOT (id_pct ≥ 95 AND aln_frac ≥ 0.80)`.

Output: `identity.tsv` + `_identity.json`.

### 7b — Pairwise identity on the consensus

Same logic, against `consensus_alleles.fasta` (skipped if < 2 records). Output: `identity_consensus.tsv` + `_identity_consensus.json`.

### 7c — Consensus completeness re-check

**Code**: `matdetangler/consensus_qc.py:_qc_one_consensus()` (lines 41–77).

For each consensus record:
1. `tblastn variable_proteins.fasta → consensus` (pid ≥ 30, aln ≥ 50 aa). Track which HD genes hit.
2. `blastn flankL.fasta → consensus` (pid ≥ 85, aln ≥ 100 bp). Track flankL presence.
3. `blastn flankR.fasta → consensus`. Track flankR presence.
4. `complete` = `(all expected HD genes hit) AND (flankL present) AND (flankR present)`.

Output: `consensus_qc.tsv` — per-allele complete flag + per-gene/flank presence.

---

## Step 8 — `summary_table` → `summary.tsv`

**Code**: `matdetangler/summary_table.py:write_row()` (lines 49–120).

Reads everything from steps 4–7 and emits one wide-format TSV row per sample.

**Columns** (17 columns):
`sample, used_k, bubble_type, genome_coverage, allele1_complete, allele2_complete, allele1_coverage, allele2_coverage, allele1_vs_allele2_id_pct, allele1_vs_allele2_aln_frac, cons_allele1_complete, cons_allele2_complete, cons_allele1_vs_allele2_id_pct, cons_allele1_vs_allele2_aln_frac, allele1_path, allele2_path` (cons_* columns are `"-"` when `--make-consensus` is OFF).

Written twice:
- `<outdir>/<sample>/summary.tsv` — per-sample (single row)
- `<outdir>/summary.tsv` — aggregate (appended)

---

## Step 9 — `summarize` → `summary.json` (cross-run JSON)

**Code**: `matdetangler/summarize.py`. Two CLI modes.

### Per-sample mode (wrapper Step 9)

Invoked as `python -m matdetangler.summarize --sample-dir <SDIR> --out <SDIR>/summary.json` after Step 8.

Walks `<SDIR>` and parses:
- `picks_summary.tsv` (sample-level: k_chosen, bubble_type, complete_var/locus, etc.)
- `summary.tsv` (pairwise identity, allele_path strings)
- `picks.tsv` (per-allele segments, cov, has_both_flanks, etc.)
- `logs/03_run_per_k_k<K>.log` — per-iteration BFS trace via regex (`nhop, seeds, cov_pass, nhood, n_var, cls, n_arms`), `finished_nhop` from `[accept]` lines, `per_k_pick` from `[pick]` lines.

Output: a single JSON with `{args, version, sample}` keys. See `test/Pcub40/known_results.json` for the schema in action.

### Cross-sample aggregate mode (test harness)

Invoked as `python -m matdetangler.summarize --results-dir <outdir> --out aggregate.json`. Walks each `<outdir>/<sample>/` and packs all per-sample entries under `samples: {…}`.

Used by `test/Pcub40/installation_run_test.sh` and `test/Pcub40/analysis_run_test.sh` to diff a fresh run against `test/Pcub40/known_results.json`.

---

## Cross-sample clustering — `MATdetangler-cli cluster`

Standalone, post-`batch`. Two subcommands.

### `cluster align`

**Code**: `matdetangler/cluster.py:align()` (lines 52–96).

1. Walk `--results-dir`, pull each sample's `primary_alleles.fasta`.
2. For each allele, identify the HD-core via `tblastn variable_proteins → allele` and extract a window of `core ± 50 bp` (configurable via `--core-pad`).
3. Run **one** MAFFT alignment over the union of all sample cores (`cores.fasta` → `cores.aln.fasta`).
4. Compute pairwise similarity from the alignment: `sim_ij = # identical aligned columns between i and j / min(core_len_i, core_len_j)`.
5. Emit `cores_pairs.tsv` (every pair: `a, b, sim, alnid`), `cores_meta.tsv` (per-allele core span + sample).

### `cluster cut`

**Code**: `matdetangler/cluster.py:cut()`.

Reads cached `cores_pairs.tsv` (no realignment) and applies **single-linkage clustering** at `--thresh` (default 0.90): two alleles join a cluster if `sim ≥ thresh`. Clusters are numbered by size (1 = largest); singletons get their own ID.

Outputs:
- `allele_classification.tsv` — `allele, sample, cluster`
- `allele_distance_matrix.tsv` — symmetric `1 - sim` matrix

`cut` is cheap → re-run at any `--thresh` value without redoing `align`.

---

## Determinism + reproducibility

- **`PYTHONHASHSEED=0`** is exported by `MATdetangler-cli` (line 21). Locks Python set/dict iteration order. Two independent runs of the Pcub40 demo on the same install are byte-identical at the level the `installation_run_test.sh` test checks (including the exact graph segments per allele).
- **`test/Pcub40/`** has the canonical reproducibility harness:
  - `installation_run_test.sh` — strict tolerance, exits 1 on any divergence; verifies "same code + env + seed = same output".
  - `analysis_run_test.sh` — biology-focused report, always exits 0; structured per-sample + per-category change log; used when changing the implementation to assess which Pcub40 samples got called differently.
  - `known_results.json` — 32-sample baseline regenerated from the canonical `MATdetangler` conda env (BLAST 2.17.0, MAFFT v7.525).
  - `run_args.json` — canonical run configuration.

---

## Defaults table

| param | default | module | role |
|---|---|---|---|
| `init_nhop` | 3 | per_k_caller | BFS starting hop count |
| `max_nhop` | 8 | per_k_caller | BFS max hop count (reduced from 10) |
| `lo_mult` | 0.2 | per_k_caller | cov filter lower band × D_k |
| `hi_mult` | 2.0 | per_k_caller | cov filter upper band × D_k |
| `divergence_threshold` | 0.01 | per_k_caller | RC-aware HW edit-distance dedup threshold |
| `max_paths` | 1000 | bubble_classifier | cap on simple paths per anchor |
| `max_path_length` | 50 | bubble_classifier | cap on per-path node count |
| `max_bp_since_var` | 5000 | bubble_classifier | bp-aware path enumeration cap |
| `min_allele_bp` | 0 | per_k_caller | hard min-bp floor (0 disables) |
| `locus_padding` | 4000 | per_k_caller | tblastn boundary padding |
| `contig_depth_size_cut` | 5000 | input_process | min seg length for genome_cov median |
| `seeds_mode` | "both" | per_k_caller | BFS seed source: var, flank, or both |
| `cov_filter` | "on" | per_k_caller | cov-OFF main + cov-ON fallback |
| `tblastn_pid` / `tblastn_aa` | 30 / 50 | input_process | HD protein BLAST thresholds |
| `flank_pid` / `flank_minlen` | 85 / 100 | input_process | flank BLAST thresholds |
| `ks` | "k45,k53" | MATdetangler-cli | k subdirs to try |
| `expected_count` | 2 | MATdetangler-cli | 1 = haploid, 2 = dikaryon |
| `SKIP_PICK` | 1 (ON) | MATdetangler-cli | step 4 skipped by default; opt in via `--no-skip-pick` |
| `SKIP_CONSENSUS` | 1 (ON) | MATdetangler-cli | step 6 skipped by default; opt in via `--make-consensus` |
| cluster `--thresh` | 0.90 | cluster.py | single-linkage join threshold |
