# MATdetangler — Methods

Academic-style pseudocode + rationale for each step. Mirror of the production code:
modules referenced here are the actual files under `matdetangler/`.

## Notation

| symbol | meaning |
|---|---|
| `locus.fasta` | single-record DNA reference spanning HD region + both flanks |
| `proteins.fasta` | user-curated protein fasta (HD1, HD2, …); intron-free |
| `contigs.fasta(k)` | SPAdes per-k contigs |
| `gfa(k)` | `assembly_graph_after_simplification.gfa` for that k — segments (`S`) + links (`L`) |
| `D_k` | genome diploid coverage at K (auto-estimated per-k as the median of SPAdes `contigs.fasta` `cov_` headers across contigs ≥ `--contig-depth-size-cut` bp; `--genome-coverage` overrides for all k's) |
| `dp(s)` | per-segment depth tag (`DP:f:`, or `KC:i:` / length) |

---

## 1. Input processing → query construction  (`input_process.py`)

```
inputs:
    locus_fa     : single-record DNA reference
    proteins_fa  : curated protein fasta (HD1, HD2, ...)
    max_flank_len    = 2000          # bp per side, tunable (--max-flank-len)
    envelope_padding = 500           # bp per side, tunable (--envelope-padding)
    min_flank_len    = 100           # warn if trimmed flank is shorter (--min-flank-len)
    tblastn_pid      = 30.0          # HSP filter, fixed (function default)
    min_aa           = 50            # HSP filter, fixed (function default)

# 1.1 sanitize protein IDs — _sanitize_protein_ids() — first-token, illegal-chars→'_',
#      and numeric suffix on collisions (blast requires unique sseqids)
proteins_sanitized = _sanitize_protein_ids(proteins_fa)

# 1.2 tblastn(proteins -> locus), keep HSPs above (pid >= 30, aln >= 50 aa).
#      tblastn_locate() runs the search + the confidence-ordered greedy below.
all_hsps = tblastn(proteins_sanitized -> locus, min_pid=30.0, min_aa=50)

# 1.3 CONFIDENCE-ORDERED greedy with cross-protein non-overlap. HD1 and HD2 are
#      paralogs: each protein query has real strong HSPs in its own locus region AND
#      weaker partial HSPs in the sister's region. Walk all HSPs in confidence order
#      (score = pid * aln_aa, descending); accept each HSP UNLESS it overlaps an
#      already-accepted HSP from a DIFFERENT protein by > 0.5 x shorter-HSP length.
#      Same-protein HSPs never conflict — divergent N/C-terminus HSPs from the SAME
#      protein still extend that gene's span.
accepted = []
for (qid, h) in sorted(all_hsps, key=lambda x: -x.score):
    if not any(k.qid != qid and overlap(h, k) / shorter(h, k) > 0.5 for k in accepted):
        accepted.append((qid, h))

# 1.4 per-gene span = union of its accepted HSPs; annotation (pid, aln_aa, strand)
#      from the BEST accepted HSP of that gene (NOT summed, so aln_aa never exceeds
#      one HSP's contribution — the 643 aa HD1 protein no longer reports aln_aa=720).
for each gene g:
    g.start, g.end = min(h.s), max(h.e) for h in accepted with qid == g
    g.pid, g.aln_aa, g.strand = best accepted HSP of g (by score)

# 1.5 intergenic regions between adjacent gene spans are reported AS-IS — not split
#      into either gene's territory. They live inside the HD envelope by virtue of
#      being between gene starts/ends.
intergenic = [(a.qid, a.end, b.qid, b.start, max(0, b.start - a.end - 1))
              for adjacent (a, b) in sorted_by_start(genes)]

# 1.6 HD envelope = union of all gene spans, padded
env_start = max(1,        min(g.start for g in genes) - envelope_padding)
env_end   = min(|locus|,  max(g.end   for g in genes) + envelope_padding)

# 1.7 TRIM flanks
flankL = locus[max(0, env_start - 1 - max_flank_len) : env_start - 1]
flankR = locus[env_end : min(|locus|, env_end + max_flank_len)]

# 1.8 derive the path-enumerator bp budget
derived_max_locus_len = envelope_size + |flankL| + |flankR|

# 1.9 emit
write variable_proteins.fasta     (sanitized copy)
write variable_nt.fasta            (per-gene NT span on locus, strand-corrected)
write flankL.fasta, flankR.fasta   (trimmed)
write manifest.json                (chrom, chrom_len, envelope_start/end/size,
                                    derived_max_locus_len, variable_genes (per-gene
                                    start/end/strand/pid/aln_aa/n_hsps_accepted),
                                    intergenic_ranges, flankL/flankR locus coords +
                                    len, fasta paths, and the input args
                                    max_flank_len / envelope_padding / min_flank_len)
```

### Rationale

- **Confidence-ordered greedy, NOT naive union.** Earlier versions took the union of all
  HSPs per protein and got bitten twice: (a) cross-paralog leakage — HD1's weak HSPs in
  the HD2 region pushed HD1's reported span into HD2's territory; (b) intra-protein
  HSP overlap — naively summing aln_aa across overlapping same-protein HSPs inflated
  HD1's reported `aln_aa` from its true 643 aa to 720 aa. Walking HSPs in confidence
  order makes the high-quality real hits land their territorial claims first; the
  low-quality paralog-leakage HSPs arrive later and are rejected because the range is
  already claimed by a different protein. `aln_aa` is now reported as the BEST accepted
  HSP's value rather than a sum, so it can never exceed the protein's true length.
- **Intergenic regions are kept as intergenic.** We do NOT split the gap between adjacent
  genes into either's territory. The space between HD1.end and HD2.start (e.g. 405 bp on
  Pcub) is real intergenic sequence — it stays inside the envelope by virtue of being
  bracketed by gene spans, and is reported explicitly as `intergenic_ranges` in the
  manifest.
- **Envelope padding (default 500 bp).** Even after the cleaned span, the very ends of a
  CDS (divergent N/C-terminus) can fall below the pid cutoff and never appear as HSPs.
  500 bp padding each side compensates. Because the cleaning step rejects cross-protein
  overlap BEFORE padding, the 500 bp padding cannot cause HD1 and HD2 envelopes to
  collide.
- **Trimmed flanks.** A long flank is mostly wasted: 2000 bp on each side is enough for
  reliable anchoring during BFS and blastn. Smaller `max_flank_len` ⇒ tighter
  `derived_max_locus_len` ⇒ tighter path search.
- **Manifest.** All downstream modules read derived values (`derived_max_locus_len`,
  envelope coords, flank coords, intergenic ranges) from `<queries_dir>/manifest.json`
  so the pipeline doesn't carry duplicate CLI state.

---

## 2. Genome-coverage estimate  (MATdetangler wrapper, PER k inside the per-k loop)

Run inside the per-k loop just before the per-K caller. Fires only when the user did
NOT pass `--genome-coverage` explicitly (the user value, if supplied, overrides all k's).
The estimate `D_k` is consumed by `matdetangler.run_per_k` (depth-filter band +
per-allele depth) and by `matdetangler.pick_k` (diploid-balance tie-break).

```
# 2.05  Median of SPAdes contigs.fasta `cov_` field, restricted to contigs >= L bp
#       where L = --contig-depth-size-cut (default 5000). Reads-free.
contigs_fa = "$SPADES_DIR/k${K}/contigs.fasta"
covs = [float(cov_) for h in headers(contigs_fa)
        for (length, cov_) in [parse('length_(\d+)_cov_(\d+(?:\.\d+)?)', h)]
        if length >= --contig-depth-size-cut]
D_k = median(covs) if covs else 0.0
write_to("genome_cov_spades_k${K}.txt", D_k)
log("[k${K}] genome_cov = ${D_k} (median cov_ across contigs >= ${L} bp)")
```

### Rationale

- **Reads-free.** SPAdes already computed contig coverage during assembly; there's no
  reason to re-map reads with bowtie2 just to recover it. The pipeline is reads-free
  through step 5 as a result; step 6 (`map_consensus`) is the first step that needs reads.
- **Same units as GFA segment depths.** Both SPAdes contigs.fasta `cov_` and GFA `DP:f:`
  segment depths are in bp-equivalent units (k-mer count divided by `(length - k + 1)`,
  scaled by SPAdes to per-base coverage). So `D_k` is directly comparable to
  per-segment depth in the per-K caller's depth filter (`[lo_mult × D_k, hi_mult × D_k]`).
- **Median, not mean.** Contig depths have a heavy right tail from collapsed-repeat
  short contigs (`cov_` up to ~10⁷ on 46 bp tips in real fungal data). Median is robust
  to those.
- **Length filter at 5 kb.** Stays well within the noise-free regime; smaller contigs
  contribute too much variance from short-tip assembly artifacts. Tunable via
  `--contig-depth-size-cut` if your assembly is unusually fragmented.
- **PER k, not once per sample.** SPAdes coverage values are in k-mer-multiplicity ×
  `R/(R-k+1)` units, so they DIFFER across k's for the same sample. The k=33 estimate is
  ~25% higher than k=53 for R=150 bp reads. Computing per-k means the per-K caller's
  filter band is automatically scaled to that k's depth distribution.
- **Cached as `genome_cov_spades_k<K>.txt`.** `--continue` reuses the cached value.
- **What if no contig is long enough?** `D_k = 0.0`. The caller treats this as
  "disable coverage-based filter" — BFS still runs, the depth filter just never fires.
  The picker's diploid-balance tie-break becomes a no-op.

---

## 3. Per-K allele caller  (`matdetangler.run_per_k`)

For each K in `--ks`, a single self-contained driver pulls the entire per-(sample, K)
result from `gfa(K)` + queries. Replaces the legacy `graph_path_search.py` +
`pick_alleles.py` chain.

The driver does three things in sequence: BLAST the GFA segments against the curated
queries (§3.1), aggregate the hits into a per-segment label TSV (§3.2 labeler), then
run the classifier + emission chain (§3.3 onward).

### 3.1 Pipeline (driver)

```
for each K in --ks:

    # 3.1.1 Extract S-line segments from gfa(K) → segs.fa (one record per segment)
    write_each_S_line_as_fasta(gfa(K)) → segs.fa
    makeblastdb -in segs.fa -dbtype nucl -out segs_db

    # 3.1.2 Full outfmt-6 BLAST against GFA segments (12 cols incl. sstart/send)
    blastn  -query flankL  -db segs_db -dust no -evalue 1e-5 -outfmt 6 → flankL_blastn.tsv
    blastn  -query flankR  -db segs_db -dust no -evalue 1e-5 -outfmt 6 → flankR_blastn.tsv
    tblastn -query HDs     -db segs_db          -evalue 1e-5 -outfmt 6 → HD_tblastn.tsv

    # 3.1.3 Labeler — aggregate BLAST TSVs into per-segment label rows (§3.2)
    seg_label_hits.tsv = emit_seg_label_hits(gfa(K),
                          [(flankL_blastn.tsv, "flank"),
                           (flankR_blastn.tsv, "flank"),
                           (HD_tblastn.tsv,    "var")],
                          min_alnlen=100)

    # 3.1.4 Caller — full algorithm in §3.3 onward
    res = find_alleles(seg_label_hits.tsv, gfa(K),
                       genome_cov=D_k,
                       init_nhop=3, max_nhop=10,
                       lo_mult=0.25, hi_mult=2.0,
                       var_proteins_ref=HDs.fasta,
                       expected_var_tags={HD1, HD2, ...},
                       locus_padding=1500,
                       divergence_threshold=0.01)

    write_fasta(res.alleles, <outdir>/<K>/alleles.fasta)
    write_result_tsv(res,    <outdir>/<K>/result.tsv)
```

Per-K outputs under `<outdir>/<K>/`:

| file | what |
|---|---|
| `flankL_blastn.tsv`, `flankR_blastn.tsv`, `HD_tblastn.tsv` | full 12-col outfmt-6 BLAST hits |
| `seg_label_hits.tsv` | labeler output (one row per (seg, tag, region); §3.2) |
| `alleles.fasta` | emitted allele/chimera records |
| `result.tsv` | one row, 21 columns (§3.11) |

**Rationale.** One BLAST DB per (sample, K) — the new caller blasts directly against
GFA segments; the legacy anchor-search step is no longer in the default flow (its
outputs were consumed by the retired picker chain). Full outfmt-6 is required because
the labeler reads `sstart`/`send` to know where each hit lands on its segment.

### 3.2 Labeler — `seg_label_hits.tsv` format

`seg_label_hits.tsv` schema — one row per `(seg_id, sub-region, tag)`:

```
seg_id   seg_length   tag      kind    start   end    strand
NODE_3   12450        flankL   flank   0       4310   +
NODE_7   8200         HD1      var     0       1020   +
NODE_7   8200         HD2      var     1450    2300   +
NODE_9   3100         flankR   flank   200     1800   +
```

The unpacked format means a single segment with two var hits has two rows; the
segment processor (§3.3) sorts by `start` and partitions the segment at midpoints.
Multi-tag co-located hits (e.g. an HD1/HD2 fusion segment) are joined with `+` and
emitted in **one** row whose `kind` is `var`.

The labeler runs three merge stages before emission:

- **Stage A — intra-tag merge**: overlapping/adjacent hits with the same tag on the
  same segment collapse to their union.
- **Stage B — sweepline same-kind merge**: among hits of the same `kind` (flank vs
  var), contiguous regions are emitted as single intervals; a region covered by
  multiple tags gets a `+`-joined label.
- **Stage C — var-priority clip**: flank regions are clipped against var regions
  (var wins). A segment carrying both `HD1` and `flankR` keeps the HD1 span pure-var
  and only the residual non-overlap as flank.

### 3.3 Segment preprocessing — P1 directional split + P2 bubble BFS

The classifier consumes a graph with three node categories: **pure-flank** (label
has flankL or flankR, no var-gene tag), **pure-var** (var-gene tag, no flank tag),
and **unlabeled connector** (neither). The input may also contain **composite**
nodes (a single GFA segment with BOTH flank and var tags, e.g. `HD1+flankL`). The
preprocessing removes composites and identifies the bubble:

**P1. Directional split — unique-label rule (one sub-segment per unique-tag run)**

Each multi-hit segment becomes a chain of sub-segments, one per
**unique-tag run** (a maximal consecutive group of hits sharing the same
tag), ordered along the parent segment's stored strand:

```
n_runs = 0   →  1 piece (parent segment kept whole, no label)
n_runs = 1   →  1 piece (parent kept whole, labeled with the run's tag)
n_runs ≥ 2   →  one sub-segment per run, ordered by run start coordinate.
                 Boundaries between adjacent sub-segments are placed at
                 the midpoint of the inter-run gap (between the previous
                 run's max end and the next run's min start). Each
                 sub-segment carries exactly the run's tag.
```

A "run" merges consecutive same-tag hits into one labeled region. So a
segment with hits `[flankL, HD1, HD1, HD2, HD2, flankR]` (6 raw hits)
becomes 4 sub-segments `[flankL]#1 / [HD1]#2 / [HD2]#3 / [flankR]#4`,
**not** 6 — there are no unlabeled "linker" sub-segments between adjacent
fragments of the same gene. (Runs only merge consecutive hits with the
same tag; a same-tag pair separated by a different-tag hit stays as two
distinct sub-segments, preserving tandem-arrangement signal.)

Sub-node IDs are coordinate-free: `{parent_seg}#1`, `{parent_seg}#2`, …
(1-based, ordered by run start). The provenance dict
`{sub_id: (parent_seg, start, end, strand)}` carries the actual
coordinates for sequence materialization downstream — IDs stay clean for
display in `bubble.txt` / `bubble.gfa` / `bubble.png`.

Example (the IBUG-14475 k45 closed bubble, parent 18741825 with 6 raw hits):

```
input:  ─── flankL + 2×HD1 + 2×HD2 + flankR composite (6148 bp) ───
              ↑ flankL hit at [0-421]
              ↑ HD1    hits at [993-1986], [2404-3067]   (merged into one run)
              ↑ HD2    hits at [3615-4284], [4498-5518]  (merged into one run)
              ↑ flankR hit at [6024-6148]
P1 out: ─── [flankL]#1 ─── [HD1]#2 ─── [HD2]#3 ─── [flankR]#4 ───
            (4 pieces, one per unique-label run; boundaries at midpoints)
```

After P1, every node carries exactly one label (one tag per sub-segment).
No composites remain. Each neighbor edge of the original segment attaches
to the sub-segment whose run covers its connection end (GFA L-line
orientation gives this).

History note: earlier rules tried (a) `n_label`-based 3-piece-max
(collapsed same-tag hits and absorbed middle labels into one piece) and
(b) one piece per raw BLAST hit (the all-split rule). The current
unique-label-run rule is the middle ground: each gene/flank contributes
exactly one labeled sub-segment regardless of how many BLAST fragments it
produced, but multiple genes stay separate. This keeps bubble visualization
clean (no unlabeled linker pieces between adjacent same-tag fragments)
while still letting the BFS / classifier distinguish multi-gene composites.

**P2. Bubble BFS** — starting from every pure-var node, BFS through unlabeled-only
neighbors. The set of nodes reached forms the **bubble**; everything else is
**flank-region**:

```
Bubble       = pure-var nodes ∪ unlabeled connectors reachable from var
                                  via unlabeled-only paths
Flank-region = nodes - Bubble
```

P2 is one BFS — no iteration over rounds, no parameters. A long unlabeled spacer
between var and a flank gets pulled INTO the bubble (it becomes part of the bubble
interior); a noise chain dangling off a pure flank with no var-reachable path stays
in flank-region. Flank-region naturally contains all pure flanks plus any "shoulder"
connectors that only reach flank, not var.

After preprocessing the bubble is a connected subgraph that may be: a linear chain
(synthetic clean closed_bubble), a Y-fork at a Uvar joint (AG10), a cycle (AG17), a
multi-fork hub (where one var node branches to 3+), or multiple disjoint components
(when var content lives in separate locations). The classifier does NOT make
structural assumptions about the bubble's internal shape. It just counts paths
through it.

### 3.4 Anchors (universal-leaf rule)

An **anchor** is a bubble node that meets the outside world, defined purely by
graph topology:

```
anchor = bubble node n such that
            n has ≥1 neighbor OUTSIDE the bubble  (any label)
        OR  n has bubble-degree ≤ 1               (dead-end leaf)
```

This generalizes the prior "flank-adjacent only" rule: flank-adjacent nodes still
qualify (they have a flank-labeled neighbor outside the bubble), and Uvar dead-end
leaves also qualify (they have ≤ 1 neighbor inside the bubble — bubble truncated by
BFS hop limit, fragmented assembly, or chromosome end).

L/R distinction is decorative — derived from flank labels when present — but does
NOT gate arm enumeration:

```
flank_L_adj = anchors with a flankL-labeled outside neighbor
flank_R_adj = anchors with a flankR-labeled outside neighbor
```

Anchors that are neither L- nor R-decorated are "bare leaves" — they get used as
path endpoints with a different downstream label (dangling rather than closed; §3.5).

### 3.5 Arms

An **arm** is a var-bearing simple path between two distinct anchors, classified by
its endpoint flank-status:

```
closed   = both anchors are flank-decorated, one on L side and
           one on R side (a proper L→R locus walk)
dangling = exactly one endpoint is flank-decorated (the other is a
           bare leaf — one side anchored, other side dangling)
ignored  = same-side (L→L, R→R) or leaf-only (no flank-adjacent
           endpoint) — not biologically meaningful as arms
```

Single-node bridge: ONE var node that is adjacent to BOTH a flankL node AND a flankR
node (a composite segment spanning the whole locus) counts as a closed arm of length 1.

"Simple path" has its standard graph-theory meaning: no node revisit. This guarantees
finite enumeration on cyclic bubbles.

**No further dedup at the classifier.** If two arms route through shared Uvar hubs in
substantially different combinations, those ARE distinct arms — that's structural
complexity and it correctly drives the verdict toward 'complexed'. (Sequence-level
dedup happens later in `_emit_result`; §3.8.)

### 3.6 Verdict

Precedence rule applied first:

- Var nodes split across 2+ connected components of the **full graph** ⇒ `separate`

Then by arm count:

| arms | configuration | verdict |
|------|---------------|---------|
| 0    | no var-bearing route                    | complexed |
| 1    | one arm                                  | single |
| 2    | both closed                              | closed_bubble |
| 2    | mixed (≥1 dangling)                      | open_bubble |
| ≥3   | three or more                            | complexed |

**Why the spec is deductive.** Every step is a single graph operation with no
parameters: P1 rewrite (split nodes by label coordinate), P2 BFS (var through
unlabeled), anchors (adjacency / degree test), arms (simple-path enumeration with
endpoint-flank-status filter), verdict (count → table). No `MAX_LINKER_PADDING`
parameter, no validity horizon, no clean-chain check, no var-series/flank-series
special handling, no minimal-var-set dedup, no arm-membership consolidation. The
composite mess is handled upstream by P1; the long-spacer mess is handled by P2;
the rest is a single counting step.

### 3.7 `find_alleles` BFS loop

`find_alleles(seg_label_hits_tsv, gfa_path, genome_cov, init_nhop=3, max_nhop=10, var_proteins_ref=…, expected_var_tags=…, locus_padding=4000, lo_mult=0.2, hi_mult=2.0, divergence_threshold=0.01, queries_dir=…, out_candidate_fa=…, seeds_mode="both", cov_filter=True, max_paths=50, max_path_length=15, min_allele_bp=3000)`

The orchestrator runs the BFS at increasing hop counts, collects per-network
candidates from every iteration, optionally short-circuits when a fully
complete candidate appears, then picks across the pool. Loop axes
**seeds_mode** and **cov_filter** are caller-time configuration (not loop
variants), so the loop is just over `nhop`.

```python
# Seeds — set ONCE from seeds_mode
seeds = var_segs              if seeds_mode == "var"
      | flank_segs            if seeds_mode == "flank"
      | (var_segs|flank_segs) if seeds_mode == "both"          # default

candidates = []
short_circuit = False

for nhop in init_nhop..max_nhop:                  # default 3..10
    res = _try_one_pass(seeds, nhop, cov_filter,  # cov_filter = True by default
                        max_paths, max_path_length)
    log(nhop, |nhood|, n_var, class, n_arms,
        bfs_limits_hit)                            # logs cap hits when fired

    if res.class == "no_var" or res.n_var == 0:
        continue

    pools = _emit_result(res, ..., return_pools=True)
    for i, pool in enumerate(pools):
        if not pool.alleles: continue
        c = {
            iter_id, net_in_iter=(i+1 if N>1 else 0),
            verdict=pool.sub_verdict,             # NEVER 'separate' at pool level
            complete_var,                          # tri-state 0/1/2 — per-allele MIN
            complete_locus,                        # tri-state 0/1/2 — per-allele MIN
            alleles, allele_cov, basepair, diploid_dist,
            ...
        }
        candidates.append(c)

        # Acceptance gate (α — hard short-circuit): tri-state == 2 required
        if (c.verdict == "closed_bubble"
                and c.n_dedup >= 2
                and c.complete_locus >= 2
                and c.complete_var >= 2):
            short_circuit = True

    if short_circuit: break

write_candidate_fasta(candidates, out_candidate_fa)

# === Unified 4-tier rank (used at all 4 picker/dedup sites) =======
best = min(candidates, key=lambda c: (
    bubble_priority(c.verdict, c.n_dedup),     # 0. K-picker priority order
    -c.complete_locus,                          # 1. tri-state DESC
    -c.complete_var,                            # 2. tri-state DESC
    c.diploid_dist,                             # 3. |allele_cov/D_k − ½|  ASC
))

# === Cross-network sibling dedup (same iteration only) ============
siblings = [c for c in candidates if c.iter_id == best.iter_id]
final = dedup_ranked_cross_network(siblings, divergence_threshold)
                                              # same 4-tier rank, RC-aware edlib HW

# === Sample verdict ================================================
surviving_nets = {c.net for c in final}
sample_verdict = ("separate" if len(surviving_nets) >= 2
                  else siblings[0].verdict)
```

#### Loop axes that were dropped

Earlier the loop had `phase × nhop × cov_filter` = 24 iterations. Empirical
results on a 148-sample whitelist showed:

| Loop axis | Wins on which samples | Decision |
|---|---|---|
| Phase 2 (`flank_fallback`) | 0 / 144 picks | dropped — moved to `seeds_mode` config |
| `cov_filter=False` pass | 1 / 292 per-k results, never wins K-pick | dropped — moved to `cov_filter` config |

Iteration count is now `max_nhop − init_nhop + 1` = up to 8 (default 3..10).

#### Key invariants

- **Acceptance (α) is a hard short-circuit**, not a fallback. Fires the
  moment ANY candidate is `closed_bubble`, `n_dedup ≥ 2`, AND both tri-state
  completeness scores at level 2 (all expected var tags AND both flanks
  present in EVERY emitted allele of the candidate's pool — the per-allele
  MIN aggregation prevents the "two single-HD fragments union to look
  complete" pathology).
- **No candidate ever carries `verdict="separate"`.** When the classifier
  returns `class="separate"`, each network's sub-classification
  (`closed_bubble`/`open_bubble`/`single`/`complexed`) becomes its
  candidate's verdict. The `separate` label is applied AT THE SAMPLE LEVEL
  only — after dedup, if ≥ 2 distinct networks survived.
- **Cross-iteration mixing is disallowed**: the picker chooses ONE
  candidate, and only its same-`iter_id` siblings join it for the
  cross-network dedup.
- **`candidate_allele.fasta`** captures every emitted sequence across the
  loop, with headers tagging `(iter_id, network, cov, verdict)` — purely
  for audit / debugging.

#### `_try_one_pass` (single iteration)

Each `_try_one_pass` BFS-expands seeds by `nhop` hops, optionally applies a
depth filter `[lo_mult × D_k, hi_mult × D_k]` (default `[0.2, 2.0]`),
restricts edges, runs P1 directional split (§3.3), runs the classifier
(§3.4–3.6), and stashes `_provenance` + `_var_per_node` + `_bfs_limits` on
the result. Self-loop GFA edges (where source == target) are skipped.

The classifier's `class="separate"` recursion (per-network sub-classification)
also happens here — sub_results are attached to `res` and consumed by
`_emit_result(return_pools=True)` to produce per-network pools.

#### BFS path enumeration — shortest-first, with hard caps

The classifier's arm enumeration uses **BFS by path length** (not DFS). At
each iteration of the outer BFS, all simple paths of length L are emitted
before any path of length L+1 starts. So when `max_paths` fires, the
SHORTEST `max_paths` paths survive — the real diploid pair (typically 2–6
nodes) always makes it into the cap even on fork-explosion samples where
DFS would plunge deep on one branch and miss the other arms.

Two hard caps:
- `max_paths = 50` — total number of simple paths per starting anchor
- `max_path_length = 15` — drop any path whose node count exceeds this

Each iteration's log line shows `⚠ limits: max_paths_hit=N max_path_length_hit=M`
when either cap fired (suppressed when both = 0).

#### Completeness — graph-level flank check

Flank presence per allele is detected by a single blastn pass over the
candidate FASTA (`queries/flankL.fasta` + `queries/flankR.fasta`, pid ≥ 85%,
aln ≥ 100 bp). Per-allele `cl_level = int(has_flankL) + int(has_flankR)` ∈
{0, 1, 2}. The pool-level `complete_locus` is the MIN across surviving
alleles — see §3.8.

### 3.8 Emission rule (trim → dedup → emit)

`_emit_result(return_pools=False)` produces ONE merged sample-level emission
(legacy path, used when the orchestrator finalizes a single picked iteration).
`_emit_result(return_pools=True)` returns the raw per-pool list (new path, used
by `find_alleles` to collect candidates across iterations — §3.7).

Per pool (= per network when classifier is "separate", or single global pool
otherwise), the pipeline is:

```
candidate paths from classifier
    │  closed_arms + dangling_arms + var_components
    ▼
[node-level path-trim]              ← always (uses var_per_node, no BLAST)
    │  ordered paths    → subpath [first_var, last_var]
    │  var_components   → only var-labeled nodes
    │
    │  (pre-trim path kept SEPARATELY — used by the emission step below
    │   and by graph-level flank-presence check — §3.7 "Completeness")
    ▼
build var-trimmed sequences via provenance (orig_seg, start, end, strand)
    │  → used for DEDUP RANKING + HD-only divergence compare ONLY
    │
    ▼
[tblastn locus-trim] for dedup compare   ← when var_proteins_ref provided
    │  one tblastn pass: var_proteins → all candidates as multi-fasta
    │  per candidate: trim to [min(sstart)−padding, max(send)+padding]
    │  drop candidates with no HD hits
    │  collect found_var_tags per surviving candidate
    ▼
[completeness-first dedup]              ← runs on var-trimmed → HD-only slice
    │  RANK candidates by (smaller = better):
    │    (−cl_level, −cv_level, diploid_dist)
    │  where cl_level/cv_level are PER-CANDIDATE tri-states (0/1/2):
    │    cv_level = how many expected_var_tags found in THIS candidate
    │               (0=none, 1=some, 2=all)
    │    cl_level = how many of {flankL, flankR} present in THIS candidate
    │  WALK in rank order:
    │    keep s iff is_divergent(s, k, threshold) for every already-kept k
    │    where is_divergent uses RC-aware edlib HW edit-distance:
    │      ed_fwd = edlib.align(q, t,         mode=HW, task=distance)
    │      ed_rc  = edlib.align(q, RC(t),     mode=HW, task=distance)
    │      identity = 1 − min(ed_fwd, ed_rc) / |q|
    │    threshold default = 1% (collapse if identity ≥ 99%)
    ▼
[EMIT-time re-trim: non-joint slice + tblastn-trim]
    │  For each dedup-surviving candidate (separately from the dedup compare):
    │    1. "joints" = nodes that appear in ≥ 2 of the pool's pre-trim paths
    │       (the bubble's anchor/fork nodes shared by sibling arms; for a
    │        single-arm pool joints = ∅).
    │    2. Slice the pre-trim path from FIRST non-joint node to LAST
    │       non-joint node (end-joints stripped; internal joints preserved).
    │    3. Build the emitted sequence from the sliced path via provenance.
    │    4. Apply tblastn locus-trim (same padding as the dedup-compare trim).
    │  Rationale: dedup decides WHICH alleles survive (on the var-window
    │  content where the biology lives); emission decides WHAT we write to
    │  primary_alleles.fasta (the arm-unique content of each surviving
    │  allele, with flank-adjacent joint nodes excluded but the
    │  HD-flanking padding preserved). Falls back to the raw non-joint
    │  sequence if tblastn drops the survivor.
    ▼
emit: allele1 / allele1+allele2 / chimera1..N (by post-dedup count)
    │
    ▼
[per-allele depth + tri-state completeness + extend_bounds]
    allele_cov = length-weighted DP:f: mean per allele path
    surviving_tags = union of found_var_tags across survivors (info-only)
    per-allele cv_level = 0/1/2  (this allele's tblastn vs expected_var_tags)
    per-allele cl_level = 0/1/2  (this allele's blastn vs flankL+flankR)
    pool-level complete_var  = MIN(per-allele cv_level across surviving)
    pool-level complete_locus = MIN(per-allele cl_level across surviving)
    diploid_dist = |Σ(allele_cov)/D_k − 1| ... actually we use
                   |mean(allele_cov)/D_k − ½| per pool
    extend_bounds = per allele, (innermost_flankL_node,
                    innermost_flankR_node)
```

#### Why per-allele MIN aggregation (not union)

A naive union — "pool is complete iff the surviving alleles' tag sets UNION
to expected_var_tags" — has a pathology: two single-HD fragments (one allele
with HD1 only, another with HD2 only) union to look fully complete, even
when neither allele individually covers the locus. This was the KYH069 n1
case: two 308 bp single-HD fragments in a disjoint network scored
`complete_var = 2` under union aggregation.

Per-allele MIN aggregation requires EVERY surviving allele to score at the
same level — a pool scores `complete_var = 2` only when every emitted
allele individually has all expected var tags. KYH069 n1's per-allele
levels become {1, 1} → pool MIN = 1, demoted below true diploid networks.

#### Pool-level outputs the orchestrator consumes

Each pool returns (used both by the legacy single-merge path and by the new
per-iteration candidate collection):

```python
{
  alleles, allele_cov, extend_bounds, allele_segments,
  n_raw, n_dedup, found_tags_surviving,
  verdict,                                  # per-pool sub-verdict (NOT "separate")
  complete_var, complete_locus,
  basepair, diploid_dist,                   # for cross-iteration ranking
  raw_dedup,                                # surviving (name, seq) pre-merge
}
```

| post-dedup `n` | emitted records           | naming inside pool |
|----------------|---------------------------|--------------------|
| 1              | `allele1`                 | (prefix + `allele1`) |
| 2              | `allele1`, `allele2`      | (prefix + `allele{i}`) |
| ≥ 3            | `chimera1`, `chimera2`, … | (prefix + `chimera{i}`) |

(prefix is `""` for non-separate, `n{i}_` for separate's i-th network — so
output is e.g. `n2_allele1`, `n2_allele2`, …)

A `closed_bubble` whose two arms collapse to one (dedup→1) is reported as
`verdict=closed_bubble` with one `allele1`, NOT downgraded to `single`.
`n_after_dedup` tells the consumer how many distinct sequences came out.

#### Why completeness-first ranking, not length-first

The old rule kept the LONGEST representative of each dedup equivalence class.
Problem: in transitive clusters where A↔B and B↔C are within 5% but A↔C is
above 5%, length-first kept A and C (both endpoints), giving 2 alleles.
Completeness-first picks B (the central more-complete representative), whose
5% radius now covers BOTH A and C → 1 allele.

Empirical: switching from length-first to completeness-first reduced total
allele count from 407 → 395 across the 148-sample whitelist, with the
biggest reductions in the 5–7 kb and ≥ 7 kb buckets (long chimeras that
were keeping their length-tied tail variants). The 4–5 kb diploid bucket
was untouched.

#### RC-aware divergence — why both strands

`_identity()` aligns the shorter sequence (as query) against BOTH the target
and `RC(target)`, taking the better edit-distance. This is needed because
graph walks through the same bubble can be emitted in either orientation
depending on which anchor BFS started from. Without RC-awareness, two arms
that are reverse-complements of the same allele would survive dedup as
distinct (forward HW edit-distance between them would be near-zero on the
shared part, but the orientation flip looks like total dissimilarity to
infix-anchored alignment if the target is shorter than the query).

#### HD-only divergence — comparing the var region, not the flanks

The dedup comparison uses only the **HD-only span** of each candidate, not
the full padded sequence. Two alleles whose HD cores diverge 30% but whose
flanks are 100% identical would look ~85% identical at the whole-allele
level — the flanks dominate the length. Restricting the comparison to the
var region surfaces the biologically-meaningful divergence signal.

Mechanics:
- `_tblastn_trim_each(seqs, ref, padding)` returns `(padded_seqs,
  found_tags, hd_span)`. The third value records the HD-only coords
  `(lo, hi)` *within each padded slice*, so the caller can do `seq[lo:hi]`
  to get the un-padded HD region.
- `_dedup_named_ranked(seqs, threshold, key_fn, compare_seqs=…)` takes an
  optional `compare_seqs={name: hd_only_seq}` dict. Sort order still uses
  the full rank, but `is_divergent()` is computed on the HD-only slice.
- Emitted alleles remain the full padded ~9 kb sequences — only the
  divergence metric narrows to the var region.

#### Hard minimum-length floor (`min_allele_bp`, default 3000)

Three-stage filter that drops sub-3kb fragment candidates:
1. **Pre-locus-trim** — applied to raw walk content before tblastn, so
   fragment-network walks (e.g. KYH069 n1's 308 bp single-HD segments)
   never reach the tblastn step.
2. **Post-locus-trim** — second pass after the trim window narrows
   sequences, catches any survivor that fell below the floor.
3. **Cross-network sibling dedup** — same floor applied when gathering
   per-iteration siblings, catches fragments that came in via a sibling
   network's emission.

Clean closed_bubble samples (real 4.5–5 kb diploid pairs) are unaffected;
fragment-only "n1 308 bp" pools are dropped entirely.

The current `_blastn_flank_presence()` helper still exists for fallback /
diagnostic use but is NOT consulted during normal emission — the
graph-level path check (§3.7) replaces it. Calling code that wants
sequence-level confirmation can still invoke it.

#### Why two layers of trim

Node-level path-trim is cheap (no BLAST, uses graph-level `var_per_node`
membership) and strips most of the conserved flank. The tblastn locus-trim
then catches cases where the labeler missed a var hit or the node-trim
left extra context; it is the content-based safety net.

### 3.9 Arm-sequence reconstruction

Each arm is a path of post-P1 node IDs. The seg-processor supplies
`provenance[node] = (orig_seg, start, end, strand)` for each sub-segment. The
reconstructor walks the path, slices the stored-strand DNA at `[start:end]` per
node, applies reverse-complement where `strand == "-"`, and concatenates.
`arm_sequences.py` exposes a joint-skipping variant (for diagnostics);
`_emit_result` uses full-path reconstruction since its trim chain owns boundary logic.

### 3.10 Module layout

The classifier + driver live under `matdetangler/graph_classifier/` and
`matdetangler/`. Each module has a single responsibility:

| module | role | key function |
|---|---|---|
| `graph_classifier/labeler.py` | Aggregator: BLAST TSVs + GFA → `seg_label_hits.tsv` (unpacked, one row per tag) | `emit_seg_label_hits(gfa_path, blast_tsvs, out_tsv)` |
| `graph_classifier/seg_processor.py` | P1 directional split; emits provenance `{sub_id: (orig_seg, start, end, strand)}` | `directional_split(seg_labels, seg_length, edges, edge_endpoints)` |
| `graph_classifier/bubble_bfs.py` | P2 Bubble BFS from var through unlabeled | `bubble_bfs(adj, var_nodes, unlabeled)` |
| `graph_classifier/bubble_classifier.py` | Anchors + arms + verdict over the post-P1 graph | `classify(nodes, edges, label_per_node, var_per_node)` |
| `graph_classifier/arm_sequences.py` | Reconstruct allele DNA per arm; skip joint/fork nodes (shared between arms) | `reconstruct_arms(arms, provenance, gfa_seqs, skip_joints=True)` |
| `graph_classifier/per_k_caller.py` | End-to-end caller: BFS loop + trim → dedup → emit | `find_alleles(seg_label_hits_tsv, gfa_path, genome_cov, …)` |
| `run_per_k.py` | Per-(sample, K) driver: extract segs → BLAST → labeler → find_alleles | `main()` (CLI) |

Test suite: `test/graph_classifier/directional_test.py` covers the classifier with
8 hand-built synthetic cases — composites, long unlabeled spacers, Y-fork bubbles,
dangling arms, multi-arm complexity, disjoint var subgraphs, and the AG10 / AG17 /
AJB36 / SA93 / AG5 shapes. All 8 pass. AG17 real-data run emits 2 divergent alleles
(5217 bp + 5045 bp, verdict `closed_bubble`) end-to-end through
`per_k_caller.find_alleles`.

### 3.11 Output: result.tsv + caller dict

Each `find_alleles()` call produces the **best allele set for one GFA at one k**.
The output dict:

```python
{
    # Classification
    "verdict":         "closed_bubble" | "open_bubble" | "single" |
                       "separate" | "complexed",
    "bubble_type":     alias for `verdict` (clarity in TSV columns),
    "k":               the K label passed in (e.g. "k53"),

    # Allele records (post-trim, post-dedup)
    "alleles":         [("allele1", "ACGT..."), ("allele2", "ACGT..."), ...]
                       OR
                       [("chimera1", "ACGT..."), ("chimera2", "ACGT..."), ...],
    "component_list":  list of emitted names (parallel to alleles),
    "allele_cov":      list[float],  # length-weighted DP:f: mean per allele
    "extend_bounds":   list[(str|None, str|None)],
                       # per allele, (innermost_flankL_node,
                       # innermost_flankR_node). Either side is None if
                       # no flank-labeled neighbor was found on that side.
    "basepair":        sum of emitted allele lengths,

    # Completeness (post-trim tblastn over surviving candidates,
    #               graph-level flank-presence check — §3.7)
    "complete_var":    bool | None,   # all expected_var_tags hit?
    "complete_locus":  bool | None,   # at sample level: aliases complete_var
                                       # (per-candidate level: complete_var AND
                                       #  has_flankL AND has_flankR)
    "locus_coverage":  float | None,  # fraction of expected_var_tags found
    "found_var_tags":  sorted list[str],

    # Diagnostics
    "n_candidates":    int,           # before dedup (post-trim)
    "n_after_dedup":   int,           # final emitted count
    "divergent":       bool,          # n_after_dedup >= 2
    "n_hops_used":     int,           # BFS hop count at acceptance/exhaustion
    "phase":           "main" | "flank_fallback",
    "cov_filter_used": bool,          # was the depth filter on at this iteration

    # Segment provenance
    "segments":        sorted list[str],   # union of orig-seg IDs across alleles
    "segments_labeled": list[(seg_id, "tag1+tag2+...")],
    "allele_segments": list[list[str]],    # PER allele, sub-node IDs along the
                                            # FULL anchor-to-anchor walk: the
                                            # closed_arm path PLUS the
                                            # flank-bearing anchor neighbors on
                                            # either end. ({parent}#N format —
                                            # see §3.3.) This is the "fullwalk"
                                            # view used by bubble.txt / bubble.gfa
                                            # / bubble.png; it is INTENTIONALLY
                                            # longer than the emitted FASTA
                                            # sequence (which is the non-joint
                                            # slice + tblastn-trim — see §3.8).
    "subnode_seqs":    {sub_id: dna_seq},  # materialized sub-node sequences
                                            # for IDs referenced by emitted
                                            # alleles; run_per_k writes to
                                            # <k>/subnode_seqs.fasta
    "genome_cov":      float,         # echoed from input

    "info":            { ...classifier stats (n_arms, n_closed, anchors, etc.) },
}
```

`run_per_k.py` flattens this dict into a 22-column `result.tsv`:

```
sample  k  bubble_type  components  complete_var  complete_locus  locus_coverage
basepair  genome_cov  allele_cov  n_cand  n_dedup  divergent  n_hops_used  phase
cov_filter_used  allele_lens  segments  segments_labeled  found_var_tags
extend_bounds  allele_segments
```

`extend_bounds` is encoded as `L:R;L:R;…` (one `L:R` pair per emitted allele,
`;`-joined). `-` placeholder where no flank was reached on that side.
`allele_segments` is `;`-joined per allele, `,`-joined within an allele.

#### Side files written by `run_per_k.py` per (sample, k):

| file | content |
|---|---|
| `alleles.fasta` | post-dedup picked alleles, FASTA. Post-`min_allele_bp` (default 3000 bp). |
| `longest_alleles.fasta` | length-first RC-aware dedup over the full candidate pool — keeps the LONGEST representative of each edit-distance equivalence class (HD-only divergence, 5%). Wider net than `alleles.fasta`; the bash wrapper unions per-k versions into a sample-level `<sample>/longest_alleles.fasta` with `k{NN}_` prefixed headers. |
| `result.tsv` | the 22-column row above. `complete_var` and `complete_locus` are now tri-state integers (0=none, 1=some, 2=all). |
| `seg_label_hits.tsv` | labeler output (§3.2) |
| `subnode_seqs.fasta` | one record per unique `{parent}#N` sub-node ID, with the materialized sub-region DNA (RC'd if post-P1 strand was `-`). Consumed by `graph_paths.py` to draw `bubble.gfa` / `bubble.png` with coord-free IDs. |
| `candidate_allele.fasta` | forensic record — every emitted sequence across all BFS iterations (post `min_allele_bp` filter). Headers carry `cand{id}_h{nhop}_n{net}_c{cov}_{verdict}_{allele_name}`. Not consumed by downstream stages. |
| `flankL_blastn.tsv`, `flankR_blastn.tsv`, `HD_tblastn.tsv` | raw BLAST inputs to the labeler |

Cross-K integration (picking the best K across the per-k results) lives in
`matdetangler.pick_k`, documented in §4. `per_k_caller`'s contract is **one k in,
one verdict + allele set out**.

### 3.12 Upstream plumbing notes

P1 (§3.3) needs per-end label coordinates on each GFA segment. The labeler
(§3.2) produces these from BLAST output (`sstart`/`send` of each hit on a segment)
and writes them to `seg_label_hits.tsv`.

If per-end coordinates are absent (e.g. when a third-party tool feeds this
classifier without sstart/send), the fallback is to retain the composite node
as-is and use a simpler anchor rule (composite carrying a flank tag is also an
anchor). This is less crisp but still correct for common shapes; the production
classifier should use directional input.

---

## 4. Cross-K consolidation  (`matdetangler.pick_k`)

After the per-K caller has run on every K, `pick_k` reads each sample's per-K result
files and picks the best K.

```
for each sample dir <outdir>/:
    candidates = []
    for each <outdir>/<K>/result.tsv:
        row = parse(result.tsv)
        candidates.append(row)
    if not candidates: emit NO_RESULT row; continue

    candidates.sort(key=_score)
    winner = candidates[0]

    # emit:
    write_legacy_picks_row(winner)              → picks.tsv          (one row per emitted allele)
    write_summary_row(winner)                   → picks_summary.tsv   (one row per sample)
    cp <outdir>/<K_winner>/alleles.fasta        → primary_alleles.fasta
```

### Priority order (preferred K first)

```
1. closed_bubble                                              (the ideal — 2 closed arms)
2. open_bubble                                                (1 closed + 1 dangling)
3. separate            (only if n_dedup == 2)                 (disjoint var components, but only 2 alleles distinct)
4. single                                                     (haploid / homozygous: 1 emitted allele)
5. complexed           (only if n_dedup == 2)                 (diploid behind a complex bubble)
6. separate            (n_dedup != 2 — chimeric emission)
7. complexed           (n_dedup != 2)
8. no_var
9. error / no result
```

### Tie-break within the same priority class (lower is better)

1. `complete_var` (True > False > None)
2. higher `locus_coverage`
3. closer to ½ × `D_k` (diploid signature: `|mean(allele_cov) / D_k − 0.5|`)
4. higher total `basepair`
5. lower `n_candidates`

### Outputs

| file | what |
|---|---|
| `picks.tsv` | legacy 12-column schema, one row per emitted allele (compatible with downstream `graph_paths` + `summary_table`) |
| `picks_summary.tsv` | sample-level new schema, one row per sample — verdict, n_dedup, allele_cov, extend_bounds, k_chosen, all_k_tried |
| `primary_alleles.fasta` | concatenated allele records from the chosen K per sample; headers `<sample>_<k>_<allele_name>` |

`picks.tsv` schema (12 cols):

```
sample  allele  origin  k  type  len  from_contig  segments
cov     n_variable_genes  has_both_flanks  is_degHD
```

`picks_summary.tsv` schema (14 cols):

```
sample  k_chosen  bubble_type  n_dedup  complete_var  complete_locus
locus_coverage  basepair  genome_cov  allele_cov  n_cand
extend_bounds  components  all_k_tried
```

### Rationale

- **Priority by topology + dedup count, not by sequence score.** The new caller's
  classifier already states what the bubble looks like; `n_dedup` says how many
  distinct sequences fell out of the trim+dedup chain. The cross-K winner is the K
  where the bubble is cleanest (closed > open > separate > single > complexed) AND
  the dedup landed on the expected ploidy (2 for diploid).
- **A diploid `complexed` with `n_dedup == 2` still wins over a tangled `separate`.**
  Complexed-with-2-survivors means the bubble was structurally tangled but the
  sequence content collapsed to two — that's a real allele pair behind graph noise.
  We prefer it over a `separate` call with > 2 chimeras.
- **Legacy-compatible picks.tsv schema.** Downstream stages (`graph_paths`,
  `summary_table.py`) read picks.tsv with the column order they always expected;
  the new pipeline preserves that column order so no downstream code needs to change.
  Some legacy columns (`from_contig`, `is_degHD`) are placeholder values.
- **`picks_summary.tsv` carries the new fields.** Anything specific to the per-K
  caller (extend_bounds, n_dedup, k_chosen, all_k_tried) goes here so the legacy
  format isn't bloated.

---

## 5. Annotated allele walks  (`graph_paths.py`)

Each picked allele has a recorded segment walk in `picks.tsv` (col 8, `segments`).
`graph_paths` reads the walk verbatim, labels each segment by content (tblastn proteins
+ blastn flanks), and emits per-sample plots + walks.

Inputs the wrapper supplies: `--primary-alleles`, `--picks-tsv`, `--spades-dir`. (The
legacy `--cand-ann` input is no longer required — the new pipeline doesn't write
`bubble_alleles.ann.tsv`.)

### Performance: labels only the segments in the recorded paths

Earlier versions ran `label_nodes` on **every segment of the whole GFA** (100k–
300k segs per call) for each allele in `primary_alleles.fasta`. With 24 candidate
alleles under `--skip-pick`, that was 24× redundant whole-GFA labeling — step 5
could take 30–60 minutes per sample.

Now `graph_paths.run` (via the `_per_allele_gfas_and_paths` helper + a per-GFA cache):
- **per-GFA cache**: parses each unique GFA only once, even when multiple alleles
  share it (the `gfa_cache` dict inside `run`).
- **recorded-path subset**: the union of all alleles' recorded segments (from
  picks.tsv) is the only set of segments labeled — ~tens of segments instead of
  100k–300k. blast queries run on tiny temp DBs.

Combined effect: step 5 now completes in **~0.5–2 seconds** (zero new blast work
when the cache is warm; all hits are USE).

| output | what |
|---|---|
| `bubble.txt` | `allele1: <walk>` / `allele2: <walk>` (or one row per candidate under `--skip-pick`) |
| `bubble.gfa` | sub-GFA of just the walks' segments + L-links (loads in Bandage) |
| `bubble.dot` | Graphviz DOT |
| `bubble.tsv` | edge list |
| `bubble.png` | matplotlib: one row per allele. **Node coloring**: yellow = carries any variable gene (e.g. HD1/HD2; flank tag, if any, ignored for color); blue = flank-only (label is purely flankL/flankR, no variable gene); white = pure-number / unlabeled. **Dashed gray cross-arm lines**: (a) one per GFA segment ID shared between the two arms — same node = definite homology; (b) fallback only when no flank-only segment ID is shared for a given flank type — one extra line connects the outermost flank-only node of each arm (leftmost flankL = arm entry, rightmost flankR = arm exit). Capped at 4 rows displayed (extras dropped from the picture; the title notes `[showing 4 of N]`). |

When both alleles came from the same K, node IDs are bare. When they came from
different K's (cross-K pair), node IDs are prefixed `K33:` / `K55:` etc. so the user
can tell at a glance which graph each allele lives in.

**Row orientation — canonical to the protein query order.** Each row's
var-tag sequence (e.g. `[HD2, HD1]` or `[HD1, HD2]`) is compared to the
order tags appear in `queries_dir/variable_proteins.fasta` (or
`Suilu4_HDs.fasta`). Rows whose order is the exact reverse of the
reference order get walked backwards before drawing, so HD2/HD1 ends align
consistently L→R across rows AND across samples. The title shows
`{sample}: {bubble_type} (N alleles)`.

---

## 6. Mapping + read-derived consensus  (`map_consensus.sh`)

```
bowtie2-build idx ← primary_alleles.fasta
bowtie2 --end-to-end --very-sensitive -k 1 --no-unal -1 R1 -2 R2 | samtools sort → reads.sorted.bam
samtools index → reads.sorted.bam.bai
# one samtools consensus call per allele record, records rewritten to allele_1, allele_2, ...
# and concatenated into one multi-record FASTA:
samtools consensus -r <allele_i> -f fasta -o tmp reads.sorted.bam → consensus_alleles.fasta

# coverage_core.py: per-allele depth restricted to HD-core
hd_iv(a)          = positions in any variable-gene tblastn hit on a
                    (tblastn variable_proteins → allele, pid >= 25, aa >= 50)
core_meandepth(a) = mean depth over hd_iv(a)
```

Outputs: `reads.sorted.bam`, `reads.sorted.bam.bai`, `consensus_alleles.fasta`,
`coverage.tsv`.

`coverage.tsv` columns: `whole_meandepth` (raw) **and** `core_meandepth` (HD-core only).
Summary reports the core depth — it's what tells you whether a real dikaryon allele is
at ≈ D/2.

---

## 7. Consensus QC  (`pairwise_identity.py` re-run + `consensus_qc.py`)

The pick-level numbers from steps 3–4 are based on the graph-derived allele sequences.
The read-derived consensus can drift (low-coverage stretches become N's, etc.), so we
**re-check** divergence and completeness on the consensus. The wrapper runs three
sub-steps:

- **7a — pairwise identity on PICKS.** `pairwise_identity` on
  `primary_alleles.fasta`. HD-core MAFFT protocol via `_align_cores`: trim each
  candidate + locus reference to HD-core ± pad, ONE MAFFT alignment, score identity
  only over HD-core columns projected from the reference's ungapped HD-protein span.
  Writes `identity.tsv`.
- **7b — pairwise identity on CONSENSUS.** Same `pairwise_identity` module
  with the same HD-core MAFFT protocol, fed `consensus_alleles.fasta` instead.
  Writes `identity_consensus.tsv`. Skipped when the consensus has < 2 records.
- **7c — consensus completeness re-check.** `consensus_qc.py`: per-allele
  tblastn(variable_proteins) + blastn(flankL/R) on every consensus record.

```
# 7c body — for each record in consensus_alleles.fasta:
vars_hit     = { gene : tblastn(variable_proteins → consensus, pid >= 30, aa >= 50) hit on gene }
aa_cov       = sum of length(aa) across qualifying tblastn hits on this record
has_flankL   = blastn(flankL → consensus, pid >= 85, len >= 100) has any hit
has_flankR   = blastn(flankR → consensus, pid >= 85, len >= 100) has any hit
has_all_vars = (|vars_hit| == nvar_total)
complete     = has_all_vars ∧ has_flankL ∧ has_flankR
```

Output: `consensus_qc.tsv` with columns
`allele, len, vars_hit, vars_total, has_all_vars, aa_cov, has_flankL, has_flankR, complete`.

Both numbers land in `summary.tsv` alongside the pick-level versions:
`cons_allele{1,2}_complete`, `cons_allele1_vs_allele2_id_pct`, `cons_allele1_vs_allele2_aln_frac`.

---

## 8. Cross-sample mating-type clustering  (`MATdetangler cluster`)

Run after `batch`. Phase A is expensive (one MAFFT alignment of all picked alleles);
phase B is cheap (single-linkage cut on the saved matrix) and re-runnable at any threshold.

```
# A: align (matdetangler/cluster.py :: align)
#    - pool every allele from every results/<sample>/primary_alleles.fasta
#    - for each allele, tblastn the curated variable proteins to find the HD-core
#      span (leftmost-to-rightmost HD-protein hit) and pad by +/- 50 bp
#    - ONE MAFFT --auto run on all cores
#    - emit:
#       cores.fasta        per-allele HD-core +/- 50 bp
#       cores.aln.fasta    MAFFT MSA
#       cores_pairs.tsv    pairwise (a, b, sim, alnid); sim = matches / shorter-core length
#       cores_meta.tsv     (allele, core_len) per record

# B: cut (matdetangler/cluster.py :: cut)
#    - load cores_pairs.tsv + cores_meta.tsv (no re-alignment -- cheap, re-runnable
#      with different --thresh values)
#    - single-linkage union-find: join (a, b) when sim >= --thresh (default 0.90)
#    - clusters numbered 1, 2, ... by size (largest first); singletons keep their own id
#    - emit:
#       allele_classification.tsv  columns: allele, sample, cluster
#       allele_distance_matrix.tsv symmetric matrix of (1 - sim)
```

---

## Performance levers

| concern | knob | effect |
|---|---|---|
| HD envelope too tight (left or right) | `--envelope-padding` (default 500 bp) | widens envelope each side |
| Flank too short for downstream anchoring | `--max-flank-len` (default 2000) | also raises `derived_max_locus_len` |
| Near-identical paths surviving as distinct emitted candidates | `divergence_threshold` (default 0.05, edlib HW edit-distance in `find_alleles`) | raise → keep more variants as distinct; lower → merge more |
| BFS blowing up the bubble at a high-copy hub | `lo_mult` / `hi_mult` (defaults 0.25 / 2.0 in `find_alleles`) | per-pass coverage filter band `[lo × D_k, hi × D_k]`; tighter band → fewer noisy segments enter the bubble |
| Caller wall-time blowing up on tangled k=33 graphs | `max_nhop` (default 10 in `find_alleles`) | tighter cap; BPL1195-shape combinatorial blow-ups still possible — see TODO `_enum_paths` node-visit budget |
| Locus-trim too aggressive or too loose | `locus_padding` (default 1500 bp in `find_alleles`) | bp pad each side of the tblastn HD-hit envelope before sequence emission |
| Cross-K picks the wrong K | `pick_k` priority order + tie-break vector | adjustable (see §4); typically k53 wins for Suilu, k77 for Pcub |
| BLAST wall-time on iterative dev | per-K cached BLAST TSVs in `<outdir>/<K>/` | first run blasts; downstream stages re-read; wipe the per-K dir to force re-blast |

---

## Failure modes the algorithm protects against

| failure | guard |
|---|---|
| Intronic stops in user proteins | rejected at input by `_validate_fasta` (expected to be protein, not NT) |
| Duplicate protein IDs (blast unique-sseqid constraint) | sanitizer suffixes them `_1`, `_2` |
| HD envelope shrunken because tblastn missed divergent ends | union over ALL accepted HSPs + padding |
| Flank too long → search blowup | `--max-flank-len` trim |
| Decoration-variant paths surviving as separate candidates | edlib HW edit-distance dedup at 5% inside `find_alleles` (sequences ≥ 95% identical collapse) |
| GFA self-loop edges (`L\tA\t+\tA\t+\t...`) crashing adjacency build | skipped in `bubble_bfs.build_adj` and `seg_processor.directional_split` (1-element frozenset case) |
| Two GFA components disjoint but both var-bearing | classifier returns `separate`; cross-K picker may pick a K where they connect |
| No flank labels in the BFS-expanded subgraph | universal-leaf anchor rule (§3.4): bubble dead-end Uvars qualify as anchors |
| Long conserved flank diluting k-mer similarity | replaced with edlib HW edit-distance after tblastn locus-trim |
| Picks not actually divergent | `is_divergent` check in `_emit_result` |
| Consensus drift from pick after read mapping | `consensus_qc.py` re-check + consensus divergence pass |
| Read-derived gap in HD body | `core_meandepth` at the HD-core positions |
| Cross-k pair makes a confusing bubble | per-K node-ID prefix in `bubble.png` |


## See also

- `TODO.md` — open algorithmic improvements (per-iteration timeout,
  depth-weighted decomposition for n>2 cases, Uvar-leaf-only emission,
  label/text gene order)
- `METHODS.legacy.md` — preserved methods for the retired
  `graph_path_search.py` + `pick_alleles.py` chain and the synthetic
  test-bench classifier (Appendix A of the legacy file)
- `README.md` — installation + quickstart
