# MATdetangler

Recover the two alleles of a tetrapolar mating-type-like locus from paired-end short reads,
using a **pre-built SPAdes assembly graph** as the substrate. Originally designed for the
fungal HD locus (HD1 / HD2 + conserved flanks); generalizes to any analogous locus with
conserved flanks and one or more variable genes (P/R, PR, idiomorph-like loci).

> Methods + rationale + pseudocode: **[METHODS.md](METHODS.md)**.

---

## What it actually does

Steps 1–5 are **reads-free** — only the SPAdes assembly + locus reference are needed.
Reads (`--reads-r1`/`--reads-r2`) are required ONLY when you opt in to the read-derived
consensus path (steps 6 and 7) via `--make-consensus`.

1. **Input processing → query construction** — tblastn user-supplied HD proteins against
   the locus reference. HSPs are resolved in confidence order with cross-protein
   non-overlap (so HD1's weak partial HSPs in HD2's region don't leak into HD1's
   reported span), per-gene spans come from the union of accepted same-protein HSPs,
   intergenic regions are reported as-is, and the HD envelope is the union of gene
   spans padded by 500 bp each side (module default). Flanks are then auto-derived as
   the locus sequence outside the envelope (trimmed to `--max-flank-len` per side).
2. **Anchor search** — per k, tiered tblastn / blastn / tblastx of the curated queries
   against (2.1) `contigs.fasta` and, when 2.1 comes up short, (2.2) the GFA segments.
   Step 2.05 estimates `genome_cov` per k from SPAdes `contigs.fasta` `cov_` headers
   (median across contigs >= `--contig-depth-size-cut` bp, default 5000). Bp-equivalent
   units; same scale as the GFA `DP:f:` segment depths used downstream. Cached as
   `genome_cov_spades_k<K>.txt`.
3. **GFA deep search — anchor + BFS** — anchor on step-2 hits, BFS-expand the GFA
   neighborhood, enumerate paths from any flankL-bearing to any flankR-bearing
   segment, drop RC-mirror duplicates, and widen the BFS until MAFFT-core clustering
   reports `--expected-count` truly distinct alleles. Fires under `DEEP_SEARCH=always`
   or when step 2 produced fewer than `--expected-count` HD-bearing contig anchors.
4. **Allele selection (pick_alleles)** — two-level pair-based selection. **Level 1 (per K)**:
   enumerate every C(N,2) pair; drop duplicates (HD-core MAFFT id and aln_frac above
   module-level `dup_id` / `dup_frac`); rank surviving pairs by a 7-element lexicographic
   tuple **(var_compl_rank, topology_rank, flank_compl_rank, -HD_core_id,
   joint_HD_aa_cov, joint_flank_bp_cov, -|len_diff|)**. Completeness is primary —
   `var_compl_rank` encodes how many of the variable genes each allele carries
   (lo*10+hi, so 22 > 12 > 11 > 2 > 1 > 0, balanced pairs beat lopsided ones).
   Topology comes next (closed_bubble > open_bubble > detached > other from
   set-intersection of segment IDs). Divergence is a tiebreaker within the same
   completeness/topology layer. MAFFT MSA is computed ONCE per completeness layer,
   only on candidates participating in the leading layer's pairs. **Level 2 (cross K)**:
   each K's winning pair is re-scored with the same tuple, and the K with the
   highest pair wins. Skipped by `--skip-pick` (default ON).
5. **Annotated allele walks (graph_paths)** — render each picked (or candidate) allele
   back on the GFA: ASCII walk, sub-GFA, DOT, PNG, edge-list TSV. Reads-free.
6. **Mapping + read-derived consensus (map_consensus.sh)** — `--make-consensus` opts in.
   Competitive bowtie2 `--end-to-end`, `samtools consensus` per allele, HD-core depth.
7. **Consensus QC (consensus_qc + pairwise_identity rerun)** — `--make-consensus` opts in.
   tblastn / blastn re-check completeness on the consensus; MAFFT divergence re-run on the
   consensus pair. Keep both pick-level AND consensus-level numbers side-by-side in
   `summary.tsv`.
8. **Cross-sample mating-type clustering (`MATdetangler cluster`)** — separate post-batch
   command. `align` extracts each picked allele's HD-core (variable-gene span +/- 50 bp),
   runs ONE MAFFT, and emits the pairwise similarity matrix; `cut` does single-linkage
   on the cached matrix at `--thresh` (default 0.90) and writes `allele_classification.tsv`
   + `allele_distance_matrix.tsv`. Cheap to re-run `cut` at any threshold.

## Installation

One-time conda environment (covers every external binary the pipeline needs:
SPAdes, BLAST+, bowtie2, samtools, MAFFT, plus Python + matplotlib):

```bash
conda env create -f install/env.yml      # ~5 min on first run
conda activate matdetangler
```

The pipeline is pure Python stdlib + a handful of subprocess calls — no biopython,
numpy, or other heavy Python deps. If you already have all the external binaries on
`$PATH` (`spades.py`, `makeblastdb`, `tblastn`, `blastn`, `bowtie2`, `samtools`,
`mafft`) you can skip the conda env entirely.

## Quick start

```bash
# Step 0 — build per-k SPAdes assemblies (skip if you already have them)
MATdetangler-spades \
  --sample Tu127439 \
  --reads-r1 reads/R1.fq.gz --reads-r2 reads/R2.fq.gz \
  --outdir examples/Tu127439_spades/ \
  --ks 33,45

# Steps 1-5 (reads-free, default) — recover the two alleles
MATdetangler run \
  --sample Tu127439 \
  --spades-dir examples/Tu127439_spades/ \
  --locus-ref examples/Suilu_locus/Suilu4_MATA.fasta \
  --proteins  examples/Suilu_locus/Suilu4_HDs.fasta \
  --outdir results/ \
  --ks k33,k45

# Steps 1-8 with read-derived consensus + QC (opt in)
MATdetangler run \
  --sample Tu127439 \
  --reads-r1 reads/R1.fq.gz  --reads-r2 reads/R2.fq.gz \
  --spades-dir examples/Tu127439_spades/ \
  --locus-ref examples/Suilu_locus/Suilu4_MATA.fasta \
  --proteins  examples/Suilu_locus/Suilu4_HDs.fasta \
  --outdir results/ \
  --ks k33,k45 \
  --make-consensus

# Batch (4-col TSV: sample r1 r2 spades_dir; r1/r2 may be "-" when --make-consensus is off)
MATdetangler batch --samplesheet samples.tsv \
  --locus-ref examples/Suilu_locus/Suilu4_MATA.fasta \
  --proteins  examples/Suilu_locus/Suilu4_HDs.fasta \
  --outdir results/ --threads 8 --ks k33,k45

# Only steps 6+7 on a sample already processed (primary_alleles.fasta on disk)
MATdetangler consensus \
  --sample Tu127439 \
  --reads-r1 reads/R1.fq.gz --reads-r2 reads/R2.fq.gz \
  --outdir results/ \
  --locus-ref examples/Suilu_locus/Suilu4_MATA.fasta

# Cross-sample mating-type clustering (step 8): one `align` pass + as many `cut`s as you want
MATdetangler cluster align \
  --results-dir results/ \
  --queries-dir results/Tu127439/queries \
  --out-dir results/clusters/

MATdetangler cluster cut \
  --align-dir results/clusters/ \
  --results-dir results/ \
  --out-dir results/clusters/ \
  --thresh 0.90
```

`MATdetangler cluster align` builds `cores.fasta` (HD-core span ±50 bp per allele),
`cores.aln.fasta` (one MAFFT of all cores), `cores_pairs.tsv` (pairwise HD-core
identities), and `cores_meta.tsv`. `cluster cut` reads the cached pairs/meta —
no re-alignment — and writes `allele_classification.tsv` (columns: allele, sample,
cluster) + `allele_distance_matrix.tsv` (symmetric `1 - sim`). Re-run `cut` with
different `--thresh` values to explore the mating-type partition.

## Inputs

### Required

| flag | what |
|---|---|
| `--locus-ref FASTA` | single-record DNA reference spanning the HD region + both flanks. ~10–15 kb is plenty. |
| `--proteins FASTA` | **curated** protein fasta (HD1, HD2, …). Real proteins, intron-free. Records with duplicate IDs are auto-suffixed (`_1`, `_2`, …). |
| `--sample`, `--spades-dir`, `--outdir` | per-sample inputs |
| `--reads-r1`, `--reads-r2` | **OPTIONAL** — only required when `--make-consensus` is set. The pipeline is reads-free through step 5 by default. |

### Optional

| flag | default | what |
|---|---|---|
| `--max-flank-len` | 2000 | trim each flank to at most this many bp; also sets the path-enumerator bp budget |
| `--ks LIST` | `k45` | which per-k subdirs of `--spades-dir` to consider (single k by default for now; comma-separated to sweep multiple) |
| `--make-consensus` | off | Run step 6 (`map_consensus.sh`) and the consensus-derived parts of step 7 (`identity_consensus`, `consensus_qc`). Default OFF — the pipeline is reads-free through step 5. `--make-consensus` opts in; `--reads-r1`/`--reads-r2` become required. |
| `--skip-pick` / `--no-skip-pick` | skip ON | QC/dev mode: skip step 4 (pick_alleles) and step 8 (summary row); alias `bubble_alleles.fasta` as `primary_alleles.fasta` so steps 5, 6, 7 run on the full candidate pool. `--no-skip-pick` runs the full pipeline. |
| `--continue` | off | Re-run with additional k's, reusing existing sample-level outputs: step 1 (`input_process`) is skipped if `queries/manifest.json` exists; per-k genome coverage reuses cached `genome_cov_spades_k<K>.txt`; per-k steps 2.1/2.2/3 are skipped for any k whose `bubble_alleles_<k>.fasta` is already on disk (new k's run). Implies `--no-skip-pick` so the picker + downstream run on the combined pool. Pass `--skip-pick` AFTER `--continue` to QC the combined pool instead. |
| `--re-blast` | off | Wipe the per-sample cached blast result TSVs (`<outdir>/<sample>/blast_*.tsv`) before running. Invoke whenever you've changed blast parameters (pid/minlen/threshold cutoffs) — otherwise the named result files on disk are reused as-is (file-existence check only; no input hashing). |
| `--threads` | 4 | blast / bowtie2 / mafft thread count |
| `--expected-count {1,2}` | 2 | 1 = haploid, 2 = dikaryon. Also stops the per-k BFS widen loop as soon as MAFFT-core clustering finds this many distinct alleles. |
| `--genome-coverage FLOAT` | auto-estimated per k (median of SPAdes contigs.fasta `cov_` across contigs ≥ `--contig-depth-size-cut` bp) | mean genomic depth, bp-equivalent units (same scale as GFA `DP:f:` segment depths). Drives the coverage-tier tiebreak in `pick_alleles` AND the repeat-detection cutoff in `graph_path_search` (segments at depth ≥ `cov_repeat_factor × this` are flagged). When supplied, overrides the per-k estimate for ALL k's. |
| `--contig-depth-size-cut INT` | 5000 | minimum contig length (bp) used by the per-k genome-cov estimator. Lower → more contigs but more short-tip noise. |
| `--asymmetric-bfs` | off | enable repeat-aware BFS: repeats (segments at depth > 2.0 x genome_cov) get absorbed into the neighborhood but don't expand from. Used to dampen combinatorial blow-up at high-copy hubs. |
| `--blastn-minlen` | 100 | minimum blastn hit length used for variable-gene classification of contigs/segments |
| `--max-locus-len` | manifest's `derived_max_locus_len` (= envelope + 2x flank; ~8-15 kb on Pcub-sized loci) | **MAFFT-trim threshold for step 4 picker only.** Candidates longer than this get tblastn-trimmed to HD-envelope ± 2500 bp pad BEFORE pairwise MAFFT. Keeps step 4 fast (~8 kb seqs) instead of MAFFT on raw 50-250 kb whole-contig anchors. Pass `--max-locus-len 0` to disable the trim entirely. (Decoupled from step 3's DFS bp cap as of 2026-05-30 — see `--max-walk-bp`.) |
| `--max-walk-bp` | 0 (no cap) | bp budget for step 3 DFS path enumeration. The simple-path rule + module-level `max_nodes=15` already bound DFS on real GFAs; this is a safety net. Independent of `--max-locus-len`. |
| `--init-hops` | (module default 5) | initial radius of the per-k BFS widen loop. |
| `--max-hops` | (module default 10) | step 3 BFS widen-loop max. The loop early-breaks when `--expected-count` distinct alleles have been seen AND the neighborhood topology has settled into `closed_bubble` or `complexed`. |

### Blast cache (always on; live in main results)

Heavy "search contigs / GFA" blasts (step 2 anchor search, step 3 segment labelling) are
cached as descriptive TSVs ALONGSIDE the per-sample outputs. There is no hidden cache
directory: cache entries ARE the result files.

```
results/<sample>/
  blast_hd_tblastn_contigs_k45.tsv         ← step 2 anchor (contigs)
  blast_hd_blastn_contigs_k45.tsv          ← (tier-2 fallback if fired)
  blast_hd_tblastx_contigs_k45.tsv         ← (tier-3 fallback if fired)
  blast_flankL_blastn_contigs_k45.tsv      ← --include-flanks debug only
  blast_hd_tblastn_segments_k45.tsv        ← step 2.2 anchor (segments)
  blast_gfa_hd_tblastn_nhood_k45_hops5.tsv ← step 3 BFS-neighborhood labelling
  blast_gfa_flankL_blastn_nhood_k45_hops5.tsv
  blast_gfa_flankR_blastn_nhood_k45_hops5.tsv
  ...
```

**Freshness check**: file existence ONLY. There is no hash, no `query.fa` byte-verify,
no mtime check. If the file exists, it's reused; if it doesn't, blast runs and writes
it. Filters (e.g. `min_pid`) are re-applied at read time, so tightening a threshold
without re-blasting subsets the cached file correctly.

| flag | behavior |
|---|---|
| (default — always on) | The wrapper passes `--blast-out-dir <outdir>/<sample>/` to every step. Subsequent runs skip blast for any cached TSV. |
| `--re-blast` | Wipe `<outdir>/<sample>/blast_*.tsv` before the run. **Use this whenever you change blast parameters** (pid/minlen/threshold). The file-existence check has no way to know the cached result was generated with different cutoffs. |

Detail blasts (`pairwise_identity.detect_core_span`, completeness checks, consensus
QC, multi-sample cluster) are NOT cached on disk — they're small, sample-internal,
and re-run each time. `detect_core_span` is memoized in-process so the same
candidate is tblastn'd at most once per run.

### Deep search (env var)

`DEEP_SEARCH=auto|always|never` (env). `auto` (default) fires GFA deep search only when
the contig pick produced fewer than `--expected-count` complete picks.

## Preparing the SPAdes graphs (step 0)

MATdetangler needs `contigs.fasta` + `assembly_graph_after_simplification.gfa` for
**each k** in the sweep. The included wrapper `MATdetangler-spades` does this in the
exact layout downstream steps expect:

```bash
# Single sample (submits one SLURM job per k; `run` subcommand is optional)
MATdetangler-spades \
  --sample Tu127439 \
  --reads-r1 R1.fq.gz --reads-r2 R2.fq.gz \
  --outdir examples/Tu127439_spades/ \
  --ks 33,45,53

# Batch (3-col TSV: sample r1 r2) — submits a SLURM array
MATdetangler-spades batch \
  --samplesheet samples.tsv \
  --outdir examples/Tu127439_spades/ \
  --ks 33,45,53 \
  --threads 8 --mem 24 --time 2:30:00 \
  --account my-acct --partition standard

# Laptop / no SLURM
MATdetangler-spades batch --samplesheet samples.tsv --outdir examples/Tu127439_spades/ --local
```

Defaults:

- `--ks 33,45,53` (override with any comma-separated list, e.g. `21,33,55`)
- `--threads 8`, `--mem 24` (GB), `--time 2:30:00`
- **No BayesHammer** (uses `--only-assembler`), assuming the reads were already adapter-
  and quality-trimmed. Add `--bayes-hammer` if you want SPAdes to do the read correction.
- Layout written: `<outdir>/<sample>/k<K>/{contigs.fasta, assembly_graph_after_simplification.gfa}`
- **Skip-safe**: any (sample, k) whose two output files already exist is skipped, so
  re-running the wrapper only fills in the missing tasks.

Why one SPAdes job per k. SPAdes' multi-k mode (`-k 33,45,53` in one call) keeps only
the FINAL-k GFA — the intermediate k GFAs are stripped. MATdetangler's `graph_path_search`
benefits from the full sweep, so we run one assembly per (sample, k) pair instead.
`matdetangler/paths.py` resolves either `k33/` (lowercase) or `K33/` (uppercase) layouts,
and falls back to a top-level final-k GFA + contigs if no per-k subdir is found.

## Outputs (per sample, in `results/<sample>/`)

![bubble schema](results/AU340/bubble.png)

| file | what |
|---|---|
| **`primary_alleles.fasta`** | **The picked pair (canonical output)** — straight from the GFA / contig walk. 0 to `--expected-count` records. |
| `consensus_alleles.fasta` | `samtools consensus` per allele, all records in one multi-record FASTA (IDs: `allele_1`, `allele_2`, …). Only written when `--make-consensus`. |
| `picks.tsv` | per-pick: k, type, length, source contig + GFA segments, depth, `n_variable_genes`, `has_both_flanks`. The upstream `bubble_alleles.ann.tsv` records `flankL` and `flankR` as separate booleans so partial / open-bubble walks remain distinguishable. |
| `anchor_contig.fasta` (+ `_k<K>.fasta`) | step 2.1 — whole-contig anchors that hit ≥1 HD or flank query. Record names: `<sample>__bubble_<k>_<contig_id>` (the "bubble" tag is a legacy of when contig-anchored records were called bubble candidates). Companion `anchor_contig.ann.tsv` is 6 cols: name, len, k, kind, hd_genes, flanks. |
| `anchor_segments.fasta` (+ `_k<K>.fasta`) | step 2.2 (conditional) — GFA-segment-level anchors. Record names: `<sample>__seg_<k>_<segment_id>`. Same 6-col ann.tsv schema as `anchor_contig.ann.tsv`. Only written when step 2.1 yielded too few HD-bearing anchors. |
| `bubble_alleles.fasta` (+ `_k<K>.fasta`), `bubble_alleles.ann.tsv` | step 3 — flank-to-flank path walks through the GFA bubble, RC-dedup'd + MAFFT-clustered. `bubble_alleles.ann.tsv` columns: name, len, k, segments, variable_genes, **flankL**, **flankR**, cov. |
| `picker_candidates.fasta`, `picker_candidates.ann.tsv` | step 4 input pool — emitted at end of step 3 by `_emit_picker_candidates`: filtered anchor contigs (drop those with 0 variable genes AND <2 flanks) ∪ all bubble paths. `cand_combined.fasta` + `cand_combined.ann.tsv` are kept as back-compat aliases of the same files. The picker reads either name. |
| `reads.sorted.bam` (+ `.bai`) | competitive end-to-end bowtie2 mapping reads → picks |
| `coverage.tsv` | per-allele depth: `whole_meandepth` AND `core_meandepth` (mean depth across HD-core positions only) |
| `identity.tsv` | MAFFT id_pct + aln_frac on the **picks** |
| `identity_consensus.tsv` | MAFFT id_pct + aln_frac on the **read-derived consensus pair** |
| `consensus_qc.tsv` | tblastn(proteins → consensus) + blastn(flanks → consensus); per-allele `complete` flag |
| `bubble.txt` | two-line labeled walks (one per allele); matches `summary.tsv` path columns |
| `bubble.png` | matplotlib: each allele on its own row; flank-bearing nodes linked by dashed gray |
| `bubble.gfa` / `.dot` / `.tsv` | sub-GFA + Graphviz + edge list for any network library |
| `summary.tsv` (per-sample) | one row, same schema as the aggregate's row for this sample |
| `queries/` | the auto-derived `variable_proteins.fasta`, `flankL.fasta`, `flankR.fasta`, `variable_nt.fasta`, `manifest.json` |
| `genome_cov_spades_k<K>.txt` | per-k genome coverage (median of contigs.fasta `cov_` ≥ `--contig-depth-size-cut`). Bp-equivalent units. Cached for `--continue`. |
| `blast_*.tsv` | cached blast hit tables — see [Blast cache](#blast-cache-always-on-live-in-main-results). Wipe with `--re-blast`. |
| `logs/` | per-step stdout/stderr |

Plus one aggregate: `results/summary.tsv` (one row per sample). Schema:

```
sample  used_k  bubble_type  genome_coverage
allele1_complete          allele2_complete          ← pick-level (path completeness)
allele1_coverage          allele2_coverage          ← consensus-mapped HD-core depth
allele1_vs_allele2_id_pct allele1_vs_allele2_aln_frac ← pick-level divergence
cons_allele1_complete     cons_allele2_complete     ← consensus tblastn re-check
cons_allele1_vs_allele2_id_pct cons_allele1_vs_allele2_aln_frac ← consensus divergence
allele1_path  allele2_path
```

## Dependencies

All wrapped by the conda env in `install/env.yml`. External tools (on `$PATH`):
- `spades.py` ≥ 3.15 (only needed for step 0; `MATdetangler-spades`)
- BLAST+ (`makeblastdb`, `blastn`, `tblastn`) ≥ 2.10
- `bowtie2` ≥ 2.4 (and `bowtie2-build`)
- `samtools` ≥ 1.15
- `mafft`

Python ≥ 3.9 (uses PEP 604 union types). The pipeline imports only the standard
library + `matplotlib` (for `bubble.png`). No biopython / numpy required.

## License & citation

TBD. If you use MATdetangler in a paper, please cite this repository.
