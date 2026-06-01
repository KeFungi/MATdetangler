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
| `D` | genome diploid coverage (auto-estimated per-k as the median of SPAdes `contigs.fasta` `cov_` headers across contigs ≥ `--contig-depth-size-cut` bp; `--genome-coverage` overrides for all k's) |
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

## 2. Anchor search  (`anchor_search.py`)

Same module runs in two source modes per k. The wrapper schedules them per-k, with
2.2 firing only when 2.1 came up short.

### 2.05 Genome-coverage estimate  (MATdetangler wrapper, PER k inside the per-k loop)

Run inside the per-k loop just before step 2.1 anchor search. Fires only when the user
did NOT pass `--genome-coverage` explicitly (the user value, if supplied, overrides
all k's). The estimate `D_k` is consumed by `graph_path_search` (repeat-by-coverage
rule, see §3) and by `pick_alleles` (coverage-tier tiebreak, see §4).

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
  scaled by SPAdes to per-base coverage). So `genome_cov` is directly comparable to
  per-segment depth in step 3's repeat cutoff (`depth > cov_repeat_factor × genome_cov`).
- **Median, not mean.** Contig depths have a heavy right tail from collapsed-repeat
  short contigs (`cov_` up to ~10⁷ on 46 bp tips in real fungal data). Median is robust
  to those.
- **Length filter at 5 kb.** Stays well within the noise-free regime; smaller contigs
  contribute too much variance from short-tip assembly artifacts. Tunable via
  `--contig-depth-size-cut` if your assembly is unusually fragmented.
- **PER k, not once per sample.** SPAdes coverage values are in k-mer-multiplicity ×
  `R/(R-k+1)` units, so they DIFFER across k's for the same sample. The k=33 estimate is
  ~25% higher than k=53 for R=150 bp reads. Computing per-k means the repeat cutoff used
  by step 3 at k=K is automatically scaled to that k's depth distribution.
- **Cached as `genome_cov_spades_k<K>.txt`.** `--continue` reuses the cached value.
- **What if no contig is long enough?** `D_k = 0.0`. The pipeline treats this as
  "disable coverage-based repeat detection" — step 3 still runs, the repeat cutoff just
  never fires. The picker's coverage-tier tiebreak becomes a no-op.
- **(Historical)** Before the reads-free refactor (2026-05-30), the wrapper estimated
  genome cov via bowtie2 `--end-to-end --very-sensitive` to the whole locus, sampling
  depth from the interior of the flank windows with a 150 bp edge buffer. That value was
  in bp-coverage units, mismatched against k-mer-units segment depths — the new SPAdes
  source is units-consistent AND cheaper.

### 2.1 Contig anchors  (`anchor_search.py --source contigs`)

Per k, fast tiered anchor-search on `contigs.fasta`:

```
for k in --ks:
    contigs = read contigs.fasta(k)
    gene_per_contig = {}

    # TIER 1 — tblastn variable_proteins (ALWAYS)
    tier1 = tblastn(variable_proteins → contigs, pid ≥ 30, aln ≥ 50 aa)
    n_complete_1 = |{c : gene_per_contig[c] == ALL variable genes AND len(c) ≥ min_allele_len}|

    # TIER 2 — blastn variable_nt (FALLBACK 1)
    if n_complete_1 < --min-complete-tblastn (default 2):
        tier2 = blastn(variable_nt → contigs, pid ≥ 80, len ≥ 100)

    # TIER 3 — tblastx variable_nt (FALLBACK 2, expensive)
    if n_complete_2 < --min-complete-tblastn:
        tier3 = tblastx(variable_nt → contigs, pid ≥ 30, aln ≥ 50 aa)

    # emit anchors (record name "<sample>__bubble_<k>_<contig_id>")
    cands_k = { c ∈ (tier1 ∪ tier2 ∪ tier3) : len(c) ≥ --min-allele-len }
    write anchor_contig.fasta + anchor_contig.ann.tsv
        # 6 cols: name, len, k, kind (= "bubble"), hd_genes (comma-sep or "-"),
        #         flanks (comma-sep or "-")
    # flank blastn (--include-flanks) is debug-only; by default only HD-bearing anchors emit.
    # The wrapper renames to anchor_contig_<k>.{fasta,ann.tsv} after each k.
```

### 2.2 GFA-segment anchors  (`anchor_search.py --source segments`)

Fires per-k only when 2.1 returned 0 HD-bearing anchors OR when step 3 came back short.
Same tiered logic, but the database is the GFA segments (S-lines from
`assembly_graph_after_simplification.gfa`) instead of `contigs.fasta`. Useful when the
locus is fragmented across many short segments that didn't get assembled into a contig.

```
gfa_segs   = parse_S_lines(spades_dir/k/assembly_graph_after_simplification.gfa)
write each segment as its own fasta record (seg_id) → segments_<k>.fa
... same tier 1/2/3 chain as 2.1, db = segments_<k>.fa ...
# wrapper hard-codes --min-len 200 for segments (vs --min-allele-len=2000 for contigs),
# since GFA segments are typically much shorter than contigs.
write anchor_segments.fasta + anchor_segments.ann.tsv (same 6 cols as 2.1,
                                                       kind = "seg",
                                                       record name "<sample>__seg_<k>_<seg_id>")
# The wrapper renames to anchor_segments_<k>.{fasta,ann.tsv} after each k.
```

### Rationale

- **Three escalating tiers.** Each lower tier is more divergence-tolerant and slower
  than the one above, and only fires when the previous tier(s) didn't already surface
  `--min-complete-tblastn` (default 2) candidates that carry every variable gene.
  - **Tier 1 (tblastn)** at 30 % pid finds any sequence with detectable HD protein
    homology. Usually sufficient.
  - **Tier 2 (blastn)** at 80 % pid, 100 bp rescues sequences whose tblastn pid was
    just under the cutoff but whose NT sequence is well-conserved.
  - **Tier 3 (tblastx)** translates BOTH the variable-gene NT panel AND the subject
    in 6 frames before alignment. ~36× slower than tblastn; reserved for the rare
    deeply diverged case where tiers 1+2 didn't reach the threshold.
- **Two sources (contigs vs segments)** because SPAdes sometimes assembles the locus
  into one contig per allele (2.1 hits) and sometimes leaves it fragmented across
  many GFA segments (2.2 rescues those).
- **HD-only emission by default.** Flank-only anchors get filtered out so the seeds
  to step 3's BFS are HD-bearing. `--include-flanks` (debug) re-enables flank-only
  emission.

---

## 3. GFA deep search — anchor + BFS  (`graph_path_search.py`)

Runs only when (a) `DEEP_SEARCH=always` or (b) `DEEP_SEARCH=auto` and step 2 produced
fewer than `--expected-count` HD-bearing contig anchors.

The wrapper iterates k by k (default `k45` — single k). Per-K early break was REMOVED
(2026-05-30): the wrapper always runs every k in `--ks` so the Level-2 cross-K picker
in step 4 sees per-k pairs from every k. Inside `graph_path_search.run` itself, the
in-Python per-k loop DOES still short-circuit on `n_distinct >= expected_count`
(MAFFT-core distinct) — that's a within-call optimization, immaterial when the wrapper
passes a single k.

Default constants (function signature of `graph_path_search.run`):
`init_hops=5`, `max_hops=10`, `expected_count=2`, `dup_id=0.95`, `dup_frac=0.80`,
`max_nodes=15`, `max_paths=50000`, `max_walk_bp=0` (no DFS bp cap),
`cov_repeat_factor=2.0`, `min_allele_len=2000`.

The widen loop **breaks** when `n_distinct ≥ expected_count` AND the
neighborhood topology has settled into `closed_bubble` OR `complexed`. Any
other topology (`no_main`, `open_bubble`, `detached`, `single`) keeps widening
until `max_hops`. The widen-on-non-bubble rule prevents premature exit on
half-formed bubbles whose two arms have been seen but haven't connected yet.

```
for k in --ks:

    # 3.1 Anchor discovery (already done in step 2.1 / 2.2; this stage just LOADS)
    #     - Step 2.1 wrote anchor_contig_<k>.fasta + .ann.tsv  (contig anchors)
    #     - Step 2.2 wrote anchor_segments_<k>.fasta + .ann.tsv (direct GFA segment anchors)
    #     `collect_anchors_per_k` partitions records by k from the names and,
    #     by default (`--seeds-from hd`), filters to records whose ann.tsv
    #     hd_genes column is non-"-" (carry a variable gene). `--seeds-from all`
    #     keeps flank-only seeds too.

    # 3.2 Bridge anchor contigs → GFA segments via contigs.paths (authoritative oriented walks)
    paths_file = find_contigs_paths_file(spades_dir/k/)          # contigs.paths / final_contigs.paths / scaffolds.paths
    contig_to_segs = parse_contigs_paths(paths_file)
    # Contigs absent from contigs.paths are dropped with a warning.
    anchor_segs = ∪ { contig_to_segs[c] : c ∈ anchor_contigs }  ∪  direct_segment_anchors

    # 3.3 High-copy hub detection (used only when --asymmetric-bfs is on)
    # Coverage-only — segments with depth > cov_repeat_factor × D are flagged so
    # they can be VISITED but not EXPANDED in the next BFS layer.
    repeat_segs  = { s : dp(s) > cov_repeat_factor × D }          # cov_repeat_factor = 2.0 default

    # 3.4 Widen-BFS loop, stop when MAFFT-core clustering says enough distinct alleles
    hops = INIT_HOPS                                              # default 5
    while hops ≤ MAX_HOPS:                                        # default 10
        neighborhood = bfs_expand(L-line undirected adj, anchor_segs, hops,
                                  repeat_segs if --asymmetric-bfs else None)
        labels, var_per = label_segments(neighborhood)            # tblastn proteins + blastn flanks
        starts   = { s ∈ neighborhood : "flankL" ∈ labels[s] }
        ends     = { s ∈ neighborhood : "flankR" ∈ labels[s] }
        var_segs = { s ∈ neighborhood : var_per[s] ≠ ∅ }
        if starts == ∅ or ends == ∅ or var_segs == ∅:
            hops += 1; continue                                   # too tight; widen

        adj_sub  = restrict_adj_to_subgraph(adj_dir, neighborhood)
        paths    = enumerate_paths(adj_sub, segs, starts, ends, var_segs,
                                   max_bp, max_nodes, max_paths)
        # enumerate_paths drops RC-mirror duplicates at emit time
        # (`mirror_path(P) ∈ emitted` ⇒ skip) — see §3-helpers.

        complete = [ p : is_complete_path(p) ]

        # cluster_alleles — single-shot MSA pipeline:
        #   a. tblastn-trim each candidate + reference to HD-core ± 1000 bp pad
        #      (detect_core_span is memoized in-process so the same seq is
        #       blasted at most ONCE per run, no matter how many widen iterations).
        #   b. ONE MAFFT --adjustdirection --auto MSA of trimmed reference + all
        #      trimmed candidates (the reference is the coordinate ruler).
        #   c. HD-core columns = MSA cols where the (trimmed) reference's
        #      un-gapped position falls in its HD-core span.
        #   d. Pairwise text identity on the HD-core columns only; aln_frac =
        #      aligned-bp / |HD-core cols|.
        #   e. Greedy single-link cluster at id >= dup_id AND frac >= dup_frac.
        recs_for_cluster = [ (label, reconstruct(p, segs)) for p in complete ]
        clusters = cluster_alleles(recs_for_cluster, queries_dir, locus_ref_fa,
                                   id_thresh=dup_id, frac_thresh=dup_frac)
        topology = classify_neighborhood_topology(nhood, adj_und, segs, labels,
                                                    var_per, min_core_len=1000)
        # Break when `n_distinct >= expected_count` AND topology has settled
        # into a stable shape — either `closed_bubble` (the ideal) or
        # `complexed` (main segments don't share flank anchors and widening
        # won't change that). For unstable topologies (no_main, etc.) keep
        # widening until max_hops.
        if len(clusters) ≥ expected_count and topology.type in {"closed_bubble", "complexed"}: break
        hops += 1                                                 # too few distinct OR topology still unstable — widen

    # 3.5 Emit per-k candidates with the SAME dedup pipeline
    #     (so bubble_alleles_<k>.fasta has ONE record per MAFFT cluster,
    #      not all canonical-distinct paths). The cluster representative is
    #      the LONGEST canonical sequence in the cluster.
    canon = { min(reconstruct(p), rc(reconstruct(p))) → longest_p for p in complete }
    clusters = cluster_alleles(canon, queries_dir, locus_ref_fa,
                               id_thresh=dup_id, frac_thresh=dup_frac)
    representatives = [ longest(cluster) for cluster in clusters ]
    for p in representatives:
        seq = reconstruct(p, segs); if |seq| < min_allele_len: continue
        emit (seq, segments, var_hits, flankL, flankR, depth)
        →  bubble_alleles_<k>.{fasta, ann.tsv}
        # bubble_alleles_<k>.ann.tsv: no header — name, len, k, seg_str
        # (e.g. "12+,34-,56+"), variable_genes (comma-list or "-"), has_flankL (T/F),
        # has_flankR (T/F), depth (length-weighted mean), plus two back-compat
        # placeholder columns retained for column-index stability.

    # 3.6 If MAFFT-core distinct ≥ expected_count INSIDE graph_path_search.run, the
    #     in-Python per-k loop early-breaks. The OUTER wrapper loop, however, no
    #     longer early-breaks across k's (2026-05-30): every k in --ks runs to feed
    #     the Level-2 cross-K picker.

# After the per-k loop:
#   cumulative bubble_alleles.fasta is REBUILT from all per-k files
#   (each step 3 call truncates bubble_alleles.fasta at start, so the cumulative
#    cannot be appended in-loop; it's reconstructed once at the end).
#
# Per-k side outputs:
#   - seg_labels_<k>.tsv      : seg_id<TAB>label_str  (per-segment "+"-joined feature
#                                tags from the final widen-loop iteration; consumed by
#                                pick_alleles for pair-level topology classification)
#   - bubble_topology.tsv     : sample, k, type, n_main_seg, n_path, n_shared_anchor,
#                                n_shared_flank_anchor, hd_seg_lens — one row per k
#
# End-of-step-3:
#   _emit_picker_candidates() writes
#     - picker_candidates.fasta + picker_candidates.ann.tsv
#       = filtered anchor contigs ∪ all bubble paths (filter: drop any anchor contig
#         whose anchor_contig_<k>.ann.tsv row has 0 variable genes AND lacks at least
#         one of flankL/flankR). See §4 for picker consumption.
```

### Subroutines

**`label_segments(neighborhood, ...)`** — assigns content tags per segment via blast against the queries.

Each segment can carry MULTIPLE labels because `feats[sid]` is a list that gets appended to
across the blast passes:

| pass | query | tag(s) appended on hit |
|---|---|---|
| tblastn | variable_proteins (HD1 / HD2 / ...) | the variable gene name (`HD1`, `HD2`, ...) |
| blastn (fallback)  | variable_nt — only if tblastn missed a gene globally | same |
| blastn | flankL.fasta | `flankL` |
| blastn | flankR.fasta | `flankR` |

The final label string is `"+".join(feats[sid])` after dedup. A segment can therefore
legitimately carry, e.g.:

- `HD1+HD2`              — both variable genes hit the same long segment
- `HD1+flankL`           — HD1 on the same segment that extends into flankL
- `flankL+flankR`        — a tiny segment between the two flanks where both blasts hit

**Per-variable-gene tracking** is independent of the label string. `var_per` is a `{seg_id: set}`
mapping, so a segment with `var_per[sid] = {"HD1", "HD2"}` contributes both genes to a path's
completeness count.

**Downstream membership tests** all use `labels.get(sid, "").split("+")`:

```python
seen_flankL = "flankL" in labels[seg].split("+")    # is_complete_path
seen_flankR = "flankR" in labels[seg].split("+")
flank_segs  = { s : "flankL" in toks  or  "flankR" in toks }   # classify_neighborhood_topology
```

So an `HD1+flankL` segment correctly counts as BOTH an HD1 carrier (via `var_per`) AND a
flankL-bearing anchor (via the label split).

**`bfs_expand(adj, seeds, hops, repeat_segs=None)`** — undirected BFS up to `hops` levels. With
`repeat_segs` supplied (passed through `--asymmetric-bfs`), a repeat segment in the frontier
gets visited but does NOT expand its own neighbors — paths can still REACH a repeat, but the
repeat doesn't drag in its many adjacencies and blow up the search.

**`enumerate_paths(adj_dir, segs, starts, ends, must_visit_any, max_bp, max_nodes, max_paths)`**
— iterative DFS with state `(seg, orient, path, visited, total_bp)`. `max_bp <= 0` is the
no-bp-limit sentinel (default in the new wrapper). RC-mirror dedup is built in at emit time:

```
emitted ← ∅                                                # set of path-tuples already emitted

for s0 ∈ sorted(starts), o0 ∈ {+, -}:
    stack ← [(s0, o0, [(s0, o0, 0)], {s0}, |segs[s0]|)]
    while stack:
        if |paths| ≥ max_paths: return paths
        seg, orient, path, visited, bp = stack.pop()
        if seg ∈ ends and visited ∩ must_visit_any ≠ ∅:
            # SAFE-FORM RC-mirror dedup: only drop if the mirror was already emitted.
            # Naive form (drop when tuple(path) > tuple(mirror)) would lose walks whose
            # mirror is never enumerated (the typical HD case with disjoint flankL/flankR).
            if mirror_path(path) ∉ emitted:
                emitted.add(tuple(path)); paths.append(path)
            continue
        if |path| ≥ max_nodes: continue
        if max_bp > 0 and bp ≥ max_bp: continue            # bp budget; 0 = no cap
        for (next_seg, next_orient, overlap) ∈ adj_dir[(seg, orient)]:
            if next_seg ∈ visited: continue                 # simple path
            new_bp = bp + |segs[next_seg]| − overlap
            if max_bp > 0 and new_bp > max_bp: continue     # bp prune; 0 = no cap
            stack.push((next_seg, next_orient,
                         path + [(next_seg, next_orient, overlap)],
                         visited ∪ {next_seg}, new_bp))
```

**`mirror_path(P)`** — same physical walk traversed in opposite direction; computed purely
from path topology (no sequence reconstruction):

```
mirror[0]   = (P[-1].seg, flip(P[-1].o), 0)
mirror[i>0] = (P[-1-i].seg, flip(P[-1-i].o), P[-i].ov)     # overlaps shift one slot
```

`mirror(mirror(P)) == P`. In the typical HD topology (starts = flankL, ends = flankR,
disjoint), mirror walks would have to seed from a flankR segment, which `enumerate_paths`
never does — so no paths get dropped. When starts ∩ ends ≠ ∅ (or some other symmetric
case), the safe-form dedup collapses each mirror pair to one emission.

**`is_complete_path(path)`** — coverage-only completeness check:

```python
seen_genes = set()
sL = sR = False
for sid, _, _ in path:
    seen_genes |= var_per.get(sid, set())          # ALL variable genes on this seg
    toks = labels.get(sid, "").split("+")          # tokenize the multi-tag label
    if "flankL" in toks: sL = True
    if "flankR" in toks: sR = True
return (len(seen_genes) == nvar_total) and sL and sR
```

**Multi-tag labels are tokenized on `"+"` before membership testing**, so a segment whose
label is `"HD1+flankL"` correctly counts as BOTH an HD1 carrier and a flankL-bearing
anchor — there is no string-match pitfall where `"flankL"` would fail because the raw
label is `"HD1+flankL"`. Worked examples:

| label | `var_per[sid]` | tokens | adds to `seen_genes` | sets sL | sets sR |
|---|---|---|---|---|---|
| `HD1+HD2` | `{HD1, HD2}` | `[HD1, HD2]` | `{HD1, HD2}` | — | — |
| `HD1+flankL` | `{HD1}` | `[HD1, flankL]` | `{HD1}` | ✓ | — |
| `flankL+flankR` | `∅` | `[flankL, flankR]` | (nothing) | ✓ | ✓ |

No length condition. A 3-segment compact walk and a 30-segment wandering walk that
both cover the same content are equally complete.

### Rationale

- **Anchor-then-BFS** instead of labeling every GFA segment: a fungal GFA has 100k+
  segments, most of them irrelevant noise. By starting from HD-anchor segments and
  BFS-expanding, we restrict labeling and DFS to a neighborhood of ~hundreds of segments.
  Wall-time per k drops from ~30 min to a few minutes.
- **HD-only seeds by default.** `--seeds-from hd` filters anchor records to those that
  carry a variable gene (ann.tsv `hd_genes` column non-`-`). Flank-only seeds typically
  pull in noise from the conserved flanks shared across distant loci, and the BFS reaches
  flank-bearing segments from HD seeds via L-line adjacency anyway. `--seeds-from all`
  is available for debug.
- **Widen iteratively until MAFFT says enough distinct alleles AND the neighborhood is
  a closed bubble**, not until the first complete walk: a single complete walk often has
  many decoration-variant siblings (different repeat-arm choices, overlap-variant L-lines).
  Stopping at the first complete walk would almost guarantee missing the second allele.
  The widen loop uses MAFFT-core clustering (id ≥ `--dup-id`, frac ≥ `--dup-frac`) plus
  `classify_neighborhood_topology` returning `closed_bubble` as the joint stopping
  criterion; the loop only breaks when BOTH the cluster count reaches `expected_count`
  AND the BFS neighborhood has converged to a closed bubble (else widen).
- **MSA-with-reference for dedup, batched on large pools.** `cluster_alleles` aligns
  trimmed candidates together with the trimmed locus reference (not pairwise). For pools
  of size > `batch_size` (default 10), processing is BATCHED: each round runs ONE MAFFT
  MSA on (carryover cluster-rep survivors + next `batch_size` new candidates), and
  survivors propagate into the next round. Below the batch threshold the function falls
  through to a single MAFFT call. The HD-core columns are defined by the reference's
  HD-protein span projected to MSA columns — a single fixed set within each batch's MSA,
  shared across all candidate pairs in that batch. Identity is a simple character
  comparison on those columns. The flanks anchor the alignment frame but don't contribute
  to identity, so two distinct HD alleles at ~70-85 % HD-core identity get correctly
  separated even though their full-locus identity is ~88-90 %.
- **detect_core_span is memoized in-process**, so the same candidate sequence is
  tblastn'd at most once per run regardless of how many cluster_alleles invocations
  happen across widen-loop hops.
- **Coverage-only high-copy detection** (`--asymmetric-bfs` mode). Segments whose
  depth exceeds `cov_repeat_factor × D` are flagged. In asymmetric mode the BFS
  may VISIT them but not expand THROUGH them, dampening combinatorial blow-up at
  high-copy hubs. Symmetric BFS is the default.
- **RC-mirror dedup inside enumerate_paths.** A walk and its reverse-traversal mirror
  represent the same physical sequence. Dropping the duplicate at emit time (safe-form:
  only when both members of a mirror pair are actually enumerated) keeps the DFS output
  free of trivial doubles in symmetric-graph cases.
- **No length condition in completeness.** Length matters in the picker as a tiebreaker,
  not in declaring whether a candidate "is an allele".

---

## 4. Allele selection  (`pick_alleles.py`, skipped under `--skip-pick`)

### Picker-input pool (emitted by step 3, `_emit_picker_candidates`)

At the end of `graph_path_search.run` the module writes
`picker_candidates.fasta` + `picker_candidates.ann.tsv` =
filtered anchor contigs ∪ all bubble paths. The filter (in
`matdetangler/graph_path_search.py::_emit_picker_candidates`) drops any anchor
contig whose per-k `anchor_contig_<k>.ann.tsv` row shows 0 variable genes AND
lacks at least one of (flankL, flankR). Those flankL-only / flankR-only /
pure-junk anchors are needed by step 3 for BFS seeding but would only form
`(var_compl=0, flank_compl≤1)` bottom-of-tuple pairs at the picker. Bubble
paths are end-to-end walks by construction and pass through unconditionally.
The wrapper just copies the two files to `cand_combined.{fasta,ann.tsv}` for
back-compat naming and hands them to the picker.

Note: when `--skip-pick` is on (the current default), this step and the
per-sample summary row (wrapper step 8) are skipped. Wrapper steps 5, 6, 7 still
run on `bubble_alleles.fasta` aliased as `primary_alleles.fasta` so you get
per-candidate consensus, identity, and QC without the picker collapsing to a
pair. `--no-skip-pick` re-enables the full pipeline.



**Two-level pair-based selection (2026-05-30 rewrite; completeness-first
re-ranking 2026-05-31).** Level 1 picks the best divergent pair within each k;
Level 2 picks the K whose pair wins on the same score tuple. Completeness ranks
ABOVE topology so that a pair where each allele carries the variable genes
beats a topologically-cleaner pair where one allele is missing them.

```
# 4.1 resolve every candidate's GFA segment list
#     - path-origin (bubble_alleles.fasta): segments already in ann.tsv col 3
#     - bubble-origin (anchor_contig.fasta): look up the contig name in
#       <spades_dir>/<k>/contigs.paths to get its SPAdes-paired-end-resolved walk

# 4.2 load per-segment labels (written by step 3 to seg_labels_<k>.tsv)
#     labels[seg_id] = "+"-joined feature tags ("flankL", "flankR+HD1", ...)

# 4.3 pair-level scoring — sequential / lexicographic tuple, highest beats lower
def pair_score(a, b):
    ids_a, ids_b = seg_ids_of(a), seg_ids_of(b)
    shared = ids_a ∩ ids_b
    flankL_shared = any(label[s] contains "flankL"  for s ∈ shared)
    flankR_shared = any(label[s] contains "flankR"  for s ∈ shared)
    topology = (3 if flankL_shared ∧ flankR_shared   else  # closed_bubble
                2 if flankL_shared ∨ flankR_shared   else  # open_bubble
                1 if shared                          else  # detached
                0)                                          # other (no shared)

    # 4.3a per-allele completeness ranks — encoded as lo*10 + hi so that
    # balanced pairs (each allele has 1 gene/flank) outrank lopsided pairs
    # (one has 2, other has 0). With nvar_total=2 the rank table is:
    #   (2,2)=22, (2,1)=12, (1,1)=11, (2,0)=2, (1,0)=1, (0,0)=0
    nvars_a, nvars_b = #variable_genes_in(a), #variable_genes_in(b)
    var_compl = sorted((nvars_a, nvars_b))[0]*10 + sorted((nvars_a, nvars_b))[1]

    flanks_a = (a in hL) + (a in hR);  flanks_b = (b in hL) + (b in hR)
    flank_compl = sorted((flanks_a, flanks_b))[0]*10 + sorted((flanks_a, flanks_b))[1]

    # one MAFFT MSA per completeness layer (see 4.4) — gives HD-core columns
    # projected from the locus-ref's ungapped HD-protein-hit span. id and
    # aln_frac scored ONLY over those columns.
    alnid, aln_frac = hd_core_pair_identity_from_MSA(a, b)

    return ( var_compl,                                 # PRIMARY — each allele's var-gene completeness
             topology,                                  # closed_bubble > open > detached > other
             flank_compl,                               # each allele's flank completeness
             -alnid,                                    # MORE divergent HD-core > less (tiebreak)
             aa_cov(a) + aa_cov(b),                     # joint HD aa-coverage
             flank_bp_cov(a) + flank_bp_cov(b),         # joint flank bp-coverage
             -abs(len(a) - len(b)) )                    # similar-length tiebreak

# 4.4 best pair within k — completeness-first optimization
def best_pair_in_k(cands):
    # cheap: compute var_compl rank for every C(N,2) pair (just counting)
    pair_vc = [(a, b, var_completeness(a, b)) for a, b ∈ pairs(cands)]

    # descend var_compl layers from highest. Only run the EXPENSIVE HD-core MAFFT
    # MSA on candidates in the leading layer. If all those pairs are _is_dup,
    # descend to the next layer and repeat. (Was topology-layered before
    # 2026-05-31; with completeness promoted above topology, layering must
    # follow the new primary key.)
    for layer ∈ sorted(unique_var_compl_ranks, descending):
        pairs_in_layer = [(a, b) for a, b, v ∈ pair_vc if v == layer]
        subset = union of candidates in pairs_in_layer
        hd_pair_id = ONE_MAFFT_MSA(subset, locus_ref, proteins)  # _align_cores

        for (a, b) ∈ pairs_in_layer:
            alnid, aln_frac = hd_pair_id[(a, b)]
            if alnid ≥ DUP_PID ∧ aln_frac ≥ DUP_AF: continue       # dup — drop
            score = pair_score(a, b)
            best ← argmax score
        if best: return best
    return None    # no divergent pair at any completeness layer

# 4.5 Level 2: cross-K winner
pair_by_k = { k: best_pair_in_k(by_k[k]) for k ∈ ks if by_k[k] }
best_k = argmax_k pair_by_k[k].score
picks  = pair_by_k[best_k]

# 4.6 Fallback: no valid divergent pair in any K
#     -> return 1 allele (highest aa_cov + flank_bp_cov single)

# 4.7 Outputs
#     primary_alleles.fasta — picks (1 or 2 records); headers carry
#       sample/allele/origin/k/type/len/from/cov tags.
#     picks.tsv — 12-col, with header:
#       sample  allele  origin  k  type  len  from_contig  segments  cov
#       n_variable_genes  has_both_flanks  is_degHD
#       (`is_degHD` is a back-compat placeholder column.)
#       Constants: DUP_PID = 0.95, DUP_AF = 0.80.
```

### Rationale

- **Pair-level instead of per-candidate.** Earlier picker scored each candidate
  independently then tried (top1, top2) etc. for divergence. The new picker
  evaluates every pair directly — pair_topology, joint flank completeness, and
  pair divergence are inherently pair properties.
- **Completeness above topology (2026-05-31).** A closed-bubble pair where one
  arm is just a flank fragment (e.g. AJB36's old pick: a HD1-only contig + a
  flankL-only fragment) shouldn't beat an "other" topology pair where both
  alleles carry the variable genes individually. Promoting `var_compl` above
  `topology` fixed cases like AJB36 (now picks a 2/2-vars, both-flanks contig
  at single-allele depth instead of two collapsed 2× fragments) without
  changing the picks of samples whose top topology was already CB with
  fully-complete alleles.
- **HD-core MAFFT identity, not whole-allele.** Whole-allele identity is
  dominated by shared flanks (alleles in a dikaryon ARE the same flanks ± HD
  divergence). HD-core columns isolate the locus-specific divergence. Matches
  step 3's `cluster_alleles` — one consistent metric across the pipeline.
- **Completeness-first optimization.** `var_compl` is just counting (essentially
  free) and is the primary score key. Compute it for all pairs, identify the
  completeness-leading layer, run MAFFT MSA only on candidates in that layer.
  Most samples have a clear leader with 2-4 candidates, so the MSA is tiny.
  Falls back to lower completeness layers only if the top layer is all dups.
- **Cross-K independent of step 3 early-break.** The wrapper's per-K outer loop
  now runs EVERY k (no early break) so every k's candidates land in the
  combined pool and Level 2 has real choices.
- **Fallback to 1 allele.** When no valid divergent pair exists in any K
  (e.g. haploid sample), return the single allele with the highest aa+flank
  coverage. Prevents zero-record `primary_alleles.fasta`.

### 4.5 Post-pick orientation + outer-flank trim

After the pair (or singleton) is picked, two passes touch each record's walk
before `primary_alleles.fasta` is written:

```
# 4.5a normalize orientation (pick_alleles._normalize_orientation)
#     for each picked allele's segment walk:
#       L_pos  = mean walk index of segments whose label contains "flankL"
#       R_pos  = mean walk index of segments whose label contains "flankR"
#       if L_pos > R_pos:  walk = reverse_complement(walk); flip the sequence too
#     composite labels ("HD1+flankL") count toward both flank classes.
#     guarantees all alleles point flankL → variable genes → flankR consistently
#     before any downstream tool reads them.

# 4.5b outer-flank trim (pick_alleles._trim_outer_flanks)
def is_flank_only(seg):
    toks = labels[seg].split("+")
    return toks and all(t.startswith("flank") for t in toks)   # no var gene

def has_var(seg):
    return any(not t.startswith("flank") for t in labels[seg].split("+"))

# PAIR-pick mode (two picks present): trim only outer tails so a flanks-only-
# no-HD allele can never lose its content. For each arm independently:
#     walk = picked allele's segment walk
#     first_var, last_var = first / last indices of has_var(walk[i])
#     if neither arm has a var-gene segment: skip trim for this arm
#     for the LEFT tail walk[:first_var]:
#         find innermost flank-only seg whose ID also appears in the OTHER arm's
#         walk (shared = definite homology anchor)
#         if found at index i:  trim_left = i
#     for the RIGHT tail walk[last_var+1:]:
#         find innermost shared flank-only seg, set trim_right = i
#     new_walk = walk[trim_left : trim_right+1]

# SINGLETON-pick mode (one pick): aggressive — collapse to the span between
# the leftmost and rightmost FLANK-tagged segments. Drops everything beyond
# the outermost flank anchors even if var genes sit in the tails. Picks for
# the C1-1 / C1-4 / C3-1 haploid samples exercise this branch.
#     first_flank = min i where any token in labels[walk[i]] starts with "flank"
#     last_flank  = max i where same
#     new_walk = walk[first_flank : last_flank+1]

# apply_trim (shared): reconstruct nucleotide from GFA segments + L-line overlaps.
# SANITY GUARD: if the reconstructed slice length > original contig length,
# the SPAdes contig spans only a subset of large GFA segments — ABORT trim
# and keep the contig intact. (Without the guard, trimming C1-1 inflated
# 8843 bp → 156981 bp because two of the picked walk's segments are megabase
# scaffolds.) The contig is logged but not modified.
```

The trim is the only step that mutates `primary_alleles.fasta` *sequence* after
selection — every prior step writes whole assembled contigs / DFS-reconstructed
paths. picks.tsv's `segments` column always reflects the post-trim walk for
downstream `graph_paths` rendering.

---

## 5. Annotated allele walks  (`graph_paths.py`)

Each primary-alleles record has a recorded segment walk (from §3's ann.tsv).
graph_paths reads the walk verbatim, labels each segment by content (tblastn proteins
+ blastn flanks), and emits:

Inputs the wrapper supplies: `--primary-alleles`, `--cand-ann`
(`bubble_alleles.ann.tsv` from step 3), `--spades-dir`, and `--picks-tsv`
(absent under `--skip-pick`; the `--cand-ann` fallback resolves per-allele
source-k GFA + segment walk in that case via `_paths_from_cand_ann`).

### Performance: labels only the segments in the recorded paths

Earlier versions ran `label_nodes` on **every segment of the whole GFA** (100k-
300k segs per call) for each allele in `primary_alleles.fasta`. With 24 candidate
alleles under `--skip-pick`, that was 24× redundant whole-GFA labeling — step 5
could take 30-60 minutes per sample.

Now `graph_paths.run` (via the `_per_allele_gfas_and_paths` /
`_paths_from_cand_ann` helpers + a per-GFA cache):
- **per-GFA cache**: parses each unique GFA only once, even when multiple alleles
  share it (the `gfa_cache` dict inside `run`).
- **recorded-path subset**: when `--cand-ann` (step-3 ann.tsv) or `--picks-tsv`
  is supplied, the union of all alleles' recorded segments is the only set of
  segments labeled — ~tens of segments instead of 100k-300k. blast queries run
  on tiny temp DBs.

Combined effect: step 5 now completes in **~0.5-2 seconds** (zero new blast
work when the cache is warm; all hits are USE).

| output | what |
|---|---|
| `bubble.txt` | `allele1: <walk>` / `allele2: <walk>` (or one row per candidate under `--skip-pick`) |
| `bubble.gfa` | sub-GFA of just the walks' segments + L-links (loads in Bandage) |
| `bubble.dot` | Graphviz DOT |
| `bubble.tsv` | edge list |
| `bubble.png` | matplotlib: one row per allele. **Node coloring**: yellow = carries any variable gene (e.g. HD1/HD2; flank tag, if any, ignored for color); blue = flank-only (label is purely flankL/flankR, no variable gene); white = pure-number / unlabeled. **Dashed gray cross-arm lines**: (a) one per GFA segment ID shared between the two arms — same node = definite homology; (b) fallback only when no flank-only segment ID is shared for a given flank type — one extra line connects the outermost flank-only node of each arm (leftmost flankL = arm entry, rightmost flankR = arm exit). |

When both alleles came from the same K, node IDs are bare. When they came from different
K's (cross-K pair), node IDs are prefixed `K33:` / `K55:` etc. so the user can tell at a
glance which graph each allele lives in.

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

The pick-level numbers in §4-6 are based on the graph-derived allele sequences. The
read-derived consensus can drift (low-coverage stretches become N's, etc.), so we
**re-check** divergence and completeness on the consensus. The wrapper runs three
sub-steps:

- **7a — pairwise identity on PICKS.** `pairwise_identity` on
  `primary_alleles.fasta`. Same HD-core MAFFT protocol the picker uses (§4)
  via `_align_cores`: trim each candidate + locus reference to HD-core ± pad,
  ONE MAFFT alignment, score identity only over HD-core columns projected
  from the reference's ungapped HD-protein span. Writes `identity.tsv`.
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
| Decoration-variant paths cluttering the candidates | `--dup-id` (default 0.95), `--dup-frac` (default 0.80) | MAFFT-core clustering thresholds; raising them keeps more variants as distinct, lowering merges more |
| BFS blowing up at a high-copy hub | `--asymmetric-bfs` (off by default; module-level `cov_repeat_factor` defaults to 2.0) | turns on the absorbing-repeat rule: repeats stay in the neighborhood but don't expand their own neighbors |
| Deep search time blowing up | `--max-hops` (default 12), module-level `--max-paths` (default 50000) | tighter caps |
| MAFFT widen-loop quitting too early on false-distinct walks | `--expected-count` and the MAFFT-core check | the loop only breaks when MAFFT clusters >= expected_count truly distinct alleles AND the neighborhood topology converges to a closed bubble |
| Blast wall-time on iterative dev | per-sample `blast_*.tsv` cache (on by default; wipe with `--re-blast`) | first run pays, subsequent runs replay from disk; wipe whenever blast cutoffs change |

## Failure modes the algorithm protects against

| failure | guard |
|---|---|
| Intronic stops in user proteins | rejected at input by `_validate_fasta` (expected to be protein, not NT) |
| Duplicate protein IDs (blast unique-sseqid constraint) | sanitizer suffixes them `_1`, `_2` |
| HD envelope shrunken because tblastn missed divergent ends | union over ALL HSPs + padding |
| Flank too long → search blowup | `--max-flank-len` trim |
| Decoration-variant paths winning the pick | MAFFT-core clustering in widen loop (`cluster_alleles` collapses decoration-variants into one cluster) + picker tiebreaks on coverage and length |
| Reverse-complement of a walk emitted as a separate candidate | `mirror_path` RC-mirror dedup at emit time inside `enumerate_paths` (safe form: drops only when both members of mirror pair enumerated) |
| Picks not actually divergent | per-K MAFFT `_is_dup` filter |
| Consensus drift from pick after read mapping | `consensus_qc.py` re-check + consensus divergence pass |
| Read-derived gap in HD body | core_meandepth at the HD-core positions |
| Cross-k pair makes a confusing bubble | per-K pair finder + bubble PNG K-prefix when cross-k |
