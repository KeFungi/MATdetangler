#!/bin/bash
# Pcub40 ANALYSIS test — assess how the CURRENT implementation differs
# biologically from the committed baseline (test/Pcub40/known_results.json).
#
# Unlike installation_run_test.sh (which exits 1 on any divergence and is
# meant for verifying that the same code + env reproduces the baseline
# byte-for-byte), this script:
#   * focuses on biologically meaningful fields
#   * structures the diff as a REPORT (not pass/fail)
#   * ALWAYS exits 0 — divergence is INFORMATION, not error
#
# Use case: you've changed the BFS, the picker, the classifier, a default
# threshold, etc. You want to know "which Pcub40 samples got called
# differently?" without the failure noise of strict-segment comparison.
#
# Pipeline:
#   1. Decompress examples/Pcub40/<sample>/k{45,53}/*.gfa.gz in place.
#   2. Run MATdetangler with the canonical args (run_args.json).
#   3. Aggregate per-sample summary.json (Step 9) into a cross-sample JSON.
#   4. Compare to known_results.json on:
#        - VERDICT-level: bubble_type, k_chosen, complete_var, complete_locus,
#          n_dedup, locus_coverage.
#        - ALLELE-level: has_both_flanks, is_degHD, n_variable_genes, k,
#          type, origin.
#        - SEGMENT drift: set-equality of segments per allele.
#   5. Print a structured report:
#        - identical samples
#        - verdict changes
#        - completeness changes
#        - allele-structure changes
#        - segment-set drift only (likely graph-isomorphic walks)
#
# Usage:
#   bash test/Pcub40/analysis_run_test.sh                 # all 32 samples
#   bash test/Pcub40/analysis_run_test.sh AJB36 BD-1248   # subset
#   bash test/Pcub40/analysis_run_test.sh --slurm         # SLURM array, exits at submit
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
EX="$ROOT/examples/Pcub40"
KNOWN="$HERE/known_results.json"
ARGS_JSON="$HERE/run_args.json"
TMP_ROOT="${MATDETANGLER_TEST_TMP:-$ROOT/_test_tmp}/analysis_$$"
OUT_DIR="$TMP_ROOT/results"
mkdir -p "$OUT_DIR"

if [ ! -d "$EX" ]; then
  echo "ERROR: $EX missing. Did you run 'git lfs pull' to fetch the demo GFAs?" >&2
  exit 2
fi

SAMPLES=()
SLURM=0
for a in "$@"; do
  case "$a" in
    --slurm) SLURM=1 ;;
    -*) echo "unknown flag: $a" >&2; exit 2 ;;
    *) SAMPLES+=("$a") ;;
  esac
done
if [ ${#SAMPLES[@]} -eq 0 ]; then
  while IFS= read -r s; do SAMPLES+=("$s"); done < <(ls "$EX" | grep -v '^_')
fi

echo "[$(date)] Pcub40 ANALYSIS test: ${#SAMPLES[@]} samples, out=$OUT_DIR"

# Decompress only k45 + k53 (matches --ks below)
declare -a DECOMPRESSED_GFAS
for s in "${SAMPLES[@]}"; do
  for K in k45 k53; do
    kgz="$EX/$s/$K/assembly_graph_after_simplification.gfa.gz"
    [ -s "$kgz" ] || continue
    plain="${kgz%.gz}"
    if [ ! -s "$plain" ]; then
      gunzip -c "$kgz" > "$plain"
      DECOMPRESSED_GFAS+=("$plain")
    fi
  done
done

run_one() {
  local s="$1"
  "$ROOT/MATdetangler-cli" run \
    --sample "$s" --spades-dir "$EX/$s" \
    --locus-ref "$ROOT/examples/Pcub_locus/NC_062999.fasta" \
    --proteins  "$ROOT/examples/Pcub_locus/NC_062999_HDs.fasta" \
    --outdir "$OUT_DIR" \
    --ks k45,k53 --threads 4 --expected-count 2 --no-skip-pick \
    > "$OUT_DIR/$s.run.log" 2>&1
}

if [ "$SLURM" -eq 1 ]; then
  echo "[$(date)] SLURM mode: emitting array sbatch"
  cat > "$TMP_ROOT/_array.sbatch" <<EOF
#!/bin/bash
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-$((${#SAMPLES[@]} - 1))
#SBATCH --output=$TMP_ROOT/_%a.log
set -uo pipefail
SAMPLES=(${SAMPLES[*]})
S="\${SAMPLES[\$SLURM_ARRAY_TASK_ID]}"
source /home/yihongke/miniconda3/etc/profile.d/conda.sh; conda activate MATdetangler
"$ROOT/MATdetangler-cli" run --sample "\$S" --spades-dir "$EX/\$S" \\
  --locus-ref "$ROOT/examples/Pcub_locus/NC_062999.fasta" \\
  --proteins  "$ROOT/examples/Pcub_locus/NC_062999_HDs.fasta" \\
  --outdir "$OUT_DIR" --ks k45,k53 --threads 4 --expected-count 2 --no-skip-pick
EOF
  jid=$(sbatch --parsable "$TMP_ROOT/_array.sbatch")
  echo "[$(date)] submitted array $jid; rerun this script after it finishes to see the report."
  exit 0
fi

for s in "${SAMPLES[@]}"; do
  echo "[$(date)] $s ..."
  run_one "$s" || echo "  $s exit=$?"
done

# Aggregate + report
NEW_JSON="$TMP_ROOT/new_results.json"
python -m matdetangler.summarize \
  --results-dir "$OUT_DIR" --args-json "$ARGS_JSON" \
  --repo-dir "$ROOT" --out "$NEW_JSON"

echo
echo "[$(date)] structured biological diff vs $KNOWN"
python - "$NEW_JSON" "$KNOWN" <<'PY'
import json, sys
new = json.load(open(sys.argv[1])); known = json.load(open(sys.argv[2]))

# Each change-type catalogs samples it touched. A single sample can appear
# in multiple categories (e.g. bubble_type AND n_dedup AND bp can all
# shift together when an allele gets added/dropped).
CHANGE_TYPES = [
    "bubble_type",      # verdict label
    "k_chosen",         # which k won
    "n_dedup",          # allele count
    "complete_var",     # HD-tag completeness (tri-state)
    "complete_locus",   # flank completeness (tri-state)
    "locus_coverage",   # numerical locus coverage
    "allele_bp",        # per-allele or total bp shifted
    "allele_structure", # has_both_flanks / is_degHD / n_variable_genes / k / type / origin per allele
    "finished_nhop",    # accepting nhop (BFS depth at acceptance) shifted
    "search_nhood",     # |nhood| at the accept iteration shifted
    "search_arms",      # n_arms at the accept iteration shifted
    "bfs_limits",       # bfs_limits counters (max_paths / max_path_length / max_bp hits) shifted
    "segment_drift",    # same biology, different graph walk (set-equality of segments differs)
]
ALLELE_BIO_FIELDS = ["has_both_flanks", "is_degHD", "n_variable_genes",
                      "k", "type", "origin"]
# Run-level args (NOT per-sample — these are the global config; if these
# differ, EVERYTHING downstream is on a different footing).
ARG_FIELDS = ["max_paths", "max_path_length", "max_bp_since_var",
              "cov_filter", "init_nhop", "max_nhop", "lo_mult", "hi_mult",
              "divergence_threshold", "locus_padding", "min_allele_bp",
              "seeds", "expected_count", "ks"]

per_sample_changes = {}        # sample -> ordered list of (category, k_val, n_val)
change_buckets = {c: [] for c in CHANGE_TYPES}
identical = []

for s in sorted(set(new["samples"]) | set(known["samples"])):
    if s not in new["samples"]:
        per_sample_changes[s] = [("REMOVED_FROM_RUN", None, None)]
        continue
    if s not in known["samples"]:
        per_sample_changes[s] = [("NEW_SAMPLE", None, None)]
        continue
    n = new["samples"][s]; k = known["samples"][s]
    changes = []

    # Top-level categorical fields
    for f in ("bubble_type", "k_chosen"):
        if n.get(f) != k.get(f):
            changes.append((f, k.get(f), n.get(f)))
            change_buckets[f].append(s)
    # Numeric / tri-state
    for f in ("n_dedup", "complete_var", "complete_locus", "locus_coverage"):
        if n.get(f) != k.get(f):
            changes.append((f, k.get(f), n.get(f)))
            change_buckets[f].append(s)

    # bp: sample-level basepair OR any per-allele len
    n_a = n.get("alleles", []) or []
    k_a = k.get("alleles", []) or []
    bp_shift = False
    if n.get("basepair") != k.get("basepair"):
        bp_shift = True
    if len(n_a) == len(k_a):
        for an, ak in zip(n_a, k_a):
            if an.get("len") != ak.get("len"):
                bp_shift = True; break
    if bp_shift:
        n_bp = n.get("basepair", "?"); k_bp = k.get("basepair", "?")
        changes.append(("allele_bp", f"total={k_bp} per_allele={[a.get('len') for a in k_a]}",
                                       f"total={n_bp} per_allele={[a.get('len') for a in n_a]}"))
        change_buckets["allele_bp"].append(s)

    # Allele-structure: bool flags + per-allele identity
    struct_changed = False
    if len(n_a) == len(k_a):
        for an, ak in zip(n_a, k_a):
            for f in ALLELE_BIO_FIELDS:
                if an.get(f) != ak.get(f):
                    struct_changed = True; break
            if struct_changed: break
    if struct_changed:
        changes.append(("allele_structure",
                          [{f: a.get(f) for f in ALLELE_BIO_FIELDS} for a in k_a],
                          [{f: a.get(f) for f in ALLELE_BIO_FIELDS} for a in n_a]))
        change_buckets["allele_structure"].append(s)

    # Search-state diffs (BFS internals at the accepting iteration).
    if n.get("finished_nhop") != k.get("finished_nhop"):
        changes.append(("finished_nhop", k.get("finished_nhop"), n.get("finished_nhop")))
        change_buckets["finished_nhop"].append(s)
    # Compare per_k_trace at the chosen k — last iteration's nhood / n_arms
    n_kc = n.get("k_chosen"); k_kc = k.get("k_chosen")
    n_trace = (n.get("per_k_trace") or {}).get(n_kc) or []
    k_trace = (k.get("per_k_trace") or {}).get(k_kc) or []
    n_last = n_trace[-1] if n_trace else {}
    k_last = k_trace[-1] if k_trace else {}
    if n_last.get("nhood") != k_last.get("nhood"):
        changes.append(("search_nhood", k_last.get("nhood"), n_last.get("nhood")))
        change_buckets["search_nhood"].append(s)
    if n_last.get("n_arms") != k_last.get("n_arms"):
        changes.append(("search_arms", k_last.get("n_arms"), n_last.get("n_arms")))
        change_buckets["search_arms"].append(s)
    # BFS limits (if either side recorded them in the last iter)
    n_lim = n_last.get("bfs_limits") or {}
    k_lim = k_last.get("bfs_limits") or {}
    if any(n_lim.get(L, 0) != k_lim.get(L, 0)
            for L in ("max_paths_hit", "max_path_length_hit", "max_bp_hit")):
        changes.append(("bfs_limits", k_lim, n_lim))
        change_buckets["bfs_limits"].append(s)

    # Segment-set drift (only flag if nothing above changed AND segments differ)
    if not changes:
        seg_n = [sorted(a.get("segments", []) or []) for a in n_a]
        seg_k = [sorted(a.get("segments", []) or []) for a in k_a]
        if seg_n != seg_k:
            changes.append(("segment_drift", None, None))
            change_buckets["segment_drift"].append(s)

    if changes:
        per_sample_changes[s] = changes
    else:
        identical.append(s)

# ----- args-block (run-level) diff -----
n_args = new.get("args") or {}; k_args = known.get("args") or {}
args_diff = [(f, k_args.get(f), n_args.get(f))
              for f in ARG_FIELDS if k_args.get(f) != n_args.get(f)]

# ----- per-sample report -----
print(f"\n========== Pcub40 ANALYSIS REPORT ==========")
print(f"baseline:  {len(known['samples'])} samples in {sys.argv[2]}")
print(f"this run:  {len(new['samples'])} samples")
print()
print(f"--- RUN-LEVEL ARGS DIFF ({len(args_diff)}) ---")
if args_diff:
    print("  (these are GLOBAL; if any differ, every per-sample result is on a different config)")
    for f, kv, nv in args_diff:
        print(f"  {f}: {kv} -> {nv}")
else:
    print("  (identical — same args/config as baseline)")
print()
print(f"--- PER-SAMPLE CHANGES ({len(per_sample_changes)} samples differ) ---")
if not per_sample_changes:
    print("  (all identical)")
else:
    for s in sorted(per_sample_changes):
        cats = per_sample_changes[s]
        names = [c[0] for c in cats]
        print(f"  {s:25s}  changed: {', '.join(names)}")
        for cat, kv, nv in cats:
            if cat in ("REMOVED_FROM_RUN", "NEW_SAMPLE", "segment_drift"):
                continue
            sk = str(kv); sn = str(nv)
            if len(sk) > 70: sk = sk[:67] + "..."
            if len(sn) > 70: sn = sn[:67] + "..."
            print(f"      {cat}: {sk} -> {sn}")
print()

print(f"==============================================")
PY

# Always exit 0 — this is informational
rm -rf "$TMP_ROOT" 2>/dev/null || true
exit 0
