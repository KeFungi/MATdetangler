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

# Verdict-level fields — when these change, the analysis answer changed.
VERDICT_FIELDS = ["bubble_type", "k_chosen", "complete_var", "complete_locus",
                   "n_dedup", "locus_coverage"]
# Per-allele biological fields — when these change, the allele identity
# changed (different flank/HD content or different source k).
ALLELE_BIO_FIELDS = ["has_both_flanks", "is_degHD", "n_variable_genes",
                      "k", "type", "origin"]

verdict_changes = []
completeness_changes = []
allele_structure_changes = []
segment_drift = []
identical = []

for s in sorted(set(new["samples"]) | set(known["samples"])):
    if s not in new["samples"]:
        verdict_changes.append((s, "REMOVED", None, None)); continue
    if s not in known["samples"]:
        verdict_changes.append((s, "NEW", None, None)); continue
    n = new["samples"][s]; k = known["samples"][s]

    # Verdict-tier diffs
    vd = [(f, k.get(f), n.get(f)) for f in VERDICT_FIELDS if k.get(f) != n.get(f)]
    if vd:
        # Split: bubble_type / n_dedup / k_chosen → "verdict change"
        # complete_var / complete_locus / locus_coverage → "completeness change"
        bubble_diffs = [f for f, _, _ in vd if f in ("bubble_type", "k_chosen", "n_dedup")]
        compl_diffs  = [f for f, _, _ in vd if f in ("complete_var", "complete_locus", "locus_coverage")]
        if bubble_diffs:
            verdict_changes.append((s, "VERDICT", bubble_diffs, vd))
        if compl_diffs:
            completeness_changes.append((s, compl_diffs, vd))
        continue

    # Allele-structure diffs (per-allele biology)
    n_a = n.get("alleles", []) or []
    k_a = k.get("alleles", []) or []
    bio_changed = False
    if len(n_a) != len(k_a):
        bio_changed = True
    else:
        for an, ak in zip(n_a, k_a):
            for f in ALLELE_BIO_FIELDS:
                if an.get(f) != ak.get(f):
                    bio_changed = True
                    break
            if bio_changed: break
    if bio_changed:
        allele_structure_changes.append(s); continue

    # Segment-set drift (graph-isomorphic walks)
    seg_n = [sorted(a.get("segments", []) or []) for a in n_a]
    seg_k = [sorted(a.get("segments", []) or []) for a in k_a]
    if seg_n != seg_k:
        segment_drift.append(s); continue

    identical.append(s)

print(f"\n========== Pcub40 ANALYSIS REPORT ==========")
print(f"compared {len(new['samples'])} samples vs {len(known['samples'])} in baseline")
print()
print(f"identical (bit-for-bit on biology + segments): {len(identical)}")
print(f"  {', '.join(identical) if identical else '(none)'}")
print()
print(f"VERDICT changes ({len(verdict_changes)})  — bubble_type / k_chosen / n_dedup shifted:")
for s, tag, fields, vd in verdict_changes:
    print(f"  {s} [{tag}]: {fields}")
    if vd:
        for f, kv, nv in vd:
            print(f"      {f}: {kv} -> {nv}")
print()
print(f"COMPLETENESS changes ({len(completeness_changes)})  — complete_var/locus/lc shifted, same verdict:")
for s, fields, vd in completeness_changes:
    print(f"  {s}: {fields}")
    for f, kv, nv in vd:
        if f in fields: print(f"      {f}: {kv} -> {nv}")
print()
print(f"ALLELE-structure changes ({len(allele_structure_changes)})  — has_both_flanks / k / etc:")
for s in allele_structure_changes:
    print(f"  {s}")
print()
print(f"Segment-set DRIFT only ({len(segment_drift)})  — same biology, different graph walk:")
print(f"  {', '.join(segment_drift) if segment_drift else '(none)'}")
print()
print(f"==============================================")
PY

# Always exit 0 — this is informational
rm -rf "$TMP_ROOT" 2>/dev/null || true
exit 0
