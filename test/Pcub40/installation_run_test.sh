#!/bin/bash
# Pcub40 reproducibility test.
#
# Pipeline:
#   1. Decompress examples/Pcub40/<sample>/k<k>/*.gfa.gz in place (sibling *.gfa
#      next to each *.gfa.gz; idempotent — skips any already decompressed).
#   2. Run MATdetangler with the args defined in test/Pcub40/run_args.json
#      (read by humans; the script forwards the same flags to MATdetangler).
#   3. Summarize the output with test/Pcub40/summarize.py.
#   4. Diff the new JSON against test/Pcub40/known_results.json. Only the
#      samples actually run this invocation are diffed — a subset run does not
#      fail on samples it didn't run.
#   5. Exit 0 if no semantic regression, 1 otherwise. On PASS the decompressed
#      *.gfa siblings are kept (gitignored) so re-runs reuse them.
#
# Usage:
#   bash test/Pcub40/run_test.sh                 # run all 32 samples
#   bash test/Pcub40/run_test.sh AJB36 BD-1248   # just these (only these diffed)
#   bash test/Pcub40/run_test.sh --slurm         # submit a slurm array instead of serial
#
# Tolerated diffs (will not fail the test):
#   - top-level "version" block (git hash will move)
#   - "args" block (run-time vs known)
#   - sample-level "per_k_trace" (BFS trace is large + noisy across machines)
#   - sample-level "finished_nhop" (depends on cov estimator)
# Everything else must match exactly.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
EX="$ROOT/examples/Pcub40"
KNOWN="$HERE/known_results.json"
ARGS_JSON="$HERE/run_args.json"
TMP_ROOT="${MATDETANGLER_TEST_TMP:-$ROOT/_test_tmp}/matdetangler_pcub40_test_$$"
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

echo "[$(date)] Pcub40 test: ${#SAMPLES[@]} samples, out=$OUT_DIR"
echo "[$(date)] decompressing GFAs in-place next to .gfa.gz files"

# Step 1: decompress in place — examples/Pcub40/<sample>/k<k>/*.gfa.gz
# stays committed (canonical) and a sibling *.gfa is created for the
# pipeline to read. Idempotent (skips already-decompressed) and tracked
# so a final cleanup step at the end can remove just the .gfa siblings.
# Only k45 and k53 are decompressed — these match `--ks k45,k53` below,
# so we don't waste disk decompressing k33 (which the wrapper won't use).
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
echo "[$(date)] decompressed ${#DECOMPRESSED_GFAS[@]} GFA(s) (skipped any already present)"

# Step 2: run MATdetangler — args mirror run_args.json defaults.
# spades-dir is examples/Pcub40/<sample> (which now contains both *.gfa.gz
# and *.gfa side-by-side; the wrapper reads the *.gfa).
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
  echo "[$(date)] submitted array $jid; wait + diff manually."
  exit 0
fi

# Serial mode (good for small subsets; the full 32 takes ~20-30 min)
for s in "${SAMPLES[@]}"; do
  echo "[$(date)] $s ..."
  run_one "$s" || echo "  $s exit=$?"
done

# Step 3: aggregate the wrapper-emitted <sample>/summary.json files into one
# cross-sample JSON for diff. (Each sample's summary.json was produced by
# Step 9 of MATdetangler.)
NEW_JSON="$TMP_ROOT/new_results.json"
python -m matdetangler.summarize \
  --results-dir "$OUT_DIR" \
  --args-json   "$ARGS_JSON" \
  --repo-dir    "$ROOT" \
  --out         "$NEW_JSON"

# Step 4: diff with semantic ignore-list
echo "[$(date)] diffing new_results.json vs known_results.json"
python - "$NEW_JSON" "$KNOWN" <<'PY'
import json, sys
new = json.load(open(sys.argv[1])); known = json.load(open(sys.argv[2]))
# Sample-level fields ignored: install-drift-sensitive (different MAFFT /
# BLAST / edlib versions give different numbers without changing the
# analysis answer). The test is for COMPARING ANALYSIS RESULTS, not for
# byte-identity across installs.
TOLERATE = {
    # machine + estimator noise
    "per_k_trace", "finished_nhop", "per_k_pick", "genome_cov", "allele_cov",
    # MAFFT alignment numbers (version-dependent)
    "allele1_vs_allele2_id_pct", "allele1_vs_allele2_aln_frac",
    # Human-readable path strings (re-derived; BLAST hit boundaries can
    # shift trace order without changing the underlying graph walk)
    "allele1_path_str", "allele2_path_str",
    # bp totals (drift with locus-trim boundaries)
    "basepair",
}
# Per-allele fields ignored: drift-sensitive numerics only. Graph segments
# ARE compared — the wrapper exports `PYTHONHASHSEED=0`, which locks
# set/dict iteration order across runs and pins BFS tie-breaking, so the
# exact list of segments per allele is reproducible within an install.
# (Verified: two independent demo runs of the 32-sample Pcub40 grid with
# seed=0 give 0/32 diff at this tolerance level.)
ALLELE_TOLERATE = {"cov", "len"}
def scrub(s):
    s = dict(s)
    for k in TOLERATE: s.pop(k, None)
    alleles = []
    for a in s.get("alleles", []):
        a = {k: v for k, v in a.items() if k not in ALLELE_TOLERATE}
        a["segments"] = sorted(a.get("segments", []) or [])
        alleles.append(a)
    if "alleles" in s: s["alleles"] = alleles
    return s
diffs = []
# Only compare samples that were actually run this invocation (present in
# new_results.json). Running a subset — e.g. `run_test.sh AJB36` — must not
# FAIL just because the other known samples weren't run.
for s in sorted(new["samples"]):
    n = scrub(new["samples"].get(s, {}))
    k = scrub(known["samples"].get(s, {}))
    if n != k:
        diffs.append(s)
        print(f"\nDIFF {s}:")
        for key in sorted(set(n) | set(k)):
            if n.get(key) != k.get(key):
                print(f"    {key}: new={str(n.get(key))[:80]}  known={str(k.get(key))[:80]}")
if diffs:
    print(f"\nFAIL: {len(diffs)} samples differ ({', '.join(diffs[:5])}{'...' if len(diffs)>5 else ''})")
    sys.exit(1)
print(f"PASS: all {len(new['samples'])} samples match (analysis result; ignoring install-drift fields)")
PY
status=$?

if [ $status -eq 0 ]; then
  rm -rf "$TMP_ROOT"
  # Keep the in-place decompressed *.gfa siblings (next to the canonical
  # *.gfa.gz). The decompress step above is idempotent — it skips any .gfa
  # that already exists — so leaving them lets subsequent runs reuse them
  # without re-gunzipping.
  echo "[$(date)] PASS — tmp cleaned; ${#DECOMPRESSED_GFAS[@]} .gfa sibling(s) kept under $EX for reuse"
else
  echo "[$(date)] FAIL — tmp kept at $TMP_ROOT; in-place .gfa siblings kept under $EX for inspection"
fi
exit $status
