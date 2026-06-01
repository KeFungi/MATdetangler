#!/bin/bash
# Run all MATdetangler unit tests. Fast (no SLURM, no blast/MAFFT/samtools needed).
# Run from anywhere; uses absolute paths.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)

# Use a modern Python (>= 3.9 needed for the dict[str, ...] annotations).
# On Great Lakes the genomics module's python is 3.9.
if command -v python3 >/dev/null 2>&1; then
    PY_VER=$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
    if [ "$(echo "$PY_VER" | awk -F. '{print ($1*100 + $2)}')" -lt 309 ]; then
        if [ -f /scratch/tyjames_root/tyjames0/yihongke/Pcub/scripts/config.sh ]; then
            source /scratch/tyjames_root/tyjames0/yihongke/Pcub/scripts/config.sh
            load_modules 2>/dev/null
        fi
    fi
fi

cd "$ROOT"
fail=0
for t in test/test_input_process.py test/test_GFA_search.py test/test_pick_alleles.py test/test_graph_paths.py test/test_paths.py test/test_graph_path_search.py test/test_consensus_qc.py test/test_anchor_search.py test/test_topology_clustering.py test/test_pairwise_identity.py test/test_bfs_asymmetric.py test/test_anchor_name_parsing.py test/test_blast_utils_out_tsv.py test/test_genome_cov_from_contigs.py test/test_wrapper_cli_flags.py; do
    echo ">>> $t"
    if ! python3 "$t" 2>&1 | grep -E "^(Ran |FAIL|ERROR|OK)" ; then
        fail=1
    fi
done
echo ""
if [ "$fail" -eq 0 ]; then echo "[unit] all suites passed"; else echo "[unit] failures present"; exit 1; fi
