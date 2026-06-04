#!/bin/bash
# Competitive mapping of trimmed reads to the picked alleles, then samtools consensus per allele.
#
# Usage:
#   map_consensus.sh  <SAMPLE> <PRIMARY_ALLELES_FA> <R1.fq.gz> <R2.fq.gz> <OUTDIR> [THREADS] [QUERIES_DIR] [REPEATS_FA]
#
# Writes (in OUTDIR):
#   reads.sam                                     — competitive end-to-end mapping (-k 1, --no-unal)
#                                                    [BAM is an internal intermediate, deleted at the end]
#   consensus_alleles.fasta                       — samtools consensus per allele, all in one multi-record FASTA
#                                                    (record IDs preserved from the input reference's
#                                                    primary_alleles.fasta — same names as the BAM @SQ entries)
#   coverage.tsv                                  — per-allele whole + HD-core (repeat-masked) depth
set -uo pipefail
S="$1"; FA="$2"; R1="$3"; R2="$4"; OD="$5"; NT="${6:-4}"; QDIR="${7:-}"; REPEATS="${8:-}"
mkdir -p "$OD"
IDX="$OD/idx"
bowtie2-build --quiet "$FA" "$IDX"
bowtie2 --end-to-end --very-sensitive -k 1 --no-unal -q -p "$NT" -x "$IDX" \
        -1 "$R1" -2 "$R2" -S "$OD/reads.sam" 2> "$OD/bt2.log"
samtools sort -@ "$NT" -o "$OD/reads.sorted.bam" "$OD/reads.sam"
samtools index "$OD/reads.sorted.bam"
# coverage per allele: whole-allele AND HD-core (positions inside variable-gene tblastn hits,
# with --repeats positions masked out). The HD-core depth is the biologically interpretable
# number — a high-copy repeat (e.g. MITE) inside one allele inflates whole-allele depth 2-3x.
if [ -n "${QDIR:-}" ] && [ -f "${QDIR}/variable_proteins.fasta" ]; then
  python3 -m matdetangler.coverage_core --bam "$OD/reads.sorted.bam" --alleles "$FA" \
    --queries-dir "$QDIR" ${REPEATS:+--repeats "$REPEATS"} \
    --out-tsv "$OD/coverage.tsv"
else
  { samtools coverage "$OD/reads.sorted.bam" \
    | awk 'NR==1{print "#allele\tlen\tmapped_reads\twhole_breadth_pct\twhole_meandepth\tcore_bp\tcore_meandepth\trepeat_bp"} NR>1{print $1"\t"$3"\t"$4"\t"$6"\t"$7"\t-\t-\t-"}'
  } > "$OD/coverage.tsv"
fi
# per-allele consensus, all records concatenated into one consensus_alleles.fasta.
# Record IDs are PRESERVED from the input reference — same names as the BAM @SQ
# entries (e.g. SAMPLE_kNN_allele1). This makes the consensus FASTA usable as a
# drop-in replacement reference without remapping name conventions.
CONS="$OD/consensus_alleles.fasta"
: > "$CONS"
for al in $(grep '^>' "$FA" | awk '{print $1}' | sed 's/^>//'); do
  tmp="$OD/_cons.fa"
  samtools consensus -r "$al" -f fasta -o "$tmp" "$OD/reads.sorted.bam" 2>/dev/null
  cat "$tmp" >> "$CONS"
  rm -f "$tmp"
done
# Persist SAM (human-readable), drop the intermediate BAM(+index).
rm -f "$IDX"*.bt2 "$OD/reads.sorted.bam" "$OD/reads.sorted.bam.bai"
echo "[map_consensus] $S done -> $OD/"
