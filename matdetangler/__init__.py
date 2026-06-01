"""MATdetangler — recover the two alleles of a tetrapolar mating-type-like locus from short reads.

Given pre-built SPAdes graphs, a locus-only reference, and HD protein queries, pick the two alleles
per sample via contig-level anchor search + GFA path enumeration (flank-to-flank walks), validate
completeness, detect a degenerated copy if a reference is supplied, and emit allele fastas, an
annotated bubble graph (ASCII + GFA + DOT + PNG), and a per-sample summary. Optionally maps reads
back to the picks for samtools-consensus QC (--no-skip-consensus)."""
__version__ = "0.1.0"
