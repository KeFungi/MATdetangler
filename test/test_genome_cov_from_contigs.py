"""Unit tests for input_process.estimate_genome_cov_from_contigs.

The function parses SPAdes contigs.fasta headers (NODE_X_length_L_cov_F)
and returns the median cov_ value across contigs with length >= min_len.
Critical correctness guarantees:
  1. Decimal precision preserved exactly (no integer truncation).
  2. Robustly skips short contigs (length filter applied BEFORE median).
  3. Returns 0.0 cleanly when no contig is long enough.
  4. Large cov values (collapsed-repeat tips that can hit ~13e6) don't crash.
  5. Headers without the expected pattern are skipped silently.
"""
import os, sys, tempfile, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler.input_process import estimate_genome_cov_from_contigs


def _write_fasta(path: str, records: list[tuple[str, str]]) -> None:
    """records is [(header_minus_>, seq_or_empty_for_short_placeholder)]."""
    with open(path, "w") as f:
        for hdr, seq in records:
            f.write(f">{hdr}\n{seq or 'N'}\n")


class TestEstimateGenomeCov(unittest.TestCase):

    def setUp(self):
        self.td = tempfile.mkdtemp(prefix="_gcvtest_")
        self.fa = os.path.join(self.td, "contigs.fasta")

    def tearDown(self):
        import shutil; shutil.rmtree(self.td, ignore_errors=True)

    def test_decimal_precision_preserved(self):
        """`cov_68.186461` must come out as exactly 68.186461, not 68.0 or 68.19."""
        _write_fasta(self.fa, [
            ("NODE_1_length_10000_cov_68.186461", "N"),
        ])
        v = estimate_genome_cov_from_contigs(self.fa, min_len=5000)
        # Compare to many digits — full SPAdes 6-decimal precision must survive.
        self.assertAlmostEqual(v, 68.186461, places=6)
        # Defensive: not rounded to integer.
        self.assertNotEqual(v, 68.0)
        self.assertNotEqual(v, 68)

    def test_typical_dataset_returns_median(self):
        """Median across long contigs only — short ones dropped before percentile."""
        _write_fasta(self.fa, [
            ("NODE_1_length_50000_cov_50.0", "N"),
            ("NODE_2_length_40000_cov_60.0", "N"),
            ("NODE_3_length_30000_cov_70.0", "N"),
            ("NODE_4_length_500_cov_99999.0", "N"),   # noise; filtered out
            ("NODE_5_length_100_cov_42.5",  "N"),     # noise; filtered out
        ])
        # median of [50.0, 60.0, 70.0] is 60.0
        self.assertEqual(estimate_genome_cov_from_contigs(self.fa, min_len=5000), 60.0)

    def test_length_filter_applied(self):
        """With min_len=20000, only the >=20000bp contigs count."""
        _write_fasta(self.fa, [
            ("NODE_1_length_30000_cov_50.0", "N"),
            ("NODE_2_length_10000_cov_100.0", "N"),   # too short for min_len=20000
            ("NODE_3_length_25000_cov_70.0", "N"),
        ])
        # min_len=20000: median([50, 70]) = 60
        self.assertEqual(estimate_genome_cov_from_contigs(self.fa, min_len=20000), 60.0)
        # min_len=5000: median([50, 100, 70]) = 70
        self.assertEqual(estimate_genome_cov_from_contigs(self.fa, min_len=5000),  70.0)

    def test_no_long_contigs_returns_zero(self):
        """All contigs below cutoff -> 0.0 (pipeline treats this as disabled)."""
        _write_fasta(self.fa, [
            ("NODE_1_length_500_cov_100.0", "N"),
            ("NODE_2_length_100_cov_200.0", "N"),
        ])
        self.assertEqual(estimate_genome_cov_from_contigs(self.fa, min_len=5000), 0.0)

    def test_single_long_contig(self):
        """Single contig at cutoff: median = its cov."""
        _write_fasta(self.fa, [
            ("NODE_1_length_5000_cov_42.123456", "N"),
            ("NODE_2_length_300_cov_999.0", "N"),
        ])
        v = estimate_genome_cov_from_contigs(self.fa, min_len=5000)
        self.assertAlmostEqual(v, 42.123456, places=6)

    def test_very_large_cov_does_not_crash(self):
        """Real data has cov up to ~13e6 on collapsed-repeat 46bp tips. Make sure
        float() handles it cleanly (would crash only if regex truncated mid-number)."""
        _write_fasta(self.fa, [
            ("NODE_1_length_50000_cov_50.0", "N"),
            ("NODE_2_length_40000_cov_60.0", "N"),
            ("NODE_3_length_46_cov_13241919.000000", "N"),   # tip; filtered by length
        ])
        # Length filter drops the absurd-cov tip — median = 55.0
        self.assertEqual(estimate_genome_cov_from_contigs(self.fa, min_len=5000), 55.0)

    def test_very_large_cov_when_above_cutoff(self):
        """If a long contig genuinely has huge cov (rare but allowed), the value
        should propagate through the median without precision loss."""
        _write_fasta(self.fa, [
            ("NODE_1_length_100000_cov_1234567.890123", "N"),
        ])
        v = estimate_genome_cov_from_contigs(self.fa, min_len=5000)
        self.assertAlmostEqual(v, 1234567.890123, places=4)

    def test_small_fractional_cov(self):
        """Very small cov like 0.851064 must NOT round to 0 or 1."""
        _write_fasta(self.fa, [
            ("NODE_1_length_50000_cov_0.851064", "N"),
            ("NODE_2_length_40000_cov_0.851064", "N"),
        ])
        v = estimate_genome_cov_from_contigs(self.fa, min_len=5000)
        self.assertAlmostEqual(v, 0.851064, places=6)
        self.assertNotEqual(v, 0)
        self.assertNotEqual(v, 1)

    def test_unparseable_header_skipped(self):
        """A header not matching the SPAdes pattern is silently skipped."""
        _write_fasta(self.fa, [
            ("NODE_1_length_50000_cov_50.0", "N"),
            ("some_other_assembler_contig_id", "N"),   # no length_/cov_ pattern
            ("NODE_2_length_40000_cov_70.0", "N"),
        ])
        # median([50.0, 70.0]) = 60.0
        self.assertEqual(estimate_genome_cov_from_contigs(self.fa, min_len=5000), 60.0)

    def test_empty_file_returns_zero(self):
        _write_fasta(self.fa, [])
        self.assertEqual(estimate_genome_cov_from_contigs(self.fa, min_len=5000), 0.0)

    def test_default_min_len_5000(self):
        """The function's default cutoff is 5000 bp — caller may omit it."""
        _write_fasta(self.fa, [
            ("NODE_1_length_5000_cov_42.0", "N"),
            ("NODE_2_length_4999_cov_999.0", "N"),    # one byte short of cutoff
        ])
        self.assertEqual(estimate_genome_cov_from_contigs(self.fa), 42.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
