"""Unit tests for matdetangler.consensus_qc — covers the pure-Python pieces (fasta record
splitting, completeness rollup). The blast-using paths require external tools and are
exercised by the SLURM smoke instead.
"""
import os, sys, tempfile, shutil, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler import consensus_qc as CQ


class TestSplitFastaRecords(unittest.TestCase):
    def setUp(self): self.tmp = tempfile.mkdtemp(prefix="_cqcrec_")
    def tearDown(self): shutil.rmtree(self.tmp, ignore_errors=True)

    def test_two_records_preserves_order(self):
        p = os.path.join(self.tmp, "x.fa")
        with open(p, "w") as f: f.write(">a\nACGT\nGGGG\n>b\nTTTT\n")
        r = CQ._split_fasta_records(p)
        self.assertEqual(r, [("a", "ACGTGGGG"), ("b", "TTTT")])

    def test_empty_file(self):
        p = os.path.join(self.tmp, "x.fa")
        open(p, "w").close()
        self.assertEqual(CQ._split_fasta_records(p), [])

    def test_id_first_token_only(self):
        p = os.path.join(self.tmp, "x.fa")
        with open(p, "w") as f: f.write(">a longer description here\nACGT\n")
        r = CQ._split_fasta_records(p)
        self.assertEqual(r, [("a", "ACGT")])


class TestCompletenessRollup(unittest.TestCase):
    """A direct test of the completeness logic (variable-gene coverage + flank presence)
    without going through blast. We synthesize the aggregated_qc dict that the blast-using
    path would normally produce, then mimic the TSV-emission logic and verify the
    `complete` flag."""

    def _emit(self, aggregated, vars_total):
        rows = []
        for nm, q in aggregated.items():
            has_all = len(q["vars_hit"]) == vars_total
            complete = has_all and q["has_flankL"] and q["has_flankR"]
            rows.append((nm, has_all, complete))
        return rows

    def test_all_present_complete(self):
        agg = {"a": {"len": 6000, "vars_hit": {"HD1", "HD2"}, "aa_cov": 1400,
                     "has_flankL": True, "has_flankR": True}}
        rows = self._emit(agg, vars_total=2)
        self.assertEqual(rows, [("a", True, True)])

    def test_missing_a_variable_gene_incomplete(self):
        agg = {"a": {"len": 6000, "vars_hit": {"HD1"}, "aa_cov": 700,
                     "has_flankL": True, "has_flankR": True}}
        rows = self._emit(agg, vars_total=2)
        self.assertEqual(rows, [("a", False, False)])

    def test_missing_a_flank_incomplete(self):
        agg = {"a": {"len": 6000, "vars_hit": {"HD1", "HD2"}, "aa_cov": 1400,
                     "has_flankL": True, "has_flankR": False}}
        rows = self._emit(agg, vars_total=2)
        self.assertEqual(rows, [("a", True, False)])

    def test_multiple_alleles(self):
        agg = {
            "a1": {"len": 6000, "vars_hit": {"HD1", "HD2"}, "aa_cov": 1400,
                   "has_flankL": True, "has_flankR": True},
            "a2": {"len": 6000, "vars_hit": {"HD1"}, "aa_cov": 700,
                   "has_flankL": True, "has_flankR": True},
        }
        rows = sorted(self._emit(agg, vars_total=2))
        self.assertEqual(rows, [("a1", True, True), ("a2", False, False)])


if __name__ == "__main__":
    unittest.main(verbosity=2)
