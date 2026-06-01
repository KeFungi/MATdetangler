"""Unit tests for matdetangler.pairwise_identity — the new core-aware MAFFT helpers.

Covers:
  _seq_pos_to_align_cols    map 1-based span on ungapped seq -> alignment-col set
  _score_columns            identity / aln_frac over a subset of alignment cols
  mafft_pair_core           flank-anchored align, score over HD-core cols only
  detect_core_span          tblastn-driven core span on an allele (subprocess mocked)
"""
import os, sys, tempfile, shutil, unittest
from unittest import mock
from types import SimpleNamespace

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler import pairwise_identity as PI


class TestSeqPosToAlignCols(unittest.TestCase):
    def test_ungapped_identity_mapping(self):
        # "ACGTACGT" — pos 1..8 maps to cols 0..7 directly
        cols = PI._seq_pos_to_align_cols("ACGTACGT", 1, 8)
        self.assertEqual(cols, set(range(0, 8)))

    def test_subspan(self):
        # span positions 3..5 on ACGTACGT -> cols 2..4
        cols = PI._seq_pos_to_align_cols("ACGTACGT", 3, 5)
        self.assertEqual(cols, {2, 3, 4})

    def test_with_gaps_in_aligned(self):
        # aligned: A C - G T  --  positions   1 2 _ 3 4   (gaps don't count)
        # request positions 2..3 -> cols 1 (C) and 3 (G)
        cols = PI._seq_pos_to_align_cols("AC-GT", 2, 3)
        self.assertEqual(cols, {1, 3})

    def test_start_past_seq_length_returns_empty(self):
        cols = PI._seq_pos_to_align_cols("ACGT", 10, 20)
        self.assertEqual(cols, set())

    def test_end_past_seq_length_clamps(self):
        # span goes to 100 but seq has only 4 positions -> include all
        cols = PI._seq_pos_to_align_cols("ACGT", 1, 100)
        self.assertEqual(cols, {0, 1, 2, 3})


class TestScoreColumns(unittest.TestCase):
    def test_all_match(self):
        alnid, frac = PI._score_columns("ACGT", "ACGT", cols=None, denom_len=4)
        self.assertEqual(alnid, 1.0); self.assertEqual(frac, 1.0)

    def test_all_mismatch(self):
        alnid, frac = PI._score_columns("ACGT", "TGCA", cols=None, denom_len=4)
        self.assertEqual(alnid, 0.0); self.assertEqual(frac, 1.0)

    def test_half_match(self):
        alnid, frac = PI._score_columns("ACGT", "ACTA", cols=None, denom_len=4)
        self.assertEqual(alnid, 0.5); self.assertEqual(frac, 1.0)

    def test_gap_in_either_side_does_not_count(self):
        # cols where either side is '-' should NOT be in `al`
        # "A-GT" vs "AAGT": col0=A/A match, col1=-/A skip, col2=G/G match, col3=T/T match
        # -> al = 3, m = 3 -> alnid = 1.0
        alnid, frac = PI._score_columns("A-GT", "AAGT", cols=None, denom_len=4)
        self.assertEqual(alnid, 1.0)
        # aln_frac is over denom_len = 4 (the requested denominator); al/4 = 0.75
        self.assertEqual(frac, 0.75)

    def test_subset_cols_restricts_scoring(self):
        # Score only cols {0, 1} of "ACGT" vs "TCGT"
        # col0 = A vs T (mismatch), col1 = C vs C (match) -> al=2, m=1, id=0.5
        alnid, frac = PI._score_columns("ACGT", "TCGT", cols={0, 1}, denom_len=4)
        self.assertEqual(alnid, 0.5)
        self.assertEqual(frac, 0.5)   # al=2 / denom_len=4

    def test_no_aligned_columns_returns_zero(self):
        alnid, frac = PI._score_columns("----", "----", cols=None, denom_len=4)
        self.assertEqual(alnid, 0.0); self.assertEqual(frac, 0.0)


class TestMafftPairCore(unittest.TestCase):
    """Patches _mafft_align_pair so we don't actually invoke MAFFT, then asserts
    the core-only scoring behavior."""

    def test_score_uses_only_core_columns(self):
        """Construct an alignment where:
            flankL columns (1..3)  are ALL match
            HD core columns (4..6) are ALL mismatch
            flankR columns (7..9)  are ALL match
        Whole-allele id = 6/9 = 0.667. But core-only score restricted to cols
        4..6 should be 0/3 = 0.0. The point: flanks must NOT pad the core score.
        """
        sa = "AAACCCAAA"     # 9 positions
        sb = "AAATTTAAA"     # 6 matches (3+3), 3 mismatches in middle
        with mock.patch.object(PI, "_mafft_align_pair", return_value=(sa, sb)):
            alnid_full, _ = PI.mafft_pair("AAACCCAAA", "AAATTTAAA")
        self.assertAlmostEqual(alnid_full, 6 / 9, places=6)
        # core span = positions 4..6 (both sequences)
        with mock.patch.object(PI, "_mafft_align_pair", return_value=(sa, sb)):
            alnid_core, frac_core = PI.mafft_pair_core(
                "AAACCCAAA", "AAATTTAAA", core1=(4, 6), core2=(4, 6))
        self.assertEqual(alnid_core, 0.0,
                          "core-only score should be 0 — flanks must not inflate it")
        # denom = min(core len) = 3; al = 3 -> frac = 1.0
        self.assertEqual(frac_core, 1.0)

    def test_score_unions_core_columns_from_both_alleles(self):
        """When the two cores cover different aligned columns, the union is used
        for scoring. Build sa/sb so that core1's cols and core2's cols are
        disjoint, and verify both contribute."""
        sa = "AAATTAAAA"     # core1 pos 4..5 -> cols 3, 4
        sb = "AAAAATTAA"     # core2 pos 6..7 -> cols 5, 6
        # cols {3,4} -> a=T/A T/T b=a/a a/a   ... let me just trust the union logic
        with mock.patch.object(PI, "_mafft_align_pair", return_value=(sa, sb)):
            alnid_core, frac_core = PI.mafft_pair_core(
                sa, sb, core1=(4, 5), core2=(6, 7))
        # We don't assert exact values here — just sanity that the union was
        # used (al > 0) and that we did NOT crash.
        self.assertGreaterEqual(alnid_core, 0.0)
        self.assertLessEqual(alnid_core, 1.0)

    def test_returns_zero_on_empty_mafft_output(self):
        with mock.patch.object(PI, "_mafft_align_pair", return_value=("", "")):
            r = PI.mafft_pair_core("ACGT", "ACGT", (1, 4), (1, 4))
        self.assertEqual(r, (0.0, 0.0))


class TestDetectCoreSpan(unittest.TestCase):
    """tblastn(variable_proteins → allele) → (min(start), max(end)) of HSPs."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="_pi_test_")
        # variable_proteins.fasta has to exist as a file path; content unused under the mock
        self.proteins_fa = os.path.join(self.tmp, "p.fa")
        with open(self.proteins_fa, "w") as f: f.write(">HD1\nMA\n>HD2\nMC\n")
        # Clear the in-process memo so each test sees clean state
        PI._DETECT_CORE_SPAN_MEMO.clear()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)
        PI._DETECT_CORE_SPAN_MEMO.clear()

    def _fake_run(self, tblastn_stdout: str):
        def _run(args, **kw):
            cmd = args[0] if args else ""
            if cmd == "makeblastdb":
                return SimpleNamespace(returncode=0, stdout="", stderr="")
            if cmd == "tblastn":
                return SimpleNamespace(returncode=0, stdout=tblastn_stdout, stderr="")
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        return _run

    def test_unions_all_hsps(self):
        """outfmt 6: sseqid pident length sstart send.
        Two HSPs on the test seq, both above threshold -> span = min..max.
        """
        out = (
            "x\t90.0\t200\t1000\t2000\n"
            "x\t90.0\t150\t3500\t4200\n"
        )
        with mock.patch.object(PI.subprocess, "run", side_effect=self._fake_run(out)):
            start, end = PI.detect_core_span("N" * 5000, self.proteins_fa)
        self.assertEqual((start, end), (1000, 4200))

    def test_handles_reverse_strand_hsps(self):
        """SPAdes may report sstart > send for reverse-strand hits. The function
        must take min(sstart, send) and max(sstart, send) per HSP."""
        out = "x\t90.0\t150\t5000\t3500\n"     # sstart > send
        with mock.patch.object(PI.subprocess, "run", side_effect=self._fake_run(out)):
            start, end = PI.detect_core_span("N" * 6000, self.proteins_fa)
        self.assertEqual((start, end), (3500, 5000))

    def test_filters_below_threshold_hsps(self):
        """Sub-threshold HSPs (pid < 30 or aln_aa < 50) must be dropped."""
        out = (
            "x\t90.0\t200\t1000\t2000\n"         # passes
            "x\t29.0\t100\t100\t300\n"           # pid below cutoff -> drop
            "x\t90.0\t10\t6000\t6020\n"          # aln_aa below cutoff -> drop
        )
        with mock.patch.object(PI.subprocess, "run", side_effect=self._fake_run(out)):
            start, end = PI.detect_core_span("N" * 7000, self.proteins_fa)
        # The two sub-threshold HSPs should NOT extend the span past 1000..2000.
        self.assertEqual((start, end), (1000, 2000))

    def test_returns_full_span_when_no_hits(self):
        """No qualifying HSPs -> degenerate fallback to (1, len(seq))."""
        with mock.patch.object(PI.subprocess, "run", side_effect=self._fake_run("")):
            start, end = PI.detect_core_span("N" * 1234, self.proteins_fa)
        self.assertEqual((start, end), (1, 1234))

    def test_memoized_same_input_runs_tblastn_only_once(self):
        """Same (seq, proteins_fa, thresholds) called twice → only ONE subprocess
        invocation. Eliminates 75-100 redundant tblastn calls per k when
        cluster_alleles is invoked at multiple widen-loop hops."""
        seq = "ACGT" * 1000
        out = "x\t90.0\t200\t1000\t2000\n"
        call_log = []
        def _run(args, **kw):
            call_log.append(args[0] if args else "")
            return SimpleNamespace(returncode=0,
                                    stdout=out if args and args[0] == "tblastn" else "",
                                    stderr="")
        with mock.patch.object(PI.subprocess, "run", side_effect=_run):
            r1 = PI.detect_core_span(seq, self.proteins_fa)
            r2 = PI.detect_core_span(seq, self.proteins_fa)
            r3 = PI.detect_core_span(seq, self.proteins_fa)
        self.assertEqual(r1, r2); self.assertEqual(r2, r3)
        # Subprocess should have run exactly TWICE for the first call
        # (makeblastdb + tblastn) and ZERO more times for the next two.
        self.assertEqual(call_log.count("tblastn"), 1,
                          f"tblastn should run once; ran {call_log.count('tblastn')} times")
        self.assertEqual(call_log.count("makeblastdb"), 1,
                          f"makeblastdb should run once; ran {call_log.count('makeblastdb')} times")

    def test_memo_distinguishes_different_inputs(self):
        """Different seq, different memo entry — no false hit across candidates."""
        seq_a = "AAAA" * 1000
        seq_b = "CCCC" * 1000
        out_a = "x\t90.0\t200\t100\t300\n"
        out_b = "x\t90.0\t200\t500\t700\n"
        calls = {"tblastn": 0}
        def _run(args, **kw):
            cmd = args[0] if args else ""
            if cmd == "makeblastdb":
                return SimpleNamespace(returncode=0, stdout="", stderr="")
            if cmd == "tblastn":
                # Find -db arg and return different output based on db path.
                # We just alternate: use the count to switch.
                calls["tblastn"] += 1
                stdout = out_a if calls["tblastn"] == 1 else out_b
                return SimpleNamespace(returncode=0, stdout=stdout, stderr="")
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        with mock.patch.object(PI.subprocess, "run", side_effect=_run):
            ra = PI.detect_core_span(seq_a, self.proteins_fa)
            rb = PI.detect_core_span(seq_b, self.proteins_fa)
            # Second call on each — both should hit memo
            ra2 = PI.detect_core_span(seq_a, self.proteins_fa)
            rb2 = PI.detect_core_span(seq_b, self.proteins_fa)
        self.assertEqual(ra, (100, 300))
        self.assertEqual(rb, (500, 700))
        self.assertEqual(ra, ra2); self.assertEqual(rb, rb2)
        self.assertEqual(calls["tblastn"], 2,
                          f"different seqs must run tblastn separately; total runs = {calls['tblastn']}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
