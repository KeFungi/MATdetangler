"""Unit tests for matdetangler.pick_alleles — focuses on the new per-K pair selection
with MAFFT-based divergence. The MAFFT-using `_is_dup` is monkey-patched with a stub
so these tests run with no external dependencies; we test the SELECTION LOGIC, not the
allele-vs-allele identity computation (that's MAFFT, exercised by the integration smoke).

Run from repo root:
    python3 test/test_pick_alleles.py
"""
import os, sys, unittest
from unittest import mock

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler import pick_alleles as PA


class TestKOf(unittest.TestCase):
    def test_path_naming(self):
        nm = "Pcub_X__path_k33_n5_L6000_d40_p3"
        self.assertEqual(PA._kof(nm), "k33")

    def test_legacy_bubble(self):
        nm = "Pcub_X__bubble_k55_NODE_1280_length_6813_cov_27"
        self.assertEqual(PA._kof(nm), "k55")

    def test_unparseable(self):
        self.assertIsNone(PA._kof("random_name_with_no_k_tag"))


class TestCovParsing(unittest.TestCase):
    def test_path_depth_in_name(self):
        # segment_alleles names: ..._d<depth>_p<idx>
        self.assertAlmostEqual(PA._cov("Pcub_X__path_k33_n5_L6000_d42.5_p7"), 42.5)

    def test_legacy_contig_cov(self):
        self.assertAlmostEqual(PA._cov("NODE_1280_length_6813_cov_27.04"), 27.04)


class _Stub:
    """Capture monkey-patched calls during selection tests."""
    def __init__(self): self.dup_calls = []


def _is_dup_stub(known_dups):
    """Returns a stub `_is_dup` that flags a fixed set of (a_seq, b_seq) pairs as dups.
    `known_dups` is a set of frozenset({seq_a, seq_b})."""
    def stub(a, b):
        return frozenset((a, b)) in known_dups
    return stub


class TestPerKPairSelection(unittest.TestCase):
    """Exercises the per-K best-pair logic by calling pick_alleles.run with a synthetic pool.
    We monkey-patch the expensive bits (_is_dup, _vars_in, _flank_hits, _degHD_hits, _cov) so
    the test focuses on selection logic alone."""

    def _setup(self, seqs, vars_per, flankL_set, flankR_set, deg_set,
               dup_sets=None, cov_per=None, nvar_total=2, var_aa_cov=None, flank_cov=None):
        """Patches pick_alleles' module-level helpers and returns the resulting `picks` dict
        after calling `run`. seqs: {name: sequence}; vars_per: {name: set(of variable-gene names)};
        flankL_set / flankR_set / deg_set: sets of names; dup_sets: list of frozenset({a, b})
        pairs to flag as duplicates by sequence; var_aa_cov / flank_cov: {name: int}
        per-candidate tblastn aa coverage and per-candidate flank nucleotide coverage. Defaults
        scale with seq length when not given (so tests that don't care about those metrics
        still work)."""
        dup_sets = dup_sets or []
        known_dups = {frozenset((seqs[a], seqs[b])) for pair in dup_sets for a, b in [tuple(pair)]}
        if var_aa_cov is None:
            var_aa_cov = {c: len(s) // 10 for c, s in seqs.items()}
        if flank_cov is None:
            flank_cov = {c: len(s) // 2 for c, s in seqs.items()}
        with mock.patch.object(PA, "_is_dup", _is_dup_stub(known_dups)), \
             mock.patch.object(PA, "_vars_in", lambda *a, **k: vars_per), \
             mock.patch.object(PA, "_vars_aa_coverage", lambda *a, **k: var_aa_cov), \
             mock.patch.object(PA, "_flank_hits", lambda seqs_d, fa, **k: flankL_set if "flankL" in fa else flankR_set), \
             mock.patch.object(PA, "_flank_coverage_bp", lambda *a, **k: flank_cov), \
             mock.patch.object(PA, "_degHD_hits", lambda *a, **k: deg_set), \
             mock.patch.object(PA, "_cov", (lambda c: cov_per.get(c, 30.0)) if cov_per else (lambda c: 30.0)), \
             mock.patch("builtins.open", create=True, side_effect=_fake_open(seqs, nvar_total)), \
             mock.patch.object(PA, "read_fasta", lambda p: seqs):
            tmpdir = "/tmp/_pa_test"; os.makedirs(tmpdir, exist_ok=True)
            return PA.run(
                sample="S", anchor_contig_fa="cand.fa", queries_dir="q",
                outdir=tmpdir, known_degHD=None,
                expected_count=2, genome_coverage=80.0,
                min_allele_len=10,  # tiny so synthetic 50-bp seqs survive
                max_locus_len=20000,
            )


def _fake_open(seqs, nvar_total):
    """Side-effect for builtins.open so reading `variable_proteins.fasta` returns enough '>' lines
    to make pick_alleles' `nvar_total` count come out right, and writes go to tmp files."""
    import io
    def opener(path, mode="r", *a, **k):
        ms = str(mode)
        if "r" in ms and "variable_proteins.fasta" in str(path):
            # one '>' line per variable gene
            return io.StringIO("".join(f">g{i}\nMAA\n" for i in range(nvar_total)))
        if "w" in ms or "a" in ms:
            return io.StringIO()  # discard writes
        return io.StringIO("")
    return opener


class TestFlankCoverage(TestPerKPairSelection):
    """The CSUF-1723 case: K33 has shorter candidates that are TRUNCATED at the flank ends
    (low flank nucleotide coverage). K55 has slightly longer candidates that reach further
    into the flanks (high flank coverage). The picker MUST prefer the K55 pair because they
    are biologically more complete. aa coverage alone is the same between K33 and K55 (both
    hit HD1+HD2 proteins fully); flank coverage is the discriminator."""

    def test_full_flank_cov_beats_truncated(self):
        seqs = {
            # K33: short, truncated at boundaries (low flank cov)
            "S__path_k33_n1_L5658_d35_p0": "A" * 5658,
            "S__path_k33_n1_L5804_d34_p1": "T" * 5804,
            # K55: slightly longer, full flank cov
            "S__path_k55_n1_L6091_d28_p0": "G" * 6091,
            "S__path_k55_n2_L6243_d27_p1": "C" * 6243,
        }
        vars_per = {c: {"g0", "g1"} for c in seqs}
        flanks_all = set(seqs.keys())
        # SAME aa cov across K33 and K55 (both fully encode both HD proteins)
        var_aa_cov = {c: 1400 for c in seqs}
        # but DIFFERENT flank cov: K33 truncated at ~200 bp into each flank, K55 reaches ~800 bp
        flank_cov = {
            "S__path_k33_n1_L5658_d35_p0": 400,   # 200 bp each side
            "S__path_k33_n1_L5804_d34_p1": 400,
            "S__path_k55_n1_L6091_d28_p0": 1600,  # 800 bp each side
            "S__path_k55_n2_L6243_d27_p1": 1600,
        }
        result = self._setup(seqs, vars_per, flanks_all, flanks_all, deg_set=set(),
                              var_aa_cov=var_aa_cov, flank_cov=flank_cov)
        picks = result["picks_names"]
        for p in picks:
            self.assertEqual(PA._kof(p), "k55",
                              f"pick {p} should be from k55 (full flank coverage); not from k33 (truncated)")


class TestAaCoverageBeforeLength(TestPerKPairSelection):
    """Mirrors the CSUF-1723 case: K33 has SHORTER complete candidates but with smaller
    tblastn aa coverage of the variable proteins (truncated alleles). K55 has slightly
    longer complete candidates with FULL aa coverage (the actual reference alleles). The
    score must prefer K55's higher-coverage pair, not K33's shorter pair.
    """
    def test_full_aa_cov_beats_truncated_shorter(self):
        seqs = {
            # K33: short but truncated (aa coverage of variable proteins low)
            "S__path_k33_n1_L5658_d35_p0": "A" * 5658,
            "S__path_k33_n1_L5804_d34_p1": "T" * 5804,
            # K55: slightly longer but FULL aa coverage of variable proteins
            "S__path_k55_n1_L6091_d28_p0": "G" * 6091,
            "S__path_k55_n2_L6243_d27_p1": "C" * 6243,
        }
        vars_per = {c: {"g0", "g1"} for c in seqs}  # all complete (some hit per gene)
        flanks_all = set(seqs.keys())
        # K33 pair: truncated coverage (HD1 + HD2 ≈ 600 aa total instead of full 1400)
        # K55 pair: ~full coverage (HD1 735 + HD2 673 = 1408 aa)
        var_aa_cov = {
            "S__path_k33_n1_L5658_d35_p0": 600,
            "S__path_k33_n1_L5804_d34_p1": 700,
            "S__path_k55_n1_L6091_d28_p0": 1400,
            "S__path_k55_n2_L6243_d27_p1": 1400,
        }
        result = self._setup(seqs, vars_per, flanks_all, flanks_all, deg_set=set(),
                              var_aa_cov=var_aa_cov)
        picks = result["picks_names"]
        for p in picks:
            self.assertEqual(PA._kof(p), "k55", f"pick {p} should be from k55 (higher aa cov)")
        # both picks should be the K55 6091 and 6243 candidates
        lengths = sorted(int(p.split("_L")[1].split("_")[0]) for p in picks)
        self.assertEqual(lengths, [6091, 6243])


class TestLengthTiebreaker(TestPerKPairSelection):
    """Within the COMPLETE tier, shorter (cleaner walk) should win over longer (decoration variant).
    The previous score tuple rewarded longer paths within max_locus_len, which caused the picker
    to choose a 9 kb decoration variant over the actual ~6 kb biological allele on CSUF-1723."""

    def test_complete_prefers_shorter(self):
        seqs = {
            "S__path_k55_n22_L9089_d30_p1": "G" * 9089,      # decoration variant of allele1
            "S__path_k55_n10_L6243_d30_p2": "A" * 6243,      # the real allele1
            "S__path_k55_n10_L6091_d30_p3": "T" * 6091,      # a divergent allele2
        }
        vars_per = {c: {"g0", "g1"} for c in seqs}
        flanks_all = set(seqs.keys())
        # ALL THREE have identical aa coverage of HD proteins (the decoration is in the
        # middle, not at the flanks). ALSO identical flank coverage (the flanks are at the
        # endpoints; the decoration adds segments BETWEEN HD1 and HD2, not before flankL or
        # after flankR). So the early tiers are all equal — length is the actual tiebreaker.
        var_aa_cov = {c: 1400 for c in seqs}
        flank_cov  = {c: 1600 for c in seqs}
        result = self._setup(seqs, vars_per, flanks_all, flanks_all, deg_set=set(),
                              var_aa_cov=var_aa_cov, flank_cov=flank_cov)
        picks = result["picks_names"]
        # the picker should NOT choose the 9089 bp decoration variant; both picks should be the
        # short, clean walks
        self.assertEqual(len(picks), 2)
        for p in picks:
            L = int(p.split("_L")[1].split("_")[0])
            self.assertLess(L, 7000, f"pick {p} (L={L}) should be a clean walk, not a 9 kb decoration variant")

    def test_partial_pair_includes_both(self):
        # Two partials (only 1 of 2 variable genes hit), both have both flanks.
        # Under the new pair-based picker (2026-05-30), the picker takes the
        # divergent pair if available — it doesn't distinguish "longer wins".
        # Both should appear in picks; per-candidate length ordering is not
        # biologically meaningful (allele1/allele2 labels are interchangeable).
        seqs = {
            "S__path_k55_n10_L4000_d30_p1": "A" * 4000,
            "S__path_k55_n15_L5500_d30_p2": "T" * 5500,
        }
        vars_per = {c: {"g0"} for c in seqs}  # only g0 — not all variable genes -> partial
        flanks_all = set(seqs.keys())
        result = self._setup(seqs, vars_per, flanks_all, flanks_all, deg_set=set(),
                              dup_sets=[])
        picks = set(result["picks_names"])
        self.assertEqual(len(picks), 2)
        self.assertIn("S__path_k55_n10_L4000_d30_p1", picks)
        self.assertIn("S__path_k55_n15_L5500_d30_p2", picks)


class TestPerKSelectionScenarios(TestPerKPairSelection):

    def test_picks_two_divergent_from_one_k(self):
        """K33 has two divergent candidates, both complete. K55 has only one. Expect picks
        from K33: a divergent pair."""
        seqs = {
            "S__path_k33_n5_L100_d30_p1": "A" * 50,        # k33 candidate A (complete, both flanks)
            "S__path_k33_n5_L100_d30_p2": "T" * 50,        # k33 candidate B (complete, both flanks), divergent
            "S__path_k55_n3_L100_d30_p1": "G" * 50,        # k55 lone candidate (complete, both flanks)
        }
        vars_per = {c: {"g0", "g1"} for c in seqs}
        flanks_all = set(seqs.keys())  # every candidate hits both flanks (synthetic full setup)
        result = self._setup(seqs, vars_per, flanks_all, flanks_all, deg_set=set(),
                              dup_sets=[])
        picks = result["picks_names"]
        self.assertEqual(len(picks), 2, f"expected 2 picks, got {picks}")
        # both from k33 (same-K policy)
        for p in picks:
            self.assertEqual(PA._kof(p), "k33", f"pick {p} should be from k33")

    def test_rejects_dup_pair_within_k(self):
        """K33 has two candidates that are flagged as MAFFT-duplicates. K55 has a divergent pair.
        Expect picks from K55."""
        seqs = {
            "S__path_k33_n5_L100_d30_p1": "AAAA",
            "S__path_k33_n5_L100_d30_p2": "AAAB",   # marked dup with the above
            "S__path_k55_n3_L100_d30_p1": "C" * 50,
            "S__path_k55_n3_L100_d30_p2": "T" * 50,
        }
        vars_per = {c: {"g0", "g1"} for c in seqs}
        flanks_all = set(seqs.keys())
        dup_sets = [("S__path_k33_n5_L100_d30_p1", "S__path_k33_n5_L100_d30_p2")]
        result = self._setup(seqs, vars_per, flanks_all, flanks_all, deg_set=set(),
                              dup_sets=dup_sets)
        picks = result["picks_names"]
        # k55 has a divergent pair, k33 doesn't -> picks should be from k55
        self.assertEqual(len(picks), 2)
        for p in picks:
            self.assertEqual(PA._kof(p), "k55",
                              f"all picks should be from k55 (k33 pair was dup); got {p}")

    def test_falls_through_to_single_when_no_pair(self):
        """Only one valid candidate (no possible pair). The picker should still return that one."""
        seqs = {"S__path_k33_n5_L100_d30_p1": "A" * 50}
        vars_per = {"S__path_k33_n5_L100_d30_p1": {"g0", "g1"}}
        flanks_all = set(seqs.keys())
        result = self._setup(seqs, vars_per, flanks_all, flanks_all, deg_set=set(),
                              dup_sets=[])
        picks = result["picks_names"]
        self.assertEqual(picks, ["S__path_k33_n5_L100_d30_p1"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
