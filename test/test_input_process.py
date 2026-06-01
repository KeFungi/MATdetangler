"""Unit tests for matdetangler.input_process — the tblastn-driven HD locator.

These tests mock subprocess.run so they exercise the parsing + confidence-ordered
overlap-resolution logic without needing blast+ binaries on PATH. The integration
smoke (real tblastn against the Pcub locus + proteins) is the existing end-to-end
SLURM run, not a unit test.
"""
import os, sys, tempfile, json, shutil, unittest
from unittest import mock
from types import SimpleNamespace

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler import input_process as ip


def _fake_run_factory(tblastn_out: str):
    """Return a fake subprocess.run that:
       - silently succeeds on makeblastdb
       - returns `tblastn_out` as stdout for tblastn
    """
    def _fake(args, **kw):
        cmd = args[0] if args else ""
        if cmd == "makeblastdb":
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if cmd == "tblastn":
            return SimpleNamespace(returncode=0, stdout=tblastn_out, stderr="")
        return SimpleNamespace(returncode=0, stdout="", stderr="")
    return _fake


class TestTblastnLocate(unittest.TestCase):
    """tblastn outfmt 6 columns are: qseqid sseqid pident length sstart send"""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="_iptest_")
        self.proteins_fa = os.path.join(self.tmp, "p.fa")
        with open(self.proteins_fa, "w") as f:
            f.write(">HD1\nMACGT\n>HD2\nMCGTA\n")
        self.locus_fa = os.path.join(self.tmp, "locus.fa")
        with open(self.locus_fa, "w") as f:
            f.write(">HD_locus\n" + "A" * 12500 + "\n")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _locate(self, tblastn_out: str):
        with mock.patch.object(ip.subprocess, "run", side_effect=_fake_run_factory(tblastn_out)):
            return ip.tblastn_locate(self.proteins_fa, self.locus_fa, threads=1)

    def test_single_hsp_per_protein(self):
        """One clean HSP per protein, no overlap. Span = that HSP; aln_aa = that HSP."""
        out = (
            "HD1\tHD_locus\t90.0\t700\t6320\t8525\n"
            "HD2\tHD_locus\t95.0\t600\t8931\t10948\n"
        )
        hits, intergenic = self._locate(out)
        self.assertEqual(set(hits), {"HD1", "HD2"})
        self.assertEqual((hits["HD1"]["start"], hits["HD1"]["end"]), (6320, 8525))
        self.assertEqual((hits["HD2"]["start"], hits["HD2"]["end"]), (8931, 10948))
        self.assertEqual(hits["HD1"]["aln_aa"], 700)
        self.assertEqual(hits["HD2"]["aln_aa"], 600)
        self.assertEqual(hits["HD1"]["n_hsps_accepted"], 1)
        self.assertEqual(len(intergenic), 1)
        self.assertEqual(intergenic[0]["intergenic_len"], 8931 - 8525 - 1)

    def test_paralog_leakage_rejected(self):
        """HD1's weak HSP that lands in HD2's region must be rejected; HD2's high-conf
        HSP claims that region first (walked in confidence order)."""
        out = (
            "HD1\tHD_locus\t90.0\t700\t6320\t8525\n"   # HD1's real strong hit
            "HD2\tHD_locus\t95.0\t600\t8931\t10948\n"  # HD2's real strong hit
            # leakage: HD1's weak hit in the HD2 region — should be rejected
            "HD1\tHD_locus\t35.0\t150\t9100\t10800\n"
            # leakage: HD2's weak hit in the HD1 region — should be rejected
            "HD2\tHD_locus\t35.0\t150\t6500\t8200\n"
        )
        hits, _ = self._locate(out)
        # HD1's span MUST NOT extend past 8525 (no leakage into HD2's region)
        self.assertEqual(hits["HD1"]["end"], 8525,
                         "HD1's span leaked into HD2's region")
        # HD2's span MUST NOT start before 8931 (no leakage into HD1's region)
        self.assertEqual(hits["HD2"]["start"], 8931,
                         "HD2's span leaked into HD1's region")
        # aln_aa is the BEST HSP only — must equal the strong-hit value, NOT the sum.
        self.assertEqual(hits["HD1"]["aln_aa"], 700,
                         "aln_aa was summed across HSPs; should be best-HSP only")
        self.assertEqual(hits["HD2"]["aln_aa"], 600)

    def test_same_protein_hsps_extend_span(self):
        """Two HSPs from the SAME protein (e.g. divergent N-terminus + conserved core)
        should both be accepted and their union becomes the gene's span. aln_aa stays
        the BEST HSP's value (not summed)."""
        out = (
            "HD1\tHD_locus\t90.0\t700\t6320\t8525\n"   # main HSP, high conf
            "HD1\tHD_locus\t40.0\t60\t5900\t6280\n"    # divergent N-term, weak but same protein
            "HD2\tHD_locus\t95.0\t600\t8931\t10948\n"
        )
        hits, _ = self._locate(out)
        # span extends LEFT to include the divergent N-term hit
        self.assertEqual(hits["HD1"]["start"], 5900,
                         "same-protein N-term HSP was not added to HD1's span")
        self.assertEqual(hits["HD1"]["end"], 8525)
        # aln_aa stays best-HSP only
        self.assertEqual(hits["HD1"]["aln_aa"], 700)
        # both HSPs accepted
        self.assertEqual(hits["HD1"]["n_hsps_accepted"], 2)

    def test_filter_thresholds(self):
        """Below-threshold HSPs (pid < 30 OR aln_aa < 50) are dropped before resolution."""
        out = (
            "HD1\tHD_locus\t90.0\t700\t6320\t8525\n"
            "HD2\tHD_locus\t95.0\t600\t8931\t10948\n"
            "HD1\tHD_locus\t29.0\t100\t5000\t5500\n"   # pid below cutoff
            "HD1\tHD_locus\t90.0\t40\t5800\t5900\n"    # aln_aa below cutoff
        )
        hits, _ = self._locate(out)
        # the two sub-threshold HSPs should NOT extend HD1's span
        self.assertEqual(hits["HD1"]["start"], 6320)
        self.assertEqual(hits["HD1"]["n_hsps_accepted"], 1)

    def test_strand_from_best_hsp(self):
        """Strand is taken from the BEST accepted HSP, not majority-voted."""
        out = (
            # best HSP is on minus strand (sstart > send)
            "HD1\tHD_locus\t90.0\t700\t8525\t6320\n"
            "HD2\tHD_locus\t95.0\t600\t8931\t10948\n"
        )
        hits, _ = self._locate(out)
        self.assertEqual(hits["HD1"]["strand"], "-")
        self.assertEqual(hits["HD2"]["strand"], "+")
        # start/end always sorted ascending regardless of strand
        self.assertLess(hits["HD1"]["start"], hits["HD1"]["end"])

    def test_intergenic_reported_not_split(self):
        """The intergenic gap between adjacent genes is reported AS-IS — gene spans
        end at their actual HSP ends and start at their actual HSP starts; the
        intergenic bp belong to NEITHER gene."""
        out = (
            "HD1\tHD_locus\t90.0\t700\t6320\t8525\n"
            "HD2\tHD_locus\t95.0\t600\t8931\t10948\n"
        )
        hits, intergenic = self._locate(out)
        # gene ends/starts UNTOUCHED by any midpoint or boundary logic
        self.assertEqual(hits["HD1"]["end"], 8525)
        self.assertEqual(hits["HD2"]["start"], 8931)
        # intergenic reported as the actual gap
        self.assertEqual(len(intergenic), 1)
        ig = intergenic[0]
        self.assertEqual(ig["left_gene"], "HD1")
        self.assertEqual(ig["right_gene"], "HD2")
        self.assertEqual(ig["left_end"], 8525)
        self.assertEqual(ig["right_start"], 8931)
        self.assertEqual(ig["intergenic_len"], 405)

    def test_abutting_genes_zero_intergenic(self):
        """Adjacent gene spans that abut (or overlap by 1 bp) report intergenic_len = 0."""
        out = (
            "HD1\tHD_locus\t90.0\t700\t6000\t8000\n"
            "HD2\tHD_locus\t95.0\t600\t8001\t10000\n"   # abuts HD1's end
        )
        hits, intergenic = self._locate(out)
        self.assertEqual(intergenic[0]["intergenic_len"], 0)


class TestWriteQueries(unittest.TestCase):
    """Drive write_queries end-to-end with a mocked tblastn so we can assert the
    manifest schema + flank trimming + envelope padding without needing real blast."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="_iptest_wq_")
        self.proteins_fa = os.path.join(self.tmp, "p.fa")
        with open(self.proteins_fa, "w") as f:
            f.write(">HD1\nMA\n>HD2\nMC\n")
        # 12500 bp locus
        self.locus_fa = os.path.join(self.tmp, "locus.fa")
        with open(self.locus_fa, "w") as f:
            f.write(">HD_locus\n")
            seq = "ACGT" * 3125  # 12500 bp
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")
        self.outdir = os.path.join(self.tmp, "queries")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, tblastn_out: str, **kw):
        with mock.patch.object(ip.subprocess, "run", side_effect=_fake_run_factory(tblastn_out)):
            return ip.write_queries(self.locus_fa, self.proteins_fa, self.outdir,
                                     threads=1, **kw)

    def test_manifest_envelope_padding_and_intergenic(self):
        out = (
            "HD1\tHD_locus\t90.0\t700\t6320\t8525\n"
            "HD2\tHD_locus\t95.0\t600\t8931\t10948\n"
        )
        m = self._run(out, envelope_padding=500, max_flank_len=2000)
        # envelope = (6320 - 500) .. (10948 + 500)
        self.assertEqual(m["envelope_start"], 5820)
        self.assertEqual(m["envelope_end"], 11448)
        self.assertEqual(m["envelope_size"], 11448 - 5820 + 1)
        # flankL = max(0, 5820-1-2000) .. 5820-1  = 3819 .. 5819  -> 2000 bp
        self.assertEqual(m["flankL"]["locus_start"], 3820)
        self.assertEqual(m["flankL"]["locus_end"],   5819)
        self.assertEqual(m["flankL"]["len"], 2000)
        # flankR = 11448 .. min(12500, 11448+2000) = 11448 .. 12500 -> 1052 bp
        self.assertEqual(m["flankR"]["locus_start"], 11449)
        self.assertEqual(m["flankR"]["locus_end"],   12500)
        # intergenic_ranges in the manifest
        self.assertEqual(len(m["intergenic_ranges"]), 1)
        self.assertEqual(m["intergenic_ranges"][0]["intergenic_len"], 405)
        # derived budget = envelope_size + flankL.len + flankR.len
        self.assertEqual(m["derived_max_locus_len"],
                         m["envelope_size"] + m["flankL"]["len"] + m["flankR"]["len"])
        # files written
        for k in ("variable_proteins_fasta", "variable_nt_fasta",
                  "flankL_fasta", "flankR_fasta"):
            self.assertTrue(os.path.exists(m[k]), f"{k} missing")

    def test_padding_clamped_to_locus_bounds(self):
        """If envelope_padding pushes past the locus boundary, it's clamped."""
        out = (
            "HD1\tHD_locus\t90.0\t700\t100\t800\n"        # near left edge
            "HD2\tHD_locus\t95.0\t600\t11700\t12400\n"    # near right edge
        )
        m = self._run(out, envelope_padding=500, max_flank_len=2000)
        self.assertEqual(m["envelope_start"], 1)             # clamped, not -400
        self.assertEqual(m["envelope_end"],   12500)         # clamped

    def test_missing_protein_raises(self):
        """If a protein produces no qualifying HSP, write_queries must fail loudly."""
        out = "HD1\tHD_locus\t90.0\t700\t6320\t8525\n"   # only HD1 hits; HD2 missing
        with self.assertRaises(SystemExit):
            self._run(out)


if __name__ == "__main__":
    unittest.main()


class TestWriteQueriesCache(unittest.TestCase):
    """write_queries(cache_dir=...) — opt-in input/output caching."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="_iptest_cache_")
        self.proteins_fa = os.path.join(self.tmp, "p.fa")
        with open(self.proteins_fa, "w") as f:
            f.write(">HD1\nMA\n>HD2\nMC\n")
        self.locus_fa = os.path.join(self.tmp, "locus.fa")
        with open(self.locus_fa, "w") as f:
            f.write(">HD_locus\n" + ("ACGT" * 3125) + "\n")   # 12500 bp
        self.outdir = os.path.join(self.tmp, "queries")
        self.cache_dir = os.path.join(self.tmp, "_cache")
        # default tblastn rows
        self.tblastn_out = (
            "HD1\tHD_locus\t90.0\t700\t6320\t8525\n"
            "HD2\tHD_locus\t95.0\t600\t8931\t10948\n"
        )

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, **kw):
        with mock.patch.object(ip.subprocess, "run",
                                side_effect=_fake_run_factory(self.tblastn_out)):
            return ip.write_queries(self.locus_fa, self.proteins_fa, self.outdir,
                                     threads=1, **kw)

    def test_no_cache_when_dir_none(self):
        """cache_dir=None → no input_process/ dir created, every call fresh."""
        self._run()
        self.assertFalse(os.path.exists(os.path.join(self.cache_dir, "input_process")))

    def test_cache_write_creates_entry_dir(self):
        """cache_dir=DIR on first call → DIR/input_process/ has locus.fasta +
        proteins.fasta + manifest.json + the 4 derived query fastas."""
        self._run(cache_dir=self.cache_dir)
        entry = os.path.join(self.cache_dir, "input_process")
        self.assertTrue(os.path.isdir(entry))
        for fname in ("locus.fasta", "proteins.fasta", "manifest.json",
                       "variable_proteins.fasta", "variable_nt.fasta",
                       "flankL.fasta", "flankR.fasta"):
            self.assertTrue(os.path.exists(os.path.join(entry, fname)),
                              f"missing {fname} in cache entry")

    def test_cache_hit_skips_subprocess_and_copies_outputs(self):
        """Second call with identical inputs hits cache → subprocess NOT invoked."""
        self._run(cache_dir=self.cache_dir)
        # remove the outdir so we can confirm cache restore writes the files
        shutil.rmtree(self.outdir)
        # Patch subprocess to RAISE if called — proves we didn't re-blast
        with mock.patch.object(ip.subprocess, "run",
                                side_effect=AssertionError("subprocess must not run on cache hit")):
            m = ip.write_queries(self.locus_fa, self.proteins_fa, self.outdir,
                                  threads=1, cache_dir=self.cache_dir)
        # The cached outputs should now be restored to outdir
        for fname in ("manifest.json", "variable_proteins.fasta", "variable_nt.fasta",
                       "flankL.fasta", "flankR.fasta"):
            self.assertTrue(os.path.exists(os.path.join(self.outdir, fname)),
                              f"cache restore did not copy {fname}")
        # Manifest values should match what the cache had
        self.assertEqual(m["chrom"], "HD_locus")

    def test_cache_stale_on_input_byte_change(self):
        """If locus.fasta bytes change, the cache is stale → re-runs."""
        self._run(cache_dir=self.cache_dir)
        # Mutate the input locus
        with open(self.locus_fa, "w") as f:
            f.write(">HD_locus\n" + ("TGCA" * 3125) + "\n")     # different bytes
        ran_blast = {"n": 0}
        def fake_track(args, **kw):
            if args and args[0] == "tblastn": ran_blast["n"] += 1
            return _fake_run_factory(self.tblastn_out)(args, **kw)
        with mock.patch.object(ip.subprocess, "run", side_effect=fake_track):
            ip.write_queries(self.locus_fa, self.proteins_fa, self.outdir,
                              threads=1, cache_dir=self.cache_dir)
        self.assertGreaterEqual(ran_blast["n"], 1,
                                 "byte-changed locus should force tblastn re-run")

    def test_cache_stale_on_param_drift(self):
        """Identical inputs but a different envelope_padding -> cache stale, re-runs."""
        self._run(cache_dir=self.cache_dir, envelope_padding=500)
        ran_blast = {"n": 0}
        def fake_track(args, **kw):
            if args and args[0] == "tblastn": ran_blast["n"] += 1
            return _fake_run_factory(self.tblastn_out)(args, **kw)
        with mock.patch.object(ip.subprocess, "run", side_effect=fake_track):
            ip.write_queries(self.locus_fa, self.proteins_fa, self.outdir,
                              threads=1, cache_dir=self.cache_dir, envelope_padding=750)
        self.assertGreaterEqual(ran_blast["n"], 1,
                                 "envelope_padding change should force re-run")
