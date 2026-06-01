"""Unit tests for matdetangler.paths — the SPAdes layout resolver.

Run from repo root:
    python3 test/test_paths.py
"""
import os, sys, tempfile, shutil, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler.paths import spades_k_paths


class TestSpadesKPaths(unittest.TestCase):
    def setUp(self): self.tmp = tempfile.mkdtemp(prefix="_pathtest_")
    def tearDown(self): shutil.rmtree(self.tmp, ignore_errors=True)

    def _touch(self, *parts):
        path = os.path.join(self.tmp, *parts)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f: f.write("# stub\n")
        return path

    def test_lowercase_k_layout(self):
        """Recommended per-k loop produces spades_dir/k33/{contigs.fasta, assembly_graph_after_simplification.gfa}."""
        self._touch("k33", "contigs.fasta")
        self._touch("k33", "assembly_graph_after_simplification.gfa")
        ctg, gfa = spades_k_paths(self.tmp, "k33")
        self.assertTrue(ctg.endswith("/k33/contigs.fasta"))
        self.assertTrue(gfa.endswith("/k33/assembly_graph_after_simplification.gfa"))

    def test_uppercase_K_layout(self):
        """SPAdes' multi-k internal layout has K33/ (uppercase) with before_rr.fasta or final_contigs.fasta."""
        self._touch("K33", "before_rr.fasta")
        self._touch("K33", "assembly_graph_after_simplification.gfa")
        ctg, gfa = spades_k_paths(self.tmp, "k33")
        self.assertTrue(ctg.endswith("/K33/before_rr.fasta"))
        self.assertTrue(gfa.endswith("/K33/assembly_graph_after_simplification.gfa"))

    def test_uppercase_K_prefers_contigs_over_before_rr(self):
        """When both exist in K{N}/, the resolver should pick contigs.fasta first per CONTIG_NAMES order."""
        # CONTIG_NAMES = ("contigs.fasta", "final_contigs.fasta", "before_rr.fasta")
        self._touch("K55", "contigs.fasta")
        self._touch("K55", "before_rr.fasta")
        self._touch("K55", "assembly_graph_after_simplification.gfa")
        ctg, _ = spades_k_paths(self.tmp, "k55")
        self.assertTrue(ctg.endswith("/K55/contigs.fasta"))

    def test_top_level_fallback(self):
        """When no per-k subdir is present but the top-level has GFA + contigs (the final-k-only
        layout left by a single multi-k SPAdes run), return those."""
        self._touch("assembly_graph_after_simplification.gfa")
        self._touch("contigs.fasta")
        ctg, gfa = spades_k_paths(self.tmp, "k55")
        self.assertTrue(ctg.endswith("/contigs.fasta"))
        # not under a k* subdir
        self.assertNotIn("/k", os.path.relpath(ctg, self.tmp))

    def test_missing_returns_none(self):
        ctg, gfa = spades_k_paths(self.tmp, "k99")
        self.assertIsNone(ctg)
        self.assertIsNone(gfa)

    def test_subdir_without_gfa_skipped(self):
        """A k{N}/ subdir that has contigs.fasta but no GFA should NOT match — caller needs both."""
        self._touch("k33", "contigs.fasta")
        # no assembly_graph file
        ctg, gfa = spades_k_paths(self.tmp, "k33")
        # both should be None (or fall through to top-level which also doesn't exist)
        self.assertIsNone(ctg); self.assertIsNone(gfa)

    def test_k_arg_is_stripped(self):
        """The k argument can be "k33", "K33", or "33" — all should hit the same dir."""
        self._touch("k33", "contigs.fasta")
        self._touch("k33", "assembly_graph_after_simplification.gfa")
        for k in ("k33", "K33", "33"):
            ctg, gfa = spades_k_paths(self.tmp, k)
            self.assertTrue(ctg.endswith("/k33/contigs.fasta"), f"k={k} resolved {ctg}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
