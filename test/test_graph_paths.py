"""Unit tests for matdetangler.graph_paths — focuses on pure-Python pieces:
  - _astr (walk-string formatting, no empty parens for unlabeled nodes)
  - _merge_tokens (LCS alignment of two arms by node-label content)
  - _dedup_loops (collapse paths that revisit a node)
  - _per_allele_gfas_and_paths (parse picks.tsv with segments column)

Run from repo root:
    python3 test/test_graph_paths.py
"""
import os, sys, tempfile, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler import graph_paths as GP


class TestAstr(unittest.TestCase):
    def test_labeled_node(self):
        self.assertEqual(GP._astr(["A"], {"A": "flankL"}), "A (flankL)")

    def test_unlabeled_node_no_empty_parens(self):
        # When the label is "" the node must render as just the id, NOT "id ()"
        self.assertEqual(GP._astr(["A"], {"A": ""}), "A")

    def test_unlabeled_node_missing_from_dict(self):
        self.assertEqual(GP._astr(["A"], {}), "A")

    def test_mixed(self):
        labels = {"A": "flankL", "B": "", "C": "HD1+HD2+flankR"}
        out = GP._astr(["A", "B", "C"], labels)
        self.assertEqual(out, "A (flankL) <-> B <-> C (HD1+HD2+flankR)")

    def test_empty_path(self):
        self.assertEqual(GP._astr([], {}), "")


class TestDedupLoops(unittest.TestCase):
    def test_no_loop_unchanged(self):
        self.assertEqual(GP._dedup_loops(["A", "B", "C"]), ["A", "B", "C"])

    def test_collapses_simple_loop(self):
        # ["flankL", "hub", "X", "hub", "flankR"] => ["flankL", "hub", "flankR"]
        # (hub appears at index 1 and 3; keep prefix to first hub, then suffix after last hub)
        path = ["flankL", "hub", "X", "hub", "flankR"]
        self.assertEqual(GP._dedup_loops(path), ["flankL", "hub", "flankR"])

    def test_multiple_loops(self):
        # First loop: A appears twice. After first collapse: [A, D, A, ...] -> wait, need a fresh case
        path = ["A", "B", "C", "B", "D", "C", "E"]
        # First pass: B at 1 and 3 -> collapse to [A, B, D, C, E]
        # Second pass: no repeats in [A, B, D, C, E] -> stop
        # But our dedup keeps prefix to first match + suffix after last match -> [A, B[0..1+1]] + path[3+1..] = [A, B] + [D, C, E] = [A, B, D, C, E]
        self.assertEqual(GP._dedup_loops(path), ["A", "B", "D", "C", "E"])

    def test_endpoints_preserved(self):
        path = ["flankL", "hub", "main", "hub", "flankR"]
        out = GP._dedup_loops(path)
        self.assertEqual(out[0], "flankL"); self.assertEqual(out[-1], "flankR")


class TestMergeTokens(unittest.TestCase):
    """LCS-by-label alignment. Empty label is treated as wildcard so hubs don't block alignment."""

    def test_identical_arms(self):
        arm1 = ["a1", "a2", "a3"]
        arm2 = ["b1", "b2", "b3"]
        labels = {"a1": "flankL", "a2": "HD", "a3": "flankR",
                  "b1": "flankL", "b2": "HD", "b3": "flankR"}
        tokens = GP._merge_tokens(arm1, arm2, labels)
        # all three positions should pair
        self.assertEqual(len(tokens), 3)
        self.assertEqual(tokens[0], ("a1", "b1", "flankL"))
        self.assertEqual(tokens[1], ("a2", "b2", "HD"))
        self.assertEqual(tokens[2], ("a3", "b3", "flankR"))

    def test_arm2_has_extra_node(self):
        arm1 = ["a1", "a3"]              # flankL -> flankR (no middle)
        arm2 = ["b1", "b2", "b3"]        # flankL -> HD -> flankR
        labels = {"a1": "flankL", "a3": "flankR",
                  "b1": "flankL", "b2": "HD", "b3": "flankR"}
        tokens = GP._merge_tokens(arm1, arm2, labels)
        # expected: (a1,b1,flankL), (None,b2,HD), (a3,b3,flankR)
        self.assertEqual(len(tokens), 3)
        self.assertEqual(tokens[0], ("a1", "b1", "flankL"))
        self.assertEqual(tokens[1], (None, "b2", "HD"))
        self.assertEqual(tokens[2], ("a3", "b3", "flankR"))

    def test_arm2_empty(self):
        arm1 = ["a1", "a2"]
        labels = {"a1": "flankL", "a2": "flankR"}
        tokens = GP._merge_tokens(arm1, [], labels)
        # all positions one-sided arm1
        self.assertEqual(tokens[0], ("a1", None, "flankL"))
        self.assertEqual(tokens[1], ("a2", None, "flankR"))

    def test_empty_label_acts_as_wildcard(self):
        # hub (label "") in both arms should align even though "" "matches" anything
        arm1 = ["a1", "hubA", "a2"]
        arm2 = ["b1", "hubB", "b2"]
        labels = {"a1": "flankL", "hubA": "", "a2": "flankR",
                  "b1": "flankL", "hubB": "", "b2": "flankR"}
        tokens = GP._merge_tokens(arm1, arm2, labels)
        self.assertEqual(len(tokens), 3)
        # middle token must align hubA with hubB; label = "" (empty)
        self.assertEqual(tokens[1][0], "hubA")
        self.assertEqual(tokens[1][1], "hubB")


class TestPerAlleleGfasAndPaths(unittest.TestCase):
    """Parses picks.tsv and resolves per-allele GFA + segment paths."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="_gptest_")
        # synthesize a picks.tsv
        self.picks_tsv = os.path.join(self.tmp, "picks.tsv")
        # header (matches the schema written by pick_alleles.py)
        with open(self.picks_tsv, "w") as f:
            f.write("sample\tallele\torigin\tk\ttype\tlen\tfrom_contig\tsegments\tcov\tn_variable_genes\thas_both_flanks\tis_degHD\n")
            f.write("S\tallele1\tpath\tk33\tcomplete\t100\tFROM\t975425+,975421-,6230157+\t30\t2/2\tTrue\tFalse\n")
            f.write("S\tallele2\tpath\tk33\tcomplete\t100\tFROM\t975425+,975421-,6230158+\t30\t2/2\tTrue\tFalse\n")
        # synthesize primary_alleles.fasta
        self.fa = os.path.join(self.tmp, "primary_alleles.fasta")
        with open(self.fa, "w") as f:
            f.write(">Pcub_S_allele1\nACGT\n>Pcub_S_allele2\nTTTT\n")
        # synthesize a spades_dir with k33/assembly_graph_after_simplification.gfa so the
        # paths resolver finds it
        self.spades_dir = os.path.join(self.tmp, "spades")
        os.makedirs(os.path.join(self.spades_dir, "k33"))
        # GFA needs to exist + have a contigs.fasta or similar; paths.py expects both
        with open(os.path.join(self.spades_dir, "k33", "assembly_graph_after_simplification.gfa"), "w") as f:
            f.write("S\t975425\tACGT\n")
        with open(os.path.join(self.spades_dir, "k33", "contigs.fasta"), "w") as f:
            f.write(">x\nACGT\n")

    def tearDown(self):
        import shutil; shutil.rmtree(self.tmp, ignore_errors=True)

    def test_resolves_segments_per_allele(self):
        gfas, paths = GP._per_allele_gfas_and_paths(self.picks_tsv, self.spades_dir, self.fa)
        # both alleles should resolve to the k33 GFA
        self.assertIn("Pcub_S_allele1", gfas)
        self.assertIn("Pcub_S_allele2", gfas)
        self.assertTrue(gfas["Pcub_S_allele1"].endswith("/k33/assembly_graph_after_simplification.gfa"))
        # paths should be strand-stripped
        self.assertEqual(paths["Pcub_S_allele1"], ["975425", "975421", "6230157"])
        self.assertEqual(paths["Pcub_S_allele2"], ["975425", "975421", "6230158"])

    def test_missing_segments_column_omits_path(self):
        # rewrite picks.tsv without segments column entirely
        with open(self.picks_tsv, "w") as f:
            f.write("sample\tallele\torigin\tk\ttype\tlen\tfrom_contig\tcov\tn_variable_genes\thas_both_flanks\tis_degHD\n")
            f.write("S\tallele1\tpath\tk33\tcomplete\t100\tFROM\t30\t2/2\tTrue\tFalse\n")
        gfas, paths = GP._per_allele_gfas_and_paths(self.picks_tsv, self.spades_dir, self.fa)
        # GFA resolution still works
        self.assertIn("Pcub_S_allele1", gfas)
        # but no recorded path
        self.assertEqual(paths, {})

    def test_dash_segments_omits_path(self):
        # candidate has segments="-" (legacy / unknown)
        with open(self.picks_tsv, "w") as f:
            f.write("sample\tallele\torigin\tk\ttype\tlen\tfrom_contig\tsegments\tcov\tn_variable_genes\thas_both_flanks\tis_degHD\n")
            f.write("S\tallele1\tpath\tk33\tcomplete\t100\tFROM\t-\t30\t2/2\tTrue\tFalse\n")
        gfas, paths = GP._per_allele_gfas_and_paths(self.picks_tsv, self.spades_dir, self.fa)
        self.assertIn("Pcub_S_allele1", gfas)
        self.assertEqual(paths, {})


if __name__ == "__main__":
    unittest.main(verbosity=2)
