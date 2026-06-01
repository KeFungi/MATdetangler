"""Unit tests for graph_path_search.collect_anchors_per_k — the name-template
parser that splits step-2 anchor fasta records by k and by kind.

Templates (set by anchor_search.py):
    contig anchors:    "<sample>__bubble_<k>_<contig_id>"
    segment anchors:   "<sample>__seg_<k>_<segment_id>"

A SPAdes contig id looks like NODE_1515_length_6281_cov_37.78, so the parser
must NOT confuse the underscores in the contig id with the template separators.

These tests exercise the NAME PARSING only. They pass `seeds_from="all"` to
bypass the HD-only filter (which would require ann.tsv companion files that
the name-parsing tests don't bother to create). The default seeds_from="hd"
behavior is covered separately where ann.tsv fixtures exist.
"""
import os, sys, tempfile, shutil, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler.graph_path_search import collect_anchors_per_k


def _write_fasta(path: str, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")


class TestCollectAnchorsPerK(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="_anchor_parse_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_bubble_name_per_k_grouping(self):
        fa = os.path.join(self.tmp, "anchor_contig.fasta")
        _write_fasta(fa, [
            ("Pcub_X__bubble_k33_NODE_111_length_5000_cov_30", "ACGT"),
            ("Pcub_X__bubble_k33_NODE_222_length_6000_cov_35", "ACGT"),
            ("Pcub_X__bubble_k45_NODE_333_length_7000_cov_28", "ACGT"),
        ])
        ctg_per_k, seg_per_k = collect_anchors_per_k(fa, None, seeds_from="all")
        self.assertEqual(seg_per_k, {})
        self.assertEqual(set(ctg_per_k.keys()), {"k33", "k45"})
        self.assertEqual(ctg_per_k["k33"],
                          {"NODE_111_length_5000_cov_30", "NODE_222_length_6000_cov_35"})
        self.assertEqual(ctg_per_k["k45"], {"NODE_333_length_7000_cov_28"})

    def test_seg_name_per_k_grouping(self):
        fa = os.path.join(self.tmp, "anchor_segments.fasta")
        _write_fasta(fa, [
            ("Pcub_X__seg_k33_12345", "ACGT"),
            ("Pcub_X__seg_k45_99999", "ACGT"),
            ("Pcub_X__seg_k53_42",    "ACGT"),
        ])
        ctg_per_k, seg_per_k = collect_anchors_per_k(None, fa, seeds_from="all")
        self.assertEqual(ctg_per_k, {})
        self.assertEqual(seg_per_k["k33"], {"12345"})
        self.assertEqual(seg_per_k["k45"], {"99999"})
        self.assertEqual(seg_per_k["k53"], {"42"})

    def test_mixed_fastas_separate_by_kind(self):
        """Each fasta has its own template — parser must not cross-pollinate."""
        bf = os.path.join(self.tmp, "anchor_contig.fasta")
        sf = os.path.join(self.tmp, "anchor_segments.fasta")
        _write_fasta(bf, [("Pcub_X__bubble_k33_NODE_1_length_5000_cov_30", "ACGT")])
        _write_fasta(sf, [("Pcub_X__seg_k33_7777", "ACGT")])
        ctg_per_k, seg_per_k = collect_anchors_per_k(bf, sf, seeds_from="all")
        self.assertEqual(ctg_per_k["k33"], {"NODE_1_length_5000_cov_30"})
        self.assertEqual(seg_per_k["k33"], {"7777"})

    def test_empty_or_missing_inputs(self):
        empty = os.path.join(self.tmp, "empty.fasta")
        open(empty, "w").close()
        # Both None
        c, s = collect_anchors_per_k(None, None)
        self.assertEqual(c, {}); self.assertEqual(s, {})
        # Missing file path (does not exist)
        c, s = collect_anchors_per_k(os.path.join(self.tmp, "nope.fasta"),
                                       os.path.join(self.tmp, "also_nope.fasta"))
        self.assertEqual(c, {}); self.assertEqual(s, {})
        # Empty file
        c, s = collect_anchors_per_k(empty, empty)
        self.assertEqual(c, {}); self.assertEqual(s, {})

    def test_contig_id_with_many_underscores_preserved(self):
        """Regression for off-by-one: SPAdes NODE_X_length_Y_cov_Z has lots of
        underscores; the parser splits ONLY on the template separators, leaving
        the contig id intact."""
        fa = os.path.join(self.tmp, "f.fasta")
        contig_id = "NODE_1515_length_6281_cov_37.782490"
        _write_fasta(fa, [(f"Pcub_NY-1593933__bubble_k33_{contig_id}", "ACGT")])
        ctg_per_k, _ = collect_anchors_per_k(fa, None, seeds_from="all")
        self.assertEqual(ctg_per_k["k33"], {contig_id},
                          f"contig id should be preserved verbatim, got {ctg_per_k}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
