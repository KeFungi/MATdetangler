"""Unit tests for matdetangler.graph_path_search — covers pure-Python pieces:
  - parse_contigs_paths (SPAdes contigs.paths format)
  - undirected_adj_from_links (collapsing L-line directed adjacency)
  - bfs_expand (BFS through L-line adjacency with hop limit)
  - restrict_adj_to_subgraph (restricting directed adjacency to a node subset)
  - is_complete_path (completeness criterion: all variable genes + any flankL + any flankR)
"""
import os, sys, tempfile, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler import graph_path_search as TS
from matdetangler import GFA_search as SA  # for build_directed_adj


class TestParseContigsPaths(unittest.TestCase):
    def test_single_component(self):
        content = ("NODE_1_length_100_cov_50.0_1\n"
                   "42+,17-,9+\n")
        with tempfile.NamedTemporaryFile(mode="w", suffix=".paths", delete=False) as f:
            f.write(content); path = f.name
        try:
            r = TS.parse_contigs_paths(path)
            self.assertIn("NODE_1_length_100_cov_50.0", r)
            self.assertEqual(r["NODE_1_length_100_cov_50.0"],
                             [("42", "+"), ("17", "-"), ("9", "+")])
        finally:
            os.unlink(path)

    def test_multi_component_concat(self):
        content = ("NODE_1_length_100_cov_50.0_1\n"
                   "42+,17-\n"
                   "\n"
                   "NODE_1_length_100_cov_50.0_2\n"
                   "9+\n")
        with tempfile.NamedTemporaryFile(mode="w", suffix=".paths", delete=False) as f:
            f.write(content); path = f.name
        try:
            r = TS.parse_contigs_paths(path)
            # both components concat under the same base name
            self.assertEqual(r["NODE_1_length_100_cov_50.0"],
                             [("42", "+"), ("17", "-"), ("9", "+")])
        finally:
            os.unlink(path)

    def test_empty_file(self):
        self.assertEqual(TS.parse_contigs_paths("/nonexistent/path"), {})


class TestUndirectedAdj(unittest.TestCase):
    def test_from_links(self):
        # L-line tuples: (s1, o1, s2, o2, ov)
        links = [("A", "+", "B", "+", 33),
                 ("B", "+", "C", "-", 33),
                 ("C", "+", "D", "+", 33)]
        adj = TS.undirected_adj_from_links(links)
        self.assertEqual(adj["A"], {"B"})
        self.assertEqual(adj["B"], {"A", "C"})
        self.assertEqual(adj["C"], {"B", "D"})
        self.assertEqual(adj["D"], {"C"})


class TestBFSExpand(unittest.TestCase):
    def setUp(self):
        # linear chain A - B - C - D - E
        self.adj = {
            "A": {"B"}, "B": {"A", "C"}, "C": {"B", "D"},
            "D": {"C", "E"}, "E": {"D"}
        }

    def test_hops_zero(self):
        self.assertEqual(TS.bfs_expand(self.adj, {"C"}, 0), {"C"})

    def test_hops_one(self):
        self.assertEqual(TS.bfs_expand(self.adj, {"C"}, 1), {"B", "C", "D"})

    def test_hops_two(self):
        self.assertEqual(TS.bfs_expand(self.adj, {"C"}, 2), {"A", "B", "C", "D", "E"})

    def test_hops_many_caps_at_reachable(self):
        # cap shouldn't go past the actual reachable set
        self.assertEqual(TS.bfs_expand(self.adj, {"C"}, 99), {"A", "B", "C", "D", "E"})

    def test_multi_seed(self):
        # seeds {A, E}: hop=1 reaches B from A, D from E
        self.assertEqual(TS.bfs_expand(self.adj, {"A", "E"}, 1), {"A", "B", "D", "E"})


class TestRestrictAdj(unittest.TestCase):
    def test_removes_nodes_outside_keep(self):
        # build a small graph with L-lines: A+ -> B+; B+ -> C+; C+ -> D+
        links = [("A", "+", "B", "+", 0),
                 ("B", "+", "C", "+", 0),
                 ("C", "+", "D", "+", 0)]
        adj = SA.build_directed_adj(links)
        sub = TS.restrict_adj_to_subgraph(adj, keep={"A", "B", "C"})
        # (A,+) -> (B,+) survives; (B,+) -> (C,+) survives; (C,+) -> (D,+) is dropped
        self.assertIn(("A", "+"), sub)
        self.assertIn(("B", "+"), sub)
        self.assertEqual(sub.get(("C", "+"), []), [])
        # the (C,+) entry is dropped if it has no surviving outgoing edges
        self.assertNotIn(("C", "+"), sub)


class TestIsCompletePath(unittest.TestCase):
    def test_complete_path(self):
        # path: [flankL] -> [HD1] -> [HD2] -> [flankR]
        path = [("L", "+", 0), ("g1", "+", 0), ("g2", "+", 0), ("R", "+", 0)]
        labels = {"L": "flankL", "g1": "HD1", "g2": "HD2", "R": "flankR"}
        var_per = {"g1": {"HD1"}, "g2": {"HD2"}}
        self.assertTrue(TS.is_complete_path(path, labels, var_per, nvar_total=2))

    def test_missing_variable_gene_incomplete(self):
        # path covers HD1 but not HD2
        path = [("L", "+", 0), ("g1", "+", 0), ("R", "+", 0)]
        labels = {"L": "flankL", "g1": "HD1", "R": "flankR"}
        var_per = {"g1": {"HD1"}}
        self.assertFalse(TS.is_complete_path(path, labels, var_per, nvar_total=2))

    def test_missing_flank_incomplete(self):
        # path covers both variable genes but no flank
        path = [("g1", "+", 0), ("g2", "+", 0)]
        labels = {"g1": "HD1", "g2": "HD2"}
        var_per = {"g1": {"HD1"}, "g2": {"HD2"}}
        self.assertFalse(TS.is_complete_path(path, labels, var_per, nvar_total=2))

    def test_combined_label_segment_counts(self):
        # one segment carries flankL AND HD1 (composite label "HD1+flankL"); another carries HD2+flankR
        path = [("a", "+", 0), ("b", "+", 0)]
        labels = {"a": "HD1+flankL", "b": "HD2+flankR"}
        var_per = {"a": {"HD1"}, "b": {"HD2"}}
        self.assertTrue(TS.is_complete_path(path, labels, var_per, nvar_total=2))

    def test_length_does_not_matter(self):
        """User clarification: completeness is about coverage, NOT length. A 3-segment path
        and a 30-segment path that both cover the same content are equally complete."""
        short_path = [("L", "+", 0), ("g1g2", "+", 0), ("R", "+", 0)]
        long_path  = ([("L", "+", 0)] + [(f"x{i}", "+", 0) for i in range(27)]
                      + [("g1g2", "+", 0), ("R", "+", 0)])
        labels = {"L": "flankL", "g1g2": "HD1+HD2", "R": "flankR"}
        for i in range(27): labels[f"x{i}"] = ""
        var_per = {"g1g2": {"HD1", "HD2"}}
        self.assertTrue(TS.is_complete_path(short_path, labels, var_per, nvar_total=2))
        self.assertTrue(TS.is_complete_path(long_path, labels, var_per, nvar_total=2))


if __name__ == "__main__":
    unittest.main(verbosity=2)
