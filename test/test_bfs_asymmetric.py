"""Unit tests for graph_path_search.bfs_expand — the asymmetric repeat-absorbing BFS.

Pure-data tests against hand-built adjacency dicts. No subprocess, no MAFFT.

Rule under test:
    A non-repeat frontier segment expands normally (all its neighbors are added).
    A repeat frontier segment is absorbed but does NOT expand its neighbors.
    -> Paths can REACH a repeat segment, but the repeat doesn't drag its many
       neighbors into the search neighborhood.
"""
import os, sys, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler.graph_path_search import bfs_expand


def _adj(*edges):
    """Convenience: undirected adjacency dict from (a, b) tuples."""
    d: dict[str, set[str]] = {}
    for a, b in edges:
        d.setdefault(a, set()).add(b)
        d.setdefault(b, set()).add(a)
    return d


class TestSymmetricBaseline(unittest.TestCase):
    """When repeat_segs is None / empty, BFS behaves exactly as before."""

    def test_symmetric_default_one_hop(self):
        adj = _adj(("A", "B"), ("B", "C"))
        visited = bfs_expand(adj, {"A"}, hops=1)
        self.assertEqual(visited, {"A", "B"})

    def test_symmetric_default_two_hops(self):
        adj = _adj(("A", "B"), ("B", "C"), ("C", "D"))
        visited = bfs_expand(adj, {"A"}, hops=2)
        self.assertEqual(visited, {"A", "B", "C"})

    def test_symmetric_default_none_arg_is_empty_set(self):
        adj = _adj(("A", "B"))
        # explicitly None
        visited_none = bfs_expand(adj, {"A"}, hops=1, repeat_segs=None)
        # empty set
        visited_empty = bfs_expand(adj, {"A"}, hops=1, repeat_segs=set())
        self.assertEqual(visited_none, visited_empty)
        self.assertEqual(visited_none, {"A", "B"})


class TestAsymmetricAbsorbButDontExpand(unittest.TestCase):
    """The two halves of the asymmetric rule."""

    def test_absorbs_repeat_neighbor(self):
        """A non-repeat frontier with a repeat neighbor -> the repeat IS added to
        the visited set (rule: paths can enter repeats)."""
        adj = _adj(("A", "R"))
        visited = bfs_expand(adj, {"A"}, hops=1, repeat_segs={"R"})
        self.assertIn("R", visited,
                       "asymmetric BFS should absorb a repeat neighbor of a non-repeat seed")

    def test_does_not_expand_from_repeat(self):
        """A repeat in the frontier does NOT add its neighbors next hop.
        A -> R -> X.  hops=2.
            hop 1: visited={A,R}, frontier={R}
            hop 2: R is a repeat -> skip; X is NOT added.
        """
        adj = _adj(("A", "R"), ("R", "X"))
        visited = bfs_expand(adj, {"A"}, hops=2, repeat_segs={"R"})
        self.assertIn("A", visited)
        self.assertIn("R", visited)
        self.assertNotIn("X", visited,
                          "asymmetric BFS should NOT expand neighbors from a repeat (X is on the far side of R)")

    def test_repeat_as_seed_does_not_expand(self):
        """When a SEED itself is a repeat, its neighbors should not be added."""
        adj = _adj(("R", "X"), ("R", "Y"))
        visited = bfs_expand(adj, {"R"}, hops=3, repeat_segs={"R"})
        self.assertEqual(visited, {"R"},
                          "a repeat seed should not pull in its neighbors under the asymmetric rule")

    def test_chain_blocked_by_middle_repeat(self):
        """A -> R -> B: even with infinite hops, B should not be reached when R is a repeat."""
        adj = _adj(("A", "R"), ("R", "B"))
        visited = bfs_expand(adj, {"A"}, hops=10, repeat_segs={"R"})
        self.assertEqual(visited, {"A", "R"},
                          "the repeat blocks BFS traversal to the far side of the chain")

    def test_alternate_non_repeat_path_still_reaches_far(self):
        """A -> R -> B (blocked) AND A -> N -> B (open): B is reachable via the
        non-repeat side, even though the repeat side is blocked."""
        adj = _adj(("A", "R"), ("R", "B"), ("A", "N"), ("N", "B"))
        visited = bfs_expand(adj, {"A"}, hops=2, repeat_segs={"R"})
        self.assertEqual(visited, {"A", "R", "N", "B"})

    def test_hops_zero_returns_seeds_only(self):
        """hops=0 still returns the seeds verbatim regardless of repeat_segs."""
        adj = _adj(("A", "B"), ("B", "C"))
        visited = bfs_expand(adj, {"A", "B"}, hops=0, repeat_segs={"B"})
        self.assertEqual(visited, {"A", "B"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
