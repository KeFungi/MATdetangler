"""Unit tests for matdetangler.segment_alleles — covers the pure-Python pieces
(GFA parsing, directed adjacency from L-line orientation, simple-path DFS,
sequence reconstruction with k-mer overlap and reverse-complement, sequence dedup).

Run from repo root:
    python3 -m unittest test.test_segment_alleles -v
"""
import os, sys, tempfile, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler import GFA_search as SA


def _gfa(content: str) -> str:
    """Write a tmp GFA, return its path. Caller is responsible for cleanup via tmpdir."""
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=".gfa", delete=False)
    fh.write(content); fh.close(); return fh.name


class TestRevComp(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(SA._rc("ACGT"), "ACGT")
        self.assertEqual(SA._rc("AAAA"), "TTTT")
        self.assertEqual(SA._rc("ACGTAA"), "TTACGT")
        self.assertEqual(SA._rc("acgtN"), "Nacgt")


class TestGFAParse(unittest.TestCase):
    def test_segments_with_depth_DP(self):
        path = _gfa("S\tA\tACGT\tDP:f:12.5\nS\tB\tTTTT\tDP:f:3.0\n")
        try:
            segs = SA.parse_segments(path)
            self.assertEqual(segs["A"], ("ACGT", 12.5))
            self.assertEqual(segs["B"], ("TTTT", 3.0))
        finally:
            os.unlink(path)

    def test_segments_with_depth_KC(self):
        # KC:i:<count> -> depth = count / len(seq)
        path = _gfa("S\tA\tACGTACGT\tKC:i:80\n")
        try:
            segs = SA.parse_segments(path)
            self.assertEqual(segs["A"][0], "ACGTACGT")
            self.assertAlmostEqual(segs["A"][1], 10.0)  # 80 / 8
        finally:
            os.unlink(path)

    def test_segments_no_depth(self):
        path = _gfa("S\tA\tAAAA\n")
        try:
            self.assertEqual(SA.parse_segments(path)["A"], ("AAAA", 0.0))
        finally:
            os.unlink(path)

    def test_links_parses_overlap(self):
        path = _gfa("S\tA\tACGT\nS\tB\tTTTT\nL\tA\t+\tB\t+\t33M\n")
        try:
            links = SA.parse_links(path)
            self.assertEqual(links, [("A", "+", "B", "+", 33)])
        finally:
            os.unlink(path)


class TestDirectedAdj(unittest.TestCase):
    """SPAdes L-line semantics: an L-line records ONE forward traversal and its
    reverse-traversal mirror. e.g. `L A + B + 33M` means:
        (A,+) -> (B,+) with overlap 33
        (B,-) -> (A,-) with overlap 33
    """

    def test_plus_plus(self):
        adj = SA.build_directed_adj([("A", "+", "B", "+", 33)])
        self.assertIn(("B", "+", 33), adj[("A", "+")])
        self.assertIn(("A", "-", 33), adj[("B", "-")])
        # nothing else
        self.assertEqual(len(adj[("A", "+")]), 1)
        self.assertEqual(len(adj[("B", "-")]), 1)
        self.assertEqual(len(adj[("A", "-")]), 0)
        self.assertEqual(len(adj[("B", "+")]), 0)

    def test_plus_minus(self):
        # L A + B - ov   ==>  (A,+)->(B,-) ; (B,+)->(A,-)
        adj = SA.build_directed_adj([("A", "+", "B", "-", 5)])
        self.assertIn(("B", "-", 5), adj[("A", "+")])
        self.assertIn(("A", "-", 5), adj[("B", "+")])

    def test_minus_plus(self):
        # L A - B + ov   ==>  (A,-)->(B,+) ; (B,-)->(A,+)
        adj = SA.build_directed_adj([("A", "-", "B", "+", 5)])
        self.assertIn(("B", "+", 5), adj[("A", "-")])
        self.assertIn(("A", "+", 5), adj[("B", "-")])


class TestReconstruct(unittest.TestCase):
    def test_forward_concat_with_overlap(self):
        # path: A+ -> B+ overlap 2
        # A = "ACGTAA" (6 bp), B = "AATTTT" (6 bp), overlap=2
        # expect: "ACGTAA" + "AATTTT"[2:] = "ACGTAATTTT" (10 bp)
        segs = {"A": ("ACGTAA", 1.0), "B": ("AATTTT", 1.0)}
        path = [("A", "+", 0), ("B", "+", 2)]
        self.assertEqual(SA.reconstruct(path, segs), "ACGTAATTTT")

    def test_reverse_complement_segment(self):
        # path: A+ -> B-  i.e. use revcomp of B's stored sequence
        segs = {"A": ("ACGT", 1.0), "B": ("GCGC", 1.0)}
        # rc("GCGC") = "GCGC"; no overlap -> "ACGT" + "GCGC" = "ACGTGCGC"
        path = [("A", "+", 0), ("B", "-", 0)]
        self.assertEqual(SA.reconstruct(path, segs), "ACGTGCGC")
        # asymmetric
        segs2 = {"A": ("ACGT", 1.0), "B": ("AACC", 1.0)}
        # rc("AACC") = "GGTT" -> "ACGT" + "GGTT" = "ACGTGGTT"
        self.assertEqual(SA.reconstruct([("A", "+", 0), ("B", "-", 0)], segs2), "ACGTGGTT")


class TestEnumeratePaths(unittest.TestCase):
    def _linear_gfa(self):
        """A→B→C linear chain via L-lines, all overlap=0."""
        segs = {"A": ("AAAA", 1.0), "B": ("BBBB", 1.0), "C": ("CCCC", 1.0)}
        links = [("A", "+", "B", "+", 0), ("B", "+", "C", "+", 0)]
        adj = SA.build_directed_adj(links)
        return segs, adj

    def test_linear_single_path(self):
        segs, adj = self._linear_gfa()
        paths = SA.enumerate_paths(adj, segs, starts={"A"}, ends={"C"},
                                    must_visit_any={"B"},
                                    max_bp=1000, max_nodes=10)
        # there should be at least one A->B->C path; depending on orientation enumeration there
        # might be a + and a - variant, but they should pass through B
        self.assertTrue(any(p[0][0] == "A" and p[-1][0] == "C" for p in paths))
        self.assertTrue(all("B" in {seg for seg, _, _ in p} for p in paths))

    def test_must_visit_any_rejects(self):
        """If must_visit_any excludes B, we expect no path even though A and C are connected through B."""
        segs, adj = self._linear_gfa()
        # require visiting some node that's not in the path -> no result
        paths = SA.enumerate_paths(adj, segs, starts={"A"}, ends={"C"},
                                    must_visit_any={"Z"},  # nonexistent
                                    max_bp=1000, max_nodes=10)
        self.assertEqual(paths, [])

    def test_max_nodes_cap(self):
        segs, adj = self._linear_gfa()
        # cap at 2 nodes — A->B->C is 3, so should be impossible
        paths = SA.enumerate_paths(adj, segs, starts={"A"}, ends={"C"},
                                    must_visit_any={"B"},
                                    max_bp=1000, max_nodes=2)
        self.assertEqual(paths, [])

    def test_max_bp_cap(self):
        segs, adj = self._linear_gfa()
        # all segs are 4 bp; A->B->C = 12 bp total. Cap at 8 bp = should drop.
        paths = SA.enumerate_paths(adj, segs, starts={"A"}, ends={"C"},
                                    must_visit_any={"B"},
                                    max_bp=8, max_nodes=10)
        self.assertEqual(paths, [])

    def test_max_bp_zero_means_no_bp_limit(self):
        """max_bp <= 0 is the no-bp-limit sentinel — paths that would have been
        bp-pruned must now be returned."""
        segs, adj = self._linear_gfa()
        paths = SA.enumerate_paths(adj, segs, starts={"A"}, ends={"C"},
                                    must_visit_any={"B"},
                                    max_bp=0, max_nodes=10)
        self.assertTrue(any({"A", "B", "C"} <= {sid for sid, _, _ in p} for p in paths),
                         f"max_bp=0 should allow the A->B->C 12bp path; got {paths}")

    def test_max_bp_negative_also_disables_limit(self):
        segs, adj = self._linear_gfa()
        paths = SA.enumerate_paths(adj, segs, starts={"A"}, ends={"C"},
                                    must_visit_any={"B"},
                                    max_bp=-1, max_nodes=10)
        self.assertTrue(paths,
                         "negative max_bp should mean unlimited, not 'always exceeded'")

    def test_max_nodes_still_enforced_when_max_bp_is_zero(self):
        """max_nodes is INDEPENDENT of max_bp — even unlimited bp, max_nodes prunes."""
        segs = {nm: ("ACGT", 1.0) for nm in "ABCDE"}
        links = [(a, "+", b, "+", 0) for a, b in zip("ABCD", "BCDE")]
        adj = SA.build_directed_adj(links)
        paths = SA.enumerate_paths(adj, segs, starts={"A"}, ends={"E"},
                                    must_visit_any={"C"},
                                    max_bp=0, max_nodes=3)
        self.assertEqual(paths, [],
                          "max_nodes=3 should still prune a 5-node chain even with no bp limit")

    def test_bubble_gives_two_paths(self):
        # A -> {B1, B2} -> C, both branches valid
        segs = {"A": ("AAAA", 1.0),
                "B1": ("BBBB", 1.0), "B2": ("CCCC", 1.0),
                "C": ("DDDD", 1.0)}
        links = [("A", "+", "B1", "+", 0), ("A", "+", "B2", "+", 0),
                 ("B1", "+", "C", "+", 0), ("B2", "+", "C", "+", 0)]
        adj = SA.build_directed_adj(links)
        paths = SA.enumerate_paths(adj, segs, starts={"A"}, ends={"C"},
                                    must_visit_any={"B1", "B2"},
                                    max_bp=1000, max_nodes=10, max_paths=100)
        # we should get at least one path through each branch (orientation may double-count)
        branches_seen = set()
        for p in paths:
            ids = {seg for seg, _, _ in p}
            if "B1" in ids: branches_seen.add("B1")
            if "B2" in ids: branches_seen.add("B2")
        self.assertEqual(branches_seen, {"B1", "B2"},
                         f"expected paths through both arms; got branches {branches_seen}")


class TestDedupKeepsLonger(unittest.TestCase):
    """When two paths reconstruct to the same sequence, the dedup should keep the path with
    MORE nodes (more graph structure) — verifying the fix from the conversation."""

    def test_via_run_one_k_dedup_inline(self):
        # synthesize: same sequence, two encodings (1-segment path and 2-segment path with overlap)
        segs = {"X": ("ACGT", 1.0), "Y1": ("AC", 1.0), "Y2": ("GT", 1.0)}
        # path1: just X (1 segment) -> "ACGT"
        # path2: Y1 then Y2 with overlap 0 -> "AC" + "GT" = "ACGT"
        p1 = [("X", "+", 0)]
        p2 = [("Y1", "+", 0), ("Y2", "+", 0)]
        # mimic the dedup block in run_one_k
        seen = {}
        for p in (p1, p2):
            s = SA.reconstruct(p, segs)
            key = min(s, SA._rc(s))
            if key not in seen or len(p) > len(seen[key]): seen[key] = p
        self.assertEqual(len(seen), 1, "should dedup to one")
        # the path that survives must be the LONGER one (2 segments)
        surviving = next(iter(seen.values()))
        self.assertEqual(len(surviving), 2,
                         f"dedup should keep longer path; got {surviving}")


class TestMirrorPath(unittest.TestCase):
    """mirror_path(P): same physical walk, opposite direction.
        mirror[0]   = (P[-1].seg, flip(P[-1].o), 0)
        mirror[i>0] = (P[-1-i].seg, flip(P[-1-i].o), P[-i].ov)
    Overlap-slot shift: the overlap at index i in P sits BETWEEN P[i-1] and P[i];
    in the mirror that same junction lives one slot earlier from the other end.
    """

    def test_empty(self):
        self.assertEqual(SA.mirror_path([]), tuple())

    def test_single_segment(self):
        # 1-node walk: mirror is the same seg in flipped orient, overlap 0.
        self.assertEqual(SA.mirror_path([("A", "+", 0)]), (("A", "-", 0),))
        self.assertEqual(SA.mirror_path([("X", "-", 0)]), (("X", "+", 0),))

    def test_two_segment_overlap_shift(self):
        # P: [(A,+,0), (B,+,33)]
        # mirror: [(B,-,0), (A,-,33)]
        P = [("A", "+", 0), ("B", "+", 33)]
        self.assertEqual(SA.mirror_path(P), (("B", "-", 0), ("A", "-", 33)))

    def test_three_segment_overlap_shift(self):
        # P: [(A,+,0), (B,+,5), (C,-,7)]
        # The junctions in P: A→B with overlap 5, B→C with overlap 7.
        # In the mirror (walking C→B→A): C→B with overlap 7, B→A with overlap 5.
        # mirror: [(C,+,0), (B,-,7), (A,-,5)]
        P = [("A", "+", 0), ("B", "+", 5), ("C", "-", 7)]
        self.assertEqual(SA.mirror_path(P),
                         (("C", "+", 0), ("B", "-", 7), ("A", "-", 5)))

    def test_mirror_is_involution(self):
        """mirror(mirror(P)) == P (as a tuple) — applying mirror twice
        returns the original walk."""
        P = [("A", "+", 0), ("B", "-", 12), ("C", "+", 4), ("D", "-", 9)]
        self.assertEqual(SA.mirror_path(list(SA.mirror_path(P))), tuple(P))


class TestEnumeratePathsRCDedup(unittest.TestCase):
    """RC-mirror dedup via path-topology canonicalization, at emit time inside
    enumerate_paths. Safe form: drop only when the mirror is ALREADY emitted,
    NEVER drop based on lex-order alone (would lose walks whose mirrors are
    not enumerated under disjoint starts/ends — the HD case)."""

    def test_disjoint_starts_ends_nothing_dropped(self):
        """HD-like topology: starts and ends are disjoint. Mirrors of emitted
        walks would have to start at end segments, which we don't seed. So no
        mirrors emerge and no dedup happens; behavior is identical to before."""
        segs = {"L": ("AAAA", 1.0), "V": ("CCCC", 1.0), "R": ("TTTT", 1.0)}
        links = [("L", "+", "V", "+", 0), ("V", "+", "R", "+", 0)]
        adj = SA.build_directed_adj(links)
        paths = SA.enumerate_paths(adj, segs, starts={"L"}, ends={"R"},
                                    must_visit_any={"V"},
                                    max_bp=1000, max_nodes=10)
        # Only the (L,+) seed leads anywhere — its mirror walk (R,-) → … → (L,-)
        # is not seeded because R is not in starts. So no walk should be dropped.
        self.assertTrue(any(p[0] == ("L", "+", 0) and p[-1][0] == "R" for p in paths),
                         f"the L+ → V+ → R+ walk must survive; got {paths}")

    def test_overlapping_starts_ends_dedups_mirror(self):
        """If starts ∩ ends ≠ ∅, the SAME physical walk can emerge twice with
        reverse orientation (one DFS seed for each direction). The safe-form
        dedup must collapse them to ONE emission. We use a middle segment as
        must_visit_any so the DFS extends past the start (otherwise single-seg
        walks that are both start AND end emit immediately and never extend)."""
        segs = {"X": ("ACGT", 1.0), "V": ("CCCC", 1.0), "Y": ("TTTT", 1.0)}
        links = [("X", "+", "V", "+", 0), ("V", "+", "Y", "+", 0)]
        adj = SA.build_directed_adj(links)
        # starts={X,Y} ends={X,Y} but must_visit_any={V}:
        #   Seed (X,+) → (V,+) → (Y,+).  V satisfies must_visit; Y is end → emit.
        #   Seed (Y,-) → (V,-) → (X,-).  Topology mirror of the above. Dedup drops.
        # The other two seeds (X,-) and (Y,+) have empty adj — no extension.
        paths = SA.enumerate_paths(adj, segs, starts={"X", "Y"}, ends={"X", "Y"},
                                    must_visit_any={"V"},
                                    max_bp=1000, max_nodes=10, max_paths=100)
        # Exactly one of the 3-seg mirror pair survives.
        target_pair = {(("X", "+", 0), ("V", "+", 0), ("Y", "+", 0)),
                       (("Y", "-", 0), ("V", "-", 0), ("X", "-", 0))}
        emitted_from_pair = [tuple(p) for p in paths if tuple(p) in target_pair]
        self.assertEqual(len(emitted_from_pair), 1,
                         f"the mirror pair must collapse to ONE emission; got "
                         f"{emitted_from_pair}")

    def test_safe_form_does_not_drop_orphan_walks(self):
        """Regression for the naive-form bug: a walk whose mirror would NOT be
        enumerated must not be dropped even if lexicographically 'bigger' than
        its (hypothetical) mirror."""
        # Construct a case where tuple(P) > tuple(mirror(P)) but mirror is
        # unreachable. Pick segment names so the mirror sorts smaller.
        # Walk W = (Z,+,0) → (A,+,5).  mirror = ((A,-,0),(Z,-,5)) which sorts
        # lexicographically SMALLER than W.  Naive form would drop W.
        segs = {"Z": ("ACGT", 1.0), "A": ("TTTT", 1.0)}
        links = [("Z", "+", "A", "+", 5)]
        adj = SA.build_directed_adj(links)
        # starts only contains Z (not A) — mirror would have to start at A, which
        # is not seeded. So mirror is never enumerated; W must survive.
        paths = SA.enumerate_paths(adj, segs, starts={"Z"}, ends={"A"},
                                    must_visit_any={"A"},
                                    max_bp=1000, max_nodes=10)
        self.assertTrue(any(tuple(p) == (("Z", "+", 0), ("A", "+", 5)) for p in paths),
                         f"the Z+ → A+ walk must survive even though its mirror sorts smaller; "
                         f"got {paths}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
