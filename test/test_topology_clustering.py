"""Unit tests for two helpers folded into graph_path_search.py:

  1. classify_neighborhood_topology — bubble topology on the BFS neighborhood
     (replaces the old bubble_topo.py output). Pure data in / pure dict out.

  2. cluster_alleles — greedy single-link MAFFT clustering with optional
     core-only divergence scoring. Patches mafft_pair / mafft_pair_core to
     return canned (id, aln_frac) so the logic can be exercised without MAFFT.
"""
import os, sys, unittest
from unittest import mock

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler import graph_path_search as G
from matdetangler.graph_path_search import classify_neighborhood_topology, cluster_alleles


def _build_segs(*specs):
    """Build the `segs` dict expected by classify(): {sid: (seq, depth)}.
    Each spec is (sid, length_bp, depth) — we don't care about actual sequence,
    just length. depth defaults to 1.0."""
    out = {}
    for s in specs:
        sid, L = s[0], s[1]
        depth = s[2] if len(s) >= 3 else 1.0
        out[sid] = ("N" * L, depth)
    return out


def _adj(*edges):
    d: dict[str, set[str]] = {}
    for a, b in edges:
        d.setdefault(a, set()).add(b)
        d.setdefault(b, set()).add(a)
    return d


class TestClassifyNeighborhoodTopology(unittest.TestCase):
    """classify_neighborhood_topology(nhood, adj_und, segs, labels, var_per, min_core_len)"""

    MIN_CORE = 1000

    def test_no_main_when_no_hd_bearing_segment(self):
        """No segment carries a variable gene -> no_main."""
        segs = _build_segs(("seg1", 2000), ("seg2", 2000))
        nhood = {"seg1", "seg2"}
        adj_und = _adj(("seg1", "seg2"))
        labels = {"seg1": "", "seg2": ""}
        var_per = {}  # nothing hits any variable gene
        r = classify_neighborhood_topology(nhood, adj_und, segs, labels, var_per,
                                             min_core_len=self.MIN_CORE)
        self.assertEqual(r["type"], "no_main")
        self.assertEqual(r["n_main_seg"], 0)
        self.assertEqual(r["n_path"], 0)

    def test_no_main_when_hd_segments_too_short(self):
        """HD-bearing but below --min-core-len -> still no_main."""
        segs = _build_segs(("seg1", 500))  # < min_core 1000
        nhood = {"seg1"}
        adj_und = {}
        labels = {"seg1": ""}
        var_per = {"seg1": {"HD1"}}
        r = classify_neighborhood_topology(nhood, adj_und, segs, labels, var_per,
                                             min_core_len=self.MIN_CORE)
        self.assertEqual(r["type"], "no_main")

    def test_single_component(self):
        """One connected blob of HD-bearing segments -> single."""
        segs = _build_segs(("m1", 1500), ("m2", 1200))
        nhood = {"m1", "m2"}
        adj_und = _adj(("m1", "m2"))     # they're connected
        labels = {"m1": "", "m2": ""}
        var_per = {"m1": {"HD1"}, "m2": {"HD2"}}
        r = classify_neighborhood_topology(nhood, adj_und, segs, labels, var_per,
                                             min_core_len=self.MIN_CORE)
        self.assertEqual(r["type"], "single")
        self.assertEqual(r["n_main_seg"], 2)
        self.assertEqual(r["n_path"], 1)

    def test_closed_bubble_with_two_shared_flank_anchors(self):
        """Two HD blobs, both sharing flankL anchor + flankR anchor -> closed_bubble."""
        segs = _build_segs(("a1", 2000), ("a2", 2000), ("fL", 800), ("fR", 800))
        nhood = {"a1", "a2", "fL", "fR"}
        adj_und = _adj(("fL", "a1"), ("fL", "a2"), ("a1", "fR"), ("a2", "fR"))
        labels = {"a1": "", "a2": "", "fL": "flankL", "fR": "flankR"}
        var_per = {"a1": {"HD1", "HD2"}, "a2": {"HD1", "HD2"}}
        r = classify_neighborhood_topology(nhood, adj_und, segs, labels, var_per,
                                             min_core_len=self.MIN_CORE)
        self.assertEqual(r["type"], "closed_bubble")
        self.assertEqual(r["n_path"], 2)
        self.assertGreaterEqual(r["n_shared_flank_anchor"], 2)

    def test_open_bubble_with_one_shared_flank_anchor(self):
        """Two HD blobs, share only flankL (no shared flankR) -> open_bubble."""
        segs = _build_segs(("a1", 2000), ("a2", 2000), ("fL", 800))
        nhood = {"a1", "a2", "fL"}
        adj_und = _adj(("fL", "a1"), ("fL", "a2"))
        labels = {"a1": "", "a2": "", "fL": "flankL"}
        var_per = {"a1": {"HD1"}, "a2": {"HD2"}}
        r = classify_neighborhood_topology(nhood, adj_und, segs, labels, var_per,
                                             min_core_len=self.MIN_CORE)
        self.assertEqual(r["type"], "open_bubble")
        self.assertEqual(r["n_shared_flank_anchor"], 1)

    def test_detached_when_no_shared_flank_anchor(self):
        """Two HD blobs share a non-flank node but no flank anchor -> detached."""
        segs = _build_segs(("a1", 2000), ("a2", 2000), ("x", 600))
        nhood = {"a1", "a2", "x"}
        adj_und = _adj(("x", "a1"), ("x", "a2"))
        labels = {"a1": "", "a2": "", "x": ""}     # x is NOT flank
        var_per = {"a1": {"HD1"}, "a2": {"HD2"}}
        r = classify_neighborhood_topology(nhood, adj_und, segs, labels, var_per,
                                             min_core_len=self.MIN_CORE)
        self.assertEqual(r["type"], "detached")
        self.assertEqual(r["n_shared_anchor"], 1)        # one shared external anchor (x)
        self.assertEqual(r["n_shared_flank_anchor"], 0)  # but it's not a flank

    def test_complexed_more_than_two_components(self):
        segs = _build_segs(("a", 1500), ("b", 1500), ("c", 1500))
        nhood = {"a", "b", "c"}
        adj_und = {}                                  # all isolated -> 3 components
        labels = {"a": "", "b": "", "c": ""}
        var_per = {"a": {"HD1"}, "b": {"HD1"}, "c": {"HD2"}}
        r = classify_neighborhood_topology(nhood, adj_und, segs, labels, var_per,
                                             min_core_len=self.MIN_CORE)
        self.assertEqual(r["type"], "complexed")
        self.assertEqual(r["n_path"], 3)

    def test_segs_outside_nhood_are_ignored(self):
        """Even an HD-bearing seg with the right length must not be counted as
        'main' if it's not in the BFS neighborhood. Regression for the
        fold-bubble-into-step-3 change (bubble_topo used to scan the whole GFA).
        """
        segs = _build_segs(("inside", 1500), ("outside", 1500))
        nhood = {"inside"}                            # 'outside' NOT in nhood
        adj_und = _adj(("inside", "outside"))         # edge exists in the GFA
        labels = {"inside": "", "outside": ""}
        var_per = {"inside": {"HD1"}, "outside": {"HD2"}}
        r = classify_neighborhood_topology(nhood, adj_und, segs, labels, var_per,
                                             min_core_len=self.MIN_CORE)
        # only 'inside' is main_seg; the outside HD blob is invisible
        self.assertEqual(r["type"], "single")
        self.assertEqual(r["n_main_seg"], 1)
        self.assertEqual(set(r["main_ids"]), {"inside"})


class TestClusterAlleles(unittest.TestCase):
    """cluster_alleles greedy single-link clustering."""

    # Helper: build (aligned, hd_cols) for mocking. hd_cols is a SINGLE shared
    # set of MSA column indices (the HD-core region in the REFERENCE's coord
    # system, projected to MSA cols). Default: every column is HD-core (since
    # test strings are intentionally short).
    @staticmethod
    def _mk(aligned: dict[str, str], hd_cols: set[int] | None = None):
        if hd_cols is None:
            length = max(len(s) for s in aligned.values()) if aligned else 0
            hd_cols = set(range(length))
        return aligned, hd_cols

    def test_identical_collapse_to_one_cluster(self):
        """Identical MSA columns (id=1.0, frac=1.0) -> all in one cluster."""
        seqs = [("a", "ACGT"), ("b", "ACGT"), ("c", "ACGT")]
        triple = self._mk({"a": "ACGT", "b": "ACGT", "c": "ACGT"})
        with mock.patch.object(G, "_align_cores", return_value=triple):
            clusters = cluster_alleles(seqs, queries_dir=None,
                                         id_thresh=0.95, frac_thresh=0.80)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(set(clusters[0]), {"a", "b", "c"})

    def test_divergent_split_to_two_clusters(self):
        """No matching MSA columns (id=0.0) -> separate clusters."""
        seqs = [("a", "ACGT"), ("b", "TGCA")]
        triple = self._mk({"a": "ACGT", "b": "TGCA"})
        with mock.patch.object(G, "_align_cores", return_value=triple):
            clusters = cluster_alleles(seqs, queries_dir=None,
                                         id_thresh=0.95, frac_thresh=0.80)
        self.assertEqual(len(clusters), 2)
        self.assertEqual({clusters[0][0], clusters[1][0]}, {"a", "b"})

    def test_chain_a_b_c_with_c_matching_only_b(self):
        """Greedy single-link via REPRESENTATIVES.
              A vs B: 5/6 matches (~0.83)
              A vs C: 0/6 (~0.0)
              B vs C: 1/6 (~0.17)
        Threshold 0.80 → A, B cluster; C separate (compared against A as rep).
        """
        seqs = [("A", "AAAA"), ("B", "AAAB"), ("C", "BBBB")]
        triple = self._mk({"A": "AAAAAA", "B": "AAAAAB", "C": "BBBBBB"})
        with mock.patch.object(G, "_align_cores", return_value=triple):
            clusters = cluster_alleles(seqs, queries_dir=None,
                                         id_thresh=0.80, frac_thresh=0.80)
        self.assertEqual(len(clusters), 2)
        members = sorted([sorted(c) for c in clusters])
        self.assertEqual(members, [["A", "B"], ["C"]])

    def test_threshold_just_below_keeps_separate(self):
        """9/10 column matches -> 0.9 id, below 0.95 threshold -> 2 clusters."""
        seqs = [("a", "ACGTACGTAC"), ("b", "ACGTACGTAA")]
        triple = self._mk({"a": "ACGTACGTAC", "b": "ACGTACGTAA"})
        with mock.patch.object(G, "_align_cores", return_value=triple):
            clusters = cluster_alleles(seqs, queries_dir=None,
                                         id_thresh=0.95, frac_thresh=0.80)
        self.assertEqual(len(clusters), 2)

    def test_align_cores_called_with_protein_fasta_when_queries_dir_supplied(self):
        """When queries_dir is given AND variable_proteins.fasta exists there,
        _align_cores receives that path as the HD-proteins fasta arg."""
        import tempfile
        # Use BYTE-DISTINCT input sequences so the canonical-dedup pre-pass
        # doesn't collapse them before any MAFFT runs; mock _align_cores to
        # return identical aligned strings so the two records cluster together.
        seqs = [("a", "ACGT"), ("b", "ACGA")]
        with tempfile.TemporaryDirectory() as t:
            prot = os.path.join(t, "variable_proteins.fasta")
            with open(prot, "w") as f: f.write(">HD1\nMAAA\n")
            triple = self._mk({"a": "ACGT", "b": "ACGT"})
            with mock.patch.object(G, "_align_cores",
                                     return_value=triple) as ma:
                clusters = cluster_alleles(seqs, queries_dir=t,
                                             id_thresh=0.95, frac_thresh=0.80)
            self.assertTrue(ma.called)
            args, _kw = ma.call_args
            # _align_cores(seqs, locus_ref_fa, hd_proteins_fa) — hd_proteins is arg 3
            self.assertEqual(args[2], prot)
            self.assertEqual(len(clusters), 1)

    def test_identity_scored_only_on_hd_core_columns(self):
        """The whole point of including the reference in the MSA + core-only
        scoring: flanks anchor the alignment but don't contribute to identity.
        Construct an alignment where flank columns disagree entirely but HD-core
        columns agree perfectly — the two candidates should cluster together."""
        seqs = [("a", "FFFFAAAAFFFF"), ("b", "GGGGAAAAGGGG")]
        aligned = {"a": "FFFFAAAAFFFF", "b": "GGGGAAAAGGGG"}
        hd_cols = {4, 5, 6, 7}   # SHARED HD-core columns from the reference
        with mock.patch.object(G, "_align_cores",
                                 return_value=(aligned, hd_cols)):
            clusters = cluster_alleles(seqs, queries_dir=None,
                                         id_thresh=0.95, frac_thresh=0.80)
        # HD-core columns all match -> id=1.0, frac=1.0 -> single cluster
        self.assertEqual(len(clusters), 1)
        self.assertEqual(set(clusters[0]), {"a", "b"})


class TestHdColsFromAlignedRef(unittest.TestCase):
    """Tests for _hd_cols_from_aligned_ref — the column-projection helper used
    inside _align_cores. Covers the ref-flip case from MAFFT --adjustdirection.

    INVARIANT we're guarding: a flip of the reference inside MAFFT must NOT
    change which biological residues of the original ungapped ref are
    selected. The set of *original-ref nucleotides* covered by the picked
    columns has to be identical regardless of strand. The MSA column indices
    will differ, but the underlying nucleotides won't.
    """

    @staticmethod
    def _rc(s: str) -> str:
        comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "-": "-"}
        return "".join(comp[c] for c in reversed(s))

    @staticmethod
    def _nts_at_cols(aligned: str, cols: set[int]) -> list[str]:
        """Pick the (non-gap) chars at the given columns, in column order."""
        return [aligned[i] for i in sorted(cols) if aligned[i] != "-"]

    def test_no_flip_simple_no_gaps(self):
        """Ungapped ref, no flip — span maps to identical indices."""
        ref = "AAAACCCCGTTT"   # len 12; HD-core 5..8 = "CCCC"
        cols = G._hd_cols_from_aligned_ref(ref, 5, 8, flipped=False)
        self.assertEqual(cols, {4, 5, 6, 7})
        self.assertEqual(self._nts_at_cols(ref, cols), list("CCCC"))

    def test_no_flip_with_gaps_in_ref(self):
        """Gaps in aligned ref shift column indices; ungapped position drives."""
        # Original ref = "AAAACCCCGTTT"; insert gaps at positions 0 and 6
        # Aligned ref = "-AAAA-CCCCGTTT" (still ungapped seq AAAACCCCGTTT)
        ref = "-AAAA-CCCCGTTT"
        # span 5..8 in ungapped = the 4 C's, which are at columns 6,7,8,9
        # (col 0 is gap, cols 1-4 are AAAA, col 5 gap, cols 6-9 CCCC, ...)
        cols = G._hd_cols_from_aligned_ref(ref, 5, 8, flipped=False)
        self.assertEqual(cols, {6, 7, 8, 9})
        self.assertEqual(self._nts_at_cols(ref, cols), list("CCCC"))

    def test_flip_simple_no_gaps_picks_same_residues(self):
        """Critical: flipping the ref must select the SAME biological residues.

        original ref:  AAAACCCCGTTT  (HD-core 5..8 = "CCCC")
        RC'd ref:      AAACGGGGTTTT  (the 4 C's are now G's at positions 5..8 of RC)
        """
        orig = "AAAACCCCGTTT"
        rc = self._rc(orig)
        self.assertEqual(rc, "AAACGGGGTTTT")
        cols = G._hd_cols_from_aligned_ref(rc, 5, 8, flipped=True)
        # In RC frame, span 5..8 maps via (12-8+1)..(12-5+1) = 5..8
        # Picked cols 4..7 (0-based) — content "GGGG" (complement of CCCC)
        self.assertEqual(cols, {4, 5, 6, 7})
        self.assertEqual(self._nts_at_cols(rc, cols), list("GGGG"))

    def test_flip_with_gaps_in_rc_ref(self):
        """RC ref has MAFFT gaps; ungapped pos still drives the span."""
        # RC of "AAAACCCCGTTT" = "AAACGGGGTTTT" (ungapped); insert gaps
        # so column indices differ from ungapped positions.
        rc_aln = "AA-AC-GGGGTTT-T"   # ungapped chars: A A A C G G G G T T T T
        # Ungapped span 5..8 of RC corresponds to the 4 G's.
        # Walking left-to-right, non-gap positions 5..8 land at col indices...
        # cols: 0=A,1=A,2=-,3=A,4=C,5=-,6=G,7=G,8=G,9=G,10=T,11=T,12=T,13=-,14=T
        # ungap 1=col0, 2=col1, 3=col3, 4=col4, 5=col6, 6=col7, 7=col8, 8=col9
        cols = G._hd_cols_from_aligned_ref(rc_aln, 5, 8, flipped=True)
        self.assertEqual(cols, {6, 7, 8, 9})
        self.assertEqual(self._nts_at_cols(rc_aln, cols), list("GGGG"))

    def test_flip_asymmetric_span_remaps_correctly(self):
        """Asymmetric span (not centered in ref) — flip is the only way to
        catch a wrong remap; a symmetric test would pass on the wrong formula."""
        # ref = 20 bp; HD-core at 3..7 (an early window)
        ref = "TTAGGGGGAACCAATTAATT"   # span 3..7 (1-based) = positions 3,4,5,6,7 = "AGGGG"
        self.assertEqual(ref[2:7], "AGGGG")
        # No-flip selection:
        cols_fwd = G._hd_cols_from_aligned_ref(ref, 3, 7, flipped=False)
        self.assertEqual(cols_fwd, {2, 3, 4, 5, 6})
        self.assertEqual(self._nts_at_cols(ref, cols_fwd), list("AGGGG"))

        # Flip path: aligned_ref is the RC of the original.
        rc = self._rc(ref)
        # RC of "AGGGG" (HD-core) is "CCCCT" (complement reversed) — in the
        # RC sequence, those 5 bases occupy positions 14..18 (= 20-7+1 .. 20-3+1).
        cols_rev = G._hd_cols_from_aligned_ref(rc, 3, 7, flipped=True)
        self.assertEqual(cols_rev, {13, 14, 15, 16, 17})
        self.assertEqual(self._nts_at_cols(rc, cols_rev), list("CCCCT"))

        # Critical biology check: the residues are RC of each other.
        fwd_nts = "".join(self._nts_at_cols(ref, cols_fwd))
        rev_nts = "".join(self._nts_at_cols(rc, cols_rev))
        self.assertEqual(self._rc(rev_nts), fwd_nts)

    def test_span_covers_whole_ref_flip_invariant(self):
        """If the span covers the entire ref, flipped vs not should both
        select every non-gap column. (Catches off-by-one in the remap.)"""
        ref_aln = "AC-GT-A"     # ungap len 5
        # forward span 1..5
        cols_fwd = G._hd_cols_from_aligned_ref(ref_aln, 1, 5, flipped=False)
        cols_rev = G._hd_cols_from_aligned_ref(ref_aln, 1, 5, flipped=True)
        non_gap = {ci for ci, c in enumerate(ref_aln) if c != "-"}
        self.assertEqual(cols_fwd, non_gap)
        self.assertEqual(cols_rev, non_gap)

    def test_single_position_span_flip(self):
        """1-bp span: flip remap must land on the symmetric base from the
        other end. This catches the classic 'forgot the +1' off-by-one."""
        ref = "AAAAACAAAAA"   # 11 bp; span at position 6 (the lone C)
        cols_fwd = G._hd_cols_from_aligned_ref(ref, 6, 6, flipped=False)
        self.assertEqual(cols_fwd, {5})
        self.assertEqual(self._nts_at_cols(ref, cols_fwd), ["C"])

        # In RC of "AAAAACAAAAA" = "TTTTTGTTTTT", the C became G at position 6.
        # Remap: scan_lo=scan_hi=11-6+1=6 — still position 6 (palindromic-length quirk).
        rc = self._rc(ref)
        cols_rev = G._hd_cols_from_aligned_ref(rc, 6, 6, flipped=True)
        self.assertEqual(cols_rev, {5})
        self.assertEqual(self._nts_at_cols(rc, cols_rev), ["G"])

        # Asymmetric test: span at position 2 of a length-11 ref → RC pos 10.
        cols_fwd2 = G._hd_cols_from_aligned_ref(ref, 2, 2, flipped=False)
        cols_rev2 = G._hd_cols_from_aligned_ref(rc, 2, 2, flipped=True)
        self.assertEqual(cols_fwd2, {1})    # col 1 = original 'A' at pos 2
        self.assertEqual(cols_rev2, {9})    # col 9 = RC pos 10 (11-2+1)
        # Biology check: rc[9] should be complement of ref[1].
        self.assertEqual(rc[9], "T")        # complement of 'A'
        self.assertEqual(ref[1], "A")


if __name__ == "__main__":
    unittest.main(verbosity=2)
