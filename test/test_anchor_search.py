"""Unit tests for matdetangler.anchor_search.

Covers:
  - --source contigs   uses contigs.fasta as the subject; output name template
                       "<sample>__bubble_<k>_<contig>"
  - --source segments  extracts S-lines from the GFA into a temp fasta; output
                       name template "<sample>__seg_<k>_<seg_id>"
  - 3-tier HD chain    tblastn (always) -> blastn(variable_nt) only when tier 1
                       insufficient -> tblastx(variable_nt) only when tiers 1+2
                       STILL insufficient
  - flank queries      always run; flank-only contigs become anchors
  - min-len filter     anchors shorter than --min-len are dropped from output
  - ann.tsv schema     6 columns: name, len, k, kind, hd_genes, flanks

Mocks the bu.tblastn_hits / blastn_hits / tblastx_hits wrappers + spades_k_paths
so no blast subprocess runs.
"""
import os, sys, tempfile, shutil, unittest
from unittest import mock

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler import anchor_search as A


def _write_fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")


def _setup_queries_dir(tmp):
    """Write a queries-dir with the four fastas anchor_search expects.
    proteins fasta has 2 records -> nvar_total = 2.
    """
    q = os.path.join(tmp, "queries"); os.makedirs(q, exist_ok=True)
    _write_fasta(os.path.join(q, "variable_proteins.fasta"),
                  [("HD1", "MAAAAAAA"), ("HD2", "MCCCCCCC")])
    _write_fasta(os.path.join(q, "variable_nt.fasta"),
                  [("HD1", "ACGT" * 50), ("HD2", "TGCA" * 50)])
    _write_fasta(os.path.join(q, "flankL.fasta"), [("flankL", "A" * 300)])
    _write_fasta(os.path.join(q, "flankR.fasta"), [("flankR", "T" * 300)])
    return q


class _Fake:
    """Bundle of mocked blast wrappers, parameterized per test."""
    def __init__(self,
                 tblastn_rows=None, blastn_nt_rows=None, tblastx_rows=None,
                 flankL_rows=None,  flankR_rows=None):
        # Each row format documented in blast_utils.py.
        # tblastn_rows / tblastx_rows: [sseqid, qseqid, pident, length]
        # blastn_nt_rows (with extra_outfmt="qseqid"): [sseqid, pident, length, qseqid]
        # flankL/flankR_rows (no extra_outfmt):       [sseqid, pident, length]
        self.tblastn_rows = tblastn_rows or []
        self.blastn_nt_rows = blastn_nt_rows or []
        self.tblastx_rows = tblastx_rows or []
        self.flankL_rows  = flankL_rows  or []
        self.flankR_rows  = flankR_rows  or []
        self.tblastn_calls = 0
        self.blastn_calls  = 0
        self.tblastx_calls = 0

    def install(self, queries_dir):
        # Stash for the side_effect callbacks
        self._flankL_path = os.path.join(queries_dir, "flankL.fasta")
        self._flankR_path = os.path.join(queries_dir, "flankR.fasta")
        return [
            mock.patch.object(A.bu, "tblastn_hits", side_effect=self._tblastn),
            mock.patch.object(A.bu, "blastn_hits",  side_effect=self._blastn),
            mock.patch.object(A.bu, "tblastx_hits", side_effect=self._tblastx),
            mock.patch.object(A.bu, "fasta_to_db", side_effect=lambda fa, t, name="db": os.path.join(t, name)),
        ]

    def _tblastn(self, query_fa, db, **kw):
        self.tblastn_calls += 1
        return list(self.tblastn_rows)

    def _blastn(self, query_fa, db, **kw):
        self.blastn_calls += 1
        # Route by which query: flankL fa -> flankL rows; flankR -> flankR rows;
        # otherwise it's the variable_nt tier-2 fallback.
        if query_fa == self._flankL_path: return list(self.flankL_rows)
        if query_fa == self._flankR_path: return list(self.flankR_rows)
        return list(self.blastn_nt_rows)

    def _tblastx(self, query_fa, db, **kw):
        self.tblastx_calls += 1
        return list(self.tblastx_rows)


class TestAnchorSearchContigs(unittest.TestCase):
    """--source contigs path: subject = contigs.fasta."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="_anchor_test_")
        self.qdir = _setup_queries_dir(self.tmp)
        # build a fake contigs.fasta with two records of known length
        self.spades_dir = os.path.join(self.tmp, "spades"); os.makedirs(self.spades_dir)
        kdir = os.path.join(self.spades_dir, "k33"); os.makedirs(kdir)
        self.contigs_fa = os.path.join(kdir, "contigs.fasta")
        _write_fasta(self.contigs_fa,
                      [("CTG_LONG_A", "A" * 3000),     # 3 kb (>= min-len 2000)
                       ("CTG_LONG_B", "C" * 2500),
                       ("CTG_SHORT_C", "G" * 500)])   # 500 bp (< min-len 2000)
        self.outdir = os.path.join(self.tmp, "out"); os.makedirs(self.outdir)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, fake, include_flanks=True):
        """Run anchor_search.run() with spades_k_paths patched to return our fake k dir.

        Default `include_flanks=True` because most tests assert flank-blast behavior.
        Individual tests can pass `include_flanks=False` to assert the no-flank default.
        """
        patches = fake.install(self.qdir) + [
            mock.patch.object(A, "spades_k_paths",
                              return_value=(self.contigs_fa, None)),
        ]
        for p in patches: p.start()
        try:
            return A.run(sample="S", source="contigs",
                          spades_dir=self.spades_dir, queries_dir=self.qdir,
                          ks=["k33"], outdir=self.outdir,
                          min_len=2000, blastn_pid=80.0, blastn_minlen=500,
                          tblastn_pid=30.0, tblastn_aa=50, threads=1,
                          include_flanks=include_flanks)
        finally:
            for p in patches: p.stop()

    def _read_ann(self):
        path = os.path.join(self.outdir, "anchor_contig.ann.tsv")
        return [ln.rstrip("\n").split("\t") for ln in open(path) if ln.strip()]

    def test_tier1_alone_when_complete(self):
        """Two tblastn-complete contigs (each carries HD1 + HD2) -> tier 2 + 3
        do NOT fire; only tblastn was called."""
        fake = _Fake(tblastn_rows=[
            ["CTG_LONG_A", "HD1", "90.0", "200"],
            ["CTG_LONG_A", "HD2", "90.0", "200"],
            ["CTG_LONG_B", "HD1", "90.0", "200"],
            ["CTG_LONG_B", "HD2", "90.0", "200"],
        ])
        self._run(fake)
        self.assertEqual(fake.tblastn_calls, 1)
        # tier 2 blastn fired only for the FLANK queries (not for the nt panel)
        self.assertEqual(fake.blastn_calls, 2,
                          "tier 2 should NOT fire on HD; only the 2 flank blasts ran")
        self.assertEqual(fake.tblastx_calls, 0,
                          "tier 3 must not fire when tier 1 satisfied --min-complete-tblastn")

    def test_tier2_fires_and_satisfies_threshold(self):
        """Tier 1 returns 1 complete contig (< default 2) -> tier 2 fires.
        Tier 2 rescues CTG_LONG_B for BOTH HD1 and HD2 -> CTG_LONG_B becomes
        complete too -> we now have 2 complete contigs, so tier 3 does NOT fire.
        """
        fake = _Fake(
            tblastn_rows=[
                ["CTG_LONG_A", "HD1", "90.0", "200"],
                ["CTG_LONG_A", "HD2", "90.0", "200"],
                # CTG_LONG_B not hit at the protein level
            ],
            # tier 2 (variable_nt) gives CTG_LONG_B BOTH genes so it becomes complete
            blastn_nt_rows=[
                ["CTG_LONG_B", "98.0", "600", "HD1"],
                ["CTG_LONG_B", "98.0", "600", "HD2"],
            ])
        self._run(fake)
        self.assertEqual(fake.tblastn_calls, 1)
        # blastn called for: tier-2 + 2 flank blasts = 3
        self.assertEqual(fake.blastn_calls, 3)
        self.assertEqual(fake.tblastx_calls, 0,
                          "tier 3 must NOT fire when tier 2 brings us to >= min_complete")

    def test_tier3_fires_when_tier1_plus_tier2_insufficient(self):
        """Tiers 1 + 2 leave us with < 2 complete contigs -> tier 3 fires."""
        fake = _Fake(
            tblastn_rows=[],     # tier 1 finds nothing
            blastn_nt_rows=[],   # tier 2 finds nothing
            # tier 3 (tblastx) hits both contigs with both genes
            tblastx_rows=[
                ["CTG_LONG_A", "HD1", "40.0", "60"],
                ["CTG_LONG_A", "HD2", "40.0", "60"],
                ["CTG_LONG_B", "HD1", "40.0", "60"],
                ["CTG_LONG_B", "HD2", "40.0", "60"],
            ])
        self._run(fake)
        self.assertEqual(fake.tblastn_calls, 1)
        self.assertEqual(fake.blastn_calls, 3,
                          "tier 2 + 2 flanks blasts ran")
        self.assertEqual(fake.tblastx_calls, 1,
                          "tier 3 fires when tiers 1+2 didn't reach min_complete")

    def test_flanks_always_run_at_contig_source(self):
        fake = _Fake(tblastn_rows=[])
        self._run(fake)
        self.assertGreaterEqual(fake.blastn_calls, 2,
                                 "flankL + flankR blasts must always run")

    def test_flank_only_contig_becomes_anchor(self):
        """No HD hits at all. flankL hits CTG_LONG_A and flankR hits CTG_LONG_B
        -> both must show up as anchors via the flank path."""
        fake = _Fake(
            tblastn_rows=[],
            flankL_rows=[["CTG_LONG_A", "95.0", "300"]],
            flankR_rows=[["CTG_LONG_B", "95.0", "300"]],
        )
        self._run(fake)
        rows = self._read_ann()
        names = {r[0] for r in rows}
        self.assertIn("S__bubble_k33_CTG_LONG_A", names)
        self.assertIn("S__bubble_k33_CTG_LONG_B", names)
        # And the flanks column should record which flank hit each anchor
        flanks = {r[0]: r[5] for r in rows}
        self.assertIn("flankL", flanks["S__bubble_k33_CTG_LONG_A"])
        self.assertIn("flankR", flanks["S__bubble_k33_CTG_LONG_B"])

    def test_min_len_filter_drops_short_anchor(self):
        """Even if a SHORT contig hits an HD gene, it gets dropped because
        len(c) < --min-len (here 2000)."""
        fake = _Fake(tblastn_rows=[
            ["CTG_SHORT_C", "HD1", "90.0", "100"],   # short contig, hits
        ])
        self._run(fake)
        rows = self._read_ann()
        names = {r[0] for r in rows}
        self.assertNotIn("S__bubble_k33_CTG_SHORT_C", names,
                          "contigs shorter than --min-len must be dropped")

    def test_ann_tsv_has_six_columns(self):
        """Regression: every emitted row has 6 columns
            name, len, k, kind, hd_genes, flanks
        """
        fake = _Fake(
            tblastn_rows=[["CTG_LONG_A", "HD1", "90.0", "200"]],
            flankR_rows=[["CTG_LONG_A", "95.0", "300"]],
        )
        self._run(fake)
        rows = self._read_ann()
        self.assertTrue(rows, "expected at least one anchor row")
        for r in rows:
            self.assertEqual(len(r), 6,
                              f"ann.tsv row must have 6 columns, got {len(r)}: {r}")
        # And the kind cell is "bubble" for --source contigs
        self.assertTrue(all(r[3] == "bubble" for r in rows),
                         "kind column must be 'bubble' for --source contigs")


class TestAnchorSearchSegments(unittest.TestCase):
    """--source segments path: extracts S-lines into a temp fasta."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="_anchor_seg_test_")
        self.qdir = _setup_queries_dir(self.tmp)
        # Build a fake GFA with 3 segments: 2 long + 1 short, 1 has "*" empty seq
        self.spades_dir = os.path.join(self.tmp, "spades")
        kdir = os.path.join(self.spades_dir, "k33"); os.makedirs(kdir)
        self.gfa = os.path.join(kdir, "assembly_graph_after_simplification.gfa")
        with open(self.gfa, "w") as f:
            f.write("S\tSEG_A\t" + "A" * 3000 + "\n")     # 3 kb
            f.write("S\tSEG_B\t" + "C" * 250 + "\n")      # 250 bp
            f.write("S\tSEG_EMPTY\t*\n")                  # sentinel; must be skipped
            f.write("L\tSEG_A\t+\tSEG_B\t+\t0M\n")
        self.outdir = os.path.join(self.tmp, "out"); os.makedirs(self.outdir)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, fake, min_len=200, include_flanks=True):
        patches = fake.install(self.qdir) + [
            mock.patch.object(A, "spades_k_paths",
                              return_value=(None, self.gfa)),
        ]
        for p in patches: p.start()
        try:
            return A.run(sample="S", source="segments",
                          spades_dir=self.spades_dir, queries_dir=self.qdir,
                          ks=["k33"], outdir=self.outdir,
                          min_len=min_len, blastn_pid=80.0, blastn_minlen=500,
                          tblastn_pid=30.0, tblastn_aa=50, threads=1,
                          include_flanks=include_flanks)
        finally:
            for p in patches: p.stop()

    def _read_ann(self):
        path = os.path.join(self.outdir, "anchor_segments.ann.tsv")
        return [ln.rstrip("\n").split("\t") for ln in open(path) if ln.strip()]

    def test_segments_source_extracts_and_skips_empty_seq(self):
        """SEG_EMPTY had '*' as seq — must NOT appear in the candidate set even
        if mocked tblastn pretended it hit something."""
        fake = _Fake(tblastn_rows=[
            ["SEG_A",     "HD1", "90.0", "200"],
            ["SEG_EMPTY", "HD2", "90.0", "200"],   # liar: should be skipped at extract
        ])
        self._run(fake)
        rows = self._read_ann()
        names = {r[0] for r in rows}
        self.assertIn("S__seg_k33_SEG_A", names)
        self.assertNotIn("S__seg_k33_SEG_EMPTY", names,
                          "S-lines with '*' empty seq must be excluded from the fasta")

    def test_segments_source_uses_seg_name_template(self):
        fake = _Fake(tblastn_rows=[
            ["SEG_A", "HD1", "90.0", "200"],
            ["SEG_A", "HD2", "90.0", "200"],
        ])
        self._run(fake)
        rows = self._read_ann()
        self.assertTrue(any(r[0].startswith("S__seg_k33_") for r in rows),
                         "names must use __seg_<k>_<segid> template")
        # 'kind' column reads 'seg'
        self.assertTrue(all(r[3] == "seg" for r in rows))

    def test_short_segment_dropped_at_min_len_200(self):
        """SEG_B is 250 bp. With --min-len 200 it stays; with --min-len 1000 it drops."""
        fake = _Fake(tblastn_rows=[
            ["SEG_B", "HD1", "90.0", "60"],
        ])
        # min_len=200 -> SEG_B (250 bp) kept
        self._run(fake)
        names = {r[0] for r in self._read_ann()}
        self.assertIn("S__seg_k33_SEG_B", names)
        # rerun with min_len=1000 -> SEG_B (250 bp) dropped
        # clean outdir first (recursively — anchor_search now writes a blast_db/ subdir)
        import shutil as _sh
        for f in os.listdir(self.outdir):
            p = os.path.join(self.outdir, f)
            _sh.rmtree(p) if os.path.isdir(p) else os.remove(p)
        self._run(fake, min_len=1000)
        names2 = {r[0] for r in self._read_ann()}
        self.assertNotIn("S__seg_k33_SEG_B", names2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
