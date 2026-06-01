"""Unit tests for the `out_tsv` contract on blast_utils.{blastn,tblastn,tblastx}_hits.

The contract (replacing the old hash-keyed blast_cache.py):
  - If `out_tsv` is None: run blast each time, return the parsed+filtered rows,
    no disk side-effect.
  - If `out_tsv` is given and the file does NOT exist: run blast, write the
    raw rows to disk, return parsed+filtered rows.
  - If `out_tsv` is given and the file DOES exist: skip blast entirely, parse
    the file from disk, apply the filter, return rows.
  - File existence is the ONLY freshness check — there's no hash, no mtime,
    no byte-verify of the query. The wrapper's `--re-blast` flag wipes
    stale files when blast parameters change.
"""
import os, sys, tempfile, unittest
from unittest import mock

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from matdetangler import blast_utils as bu


def _fake_completed_process(stdout: str):
    """Build a CompletedProcess-like object with the given stdout."""
    cp = mock.Mock()
    cp.stdout = stdout
    cp.returncode = 0
    cp.stderr = ""
    return cp


class TestOutTsvContract(unittest.TestCase):

    def setUp(self):
        self.td = tempfile.mkdtemp(prefix="_buotsv_")

    def tearDown(self):
        import shutil; shutil.rmtree(self.td, ignore_errors=True)

    # --- 1. No out_tsv: blast runs every call, no file written ---
    def test_no_out_tsv_runs_blast_each_call_no_disk(self):
        stdout = "seg1\t99.5\t800\nseg2\t98.0\t900\n"
        with mock.patch.object(bu.subprocess, "run",
                                 return_value=_fake_completed_process(stdout)) as m:
            r1 = bu.blastn_hits("q.fa", "db", min_pid=80, min_len=500)
            r2 = bu.blastn_hits("q.fa", "db", min_pid=80, min_len=500)
        self.assertEqual(len(r1), 2)
        self.assertEqual(len(r2), 2)
        # blast was invoked TWICE — no caching when out_tsv is None
        self.assertEqual(m.call_count, 2)
        # no extraneous files in self.td
        self.assertEqual(sorted(os.listdir(self.td)), [])

    # --- 2. out_tsv missing: blast runs, file is written, second call skips ---
    def test_out_tsv_writes_then_reuses(self):
        stdout = "seg1\t99.5\t800\nseg2\t98.0\t900\n"
        out = os.path.join(self.td, "blast_hd_blastn_contigs_k45.tsv")
        with mock.patch.object(bu.subprocess, "run",
                                 return_value=_fake_completed_process(stdout)) as m:
            r1 = bu.blastn_hits("q.fa", "db", min_pid=80, min_len=500, out_tsv=out)
            self.assertEqual(m.call_count, 1)
            self.assertTrue(os.path.exists(out))
            # second call should NOT invoke blast again — file existence is the cache
            r2 = bu.blastn_hits("q.fa", "db", min_pid=80, min_len=500, out_tsv=out)
            self.assertEqual(m.call_count, 1)
        self.assertEqual(r1, r2)
        self.assertEqual(len(r1), 2)

    # --- 3. Re-filter at read-time: file content is preserved verbatim, filter
    #         is re-applied on read so tightening the threshold drops more rows
    #         without re-running blast ---
    def test_tighter_threshold_at_read_filters_more_without_rerun(self):
        # raw blast output: 3 hits with pids 99.5, 85.0, 75.0
        stdout = "seg1\t99.5\t800\nseg2\t85.0\t900\nseg3\t75.0\t1000\n"
        out = os.path.join(self.td, "blast_x.tsv")
        with mock.patch.object(bu.subprocess, "run",
                                 return_value=_fake_completed_process(stdout)) as m:
            # First call writes ALL rows at min_pid=70 (no filter loss)
            r1 = bu.blastn_hits("q.fa", "db", min_pid=70, min_len=500, out_tsv=out)
            self.assertEqual(len(r1), 3)
            self.assertEqual(m.call_count, 1)
            # Second call with TIGHTER threshold — should NOT re-blast, just
            # re-filter the cached rows in memory.
            r2 = bu.blastn_hits("q.fa", "db", min_pid=90, min_len=500, out_tsv=out)
            self.assertEqual(m.call_count, 1)  # still 1 — no re-blast
            self.assertEqual(len(r2), 1)       # only seg1@99.5 survives
            self.assertEqual(r2[0][0], "seg1")

    # --- 4. Empty file is a valid cache entry (blast found nothing) ---
    def test_empty_cached_file_returns_no_hits(self):
        out = os.path.join(self.td, "blast_empty.tsv")
        open(out, "w").close()
        with mock.patch.object(bu.subprocess, "run") as m:
            r = bu.blastn_hits("q.fa", "db", min_pid=80, min_len=500, out_tsv=out)
        self.assertEqual(r, [])
        m.assert_not_called()    # cached empty → no blast run

    # --- 5. The cached file has the exact same name across blast programs ---
    def test_tblastn_and_tblastx_same_out_tsv_contract(self):
        stdout = "seg1\tquery1\t45.0\t100\nseg2\tquery1\t40.0\t60\n"
        out_tn = os.path.join(self.td, "blast_hd_tblastn_segs_k45.tsv")
        out_tx = os.path.join(self.td, "blast_hd_tblastx_segs_k45.tsv")
        with mock.patch.object(bu.subprocess, "run",
                                 return_value=_fake_completed_process(stdout)) as m:
            r_tn = bu.tblastn_hits("p.fa", "db", min_pid=30, min_aa=50, out_tsv=out_tn)
            r_tx = bu.tblastx_hits("n.fa", "db", min_pid=30, min_aa=50, out_tsv=out_tx)
        self.assertEqual(m.call_count, 2)
        self.assertEqual(len(r_tn), 2)
        self.assertEqual(len(r_tx), 2)
        self.assertTrue(os.path.exists(out_tn))
        self.assertTrue(os.path.exists(out_tx))

    # --- 6. The wrapper's --re-blast convention: deleting the file forces
    #         a fresh blast on the next call ---
    def test_deleting_cached_file_forces_rerun(self):
        stdout = "seg1\t99.5\t800\n"
        out = os.path.join(self.td, "blast_x.tsv")
        with mock.patch.object(bu.subprocess, "run",
                                 return_value=_fake_completed_process(stdout)) as m:
            bu.blastn_hits("q.fa", "db", min_pid=80, min_len=500, out_tsv=out)
            self.assertEqual(m.call_count, 1)
            bu.blastn_hits("q.fa", "db", min_pid=80, min_len=500, out_tsv=out)
            self.assertEqual(m.call_count, 1)    # cached
            os.remove(out)
            bu.blastn_hits("q.fa", "db", min_pid=80, min_len=500, out_tsv=out)
            self.assertEqual(m.call_count, 2)    # re-ran after file removal


if __name__ == "__main__":
    unittest.main(verbosity=2)
