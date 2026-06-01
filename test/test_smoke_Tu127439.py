"""End-to-end smoke test: run the full MATdetangler pipeline (steps 1-5, reads-free)
on the shipped Tu127439 SPAdes example and assert that the picker recovers a biallelic
pair of complete alleles.

Run via: python -m unittest test.test_smoke_Tu127439 -v

Skips itself if external binaries (mafft / blastn / tblastn / makeblastdb) are
missing — so it never blocks the unit-test suite.
"""
from __future__ import annotations
import os, shutil, subprocess, tempfile, unittest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
WRAPPER   = os.path.join(REPO_ROOT, "MATdetangler")
SPADES    = os.path.join(REPO_ROOT, "examples", "Tu127439_spades")
LOCUS     = os.path.join(REPO_ROOT, "examples", "Suilu_locus", "Suilu4_MATA.fasta")
PROTEINS  = os.path.join(REPO_ROOT, "examples", "Suilu_locus", "Suilu4_HDs.fasta")


def _have(*tools: str) -> bool:
    return all(shutil.which(t) is not None for t in tools)


@unittest.skipUnless(
    _have("mafft", "blastn", "tblastn", "tblastx", "makeblastdb"),
    "external binaries missing — install install/env.yml and activate the conda env first",
)
@unittest.skipUnless(
    os.path.exists(os.path.join(SPADES, "k45", "assembly_graph_after_simplification.gfa")),
    "examples/Tu127439_spades/k45/assembly_graph_after_simplification.gfa missing",
)
class TestSmokeTu127439(unittest.TestCase):
    """Full pipeline run on the Tu127439 example. Steps 1-5 only (reads-free)."""

    def test_pipeline_recovers_biallelic_pair(self) -> None:
        with tempfile.TemporaryDirectory(prefix="matdetangler_smoke_") as outdir:
            cmd = [
                WRAPPER, "run",
                "--sample", "Tu127439",
                "--spades-dir", SPADES,
                "--locus-ref", LOCUS,
                "--proteins",  PROTEINS,
                "--outdir", outdir,
                "--ks", "k33,k45",
                "--threads", "4",
                "--expected-count", "2",
                "--no-skip-pick",
            ]
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
            self.assertEqual(r.returncode, 0,
                             f"MATdetangler exit={r.returncode}\nSTDOUT:\n{r.stdout[-4000:]}\n\nSTDERR:\n{r.stderr[-4000:]}")

            picks_tsv = os.path.join(outdir, "Tu127439", "picks.tsv")
            self.assertTrue(os.path.exists(picks_tsv), f"missing {picks_tsv}")

            rows = [ln.rstrip("\n").split("\t") for ln in open(picks_tsv)]
            header, body = rows[0], rows[1:]
            self.assertEqual(len(body), 2, f"expected 2 picks, got {len(body)} ({body})")

            ix_type    = header.index("type")
            ix_nvars   = header.index("n_variable_genes")
            ix_flanks  = header.index("has_both_flanks")
            for r_ in body:
                self.assertEqual(r_[ix_type],   "complete", f"non-complete pick: {r_}")
                self.assertEqual(r_[ix_nvars],  "2/2",      f"missing variable gene: {r_}")
                self.assertEqual(r_[ix_flanks], "True",     f"missing flank: {r_}")

            id_tsv = os.path.join(outdir, "Tu127439", "identity.tsv")
            self.assertTrue(os.path.exists(id_tsv), f"missing {id_tsv}")
            id_rows = [ln.rstrip("\n").split("\t") for ln in open(id_tsv)]
            id_row = id_rows[1]
            id_pct = float(id_row[id_rows[0].index("id_pct")])
            distinct = id_row[id_rows[0].index("distinct")]
            self.assertEqual(distinct, "True", f"alleles not distinct: id_pct={id_pct}")
            self.assertLess(id_pct, 95.0, f"alleles too similar (id_pct={id_pct}) — picker did not find a divergent pair")


if __name__ == "__main__":
    unittest.main()
