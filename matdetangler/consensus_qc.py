"""Re-check completeness on the read-derived consensus sequences (post-mapping).

The pick-level completeness in picks.tsv reflects the path's tblastn coverage AT PICK TIME
(against the graph-derived candidate sequence). After read mapping + `samtools consensus`,
the consensus can drift from the pick — low-coverage stretches become N's, repeat regions
collapse, etc. This module verifies the read-supported consensus independently.

Per allele consensus FASTA, we run:
    tblastn(variable_proteins -> consensus)   -> which variable genes are covered, total aa
    blastn(flankL  -> consensus)              -> does the consensus reach the 5' flank?
    blastn(flankR  -> consensus)              -> does the consensus reach the 3' flank?

A consensus is `complete` iff:
    (every variable gene has at least one tblastn hit ≥ min_aa) AND
    (any flankL hit exists)                                       AND
    (any flankR hit exists)
— same definition as the path-level is_complete_path test in graph_path_search.

Output: <outdir>/consensus_qc.tsv
    allele  len  vars_hit  vars_total  has_all_vars  aa_cov  has_flankL  has_flankR  complete
"""
from __future__ import annotations
import os, sys, argparse, tempfile, collections
from . import blast_utils as bu


def _split_fasta_records(path: str) -> list[tuple[str, str]]:
    """Return [(name, seq), ...] preserving order."""
    out, name, body = [], None, []
    for ln in open(path):
        if not ln.strip(): continue
        if ln.startswith(">"):
            if name is not None: out.append((name, "".join(body)))
            name = ln[1:].split()[0]; body = []
        elif name is not None:
            body.append(ln.strip())
    if name is not None: out.append((name, "".join(body)))
    return out


def _qc_one_consensus(consensus_fa: str, proteins_fa: str,
                      flankL_fa: str, flankR_fa: str,
                      tblastn_pid: float, tblastn_aa: int,
                      flank_pid: float, flank_minlen: int,
                      threads: int) -> dict:
    """Return per-record QC: {allele_name: {len, vars_hit, aa_cov, has_flankL, has_flankR}}."""
    recs = _split_fasta_records(consensus_fa)
    if not recs:
        return {}
    out: dict[str, dict] = {nm: {"len": len(s), "vars_hit": set(), "aa_cov": 0,
                                  "has_flankL": False, "has_flankR": False}
                             for nm, s in recs}
    with tempfile.TemporaryDirectory() as t:
        # build a single blast DB out of all consensus records
        sf = os.path.join(t, "consensus.fa")
        with open(sf, "w") as o:
            for nm, s in recs: o.write(f">{nm}\n{s}\n")
        db = bu.fasta_to_db(sf, t, name="cons_db")
        # tblastn variable_proteins -> per-record gene set + aa coverage
        if os.path.exists(proteins_fa):
            for r in bu.tblastn_hits(proteins_fa, db, min_pid=tblastn_pid,
                                      min_aa=tblastn_aa, threads=threads):
                sid, gene = r[0], r[1]
                if sid in out:
                    out[sid]["vars_hit"].add(gene)
                    try: out[sid]["aa_cov"] += int(r[3])
                    except (IndexError, ValueError): pass
        # blastn flankL/flankR -> per-record has_flank
        if os.path.exists(flankL_fa):
            for sid in bu.hits_ids(bu.blastn_hits(flankL_fa, db, min_pid=flank_pid,
                                                   min_len=flank_minlen, threads=threads)):
                if sid in out: out[sid]["has_flankL"] = True
        if os.path.exists(flankR_fa):
            for sid in bu.hits_ids(bu.blastn_hits(flankR_fa, db, min_pid=flank_pid,
                                                   min_len=flank_minlen, threads=threads)):
                if sid in out: out[sid]["has_flankR"] = True
    return out


def run(consensus_fastas: list[str], queries_dir: str, out_tsv: str,
        tblastn_pid: float = 30.0, tblastn_aa: int = 50,
        flank_pid: float = 85.0, flank_minlen: int = 100,
        threads: int = 1) -> dict:
    """Verify each consensus FASTA (one file per allele) and write a single combined TSV.
    Returns a dict {allele_name: qc_dict} for downstream summary.
    """
    proteins = os.path.join(queries_dir, "variable_proteins.fasta")
    flankL  = os.path.join(queries_dir, "flankL.fasta")
    flankR  = os.path.join(queries_dir, "flankR.fasta")
    if not os.path.exists(proteins):
        raise SystemExit(f"consensus_qc: queries dir missing variable_proteins.fasta: {queries_dir}")
    # count variable genes total — drives the has_all_vars test
    vars_total = sum(1 for ln in open(proteins) if ln.startswith(">"))
    aggregated: dict[str, dict] = {}
    for fa in consensus_fastas:
        if not fa or not os.path.exists(fa) or os.path.getsize(fa) == 0:
            continue
        qc = _qc_one_consensus(fa, proteins, flankL, flankR,
                                tblastn_pid=tblastn_pid, tblastn_aa=tblastn_aa,
                                flank_pid=flank_pid, flank_minlen=flank_minlen,
                                threads=threads)
        aggregated.update(qc)
    cols = ["allele", "len", "vars_hit", "vars_total", "has_all_vars",
            "aa_cov", "has_flankL", "has_flankR", "complete"]
    with open(out_tsv, "w") as o:
        o.write("\t".join(cols) + "\n")
        for nm, q in aggregated.items():
            has_all = len(q["vars_hit"]) == vars_total
            complete = has_all and q["has_flankL"] and q["has_flankR"]
            o.write("\t".join([nm, str(q["len"]),
                                ",".join(sorted(q["vars_hit"])) or "-",
                                str(vars_total),
                                "TRUE" if has_all else "FALSE",
                                str(q["aa_cov"]),
                                "TRUE" if q["has_flankL"] else "FALSE",
                                "TRUE" if q["has_flankR"] else "FALSE",
                                "TRUE" if complete else "FALSE"]) + "\n")
    print(f"[consensus_qc] wrote {out_tsv}")
    return aggregated


def _cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--consensus", required=True, action="append",
                   help="path to a per-allele consensus.fasta (pass multiple times)")
    p.add_argument("--queries-dir", required=True,
                   help="dir with variable_proteins.fasta, flankL.fasta, flankR.fasta")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--tblastn-pid", type=float, default=30.0)
    p.add_argument("--tblastn-aa", type=int, default=50)
    p.add_argument("--flank-pid", type=float, default=85.0)
    p.add_argument("--flank-minlen", type=int, default=100)
    p.add_argument("--threads", type=int, default=2)
    a = p.parse_args(argv)
    run(a.consensus, a.queries_dir, a.out_tsv,
        a.tblastn_pid, a.tblastn_aa, a.flank_pid, a.flank_minlen, a.threads)


if __name__ == "__main__":
    _cli()
