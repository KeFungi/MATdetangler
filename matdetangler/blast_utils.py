"""Thin wrappers around makeblastdb / blastn / tblastn / tblastx used across MATdetangler.

Caching policy (per the user's "named main-results" model):
- Pass `out_tsv=<path>` to skip re-blast when that file already exists. The file
  IS the result, lives next to the sample's other outputs (descriptive name like
  hd_tblastn_contigs.tsv), and is parsed back from disk on subsequent calls.
- File existence is the ONLY check (no hash / no byte-verify / no mtime). Users
  invoke `--no-cached-blast` on the wrapper to clear stale files when they
  change blast parameters.
- Detail blasts (detect_core_span, completeness checks, etc.) just omit
  out_tsv and stay uncached — only the big "search contigs/GFA" pass cares.
"""
from __future__ import annotations
import os, subprocess


def makeblastdb(fasta: str, dbprefix: str | None = None) -> str:
    if dbprefix is None:
        dbprefix = fasta + ".bdb"
    subprocess.run(["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", dbprefix],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return dbprefix


def _read_blast_tsv(path: str) -> list[list[str]]:
    rows: list[list[str]] = []
    with open(path) as f:
        for ln in f:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith("#"): continue
            rows.append(ln.split("\t"))
    return rows


def _write_blast_tsv(path: str, rows: list[list[str]]) -> None:
    d = os.path.dirname(path)
    if d: os.makedirs(d, exist_ok=True)
    with open(path, "w") as f:
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def blastn_hits(query_fa: str, db: str, min_pid: float = 80.0, min_len: int = 500,
                threads: int = 1, extra_outfmt: str = "",
                out_tsv: str | None = None) -> list[list[str]]:
    """Return list of hit rows. Default outfmt: sseqid pident length.

    If `out_tsv` is given and the file exists, parse + filter from disk
    (skip blast). If `out_tsv` is given and missing, run blast and write the
    filtered rows. If `out_tsv` is None, run blast with no disk side-effect.
    """
    if out_tsv and os.path.exists(out_tsv):
        rows = _read_blast_tsv(out_tsv)
    else:
        fmt = f"6 sseqid pident length {extra_outfmt}".strip()
        o = subprocess.run(["blastn", "-query", query_fa, "-db", db,
                            "-dust", "no", "-outfmt", fmt, "-num_threads", str(threads)],
                           capture_output=True, text=True).stdout
        rows = [ln.split("\t") for ln in o.splitlines() if ln]
        if out_tsv: _write_blast_tsv(out_tsv, rows)
    out = []
    for f in rows:
        if len(f) >= 3 and float(f[1]) >= min_pid and float(f[2]) >= min_len:
            out.append(f)
    return out


def tblastn_hits(query_fa: str, db: str, min_pid: float = 30.0, min_aa: int = 50,
                 threads: int = 1, extra_outfmt: str = "",
                 out_tsv: str | None = None) -> list[list[str]]:
    """Return list of hit rows. Default outfmt: sseqid qseqid pident length.

    See `blastn_hits` for the out_tsv caching contract.
    """
    if out_tsv and os.path.exists(out_tsv):
        rows = _read_blast_tsv(out_tsv)
    else:
        fmt = f"6 sseqid qseqid pident length {extra_outfmt}".strip()
        o = subprocess.run(["tblastn", "-query", query_fa, "-db", db,
                            "-evalue", "1e-5", "-outfmt", fmt, "-num_threads", str(threads)],
                           capture_output=True, text=True).stdout
        rows = [ln.split("\t") for ln in o.splitlines() if ln]
        if out_tsv: _write_blast_tsv(out_tsv, rows)
    out = []
    for f in rows:
        if len(f) >= 4 and float(f[2]) >= min_pid and float(f[3]) >= min_aa:
            out.append(f)
    return out


def tblastx_hits(query_fa: str, db: str, min_pid: float = 30.0, min_aa: int = 50,
                 threads: int = 1, extra_outfmt: str = "",
                 out_tsv: str | None = None) -> list[list[str]]:
    """tblastx — both query AND subject translated in 6 frames before alignment.
    Default outfmt: sseqid qseqid pident length (same columns as tblastn).

    Used as the DEEPEST contig-level fallback: ~36× slower than tblastn but catches
    cases where both the protein query (tblastn) and the NT panel (blastn) miss the
    target — e.g. a variable gene whose CDS has diverged so much at the NT level that
    blastn's 80 % cutoff filters it out, while the protein query happens to also have
    diverged enough at the aa level to miss the 30 % tblastn cutoff in a way that any
    6-frame translation of the NT span would still recognize.

    See `blastn_hits` for the out_tsv caching contract.
    """
    if out_tsv and os.path.exists(out_tsv):
        rows = _read_blast_tsv(out_tsv)
    else:
        fmt = f"6 sseqid qseqid pident length {extra_outfmt}".strip()
        o = subprocess.run(["tblastx", "-query", query_fa, "-db", db,
                            "-evalue", "1e-5", "-outfmt", fmt, "-num_threads", str(threads)],
                           capture_output=True, text=True).stdout
        rows = [ln.split("\t") for ln in o.splitlines() if ln]
        if out_tsv: _write_blast_tsv(out_tsv, rows)
    out = []
    for f in rows:
        if len(f) >= 4 and float(f[2]) >= min_pid and float(f[3]) >= min_aa:
            out.append(f)
    return out


def hits_ids(rows: list[list[str]], col: int = 0) -> set[str]:
    return {r[col] for r in rows}


def fasta_to_db(fasta: str, workdir: str | None = None, name: str = "db",
                reuse: bool = True) -> str:
    """Build (or reuse) a blast nucleotide DB for `fasta`.

    workdir: where to put the .nsq/.nhr/.nin files. If None, the DB goes in the
             SAME directory as `fasta` so it persists across runs (e.g. next to
             SPAdes' contigs.fasta in the per-k assembly dir).
    name:    DB prefix.
    reuse:   if True and `<workdir>/<name>.nsq` already exists, skip
             makeblastdb entirely. Setting False forces a rebuild.
    """
    if workdir is None:
        workdir = os.path.dirname(os.path.abspath(fasta)) or "."
    os.makedirs(workdir, exist_ok=True)
    db = os.path.join(workdir, name)
    if reuse and os.path.exists(db + ".nsq"):
        return db
    subprocess.run(["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return db
