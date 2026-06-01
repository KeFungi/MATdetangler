"""Anchor search — locate sequences that hit any variable gene (HD) or either flank.

Works on either input source, selected by --source:

    --source contigs   step 2.1: search the per-k SPAdes contigs.fasta
                       output: anchor_contig.fasta + anchor_contig.ann.tsv
                       record name: "<sample>__bubble_<k>_<contig_id>"

    --source segments  step 2.2 (fallback): search the per-k GFA segments directly.
                       Segment sequences are pulled from the S-lines of the GFA into
                       a temp fasta and the same blast logic runs against them.
                       output: anchor_segments.fasta + anchor_segments.ann.tsv
                       record name: "<sample>__seg_<k>_<segment_id>"

A record is an "anchor" if it carries either of:
  - any VARIABLE GENE (HD) via the 3-tier fallback chain (tblastn -> blastn -> tblastx).
  - either FLANK (flankL / flankR) via blastn.

HD tier chain — each lower tier is more divergence-tolerant and slower than the one
above, and only fires when the previous tier did not already surface enough complete
anchors (anchors carrying ALL variable genes):

  tier 1  tblastn(variable_proteins → records)  ALWAYS    protein query (curated)
                                                          -> 6-frame translated NT subject.
  tier 2  blastn(variable_nt → records)         FALLBACK  if tier 1 found < --min-complete
                                                          anchors carrying ALL variable
                                                          genes. Stricter NT identity.
  tier 3  tblastx(variable_nt → records)        FALLBACK  if tiers 1+2 STILL did not reach
                                                          --min-complete. Both sides 6-frame
                                                          translated; ~36x slower.

Flank queries always run — they're cheap, and a pure-flank record still anchors the
downstream graph walk in step 3.
"""
from __future__ import annotations
import os, sys, argparse, tempfile
from . import blast_utils as bu
from .paths import spades_k_paths


# ---------- shared helpers ----------

def _record_lengths(fa: str) -> dict[str, int]:
    """Per-record length from a fasta (works for contigs OR segment fasta)."""
    out, name, n = {}, None, 0
    for ln in open(fa):
        if ln.startswith(">"):
            if name: out[name] = n
            name = ln[1:].split()[0]; n = 0
        else:
            n += len(ln.strip())
    if name: out[name] = n
    return out


def _extract_record(fa: str, name: str) -> str:
    """Read one record's sequence out of a fasta without needing samtools."""
    out, take = [], False
    for ln in open(fa):
        if ln.startswith(">"):
            cur = ln[1:].split()[0]
            if take: break
            take = (cur == name)
        elif take:
            out.append(ln.strip())
    return "".join(out)


def _count_complete(gene_per_rec: dict[str, set[str]],
                     lens: dict[str, int],
                     nvar_total: int,
                     min_len: int) -> int:
    """How many records carry ALL variable genes AND meet the min length threshold."""
    return sum(1 for c, genes in gene_per_rec.items()
                  if len(genes) == nvar_total and lens.get(c, 0) >= min_len)


# ---------- GFA segments -> fasta ----------

def _segments_fasta_from_gfa(gfa_path: str, out_fa: str) -> int:
    """Extract every S-line in a GFA into a fasta. Returns the number of segments."""
    n = 0
    with open(out_fa, "w") as out:
        for ln in open(gfa_path):
            if ln[:2] != "S\t": continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3: continue
            sid, seq = f[1], f[2]
            if not seq or seq == "*": continue
            out.write(f">{sid}\n")
            for i in range(0, len(seq), 80): out.write(seq[i:i + 80] + "\n")
            n += 1
    return n


# ---------- core per-fasta search (shared between contigs + segments) ----------

def _search_one_fasta(subject_fa: str, lens: dict[str, int],
                      queries_dir: str, min_len: int,
                      blastn_pid: float, blastn_minlen: int,
                      tblastn_pid: float, tblastn_aa: int,
                      tblastx_pid: float, tblastx_aa: int,
                      flank_pid: float, flank_minlen: int,
                      min_complete_tblastn: int,
                      nvar_total: int, threads: int,
                      tmpdir: str, db_name: str,
                      include_flanks: bool = False,
                      blast_out_dir: str | None = None,
                      blast_tag: str = ""):
    """Run the HD tier chain + flank blasts against one subject fasta. Returns
    anchors as-is; repeat-aware handling is in step 3's BFS expansion (paths can
    enter repeat segments, but the BFS doesn't expand neighbors FROM a repeat
    segment - prevents combinatorial blow-up at high-copy regions).

    Returns (hd_anchors, flank_anchors, gene_per, flanks_per, tier_tag,
              n_complete_after).
    """
    nt_panel   = os.path.join(queries_dir, "variable_nt.fasta")
    prot_panel = os.path.join(queries_dir, "variable_proteins.fasta")
    flankL_fa  = os.path.join(queries_dir, "flankL.fasta")
    flankR_fa  = os.path.join(queries_dir, "flankR.fasta")
    db = bu.fasta_to_db(subject_fa, tmpdir, name=db_name)
    gene_per: dict[str, set[str]] = {}
    flanks_per: dict[str, set[str]] = {}
    # Build out_tsv paths for the 5 cached blast calls (or None to skip caching).
    # All cached blast tsvs share the "blast_" prefix so the wrapper's --re-blast
    # flag can wipe them with a single glob.
    def _out(prog: str, query: str) -> str | None:
        if not blast_out_dir or not blast_tag: return None
        return os.path.join(blast_out_dir, f"blast_{query}_{prog}_{blast_tag}.tsv")
    # tier 1 (HD): tblastn (always)
    tp_rows = bu.tblastn_hits(prot_panel, db, min_pid=tblastn_pid,
                                min_aa=tblastn_aa, threads=threads,
                                out_tsv=_out("tblastn", "hd"))
    for r in tp_rows: gene_per.setdefault(r[0], set()).add(r[1])
    hp = bu.hits_ids(tp_rows, col=0)
    n_complete_t1 = _count_complete(gene_per, lens, nvar_total, min_len)
    # tier 2 (HD): blastn(variable_nt) -- only if tier 1 short
    hn: set[str] = set()
    tier2_fired = False
    n_complete_t2 = n_complete_t1
    if n_complete_t1 < min_complete_tblastn:
        tier2_fired = True
        bn_rows = bu.blastn_hits(nt_panel, db, min_pid=blastn_pid,
                                   min_len=blastn_minlen, threads=threads,
                                   extra_outfmt="qseqid",
                                   out_tsv=_out("blastn", "hd"))
        for r in bn_rows:
            if len(r) >= 4: gene_per.setdefault(r[0], set()).add(r[3])
        hn = bu.hits_ids(bn_rows, col=0)
        n_complete_t2 = _count_complete(gene_per, lens, nvar_total, min_len)
    # tier 3 (HD): tblastx(variable_nt) -- only if tiers 1+2 STILL short
    hx: set[str] = set()
    tier3_fired = False
    n_complete_t3 = n_complete_t2
    if n_complete_t2 < min_complete_tblastn:
        tier3_fired = True
        tx_rows = bu.tblastx_hits(nt_panel, db, min_pid=tblastx_pid,
                                    min_aa=tblastx_aa, threads=threads,
                                    out_tsv=_out("tblastx", "hd"))
        for r in tx_rows: gene_per.setdefault(r[0], set()).add(r[1])
        hx = bu.hits_ids(tx_rows, col=0)
        n_complete_t3 = _count_complete(gene_per, lens, nvar_total, min_len)
    # Flank anchors — OPT IN only. Default off for production runs since the BFS
    # in step 3 reaches flank-bearing segs via L-line adjacency from HD seeds anyway,
    # and starting BFS from a flank-only big contig (e.g. 25 kb) pulls in noise.
    # Pass include_flanks=True (--include-flanks) under --debug to add them.
    hL: set[str] = set()
    hR: set[str] = set()
    if include_flanks:
        if os.path.exists(flankL_fa):
            hL = bu.hits_ids(bu.blastn_hits(flankL_fa, db, min_pid=flank_pid,
                                              min_len=flank_minlen, threads=threads,
                                              out_tsv=_out("blastn", "flankL")))
            for c in hL: flanks_per.setdefault(c, set()).add("flankL")
        if os.path.exists(flankR_fa):
            hR = bu.hits_ids(bu.blastn_hits(flankR_fa, db, min_pid=flank_pid,
                                              min_len=flank_minlen, threads=threads,
                                              out_tsv=_out("blastn", "flankR")))
            for c in hR: flanks_per.setdefault(c, set()).add("flankR")
    hd_anchors    = (hp | hn | hx) - {""}
    flank_anchors = (hL | hR) - {""}
    tag = "tblastn alone"
    if tier3_fired:    tag = f"tblastn+blastn+tblastx ({n_complete_t1}/{n_complete_t2}/{n_complete_t3} complete after each tier)"
    elif tier2_fired:  tag = f"tblastn+blastn ({n_complete_t1}/{n_complete_t2} complete after each tier)"
    return hd_anchors, flank_anchors, gene_per, flanks_per, tag, n_complete_t3


# ---------- entry point ----------

def run(sample: str, source: str, spades_dir: str, queries_dir: str,
        ks: list[str], outdir: str,
        min_len: int, blastn_pid: float, blastn_minlen: int,
        tblastn_pid: float, tblastn_aa: int, threads: int,
        tblastx_pid: float = 30.0, tblastx_aa: int = 50,
        flank_pid: float = 85.0, flank_minlen: int = 100,
        min_complete_tblastn: int = 2,
        include_flanks: bool = False,
        blast_out_dir: str | None = None) -> str:
    """Per-k anchor search on either contigs or GFA segments.

    `source` is one of:
        "contigs"  -> step 2.1: search contigs.fasta per k.
                      Output filenames: anchor_contig.fasta + anchor_contig.ann.tsv
                      Record name template: "<sample>__bubble_<k>_<contig_id>"
        "segments" -> step 2.2: search GFA segments directly.
                      Output filenames: anchor_segments.fasta + anchor_segments.ann.tsv
                      Record name template: "<sample>__seg_<k>_<segment_id>"
    """
    if source == "contigs":
        out_basename = "anchor_contig"
        name_kind = "bubble"
    elif source == "segments":
        out_basename = "anchor_segments"
        name_kind = "seg"
    else:
        raise ValueError(f"source must be 'contigs' or 'segments', got {source!r}")

    os.makedirs(outdir, exist_ok=True)
    cand_fa = os.path.join(outdir, f"{out_basename}.fasta")
    cand_an = os.path.join(outdir, f"{out_basename}.ann.tsv")
    prot_panel = os.path.join(queries_dir, "variable_proteins.fasta")
    open(cand_fa, "w").close(); open(cand_an, "w").close()
    nvar_total = sum(1 for ln in open(prot_panel) if ln.startswith(">"))

    # Persist BLAST DBs (and the GFA-extracted segments fasta) inside the
    # per-sample outdir so re-runs reuse them via fasta_to_db's `reuse=True`
    # check. Previously these lived in a tempfile.TemporaryDirectory and got
    # rebuilt every invocation — ~30-60 s per k for a fungal contigs/segments DB.
    blast_db_dir = os.path.join(outdir, "blast_db")
    os.makedirs(blast_db_dir, exist_ok=True)
    n_total = 0
    for k in ks:
        ctg, gfa = spades_k_paths(spades_dir, k)
        if source == "contigs":
            subject_fa = ctg
            if not subject_fa:
                print(f"  [{k}] no contigs.fasta found under {spades_dir}; skipping"); continue
        else:  # segments
            if not gfa:
                print(f"  [{k}] no GFA found under {spades_dir}; skipping"); continue
        lens = _record_lengths(subject_fa) if source == "contigs" else {}
        if source == "segments":
            subject_fa = os.path.join(blast_db_dir, f"segments_{k}.fa")
            if not os.path.exists(subject_fa):
                n_seg = _segments_fasta_from_gfa(gfa, subject_fa)
                print(f"  [{k}] extracted {n_seg} segments from {gfa} -> {subject_fa}")
            else:
                print(f"  [{k}] reusing extracted segments fasta {subject_fa}")
            lens = _record_lengths(subject_fa)
        hd, fl, gene_per, flanks_per, tag, _ = _search_one_fasta(
            subject_fa, lens, queries_dir, min_len,
            blastn_pid, blastn_minlen, tblastn_pid, tblastn_aa,
            tblastx_pid, tblastx_aa, flank_pid, flank_minlen,
            min_complete_tblastn, nvar_total, threads,
            tmpdir=blast_db_dir, db_name=f"{source}_{k}_db",
            include_flanks=include_flanks,
            blast_out_dir=blast_out_dir,
            blast_tag=f"{source}_{k}")
        cands = sorted((hd | fl))
        n_kept = 0
        for c in cands:
            if lens.get(c, 0) < min_len:
                continue
            seq = _extract_record(subject_fa, c)
            nm = f"{sample}__{name_kind}_{k}_{c}"
            with open(cand_fa, "a") as o:
                o.write(f">{nm}\n")
                for i in range(0, len(seq), 80): o.write(seq[i:i + 80] + "\n")
            genes  = ",".join(sorted(gene_per.get(c, set())))    or "-"
            flanks = ",".join(sorted(flanks_per.get(c, set())))  or "-"
            with open(cand_an, "a") as o:
                o.write(f"{nm}\t{lens[c]}\t{k}\t{name_kind}\t{genes}\t{flanks}\n")
            n_total += 1; n_kept += 1
        print(f"  [{k}] anchors: {n_kept}  hd={len(hd)}  flank={len(fl)}  threshold={min_complete_tblastn}  used={tag}")
    print(f"[anchor_search:{source}] {sample}: {n_total} anchors -> {cand_fa}")
    return cand_fa


def _cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--sample", required=True)
    p.add_argument("--source", choices=("contigs", "segments"), default="contigs",
                   help="contigs (step 2.1) or segments (step 2.2 GFA fallback)")
    p.add_argument("--spades-dir", required=True)
    p.add_argument("--queries-dir", required=True, help="dir written by input_process.write_queries")
    p.add_argument("--ks", default="k21,k33,k55")
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-len", type=int, default=2000,
                   help="min subject length to KEEP as an anchor. For --source segments, you "
                        "probably want a smaller value (e.g. 200) since GFA segments are often "
                        "much shorter than full contigs.")
    p.add_argument("--blastn-pid", type=float, default=80.0)
    p.add_argument("--blastn-minlen", type=int, default=100)
    p.add_argument("--tblastn-pid", type=float, default=30.0)
    p.add_argument("--tblastn-aa", type=int, default=50)
    p.add_argument("--tblastx-pid", type=float, default=30.0,
                   help="(tier-3 HD fallback) min pid for tblastx hits; default 30.0")
    p.add_argument("--tblastx-aa",  type=int,   default=50,
                   help="(tier-3 HD fallback) min aln length (aa) for tblastx hits; default 50")
    p.add_argument("--flank-pid",    type=float, default=85.0,
                   help="min pid for flankL / flankR blastn hits; default 85.0")
    p.add_argument("--flank-minlen", type=int, default=100,
                   help="min aln length (bp) for flankL / flankR blastn hits; default 100")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--min-complete-tblastn", type=int, default=2,
                   help="per k, fall through to the next HD tier only if fewer than this many "
                        "records carry ALL variable genes after the current tier (default 2).")
    p.add_argument("--blast-out-dir", default=None,
                   help="If set, the heavy contigs/segments search blasts (HD tblastn, HD blastn, "
                        "HD tblastx, flankL/R blastn) write their hit tables to this directory as "
                        "named TSVs (e.g. hd_tblastn_contigs_k45.tsv). On re-run, an existing file "
                        "is reused as-is (file-existence-only check). The user is responsible for "
                        "deleting stale files when blast parameters change — the wrapper exposes "
                        "--no-cached-blast to wipe them.")
    p.add_argument("--include-flanks", action="store_true",
                   help="Also blastn flankL.fasta and flankR.fasta against the subject and "
                        "include flank-only records as anchors. Default OFF — the BFS in step 3 "
                        "reaches flank-bearing segs via L-line adjacency from HD seeds anyway, "
                        "and pure-flank big contigs pull in noise. Use under --debug to inspect "
                        "flank coverage at the contig/segment level.")
    a = p.parse_args(argv)
    run(a.sample, a.source, a.spades_dir, a.queries_dir, a.ks.split(","), a.outdir,
        a.min_len, a.blastn_pid, a.blastn_minlen,
        a.tblastn_pid, a.tblastn_aa, a.threads,
        tblastx_pid=a.tblastx_pid, tblastx_aa=a.tblastx_aa,
        flank_pid=a.flank_pid, flank_minlen=a.flank_minlen,
        min_complete_tblastn=a.min_complete_tblastn,
        include_flanks=a.include_flanks,
        blast_out_dir=a.blast_out_dir)

if __name__ == "__main__":
    _cli()
