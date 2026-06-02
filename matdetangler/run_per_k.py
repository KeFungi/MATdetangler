#!/usr/bin/env python3
"""Production per-K driver — replaces graph_path_search + pick_alleles.

For one (sample, k):
  1. Extract S-line segments from the k-GFA → temp fasta + BLAST DB
  2. blastn flankL/R against GFA segments (full outfmt-6, 12 cols)
  3. tblastn HD proteins against GFA segments (full outfmt-6)
  4. Aggregate via labeler → seg_label_hits.tsv
  5. Run per_k_caller.find_alleles → verdict + alleles + extend_bounds
  6. Write result_k<k>.tsv + alleles_k<k>.fasta into sample outdir

Usage (CLI):
  python -m matdetangler.run_per_k \\
      --sample SAMPLE --k k53 \\
      --gfa  spades/SAMPLE/k53/assembly_graph_after_simplification.gfa \\
      --queries-dir results/SAMPLE/queries \\
      --hd-proteins examples/Suilu_locus/Suilu4_HDs.fasta \\
      --cov-file results/SAMPLE/genome_cov_spades_k53.txt \\
      --outdir results/SAMPLE \\
      --threads 4
"""
from __future__ import annotations
import argparse, os, subprocess, sys, tempfile, traceback

from matdetangler.graph_classifier.labeler import emit_seg_label_hits
from matdetangler.graph_classifier.per_k_caller import find_alleles, write_fasta


def _extract_segs_from_gfa(gfa: str, out_fa: str) -> int:
    n = 0
    with open(gfa) as fin, open(out_fa, "w") as fout:
        for ln in fin:
            if ln.startswith("S\t"):
                f = ln.rstrip("\n").split("\t", 3)
                if len(f) >= 3:
                    fout.write(f">{f[1]}\n{f[2]}\n")
                    n += 1
    return n


def _blast(query_fa: str, db: str, out_tsv: str, threads: int = 4,
            kind: str = "blastn", evalue: str = "1e-5") -> int:
    cmd = [kind, "-query", query_fa, "-db", db, "-evalue", evalue,
           "-outfmt", "6", "-num_threads", str(threads)]
    if kind == "blastn":
        cmd.extend(["-dust", "no"])
    with open(out_tsv, "w") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.stderr.write(f"[run_per_k] {kind} failed:\n{p.stderr}\n")
    return p.returncode


def _expected_var_tags(hd_proteins: str) -> set[str]:
    tags: set[str] = set()
    if not os.path.exists(hd_proteins): return tags
    with open(hd_proteins) as fh:
        for ln in fh:
            if ln.startswith(">"):
                tags.add(ln[1:].strip().split()[0])
    return tags


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Per-K allele caller (production).")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--k",      required=True, help='e.g. "k53"')
    ap.add_argument("--gfa",    required=True)
    ap.add_argument("--queries-dir", required=True,
                     help="dir with flankL.fasta + flankR.fasta")
    ap.add_argument("--hd-proteins", required=True,
                     help="HD protein fasta for tblastn (locus-trim + completeness)")
    ap.add_argument("--cov-file", required=True,
                     help="single-line genome_cov_spades_k<k>.txt")
    ap.add_argument("--outdir", required=True,
                     help="per-sample output root; writes outdir/<k>/...")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--init-nhop", type=int, default=3)
    ap.add_argument("--max-nhop",  type=int, default=10)
    ap.add_argument("--lo-mult",   type=float, default=0.25)
    ap.add_argument("--hi-mult",   type=float, default=2.0)
    ap.add_argument("--locus-padding", type=int, default=1500)
    ap.add_argument("--divergence-threshold", type=float, default=0.05)
    args = ap.parse_args(argv)

    out_k = os.path.join(args.outdir, args.k)
    os.makedirs(out_k, exist_ok=True)
    flankL = os.path.join(args.queries_dir, "flankL.fasta")
    flankR = os.path.join(args.queries_dir, "flankR.fasta")
    if not all(os.path.exists(p) for p in (args.gfa, flankL, flankR, args.cov_file)):
        sys.stderr.write(f"[{args.sample}] MISSING_INPUT — gfa/flank/cov\n")
        with open(os.path.join(out_k, "result.tsv"), "w") as fh:
            fh.write(f"{args.sample}\tMISSING_INPUT\n")
        return 0

    genome_cov = float(open(args.cov_file).read().strip().split()[0])
    expected_tags = _expected_var_tags(args.hd_proteins)

    flankL_tsv = os.path.join(out_k, "flankL_blastn.tsv")
    flankR_tsv = os.path.join(out_k, "flankR_blastn.tsv")
    hd_tsv     = os.path.join(out_k, "HD_tblastn.tsv")
    hits_tsv   = os.path.join(out_k, "seg_label_hits.tsv")
    fasta      = os.path.join(out_k, "alleles.fasta")
    result_tsv = os.path.join(out_k, "result.tsv")

    with tempfile.TemporaryDirectory() as td:
        seg_fa = os.path.join(td, "segs.fa")
        db     = os.path.join(td, "segs_db")
        n_segs = _extract_segs_from_gfa(args.gfa, seg_fa)
        sys.stderr.write(f"[{args.sample}/{args.k}] n_segs={n_segs}\n")
        subprocess.run(["makeblastdb", "-in", seg_fa, "-dbtype", "nucl",
                         "-out", db], check=True, stdout=subprocess.DEVNULL)
        _blast(flankL, db, flankL_tsv, threads=args.threads, kind="blastn")
        _blast(flankR, db, flankR_tsv, threads=args.threads, kind="blastn")
        _blast(args.hd_proteins, db, hd_tsv, threads=args.threads, kind="tblastn")

    try:
        n_rows = emit_seg_label_hits(
            args.gfa,
            [(flankL_tsv, "flank"), (flankR_tsv, "flank"), (hd_tsv, "var")],
            hits_tsv, min_alnlen=100,
        )
    except Exception as e:
        with open(result_tsv, "w") as fh:
            fh.write(f"{args.sample}\tLABELER_FAIL\t{type(e).__name__}: {e}\n")
        traceback.print_exc(); return 1

    try:
        res = find_alleles(
            hits_tsv, args.gfa,
            genome_cov=genome_cov,
            init_nhop=args.init_nhop, max_nhop=args.max_nhop,
            lo_mult=args.lo_mult, hi_mult=args.hi_mult,
            k=args.k, var_proteins_ref=args.hd_proteins,
            expected_var_tags=expected_tags or None,
            locus_padding=args.locus_padding,
            divergence_threshold=args.divergence_threshold,
        )
    except Exception as e:
        with open(result_tsv, "w") as fh:
            fh.write(f"{args.sample}\tCALLER_FAIL\t{type(e).__name__}: {e}\n")
        traceback.print_exc(); return 1

    write_fasta(res["alleles"], fasta, sample=f"{args.sample}_{args.k}")

    allele_names = ",".join(res.get("component_list", [])) or "-"
    allele_lens  = ",".join(str(len(s)) for _, s in res["alleles"]) or "-"
    seg_list     = ",".join(res.get("segments", [])) or "-"
    seg_lab_list = ",".join(f"{s}:{lab}" for s, lab in res.get("segments_labeled", [])) or "-"
    allele_cov   = ",".join(f"{c:.1f}" for c in res.get("allele_cov", [])) or "-"
    extend_bounds = ";".join(f"{L or '-'}:{R or '-'}"
                              for (L, R) in res.get("extend_bounds", [])) or "-"

    header = ["sample", "k", "bubble_type", "components", "complete_var",
              "complete_locus", "locus_coverage", "basepair", "genome_cov",
              "allele_cov", "n_cand", "n_dedup", "divergent", "n_hops_used",
              "phase", "cov_filter_used", "allele_lens", "segments",
              "segments_labeled", "found_var_tags", "extend_bounds"]
    row = [
        args.sample, str(res.get("k") or "?"),
        str(res.get("bubble_type", res.get("verdict"))),
        allele_names,
        str(res.get("complete_var")), str(res.get("complete_locus")),
        f"{res.get('locus_coverage') or 0:.3f}",
        str(res.get("basepair", 0)),
        f"{res.get('genome_cov') or 0:.2f}", allele_cov,
        str(res.get("n_candidates", 0)), str(res.get("n_after_dedup", 0)),
        str(int(res.get("divergent", False))),
        str(res.get("n_hops_used") or ""),
        str(res.get("phase") or ""),
        str(int(res.get("cov_filter_used") or False)),
        allele_lens, seg_list, seg_lab_list,
        ",".join(res.get("found_var_tags", [])) or "-",
        extend_bounds,
    ]
    with open(result_tsv, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(row) + "\n")

    print(f"[{args.sample}/{args.k}] verdict={res.get('verdict')} "
          f"n_dedup={res.get('n_after_dedup')} basepair={res.get('basepair')} "
          f"complete_var={res.get('complete_var')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
