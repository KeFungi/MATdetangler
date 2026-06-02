#!/usr/bin/env python3
"""Cross-K picker — consolidates per-K result.tsv files into a single
sample-level pick.

For each sample with `<outdir>/<k>/result.tsv` files produced by
`matdetangler.run_per_k`, this picks the best K according to:

  Priority (preferred first):
    1. closed_bubble
    2. open_bubble
    3. separate          (only if n_dedup == 2)
    4. single
    5. complexed         (only if n_dedup == 2; emits as allele pair)
    6. separate          (n_dedup != 2 — chimeric emission)
    7. complexed         (n_dedup != 2)
    8. no_var
    9. anything that errored

  Tie-break (within same priority class):
    a. Higher complete_var (True > False > None)
    b. Higher locus_coverage
    c. allele_cov closer to ½ × genome_cov (diploid signature)
       — for samples where ploidy=2 expected; uses
         abs(mean(allele_cov)/genome_cov − 0.5)
    d. Longer total basepair
    e. Lower n_cand (less graph tangle)

Writes:
  <outdir>/picks.tsv        sample-level summary (one row per sample,
                             with the K choice and its result fields)
  <outdir>/primary_alleles.fasta  concatenated FASTA of all winning
                             alleles across samples (headers carry
                             <sample>_<k>_<allele_name>)

Usage:
  python -m matdetangler.pick_k --outdir <sample_dir>
       picks across that one sample's per-k subdirs
  python -m matdetangler.pick_k --batch --results-root <results>
       picks across every <results>/<sample>/<k>/result.tsv tree
"""
from __future__ import annotations
import argparse, os, sys

PRIORITY_ORDER = {
    ("closed_bubble", "any"): 1,
    ("open_bubble",   "any"): 2,
    ("separate",      "div2"): 3,
    ("single",        "any"): 4,
    ("complexed",     "div2"): 5,
    ("separate",      "other"): 6,
    ("complexed",     "other"): 7,
    ("no_var",        "any"): 8,
}


def _priority(verdict: str, n_dedup: int) -> int:
    if verdict in ("closed_bubble", "open_bubble", "single"):
        return PRIORITY_ORDER[(verdict, "any")]
    if verdict in ("separate", "complexed"):
        cls = "div2" if n_dedup == 2 else "other"
        return PRIORITY_ORDER[(verdict, cls)]
    if verdict == "no_var":
        return PRIORITY_ORDER[(verdict, "any")]
    return 99


def _load_result(path: str) -> dict | None:
    if not os.path.exists(path): return None
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        first = fh.readline().rstrip("\n").split("\t")
    # error rows (no proper header) are 1-2 cols
    if len(header) < 4 or len(first) < 4: return None
    row = dict(zip(header, first))
    # numeric fields
    for k in ("basepair", "n_cand", "n_dedup", "n_hops_used"):
        try: row[k] = int(row.get(k, "") or 0)
        except ValueError: row[k] = 0
    for k in ("locus_coverage", "genome_cov"):
        try: row[k] = float(row.get(k, "") or 0.0)
        except ValueError: row[k] = 0.0
    row["complete_var"] = row.get("complete_var", "") == "True"
    return row


def _allele_cov_mean(row: dict) -> float:
    s = row.get("allele_cov", "")
    if not s or s == "-": return 0.0
    try:
        vals = [float(x) for x in s.split(",") if x]
    except ValueError:
        return 0.0
    return sum(vals) / len(vals) if vals else 0.0


def _score(row: dict) -> tuple:
    verdict = row.get("bubble_type", "")
    n_dedup = row.get("n_dedup", 0)
    prio = _priority(verdict, n_dedup)
    ac_mean = _allele_cov_mean(row)
    gen_cov = row.get("genome_cov", 0.0)
    # tie-break vector — lower is better
    cv_balance = abs((ac_mean / gen_cov) - 0.5) if gen_cov else 1.0
    return (
        prio,                           # 1. priority class (lower is better)
        not row.get("complete_var"),     # 2. True < False (prefer complete)
        -row.get("locus_coverage", 0),  # 3. higher locus_cov better
        cv_balance,                      # 4. ½-cov balance for diploid
        -row.get("basepair", 0),         # 5. longer total bp better
        row.get("n_cand", 0),            # 6. fewer raw candidates better
    )


def pick_for_sample(sample_dir: str) -> dict | None:
    """Pick the best K across the per-k subdirs of one sample dir."""
    candidates = []
    for k in sorted(os.listdir(sample_dir)):
        kd = os.path.join(sample_dir, k)
        if not os.path.isdir(kd): continue
        r = _load_result(os.path.join(kd, "result.tsv"))
        if r is None: continue
        r["_k"]      = k
        r["_kdir"]   = kd
        r["_fasta"]  = os.path.join(kd, "alleles.fasta")
        candidates.append(r)
    if not candidates: return None
    candidates.sort(key=_score)
    return candidates[0]


def _emit_combined(sample: str, win: dict, out_fa) -> int:
    """Append winning alleles to combined FASTA. Returns n records."""
    n = 0
    fa = win.get("_fasta")
    if not fa or not os.path.exists(fa): return 0
    with open(fa) as fh:
        for ln in fh:
            out_fa.write(ln)
            if ln.startswith(">"): n += 1
    return n


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", help="single sample dir (with k<X>/result.tsv subdirs)")
    ap.add_argument("--batch", action="store_true",
                     help="iterate over every sample under --results-root")
    ap.add_argument("--results-root",
                     help="root containing sample-named subdirs (for --batch)")
    ap.add_argument("--picks-tsv", default=None,
                     help="output picks.tsv path (default: <outdir>/picks.tsv "
                          "or <results-root>/picks.tsv for batch)")
    ap.add_argument("--primary-fasta", default=None,
                     help="output FASTA path (default: alongside picks.tsv)")
    args = ap.parse_args(argv)

    if args.batch:
        if not args.results_root:
            sys.stderr.write("--batch requires --results-root\n"); return 2
        sample_dirs = [
            (s, os.path.join(args.results_root, s))
            for s in sorted(os.listdir(args.results_root))
            if os.path.isdir(os.path.join(args.results_root, s))
        ]
        picks_tsv = args.picks_tsv or os.path.join(args.results_root, "picks.tsv")
        fasta_out = args.primary_fasta or os.path.join(args.results_root, "primary_alleles.fasta")
    else:
        if not args.outdir:
            sys.stderr.write("either --outdir or --batch+--results-root\n"); return 2
        sample = os.path.basename(args.outdir.rstrip("/"))
        sample_dirs = [(sample, args.outdir)]
        picks_tsv = args.picks_tsv or os.path.join(args.outdir, "picks.tsv")
        fasta_out = args.primary_fasta or os.path.join(args.outdir, "primary_alleles.fasta")

    # Legacy-compatible picks.tsv schema (one ROW PER ALLELE) — matches the
    # columns that downstream stages (graph_paths.py, summary_table.py)
    # expect. Sample-level info (k_chosen, bubble_type, etc.) is also written
    # to a parallel picks_summary.tsv with one row per sample.
    legacy_header = ["sample", "allele", "origin", "k", "type", "len",
                     "from_contig", "segments", "cov", "n_variable_genes",
                     "has_both_flanks", "is_degHD"]
    summary_header = ["sample", "k_chosen", "bubble_type", "n_dedup",
                      "complete_var", "complete_locus", "locus_coverage",
                      "basepair", "genome_cov", "allele_cov", "n_cand",
                      "extend_bounds", "components", "all_k_tried"]

    summary_tsv = picks_tsv.replace("picks.tsv", "picks_summary.tsv")
    if summary_tsv == picks_tsv: summary_tsv += ".summary"

    n_samples = 0; n_picked = 0
    with open(picks_tsv, "w") as ptsv, \
         open(summary_tsv, "w") as stsv, \
         open(fasta_out, "w") as fa_out:
        ptsv.write("\t".join(legacy_header) + "\n")
        stsv.write("\t".join(summary_header) + "\n")
        for sample, sdir in sample_dirs:
            win = pick_for_sample(sdir)
            tried = ",".join(
                k for k in sorted(os.listdir(sdir))
                if os.path.isdir(os.path.join(sdir, k))
                and os.path.exists(os.path.join(sdir, k, "result.tsv"))
            ) if os.path.isdir(sdir) else ""
            if win is None:
                stsv.write("\t".join([sample, "-", "NO_RESULT"] + [""] * 11) + "\n")
                continue
            n_samples += 1

            # summary row (sample-level)
            stsv.write("\t".join([
                sample, win["_k"], win.get("bubble_type", ""),
                str(win.get("n_dedup", 0)),
                str(win.get("complete_var", False)),
                str(win.get("complete_locus", "")),
                f"{win.get('locus_coverage', 0):.3f}",
                str(win.get("basepair", 0)),
                f"{win.get('genome_cov', 0):.2f}",
                win.get("allele_cov", "-"),
                str(win.get("n_cand", 0)),
                win.get("extend_bounds", "-"),
                win.get("components", "-"),
                tried,
            ]) + "\n")

            # legacy per-allele rows
            allele_names = (win.get("components", "") or "").split(",")
            allele_lens  = (win.get("allele_lens",   "") or "").split(",")
            cov_str      = (win.get("allele_cov",    "") or "").split(",")
            ext_pairs    = (win.get("extend_bounds", "") or "").split(";")
            segs_field   = win.get("segments", "-")
            origin = ("path" if win.get("bubble_type") in ("closed_bubble",
                                                            "open_bubble", "single")
                            else "graph")  # complexed / separate → "graph"
            type_lbl = ("complete" if win.get("complete_var") else "partial")
            for i, name in enumerate(allele_names):
                name = name.strip()
                if not name: continue
                try: ln = int((allele_lens[i] if i < len(allele_lens) else "0") or "0")
                except ValueError: ln = 0
                cov_v = (cov_str[i] if i < len(cov_str) else "") or "0"
                ext = (ext_pairs[i] if i < len(ext_pairs) else "-:-")
                L_ok, R_ok = "-" not in ext.split(":")[0], "-" not in ext.split(":")[-1]
                # n_variable_genes "found/expected" — sample-level approximation:
                # complete_var=True → "K/K" where K is len(found_var_tags) (≈expected).
                # Else → estimate found from comma-joined list.
                # (Approximation: every emitted allele inherits the sample-level count.)
                # Without per-allele recheck, just use sample-level found_var_tags count.
                nv_found = len((win.get("found_var_tags", "") or "").split(",")) if win.get("found_var_tags") and win.get("found_var_tags") != "-" else 0
                n_variable = f"{nv_found}/{nv_found}" if win.get("complete_var") else f"?/{nv_found}"
                has_both_flanks = "True" if (L_ok and R_ok) else "False"
                ptsv.write("\t".join([
                    sample, name, origin, win["_k"], type_lbl, str(ln),
                    "-", segs_field, cov_v, n_variable, has_both_flanks, "False",
                ]) + "\n")
                n_picked += 1

            _emit_combined(sample, win, fa_out)

    print(f"[pick_k] {n_samples} samples picked, {n_picked} allele rows → {picks_tsv}")
    print(f"[pick_k] summary → {summary_tsv}")
    print(f"[pick_k] combined FASTA → {fasta_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
