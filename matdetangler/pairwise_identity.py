"""Compute allele1-vs-allele2 identity (pairwise MAFFT alignment) for a sample.

Two flavors:
  mafft_pair       — identity / aln_frac over the WHOLE aligned sequence (flanks +
                     HD core combined). Same metric the cross-sample mating-type
                     clustering uses.
  mafft_pair_core  — aligns the WHOLE sequences with MAFFT (so the conserved flanks
                     anchor the alignment) but computes identity / aln_frac ONLY over
                     the alignment columns that fall inside the HD core span on
                     either allele. Use this when comparing alleles of the same locus
                     and you don't want the highly-conserved flanks to mask the real
                     HD-level divergence.

A pair is "distinct" iff NOT (alnid >= 0.95 AND aln >= 0.80).
"""
from __future__ import annotations
import os, subprocess, tempfile, argparse
from .input_process import read_fasta

MAFFT = "mafft"


def _mafft_align_pair(s1: str, s2: str) -> tuple[str, str]:
    """Run MAFFT --auto on two sequences; return their aligned strings (with gaps).
    Returns ("", "") on failure or empty input.
    """
    if not s1 or not s2:
        return "", ""
    with tempfile.TemporaryDirectory() as t:
        fa = os.path.join(t, "p.fa")
        with open(fa, "w") as o:
            o.write(f">a\n{s1}\n>b\n{s2}\n")
        try:
            # --adjustdirection: MAFFT picks the better of (s2, rc(s2)) per record.
            # Without it, candidates on opposite strands of the same allele score
            # as falsely distinct.
            out = subprocess.run([MAFFT, "--adjustdirection", "--auto", "--quiet", fa],
                                  capture_output=True, text=True, timeout=300).stdout
        except Exception:
            return "", ""
    seqs, name = {}, None
    for ln in out.splitlines():
        if ln.startswith(">"):
            nm = ln[1:].split()[0]
            # --adjustdirection prefixes reverse-complemented records with "_R_".
            # Strip it so the caller's lookups by "a"/"b" still work.
            if nm.startswith("_R_"): nm = nm[3:]
            name = nm; seqs[name] = []
        elif name is not None:
            seqs[name].append(ln.strip())
    if "a" not in seqs or "b" not in seqs:
        return "", ""
    return "".join(seqs["a"]), "".join(seqs["b"])


def _score_columns(sa: str, sb: str, cols: range | set | None,
                    denom_len: int) -> tuple[float, float]:
    """Given two aligned strings and a set of column indices to score (None = ALL
    columns), return (identity, aligned-fraction-over-denom).
    """
    if not sa or not sb:
        return 0.0, 0.0
    L = min(len(sa), len(sb))
    iter_cols = range(L) if cols is None else (c for c in cols if c < L)
    m = al = 0
    for col in iter_cols:
        x, y = sa[col], sb[col]
        if x != "-" and y != "-":
            al += 1
            if x.upper() == y.upper(): m += 1
    if al == 0:
        return 0.0, 0.0
    alnid = m / al
    frac  = al / denom_len if denom_len else 0.0
    return alnid, frac


def mafft_pair(s1: str, s2: str) -> tuple[float, float]:
    """Return (alnid, aln_frac_of_shorter) over the WHOLE aligned sequences."""
    sa, sb = _mafft_align_pair(s1, s2)
    return _score_columns(sa, sb, cols=None, denom_len=min(len(s1), len(s2)))


# In-process memo: detect_core_span's result depends only on (seq, proteins_fa,
# min_pid, min_aa). Within ONE pipeline run, the same candidate sequence is
# passed to detect_core_span multiple times (once per cluster_alleles call:
# widen-loop Stage B at every hop, plus final emission). Without memoization
# we'd do 75-100 redundant tblastn calls per k. Memoization makes it ONE call
# per unique sequence per run. Strings hash cheaply in CPython.
_DETECT_CORE_SPAN_MEMO: dict[tuple[int, str, float, int], tuple[int, int]] = {}


def detect_core_span(seq: str, variable_proteins_fa: str,
                     min_pid: float = 30.0, min_aa: int = 50) -> tuple[int, int]:
    """tblastn(variable_proteins → seq); return the HD core span on `seq` as
    (start, end), both 1-based inclusive. Span = [min(start), max(end)] over all
    qualifying tblastn HSPs. Returns (1, len(seq)) if no qualifying hits.

    Memoized within the process — the same (seq, proteins_fa, thresholds) input
    runs blast at most once per run.
    """
    if not seq or not os.path.exists(variable_proteins_fa):
        return (1, len(seq) if seq else 0)
    # Memo key uses id() of the seq string as a fast first check, with the
    # actual content kept in the key so different equal-content strings still
    # hit. (id() alone is unsafe across runs but fine for the dict-key tuple.)
    memo_key = (len(seq), seq, min_pid, min_aa)
    cached = _DETECT_CORE_SPAN_MEMO.get(memo_key)
    if cached is not None: return cached
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "s.fa")
        with open(sf, "w") as o: o.write(f">x\n{seq}\n")
        db = os.path.join(t, "s_db")
        subprocess.run(["makeblastdb", "-in", sf, "-dbtype", "nucl", "-out", db],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        out = subprocess.run(
            ["tblastn", "-query", variable_proteins_fa, "-db", db,
             "-evalue", "1e-5", "-outfmt", "6 sseqid pident length sstart send"],
            stdout=subprocess.PIPE, text=True).stdout
    starts, ends = [], []
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 5: continue
        pid, ln_aa = float(f[1]), int(f[2])
        if pid < min_pid or ln_aa < min_aa: continue
        ss, se = int(f[3]), int(f[4])
        starts.append(min(ss, se)); ends.append(max(ss, se))
    result = (1, len(seq)) if not starts else (min(starts), max(ends))
    _DETECT_CORE_SPAN_MEMO[memo_key] = result
    return result


def _seq_pos_to_align_cols(aligned: str, start_1b: int, end_1b: int) -> set[int]:
    """Map a 1-based inclusive sequence span (start..end on the un-gapped sequence)
    to the set of column indices it occupies in the aligned string.
    """
    cols: set[int] = set()
    seq_pos = 0
    for col, ch in enumerate(aligned):
        if ch == "-": continue
        seq_pos += 1
        if seq_pos < start_1b: continue
        if seq_pos > end_1b: break
        cols.add(col)
    return cols


def mafft_pair_core(s1: str, s2: str,
                    core1: tuple[int, int], core2: tuple[int, int]) -> tuple[float, float]:
    """Align s1 vs s2 with MAFFT, then score ONLY over the union of alignment columns
    that fall inside core1 (mapped through gaps on aligned-s1) or core2 (through
    gaps on aligned-s2). Returns (alnid, aln_frac) where aln_frac uses the SHORTER
    core length as denominator — so it directly reflects HD-core coverage.
    """
    sa, sb = _mafft_align_pair(s1, s2)
    if not sa: return 0.0, 0.0
    cols_a = _seq_pos_to_align_cols(sa, core1[0], core1[1])
    cols_b = _seq_pos_to_align_cols(sb, core2[0], core2[1])
    cols = cols_a | cols_b
    core_a_len = core1[1] - core1[0] + 1
    core_b_len = core2[1] - core2[0] + 1
    denom = min(core_a_len, core_b_len) if core_a_len and core_b_len else 0
    return _score_columns(sa, sb, cols, denom_len=denom)

def run(primary_alleles_fa: str, out_tsv: str | None = None,
        queries_dir: str | None = None,
        locus_ref_fa: str | None = None) -> dict:
    """Pairwise identity for the picks.

    When `queries_dir` + `locus_ref_fa` are supplied, uses the SAME 3-step
    protocol as step 3 `cluster_alleles` (tblastn-trim each candidate to HD-core
    ± pad, MAFFT MSA with the trimmed locus reference, identity computed only
    on HD-core columns projected from the reference's ungapped span). This
    matches the picker's and the defensive-dedup's metric.

    Falls back to whole-allele `mafft_pair` when locus_ref or queries_dir is
    missing (legacy behavior).
    """
    alleles = read_fasta(primary_alleles_fa)
    names = list(alleles)
    res = {"a1": names[0] if names else None,
           "a2": names[1] if len(names) > 1 else None,
           "len_a1": len(alleles[names[0]]) if names else 0,
           "len_a2": len(alleles[names[1]]) if len(names) > 1 else 0,
           "id_pct": None, "aln_frac": None, "distinct": None,
           "metric": None}
    if len(names) >= 2:
        alnid, aln_frac, metric = None, None, None
        # Prefer the 3-step HD-core protocol when ref+queries available — matches
        # step 3 / step 5 / step 5.5 metrics so the summary's id_pct is the same
        # value the picker used.
        if queries_dir and locus_ref_fa and os.path.exists(locus_ref_fa):
            try:
                from .graph_path_search import _align_cores
                proteins = os.path.join(queries_dir, "variable_proteins.fasta")
                named = [(n, alleles[n]) for n in names[:2]]
                aligned, hd_cols = _align_cores(named, locus_ref_fa, proteins)
                if hd_cols and names[0] in aligned and names[1] in aligned:
                    sa, sb = aligned[names[0]], aligned[names[1]]
                    hd_len_a = sum(1 for ci in hd_cols if ci < len(sa) and sa[ci] != '-')
                    hd_len_b = sum(1 for ci in hd_cols if ci < len(sb) and sb[ci] != '-')
                    denom = min(hd_len_a, hd_len_b) if (hd_len_a and hd_len_b) else 0
                    alnid, aln_frac = _score_columns(sa, sb, hd_cols, denom_len=denom)
                    metric = "HD-core"
            except Exception as e:
                print(f"[pairwise_identity] HD-core protocol failed ({e}); falling back to mafft_pair")
        if alnid is None:
            alnid, aln_frac = mafft_pair(alleles[names[0]], alleles[names[1]])
            metric = "whole-allele"
        res["id_pct"] = round(100 * alnid, 1)
        res["aln_frac"] = round(aln_frac, 2)
        res["distinct"] = not (alnid >= 0.95 and aln_frac >= 0.80)
        res["metric"] = metric
    if out_tsv:
        with open(out_tsv, "w") as o:
            o.write("a1\ta2\tlen_a1\tlen_a2\tid_pct\taln_frac\tdistinct\tmetric\n")
            o.write("\t".join(str(res[k]) if res[k] is not None else "-"
                              for k in ("a1", "a2", "len_a1", "len_a2", "id_pct", "aln_frac", "distinct", "metric")) + "\n")
    return res

def _cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--primary-alleles", required=True)
    p.add_argument("--out-tsv", default=None)
    p.add_argument("--queries-dir", default=None,
                   help="when supplied with --locus-ref, computes identity on "
                        "HD-core columns (matches step 3 / step 5 / step 5.5 metric).")
    p.add_argument("--locus-ref", default=None)
    a = p.parse_args(argv)
    r = run(a.primary_alleles, a.out_tsv,
            queries_dir=a.queries_dir, locus_ref_fa=a.locus_ref)
    for k, v in r.items(): print(f"{k}\t{v}")

if __name__ == "__main__":
    _cli()
