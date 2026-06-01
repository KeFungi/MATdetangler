"""Assemble the per-sample summary row from the artifacts produced by the other modules.

Columns:
  sample, used_k, bubble_type, genome_coverage,
  allele1_coverage, allele2_coverage, allele1_vs_allele2_identity,
  allele1_path, allele2_path
"""
from __future__ import annotations
import os, argparse, json

def _read_coverage(coverage_tsv: str) -> dict[str, dict]:
    """Parse coverage.tsv. Current schema (coverage_core.py):
        #allele  len  mapped_reads  whole_breadth_pct  whole_meandepth  core_bp  core_meandepth  repeat_bp
    `meandepth` is the headline number — set to `core_meandepth` (HD-core, repeat-masked) when
    present, falling back to `whole_meandepth`. The HD-core depth is what a dikaryon balance
    check should use; whole-allele depth is inflated 2-3x for any allele carrying a high-copy
    repeat (e.g. MITE inside a stitched allele)."""
    out = {}
    if not os.path.exists(coverage_tsv): return out
    with open(coverage_tsv) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip(): continue
            f = ln.rstrip("\n").split("\t")
            def _f(i):
                if i >= len(f) or f[i] in ("", "-"): return 0.0
                try: return float(f[i])
                except ValueError: return 0.0
            whole = _f(4)
            has_core = len(f) >= 7 and f[6] not in ("", "-")
            core = _f(6) if has_core else whole
            out[f[0]] = {"len": int(f[1]) if f[1].isdigit() else 0,
                         "mapped_reads": f[2], "breadth_pct": f[3],
                         "whole_meandepth": whole, "core_meandepth": core,
                         "meandepth": core}
    return out

def _read_consensus_qc(qc_tsv: str | None) -> dict[str, dict]:
    """Parse consensus_qc.tsv (written by consensus_qc.run). Returns {allele_name: row_dict}."""
    out: dict[str, dict] = {}
    if not qc_tsv or not os.path.exists(qc_tsv): return out
    with open(qc_tsv) as fh:
        hdr = next(fh, "").rstrip("\n").split("\t")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < len(hdr): continue
            out[f[0]] = dict(zip(hdr, f))
    return out

def write_row(sample: str, used_k: str, bubble_type: str, genome_cov: float,
              coverage_tsv: str, picks_tsv: str, identity_dict: dict,
              graph_paths: dict, out_tsv: str, append: bool = True,
              consensus_identity_dict: dict | None = None,
              consensus_qc_tsv: str | None = None) -> None:
    """Write one summary row.

    Existing pick-level columns are preserved verbatim (`allele1_complete`, …,
    `allele1_vs_allele2_id_pct`, …). New consensus-level columns are appended:
        cons_allele1_complete, cons_allele2_complete,
        cons_allele1_vs_allele2_id_pct, cons_allele1_vs_allele2_aln_frac
    so existing parsers keep working.
    """
    cov = _read_coverage(coverage_tsv)
    picks_hdr = []
    picks = []
    if os.path.exists(picks_tsv):
        with open(picks_tsv) as fh:
            picks_hdr = next(fh, "").rstrip("\n").split("\t")
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) >= 2: picks.append(f)
    def pick_col(row, col_name, default="-"):
        try: i = picks_hdr.index(col_name)
        except ValueError: return default
        return row[i] if i < len(row) else default
    # pick-level completeness (TRUE iff type == "complete"; partial / element-like -> FALSE)
    def is_complete(p): return "TRUE" if pick_col(p, "type") == "complete" else "FALSE"
    a1_complete = is_complete(picks[0]) if len(picks) >= 1 else "-"
    a2_complete = is_complete(picks[1]) if len(picks) >= 2 else "-"
    def cov_for(allele_label: str) -> float:
        for k, v in cov.items():
            if k.endswith("_" + allele_label) or k == allele_label or k.endswith(allele_label):
                return v["meandepth"]
        return 0.0
    a1_cov = cov_for("allele1")
    a2_cov = cov_for("allele2")
    idp = identity_dict.get("id_pct", "-")
    aln = identity_dict.get("aln_frac", "-")
    a1_path = graph_paths.get("arm1_str", "-")
    a2_path = graph_paths.get("arm2_str", "-")
    # consensus-level: re-checked completeness + identity on the read-derived consensus
    ci = consensus_identity_dict or {}
    cons_idp = ci.get("id_pct", "-")
    cons_aln = ci.get("aln_frac", "-")
    cqc = _read_consensus_qc(consensus_qc_tsv)
    def cons_complete(label: str) -> str:
        for k, row in cqc.items():
            if k.endswith("_" + label) or k == label or k.endswith(label):
                return row.get("complete", "-")
        return "-"
    cons_a1_complete = cons_complete("allele1")
    cons_a2_complete = cons_complete("allele2")
    cols = ["sample", "used_k", "bubble_type", "genome_coverage",
            "allele1_complete", "allele2_complete",
            "allele1_coverage", "allele2_coverage",
            "allele1_vs_allele2_id_pct", "allele1_vs_allele2_aln_frac",
            "cons_allele1_complete", "cons_allele2_complete",
            "cons_allele1_vs_allele2_id_pct", "cons_allele1_vs_allele2_aln_frac",
            "allele1_path", "allele2_path"]
    row = [sample, used_k, bubble_type, f"{genome_cov:.1f}",
           a1_complete, a2_complete,
           f"{a1_cov:.1f}", f"{a2_cov:.1f}",
           str(idp), str(aln),
           cons_a1_complete, cons_a2_complete,
           str(cons_idp), str(cons_aln),
           a1_path or "-", a2_path or "-"]
    write_header = (not append) or (not os.path.exists(out_tsv))
    with open(out_tsv, "a" if append else "w") as o:
        if write_header: o.write("\t".join(cols) + "\n")
        o.write("\t".join(row) + "\n")
    print(f"[summary] {sample} -> {out_tsv}")

def _cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--sample", required=True)
    p.add_argument("--used-k", required=True)
    p.add_argument("--bubble-type", required=True)
    p.add_argument("--genome-coverage", type=float, default=0.0)
    p.add_argument("--coverage-tsv", required=True)
    p.add_argument("--picks-tsv", required=True)
    p.add_argument("--identity-json", required=True, help="JSON with id_pct/aln_frac on PICKS (pairwise_identity output)")
    p.add_argument("--graph-paths-json", required=True, help="JSON with arm1_str/arm2_str (graph_paths output)")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--no-append", action="store_true")
    p.add_argument("--consensus-identity-json", default=None,
                   help="optional JSON with id_pct/aln_frac on the CONSENSUS pair (pairwise_identity output, second pass)")
    p.add_argument("--consensus-qc-tsv", default=None,
                   help="optional TSV from consensus_qc.py — adds cons_allele{1,2}_complete columns")
    a = p.parse_args(argv)
    idn = json.load(open(a.identity_json))
    gp  = json.load(open(a.graph_paths_json))
    cons_idn = json.load(open(a.consensus_identity_json)) if a.consensus_identity_json and os.path.exists(a.consensus_identity_json) else None
    write_row(a.sample, a.used_k, a.bubble_type, a.genome_coverage,
              a.coverage_tsv, a.picks_tsv, idn, gp, a.out_tsv, append=not a.no_append,
              consensus_identity_dict=cons_idn, consensus_qc_tsv=a.consensus_qc_tsv)

if __name__ == "__main__":
    _cli()
