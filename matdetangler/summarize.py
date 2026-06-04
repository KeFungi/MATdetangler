"""Build a comprehensive cross-run-comparable JSON of a MATdetangler results dir.

Schema (single root object):
{
  "args": {
    "spades_dir":      <str | "examples/Pcub40">,
    "ks":              [...],
    "expected_count":  1 | 2,
    "min_allele_bp":   0,
    "cov_filter":      "on" | "off",
    "max_paths":       1000,
    "max_path_length": 50,
    "max_bp_since_var": 5000,
    "init_nhop":       3,
    "max_nhop":        8,
    "divergence_threshold": 0.01,
    "lo_mult":         0.2,
    "hi_mult":         2.0,
    "locus_padding":   4000,
    "seeds":           "both",
    "make_consensus":  bool,
    ...                                   # captured from the input args.json
  },
  "version": { "git_hash": "<sha>", "describe": "<git describe>" },
  "samples": {
    "<sample>": {
      # sample-level (from picks_summary.tsv)
      "k_chosen":        "k45" | "k53" | "-",
      "bubble_type":     "closed_bubble" | ...,
      "n_dedup":         int,
      "complete_var":    0|1|2,
      "complete_locus":  0|1|2,
      "locus_coverage":  float,
      "basepair_total":  int,
      "genome_cov":      float,
      "allele_cov":      [float, ...],
      "extend_bounds":   "<L_a:R_a;L_b:R_b>",
      "components":      "allele1,allele2,...",
      "all_k_tried":     "k45,k53,...",

      # per-allele (from picks.tsv)
      "alleles": [
        { "name":..., "origin":..., "k":..., "type":..., "len":..., "from_contig":...,
          "segments": [str], "cov":..., "n_variable_genes": int, "has_both_flanks": bool,
          "is_degHD": bool }, ...
      ],

      # pairwise (from summary.tsv)
      "allele1_complete": bool, "allele2_complete": bool,
      "allele1_vs_allele2_id_pct": float, "allele1_vs_allele2_aln_frac": float,
      "allele1_path_str": str, "allele2_path_str": str,

      # per-iteration trace (from logs/03_run_per_k_k<k>.log)
      "finished_nhop": int | null,                # nhop at which accept fired (None if no accept)
      "per_k_trace": {
        "k45": [ {"nhop": int, "cov_pass": "on"|"off", "n_arms": int, "cls": str, "nhood": int, ...}, ... ],
        "k53": ...
      }
    }
  }
}
"""
from __future__ import annotations
import argparse, json, os, re, subprocess, sys


def parse_picks_summary(path):
    if not os.path.exists(path): return None
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        row = fh.readline().rstrip("\n").split("\t")
        if not row or not row[0]: return None
        return dict(zip(hdr, row))


def parse_summary_tsv(path):
    if not os.path.exists(path): return None
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        row = fh.readline().rstrip("\n").split("\t")
        if not row or not row[0]: return None
        return dict(zip(hdr, row))


def parse_picks(path):
    """Returns list of per-allele dicts from picks.tsv."""
    if not os.path.exists(path): return []
    out = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            cells = ln.rstrip("\n").split("\t")
            if len(cells) < 2 or not cells[0]: continue
            d = dict(zip(hdr, cells))
            # Normalize types where possible
            for k in ("len",):
                try: d[k] = int(d.get(k, "0"))
                except (ValueError, TypeError): pass
            for k in ("cov",):
                try: d[k] = float(d.get(k, "0"))
                except (ValueError, TypeError): pass
            segs = d.get("segments", "")
            d["segments"] = [s for s in segs.split(",") if s] if segs and segs != "-" else []
            for k in ("has_both_flanks", "is_degHD"):
                v = d.get(k, "").strip().lower()
                d[k] = (v == "true") if v in ("true", "false") else d.get(k)
            out.append(d)
    return out


_TRACE_RX = re.compile(
    r"\[nhop=(\d+)\s+seeds=(\w+)\s+cov=(on|off)\]\s+"
    r"\|nhood\|=(\d+)\s+var=(\d+)\s+cls=(\S+)\s+arms=(\d+)"
)
_ACCEPT_RX = re.compile(r"\[accept\].*?coff_h(\d+)|\[accept\].*?con_h(\d+)")
_PICK_RX = re.compile(r"\[pick\] best=(\S+)\s+net=(\S+)\s+verdict=(\S+)\s+n=(\d+)")


def parse_per_k_log(path):
    """Read a logs/03_run_per_k_k<k>.log, return {trace: [...], finished_nhop, pick}."""
    if not os.path.exists(path): return {"trace": [], "finished_nhop": None, "pick": None}
    trace = []
    finished_nhop = None
    pick = None
    with open(path) as fh:
        for ln in fh:
            m = _TRACE_RX.search(ln)
            if m:
                nhop, seeds, cov, nhood, n_var, cls, arms = m.groups()
                trace.append({
                    "nhop": int(nhop), "seeds": seeds, "cov_pass": cov,
                    "nhood": int(nhood), "n_var": int(n_var),
                    "cls": cls, "n_arms": int(arms),
                })
                continue
            m = _ACCEPT_RX.search(ln)
            if m:
                finished_nhop = int(m.group(1) or m.group(2))
                continue
            m = _PICK_RX.search(ln)
            if m:
                pick = {"iter_id": m.group(1), "net": m.group(2),
                        "verdict": m.group(3), "n_dedup": int(m.group(4))}
    return {"trace": trace, "finished_nhop": finished_nhop, "pick": pick}


def git_version(repo_dir):
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       cwd=repo_dir, stderr=subprocess.DEVNULL).strip().decode()
    except Exception:
        sha = None
    try:
        desc = subprocess.check_output(["git", "describe", "--always", "--tags", "--dirty"],
                                        cwd=repo_dir, stderr=subprocess.DEVNULL).strip().decode()
    except Exception:
        desc = None
    return {"git_hash": sha, "describe": desc}


def summarize_sample(results_dir, sample):
    sd = os.path.join(results_dir, sample)
    out = {"sample": sample}
    ps = parse_picks_summary(os.path.join(sd, "picks_summary.tsv"))
    if ps:
        # Sample-level fields, with type normalization.
        cast_int = ("n_dedup", "complete_var", "complete_locus", "basepair", "n_cand")
        cast_flt = ("locus_coverage", "genome_cov")
        for k, v in ps.items():
            if k == "sample": continue
            if k in cast_int:
                try: out[k] = int(v)
                except (ValueError, TypeError): out[k] = v
            elif k in cast_flt:
                try: out[k] = float(v)
                except (ValueError, TypeError): out[k] = v
            elif k == "allele_cov":
                out[k] = [float(x) for x in v.split(",") if x and x != "-"]
            else:
                out[k] = v
    s = parse_summary_tsv(os.path.join(sd, "summary.tsv"))
    if s:
        for k in ("allele1_complete", "allele2_complete"):
            v = s.get(k, "").strip().upper()
            out[k] = (v == "TRUE") if v in ("TRUE", "FALSE") else v
        for k in ("allele1_vs_allele2_id_pct", "allele1_vs_allele2_aln_frac"):
            try: out[k] = float(s.get(k, "0"))
            except (ValueError, TypeError): out[k] = s.get(k)
        out["allele1_path_str"] = s.get("allele1_path", "")
        out["allele2_path_str"] = s.get("allele2_path", "")
    out["alleles"] = parse_picks(os.path.join(sd, "picks.tsv"))
    out["per_k_trace"] = {}
    out["finished_nhop"] = None
    out["per_k_pick"] = {}
    for log in sorted(os.listdir(os.path.join(sd, "logs")) if os.path.isdir(os.path.join(sd, "logs")) else []):
        m = re.match(r"03_run_per_k_(k\d+)\.log", log)
        if not m: continue
        k = m.group(1)
        rec = parse_per_k_log(os.path.join(sd, "logs", log))
        out["per_k_trace"][k] = rec["trace"]
        out["per_k_pick"][k] = rec["pick"]
        if rec["finished_nhop"] is not None:
            # Set sample-level finished_nhop from the picked k (if multiple k's accepted,
            # the picked k's value wins).
            if k == out.get("k_chosen"):
                out["finished_nhop"] = rec["finished_nhop"]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--results-dir",
                      help="results dir to walk (one subdir per sample) — emits cross-sample aggregate")
    src.add_argument("--sample-dir",
                      help="a single sample's results dir — emits a one-sample summary")
    ap.add_argument("--args-json", default=None,
                     help="optional JSON with run-time args/config to embed under 'args'")
    ap.add_argument("--repo-dir", default=None,
                     help="repo root for git hash; default: parent of input dir")
    ap.add_argument("--samples", default=None,
                     help="(--results-dir only) comma-separated subset of samples to include")
    ap.add_argument("--out", required=True, help="output JSON path")
    a = ap.parse_args(argv)

    args_block = {}
    if a.args_json and os.path.exists(a.args_json):
        with open(a.args_json) as fh:
            args_block = json.load(fh)

    if a.sample_dir:
        # Single-sample mode — wrapper invokes this per sample as the final step.
        sd_abs = os.path.abspath(a.sample_dir).rstrip("/")
        sample = os.path.basename(sd_abs)
        results_dir = os.path.dirname(sd_abs)
        repo_dir = a.repo_dir or os.path.dirname(results_dir)
        out = {
            "args": args_block,
            "version": git_version(repo_dir),
            "sample": summarize_sample(results_dir, sample),
        }
    else:
        repo_dir = a.repo_dir or os.path.dirname(os.path.abspath(a.results_dir))
        samples = sorted(
            d for d in os.listdir(a.results_dir)
            if os.path.isdir(os.path.join(a.results_dir, d)) and not d.startswith("_")
        )
        if a.samples:
            keep = {x for x in a.samples.split(",") if x}
            samples = [s for s in samples if s in keep]
        out = {
            "args": args_block,
            "version": git_version(repo_dir),
            "samples": {s: summarize_sample(a.results_dir, s) for s in samples},
        }
    out_dir = os.path.dirname(a.out)
    if out_dir: os.makedirs(out_dir, exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True, default=str)
    n = 1 if a.sample_dir else len(out["samples"])
    print(f"wrote {a.out} ({n} sample{'' if n == 1 else 's'})")


if __name__ == "__main__":
    main()
