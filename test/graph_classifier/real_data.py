#!/usr/bin/env python3
"""Adapt MATdetangler per-sample outputs into the classifier's input format
and run the new walk-aware classifier.

Per sample, build a TestNetwork from:
  results/<sample>/bubble.tsv          edges + node labels (picked allele paths)
The bubble.tsv is the visualization-focused subgraph holding the two (or one)
picked allele arms with their per-segment labels. Composite labels (e.g.
"flankL+HD1") are preserved on the node side. var_per is derived from any
non-"flank*" tag.

Output:
  classifier_real.tsv  — columns:
    sample  expected_class_from_picks  n_picks  classifier_verdict  n_var_series  notes
"""
from __future__ import annotations
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from test.graph_classifier.network import TestNetwork
from test.graph_classifier.classifier import classify


def _split_label_tokens(label: str) -> list[str]:
    return [t for t in (label or "").split("+") if t]


def _is_flank_token(t: str) -> bool:
    return t.lower().startswith("flank")


def _read_bubble_tsv(path: str) -> tuple[TestNetwork | None, dict[str, set[str]]]:
    """Returns (network, arm_per_node) where arm_per_node[node] is the set of
    arm labels (from bubble.tsv 'arm' column) the node belongs to."""
    if not os.path.isfile(path): return None, {}
    net = TestNetwork(name=os.path.basename(os.path.dirname(path)), expected="unknown")
    arm_per_node: dict[str, set[str]] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            line = line.rstrip("\n")
            if not line: continue
            cells = line.split("\t")
            if len(cells) < 4: continue
            a = cells[idx["node_a"]]
            b = cells[idx["node_b"]]
            la = cells[idx["label_a"]] if "label_a" in idx else ""
            lb = cells[idx["label_b"]] if "label_b" in idx else ""
            arm = cells[idx["arm"]] if "arm" in idx else ""
            for node, lab in ((a, la), (b, lb)):
                if node not in net.nodes:
                    toks = _split_label_tokens(lab)
                    var_tokens = {t for t in toks if not _is_flank_token(t)}
                    net.add_node(node, label=lab,
                                 vars=var_tokens if var_tokens else None)
                if arm:
                    arm_per_node.setdefault(node, set()).add(arm)
            net.add_edge(a, b)
    return (net if net.nodes else None), arm_per_node


def _read_segment_depths(gfa_path: str) -> dict[str, float]:
    """Parse S-lines from a GFA and return seg_id -> DP:f: value."""
    out: dict[str, float] = {}
    if not os.path.isfile(gfa_path): return out
    with open(gfa_path) as fh:
        for line in fh:
            if not line.startswith("S\t"): continue
            cells = line.rstrip("\n").split("\t")
            if len(cells) < 3: continue
            sid = cells[1]
            for tag in cells[3:]:
                if tag.startswith("DP:f:"):
                    try: out[sid] = float(tag[5:])
                    except ValueError: pass
                    break
    return out


def _filter_by_depth(net: TestNetwork, depths: dict[str, float],
                      lo: float, hi: float) -> tuple[TestNetwork, int]:
    """Return a deep copy with nodes whose depth is outside [lo, hi] removed.
    Edges with at least one removed endpoint are dropped. Returns (filtered_net,
    n_dropped). Nodes without a known depth (not in `depths`) are KEPT — we
    don't want to throw away nodes just because we couldn't look them up."""
    keep = {n for n in net.nodes
            if (n not in depths) or (lo <= depths[n] <= hi)}
    dropped = len(net.nodes) - len(keep)
    if dropped == 0:
        return net, 0
    out = net.deep_copy(new_name=net.name + f"+depthfilt[{lo:.1f}-{hi:.1f}]")
    out.nodes = keep
    out.edges = {e for e in out.edges if all(x in keep for x in tuple(e))}
    out.labels = {k: v for k, v in out.labels.items() if k in keep}
    out.var_per = {k: v for k, v in out.var_per.items() if k in keep}
    return out, dropped


def _expected_from_picks(picks_tsv: str) -> tuple[str, int]:
    """Return (expected_class, n_alleles_picked) derived from picks.tsv.

    Heuristic ground truth derived from how MATdetangler resolved the sample:
      - 0 alleles -> 'no_var' (sample failed)
      - 1 allele  -> 'single' (singleton pick / haploid-like)
      - 2 alleles, both complete, picks tagged 'bubble' or 'path' from one k:
                  -> 'closed_bubble' if both have has_both_flanks=True
                  -> 'open_bubble'   otherwise
      - 2 alleles, k differs                  -> 'complexed' (cross-k pair)
    """
    if not os.path.isfile(picks_tsv): return ("no_var", 0)
    rows = []
    with open(picks_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            cells = line.rstrip("\n").split("\t")
            if len(cells) < len(header): continue
            if not cells[idx.get("allele", 1)].startswith("allele"): continue
            rows.append(cells)
    n = len(rows)
    if n == 0: return ("no_var", 0)
    if n == 1: return ("single", 1)
    if n == 2:
        ks = {r[idx["k"]] for r in rows}
        both_flanks = all(r[idx["has_both_flanks"]] == "True" for r in rows)
        if len(ks) > 1: return ("complexed", 2)
        return ("closed_bubble" if both_flanks else "open_bubble", 2)
    return ("complexed", n)


def _load_genome_cov(path: str) -> dict[str, float]:
    out: dict[str, float] = {}
    if not os.path.isfile(path): return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            cells = line.rstrip("\n").split("\t")
            if len(cells) < len(header): continue
            try: out[cells[idx["sample"]]] = float(cells[idx["genome_cov"]])
            except (KeyError, ValueError): pass
    return out


def _classify_one(net: TestNetwork,
                   arm_per_node: dict[str, set[str]] | None = None) -> dict:
    return classify(set(net.nodes), set(net.edges), dict(net.labels),
                    {k: set(v) for k, v in net.var_per.items()},
                    arm_per_node=arm_per_node)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--sample-list", default="",
                    help="TSV with sample -> genome_cov (column names "
                         "'sample' and 'genome_cov'). Enables depth filter.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--depth-lo", type=float, default=0.25,
                    help="lower depth multiplier of genome_cov (default 0.25)")
    ap.add_argument("--depth-hi", type=float, default=2.0,
                    help="upper depth multiplier of genome_cov (default 2.0)")
    a = ap.parse_args()

    samples = sorted(os.listdir(a.results_dir))
    genome_cov = _load_genome_cov(a.sample_list) if a.sample_list else {}
    cols = ["sample", "expected_from_picks", "n_picks", "pass",
            "classifier_verdict", "n_nodes_in", "n_dropped",
            "n_arms_closed", "n_arms_dangling", "explain"]
    n_total = 0; n_match = 0; n_pass2 = 0
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for s in samples:
            sdir = os.path.join(a.results_dir, s)
            picks = os.path.join(sdir, "picks.tsv")
            bubble_tsv = os.path.join(sdir, "bubble.tsv")
            bubble_gfa = os.path.join(sdir, "bubble.gfa")
            exp, n_picks = _expected_from_picks(picks)
            net, arm_per_node = _read_bubble_tsv(bubble_tsv)
            if net is None:
                fh.write(f"{s}\t{exp}\t{n_picks}\t-\t-\t-\t-\t-\t-\tno-bubble-tsv\n")
                continue
            # Depth filter — first pass tries to find a closed bubble in the
            # "normal-coverage" subgraph, falling back to unfiltered if no
            # close_bubble verdict.
            gc = genome_cov.get(s)
            depths = _read_segment_depths(bubble_gfa) if gc else {}
            net_pass1, dropped = (net, 0)
            if gc and depths:
                net_pass1, dropped = _filter_by_depth(net, depths,
                                                       a.depth_lo * gc,
                                                       a.depth_hi * gc)
            res = _classify_one(net_pass1, arm_per_node)
            pass_used = 1
            # Second pass — restore filtered nodes if pass 1 didn't find a closed bubble
            if res["class"] != "closed_bubble" and dropped > 0:
                res2 = _classify_one(net, arm_per_node)
                if res2["class"] == "closed_bubble":
                    res = res2; pass_used = 2; n_pass2 += 1
                else:
                    # Keep pass 1's verdict (was strict; pass 2 didn't improve)
                    pass
            n_total += 1
            verdict = res["class"]
            if verdict == exp: n_match += 1
            fh.write("\t".join([
                s, exp, str(n_picks), str(pass_used), verdict,
                str(res.get("n_nodes", "-")), str(dropped),
                str(res.get("n_closed_arms", res.get("n_var_series", "-"))),
                str(res.get("n_dangling_arms", "-")),
                res.get("explain", "").replace("\t", " "),
            ]) + "\n")
    print(f"wrote {a.out}")
    print(f"{n_match}/{n_total} agreement; {n_pass2} samples rescued by 2nd-pass (filter restored)")


if __name__ == "__main__":
    sys.exit(main())
