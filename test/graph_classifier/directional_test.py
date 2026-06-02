"""Synthetic tests for the Appendix-B (directional-split) classifier.

Each test case provides:
  - seg_labels: dict[seg_id, list[Hit]]                     (per-segment hits)
  - seg_length: dict[seg_id, int]                            (segment lengths)
  - edges:      set[frozenset({seg_a, seg_b})]              (undirected)
  - edge_endpoints (optional): dict[edge, (side_a, side_b)] (L/R per side)
  - expected:   ground-truth class

The Hit class records (tag, kind, start, end, strand) with 0-based segment
coordinates. Tests exercise:
  - composite splitting (P1)
  - long unlabeled spacers absorbed into bubble (P2)
  - Y-fork bubbles (AG10 shape)
  - single-composite bridges (AJB36 shape)
  - shared-spine multi-arm bubbles (AG5-like)
  - true open/complex/separate topologies
"""
from __future__ import annotations
from test.graph_classifier.seg_processor import Hit, directional_split
from test.graph_classifier.bubble_classifier import classify


def _net(name, expected, seg_labels, seg_length, edges, edge_endpoints=None):
    return {
        "name": name, "expected": expected,
        "seg_labels": seg_labels, "seg_length": seg_length,
        "edges": edges, "edge_endpoints": edge_endpoints,
    }


def _h_flankL(s, e):  return Hit("flankL", "flank", s, e)
def _h_flankR(s, e):  return Hit("flankR", "flank", s, e)
def _h_HD1(s, e):     return Hit("HD1", "var", s, e)
def _h_HD2(s, e):     return Hit("HD2", "var", s, e)
def _h_HD(s, e, g):   return Hit(g, "var", s, e)


# ---- TEST CASES ----

def case_synthetic_closed_bubble():
    """flankL — HD1a — HD2a — flankR ; flankL — HD1b — HD2b — flankR"""
    seg_labels = {
        "flankL": [_h_flankL(0, 5000)],
        "flankR": [_h_flankR(0, 5000)],
        "HD1a":   [_h_HD1(0, 1000)],
        "HD2a":   [_h_HD2(0, 1000)],
        "HD1b":   [_h_HD1(0, 1000)],
        "HD2b":   [_h_HD2(0, 1000)],
    }
    seg_length = {k: 5000 for k in seg_labels}
    edges = {frozenset(("flankL", "HD1a")), frozenset(("HD1a", "HD2a")),
             frozenset(("HD2a", "flankR")),
             frozenset(("flankL", "HD1b")), frozenset(("HD1b", "HD2b")),
             frozenset(("HD2b", "flankR"))}
    return _net("synthetic_closed_bubble", "closed_bubble",
                seg_labels, seg_length, edges)


def case_ajb36_single_composite():
    """One composite segment carrying HD1+HD2+flankL+flankR. P1 should split
    into [flankL]-[HD1+HD2]-[flankR]; bubble = single HD; verdict = single."""
    seg_labels = {
        "L":   [_h_flankL(0, 5000)],
        "R":   [_h_flankR(0, 5000)],
        "all": [_h_flankL(0, 500), _h_HD1(700, 1500), _h_HD2(1700, 2300),
                _h_flankR(2500, 3000)],
    }
    seg_length = {"L": 5000, "R": 5000, "all": 3000}
    edges = {frozenset(("L", "all")), frozenset(("all", "R"))}
    # Provide explicit edge endpoints (canonical sorted-tuple keys)
    edge_endpoints = {
        tuple(sorted(("L", "all"))): ("R", "L"),    # L's right end to all's left
        tuple(sorted(("all", "R"))): ("R", "L"),    # all's right end to R's left
    }
    return _net("ajb36_single_composite", "single",
                seg_labels, seg_length, edges, edge_endpoints)


def case_ag17_composite_flankR():
    """Y-fork interior with composite at flankR end.
    flankL — Jl — HD_arm1 — Jr — flankR_comp(HD2+flankR)
                 \\ HD_arm2a — HD_arm2b /
    """
    seg_labels = {
        "flankL":      [_h_flankL(0, 5000)],
        "Jl":          [],
        "HD_arm1":     [_h_HD1(0, 1000)],
        "HD_arm2a":    [_h_HD1(0, 1000)],
        "HD_arm2b":    [_h_HD2(0, 1000)],
        "Jr":          [],
        "flankR_comp": [_h_HD2(0, 800), _h_flankR(900, 1500)],
    }
    seg_length = {k: max(1000, max((h.end for h in v), default=1000))
                  for k, v in seg_labels.items()}
    seg_length["flankL"] = 5000
    edges = {
        frozenset(("flankL", "Jl")),
        frozenset(("Jl", "HD_arm1")), frozenset(("Jl", "HD_arm2a")),
        frozenset(("HD_arm1", "Jr")),
        frozenset(("HD_arm2a", "HD_arm2b")), frozenset(("HD_arm2b", "Jr")),
        frozenset(("Jr", "flankR_comp")),
    }
    edge_endpoints = {
        tuple(sorted(("Jr", "flankR_comp"))): ("L", "R"),  # composite L (HD2) ↔ Jr R
    }
    return _net("ag17_composite_flankR", "closed_bubble",
                seg_labels, seg_length, edges, edge_endpoints)


def case_long_unlabeled_spacer():
    """SA93 shape: 5 unlabeled connectors between var and flankR."""
    seg_labels = {
        "flankL":   [_h_flankL(0, 5000)],
        "flankR":   [_h_flankR(0, 5000)],
        "HD1a":     [_h_HD1(0, 1000)],
        "HD2a":     [_h_HD2(0, 1000)],
        "HD1b":     [_h_HD1(0, 1000)],
        "HD2b":     [_h_HD2(0, 1000)],
        "link_0": [], "link_1": [], "link_2": [], "link_3": [], "link_4": [],
    }
    seg_length = {k: 1000 for k in seg_labels}
    edges = {
        frozenset(("flankL", "HD1a")), frozenset(("HD1a", "HD2a")),
        frozenset(("flankL", "HD1b")), frozenset(("HD1b", "HD2b")),
        frozenset(("HD2a", "link_0")), frozenset(("HD2b", "link_0")),
        frozenset(("link_0", "link_1")), frozenset(("link_1", "link_2")),
        frozenset(("link_2", "link_3")), frozenset(("link_3", "link_4")),
        frozenset(("link_4", "flankR")),
    }
    return _net("long_unlabeled_spacer", "closed_bubble",
                seg_labels, seg_length, edges)


def case_open_bubble_dangling():
    """One arm closes at flankR, the other dangles."""
    seg_labels = {
        "flankL": [_h_flankL(0, 5000)],
        "flankR": [_h_flankR(0, 5000)],
        "HD1a":   [_h_HD1(0, 1000)], "HD2a": [_h_HD2(0, 1000)],
        "HD1b":   [_h_HD1(0, 1000)], "HD2b": [_h_HD2(0, 1000)],
    }
    seg_length = {k: 1000 for k in seg_labels}
    edges = {frozenset(("flankL", "HD1a")), frozenset(("HD1a", "HD2a")),
             frozenset(("HD2a", "flankR")),
             frozenset(("flankL", "HD1b")), frozenset(("HD1b", "HD2b"))}
    return _net("open_bubble_dangling", "open_bubble",
                seg_labels, seg_length, edges)


def case_complex_3arm():
    """3 parallel arms — complexed."""
    seg_labels = {
        "flankL": [_h_flankL(0, 5000)],
        "flankR": [_h_flankR(0, 5000)],
    }
    edges = set()
    for s in "abc":
        seg_labels[f"HD1{s}"] = [_h_HD1(0, 1000)]
        seg_labels[f"HD2{s}"] = [_h_HD2(0, 1000)]
        edges.add(frozenset(("flankL", f"HD1{s}")))
        edges.add(frozenset((f"HD1{s}", f"HD2{s}")))
        edges.add(frozenset((f"HD2{s}", "flankR")))
    seg_length = {k: 1000 for k in seg_labels}
    return _net("complex_3arm", "complexed", seg_labels, seg_length, edges)


def case_separate():
    """Two completely disjoint flank-HD-flank chains."""
    seg_labels, edges = {}, set()
    for tag in "ab":
        seg_labels[f"flankL_{tag}"] = [_h_flankL(0, 5000)]
        seg_labels[f"flankR_{tag}"] = [_h_flankR(0, 5000)]
        seg_labels[f"HD1_{tag}"] = [_h_HD1(0, 1000)]
        seg_labels[f"HD2_{tag}"] = [_h_HD2(0, 1000)]
        edges.add(frozenset((f"flankL_{tag}", f"HD1_{tag}")))
        edges.add(frozenset((f"HD1_{tag}", f"HD2_{tag}")))
        edges.add(frozenset((f"HD2_{tag}", f"flankR_{tag}")))
    seg_length = {k: 1000 for k in seg_labels}
    return _net("separate", "separate", seg_labels, seg_length, edges)


def case_ag10_composite_L_chain():
    """AG10 shape: two composites at L side (HD+flankL each) sharing a flankL
    chain; unlabeled joints at both ends of the bubble; composite at R."""
    seg_labels = {
        "flankL_chain":   [_h_flankL(0, 5000)],
        "Jl":             [],
        "HD_flankL_a":    [_h_flankL(0, 400), _h_HD1(500, 1500)],
        "HD_flankL_b":    [_h_flankL(0, 400), _h_HD1(500, 1500)],
        "Jr":             [],
        "HD_flankR":      [_h_HD2(0, 800), _h_flankR(900, 1500)],
    }
    seg_length = {"flankL_chain": 5000, "Jl": 1000, "HD_flankL_a": 1500,
                  "HD_flankL_b": 1500, "Jr": 1000, "HD_flankR": 1500}
    edges = {
        frozenset(("flankL_chain", "Jl")),
        frozenset(("Jl", "HD_flankL_a")), frozenset(("Jl", "HD_flankL_b")),
        frozenset(("HD_flankL_a", "Jr")), frozenset(("HD_flankL_b", "Jr")),
        frozenset(("Jr", "HD_flankR")),
    }
    # Canonical-tuple keys (sorted); endpoints are (side_of_smaller, side_of_larger)
    def _ek(a, b, sa, sb):
        if a < b: return (a, b), (sa, sb)
        return (b, a), (sb, sa)
    edge_endpoints = dict([
        _ek("flankL_chain", "Jl",         "R", "L"),
        _ek("Jl",           "HD_flankL_a", "R", "L"),   # composite L (flankL) ↔ Jl R
        _ek("Jl",           "HD_flankL_b", "R", "L"),
        _ek("HD_flankL_a",  "Jr",          "R", "L"),   # composite R (HD1) ↔ Jr L
        _ek("HD_flankL_b",  "Jr",          "R", "L"),
        _ek("Jr",           "HD_flankR",   "R", "L"),   # Jr R ↔ composite L (HD2)
    ])
    return _net("ag10_composite_L_chain", "closed_bubble",
                seg_labels, seg_length, edges, edge_endpoints)


ALL_CASES = [
    case_synthetic_closed_bubble,
    case_ajb36_single_composite,
    case_ag17_composite_flankR,
    case_long_unlabeled_spacer,
    case_open_bubble_dangling,
    case_complex_3arm,
    case_separate,
    case_ag10_composite_L_chain,
]


def run_all():
    n_pass = n_fail = 0
    for builder in ALL_CASES:
        c = builder()
        nodes, edges, labels, var_per, _provenance = directional_split(
            c["seg_labels"], c["seg_length"], c["edges"], c["edge_endpoints"]
        )
        res = classify(nodes, edges, labels, var_per)
        ok = res["class"] == c["expected"]
        mark = "PASS" if ok else "FAIL"
        explain = res.get("explain", "")
        print(f"  [{mark}] {c['name']:<35} expected={c['expected']:<14} "
              f"got={res['class']:<14} {explain}")
        if ok: n_pass += 1
        else:
            print(f"          nodes={sorted(nodes)}")
            print(f"          labels={labels}")
            n_fail += 1
    print(f"  {n_pass} pass, {n_fail} fail")
    return n_pass, n_fail


if __name__ == "__main__":
    run_all()
