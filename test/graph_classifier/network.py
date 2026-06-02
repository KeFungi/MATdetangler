"""TestNetwork: a minimal GFA-like graph + builders for the four topology classes.

The data shape mirrors what the BFS/path-enum stages of MATdetangler consume:
  - nodes:     set of seg-id strings
  - edges:     set of frozenset({a, b}) — undirected (GFA L-lines reduced to undirected)
  - labels:    dict node -> "+"-joined tag string (e.g. "flankL", "HD1+HD2+flankR")
  - var_per:   dict node -> set of variable-gene names ({"HD1"}, {"HD1","HD2"}, ...)

The classifier under test consumes (nodes, edges, labels, var_per).
"""
from __future__ import annotations
from dataclasses import dataclass, field
from collections import defaultdict
import copy


@dataclass
class TestNetwork:
    name: str
    expected: str                                   # ground-truth class
    nodes: set[str] = field(default_factory=set)
    edges: set[frozenset] = field(default_factory=set)
    labels: dict[str, str] = field(default_factory=dict)
    var_per: dict[str, set[str]] = field(default_factory=dict)

    def add_node(self, n: str, label: str = "", vars: set[str] | None = None) -> None:
        self.nodes.add(n)
        if label:
            self.labels[n] = label
        if vars:
            self.var_per[n] = set(vars)

    def add_edge(self, a: str, b: str) -> None:
        if a not in self.nodes: self.add_node(a)
        if b not in self.nodes: self.add_node(b)
        if a != b: self.edges.add(frozenset((a, b)))

    def adj(self) -> dict[str, set[str]]:
        a: dict[str, set[str]] = defaultdict(set)
        for e in self.edges:
            x, y = tuple(e)
            a[x].add(y); a[y].add(x)
        return dict(a)

    def deep_copy(self, new_name: str | None = None, new_expected: str | None = None) -> "TestNetwork":
        return TestNetwork(
            name=new_name or self.name,
            expected=new_expected or self.expected,
            nodes=copy.copy(self.nodes),
            edges=copy.copy(self.edges),
            labels=copy.copy(self.labels),
            var_per={k: set(v) for k, v in self.var_per.items()},
        )

    def summary(self) -> str:
        adj = self.adj()
        return (f"{self.name}  expected={self.expected}  "
                f"n_nodes={len(self.nodes)} n_edges={len(self.edges)} "
                f"var_nodes={sum(1 for n in self.var_per)} "
                f"flankL={sum('flankL' in self.labels.get(n, '') for n in self.nodes)} "
                f"flankR={sum('flankR' in self.labels.get(n, '') for n in self.nodes)}")


# ---- clean case builders ----

def closed_bubble() -> TestNetwork:
    """flankL — HD1a — HD2a — flankR
             \\ HD1b — HD2b /"""
    net = TestNetwork(name="closed_bubble", expected="closed_bubble")
    net.add_node("flankL", label="flankL")
    net.add_node("flankR", label="flankR")
    for s in "ab":
        net.add_node(f"HD1{s}", label="HD1", vars={"HD1"})
        net.add_node(f"HD2{s}", label="HD2", vars={"HD2"})
        net.add_edge("flankL", f"HD1{s}")
        net.add_edge(f"HD1{s}", f"HD2{s}")
        net.add_edge(f"HD2{s}", "flankR")
    return net


def open_bubble_case1() -> TestNetwork:
    """flankL — HD1a — HD2a — flankR
             \\ HD1b — HD2b (no flankR)"""
    net = TestNetwork(name="open_bubble_case1", expected="open_bubble")
    net.add_node("flankL", label="flankL")
    net.add_node("flankR", label="flankR")
    net.add_node("HD1a", label="HD1", vars={"HD1"})
    net.add_node("HD2a", label="HD2", vars={"HD2"})
    net.add_node("HD1b", label="HD1", vars={"HD1"})
    net.add_node("HD2b", label="HD2", vars={"HD2"})
    for e in [("flankL", "HD1a"), ("HD1a", "HD2a"), ("HD2a", "flankR"),
              ("flankL", "HD1b"), ("HD1b", "HD2b")]:
        net.add_edge(*e)
    return net


def open_bubble_case2() -> TestNetwork:
    """flankL — HD1a — HD2a (no flankR)
             \\ HD1b — HD2b"""
    net = TestNetwork(name="open_bubble_case2", expected="open_bubble")
    net.add_node("flankL", label="flankL")
    net.add_node("HD1a", label="HD1", vars={"HD1"})
    net.add_node("HD2a", label="HD2", vars={"HD2"})
    net.add_node("HD1b", label="HD1", vars={"HD1"})
    net.add_node("HD2b", label="HD2", vars={"HD2"})
    for e in [("flankL", "HD1a"), ("HD1a", "HD2a"),
              ("flankL", "HD1b"), ("HD1b", "HD2b")]:
        net.add_edge(*e)
    return net


def complex_case1() -> TestNetwork:
    """3 parallel arms, all share flankL and flankR."""
    net = TestNetwork(name="complex_case1", expected="complexed")
    net.add_node("flankL", label="flankL")
    net.add_node("flankR", label="flankR")
    for s in "abc":
        net.add_node(f"HD1{s}", label="HD1", vars={"HD1"})
        net.add_node(f"HD2{s}", label="HD2", vars={"HD2"})
        net.add_edge("flankL", f"HD1{s}")
        net.add_edge(f"HD1{s}", f"HD2{s}")
        net.add_edge(f"HD2{s}", "flankR")
    return net


def complex_case2() -> TestNetwork:
    """3 arms; one closed (a), two open from flankL (b, c). No second flankR connection."""
    net = TestNetwork(name="complex_case2", expected="complexed")
    net.add_node("flankL", label="flankL")
    net.add_node("flankR", label="flankR")
    net.add_node("HD1a", label="HD1", vars={"HD1"})
    net.add_node("HD2a", label="HD2", vars={"HD2"})
    for e in [("flankL", "HD1a"), ("HD1a", "HD2a"), ("HD2a", "flankR")]:
        net.add_edge(*e)
    for s in "bc":
        net.add_node(f"HD1{s}", label="HD1", vars={"HD1"})
        net.add_node(f"HD2{s}", label="HD2", vars={"HD2"})
        net.add_edge("flankL", f"HD1{s}")
        net.add_edge(f"HD1{s}", f"HD2{s}")
    return net


def complex_case3() -> TestNetwork:
    """3 arms, all closed, but the 3rd arm joins from the opposite side
    (asymmetric — c-arm goes flankR → HD1c → HD2c → flankL)."""
    net = TestNetwork(name="complex_case3", expected="complexed")
    net.add_node("flankL", label="flankL")
    net.add_node("flankR", label="flankR")
    for s in "ab":
        net.add_node(f"HD1{s}", label="HD1", vars={"HD1"})
        net.add_node(f"HD2{s}", label="HD2", vars={"HD2"})
        net.add_edge("flankL", f"HD1{s}")
        net.add_edge(f"HD1{s}", f"HD2{s}")
        net.add_edge(f"HD2{s}", "flankR")
    # 3rd arm
    net.add_node("HD1c", label="HD1", vars={"HD1"})
    net.add_node("HD2c", label="HD2", vars={"HD2"})
    net.add_edge("flankL", f"HD2c")
    net.add_edge("HD1c", "HD2c")
    net.add_edge("HD1c", "flankR")
    return net


def complex_case4() -> TestNetwork:
    """Main arm closed, with a dangling HD1c-HD2c branch off HD1a, plus
    a separate HD1b-HD2b chain attached to flankR only.

           HD1c — HD2c
          /
    flankL — HD1a — HD2a — flankR
                    HD1b — HD2b /"""
    net = TestNetwork(name="complex_case4", expected="complexed")
    net.add_node("flankL", label="flankL")
    net.add_node("flankR", label="flankR")
    net.add_node("HD1a", label="HD1", vars={"HD1"})
    net.add_node("HD2a", label="HD2", vars={"HD2"})
    net.add_node("HD1c", label="HD1", vars={"HD1"})
    net.add_node("HD2c", label="HD2", vars={"HD2"})
    net.add_node("HD1b", label="HD1", vars={"HD1"})
    net.add_node("HD2b", label="HD2", vars={"HD2"})
    for e in [("flankL", "HD1a"), ("HD1a", "HD2a"), ("HD2a", "flankR"),
              ("HD1a", "HD1c"), ("HD1c", "HD2c"),
              ("HD1b", "HD2b"), ("HD2b", "flankR")]:
        net.add_edge(*e)
    return net


def complex_case5() -> TestNetwork:
    """HD1a connects to BOTH HD2a and HD2b (a fork at the gene-1 layer).
    HD1c-HD2c dangles off the top.

           HD1c — HD2c
          /
    flankL — HD1a — HD2a — flankR
                  \\HD2b /"""
    net = TestNetwork(name="complex_case5", expected="complexed")
    net.add_node("flankL", label="flankL")
    net.add_node("flankR", label="flankR")
    net.add_node("HD1a", label="HD1", vars={"HD1"})
    net.add_node("HD2a", label="HD2", vars={"HD2"})
    net.add_node("HD2b", label="HD2", vars={"HD2"})
    net.add_node("HD1c", label="HD1", vars={"HD1"})
    net.add_node("HD2c", label="HD2", vars={"HD2"})
    for e in [("flankL", "HD1a"), ("HD1a", "HD2a"), ("HD1a", "HD2b"),
              ("HD2a", "flankR"), ("HD2b", "flankR"),
              ("HD1a", "HD1c"), ("HD1c", "HD2c")]:
        net.add_edge(*e)
    return net


def closed_bubble_with_connector_joints() -> TestNetwork:
    """Real-data shape (AG17-like). Bubble joints are UNLABELED connector
    segments; flanks attach to the joints from outside the cycle.

      flankL — Jl — HDa — Jr — flankR
                 \\ HDb /

    Jl (left joint) and Jr (right joint) are pure connectors shared between
    both arms. Each arm goes through its own var node (HDa or HDb)."""
    net = TestNetwork(name="closed_bubble_with_connector_joints",
                       expected="closed_bubble")
    net.add_node("flankL", label="flankL")
    net.add_node("flankR", label="flankR")
    net.add_node("Jl")
    net.add_node("Jr")
    net.add_node("HDa", label="HD", vars={"HD"})
    net.add_node("HDb", label="HD", vars={"HD"})
    for e in [("flankL", "Jl"), ("Jl", "HDa"), ("HDa", "Jr"), ("Jr", "flankR"),
              ("Jl", "HDb"), ("HDb", "Jr")]:
        net.add_edge(*e)
    return net


def closed_bubble_ag17_shape() -> TestNetwork:
    """Exact AG17 shape: arm1 has 1 var seg, arm2 has 3 var segs, joints are
    unlabeled connectors, flankR is a composite (var+flankR)."""
    net = TestNetwork(name="closed_bubble_ag17_shape", expected="closed_bubble")
    net.add_node("flankL", label="flankL")
    net.add_node("Jl")
    net.add_node("Jr")
    net.add_node("HD_arm1", label="HD1+HD2", vars={"HD1", "HD2"})
    net.add_node("HD_arm2a", label="HD1", vars={"HD1"})
    net.add_node("HD_arm2b", label="HD1+HD2", vars={"HD1", "HD2"})
    net.add_node("flankR_comp", label="HD2+flankR", vars={"HD2"})
    for e in [("flankL", "Jl"),
              ("Jl", "HD_arm1"), ("HD_arm1", "Jr"),
              ("Jl", "HD_arm2a"), ("HD_arm2a", "HD_arm2b"), ("HD_arm2b", "Jr"),
              ("Jr", "flankR_comp")]:
        net.add_edge(*e)
    return net


def closed_bubble_long_flank_chain() -> TestNetwork:
    """Long flankL chain alternating with unlabeled connectors, terminating in
    a composite var+flankL+flankR node. Two arms diverge only at the final
    composite. Mirrors 2019NZAus082.1 / AG12-shape real data.

    flankL — fL_link1 — flankL — fL_link2 — flankL — composite_a(HD+flankL+flankR)
                                                    composite_b(HD+flankL+flankR)
    """
    net = TestNetwork(name="closed_bubble_long_flank_chain",
                       expected="closed_bubble")
    net.add_node("flankL_1", label="flankL")
    net.add_node("flankL_2", label="flankL")
    net.add_node("flankL_3", label="flankL")
    net.add_node("link_1")
    net.add_node("link_2")
    net.add_node("composite_a", label="HD1+HD2+flankL+flankR",
                  vars={"HD1", "HD2"})
    net.add_node("composite_b", label="HD1+HD2+flankL+flankR",
                  vars={"HD1", "HD2"})
    for e in [("flankL_1", "link_1"), ("link_1", "flankL_2"),
              ("flankL_2", "link_2"), ("link_2", "flankL_3"),
              ("flankL_3", "composite_a"), ("flankL_3", "composite_b")]:
        net.add_edge(*e)
    return net


def closed_bubble_long_unlabeled_to_flank() -> TestNetwork:
    """Long unlabeled connector chain (= MAX_LINKER_PADDING) between var
    content and the opposite flank. Mirrors SA93-shape (4–5 unlabeled
    intermediates between flankR and the nearest var).

    flankL — HDa — HDb — link — link — link — link — link — flankR
           \\ HDa2 — HDb2 /
    """
    net = TestNetwork(name="closed_bubble_long_unlabeled_to_flank",
                       expected="closed_bubble")
    net.add_node("flankL", label="flankL")
    net.add_node("flankR", label="flankR")
    net.add_node("HDa", label="HD1", vars={"HD1"})
    net.add_node("HDb", label="HD2", vars={"HD2"})
    net.add_node("HDa2", label="HD1", vars={"HD1"})
    net.add_node("HDb2", label="HD2", vars={"HD2"})
    for i in range(5):
        net.add_node(f"link_{i}")
    for e in [("flankL", "HDa"), ("HDa", "HDb"),
              ("flankL", "HDa2"), ("HDa2", "HDb2"),
              ("HDb", "link_0"), ("HDb2", "link_0"),
              ("link_0", "link_1"), ("link_1", "link_2"),
              ("link_2", "link_3"), ("link_3", "link_4"),
              ("link_4", "flankR")]:
        net.add_edge(*e)
    return net


def closed_bubble_var_spine_shared() -> TestNetwork:
    """Two arms share part of the var "spine" — the HD2 and HD2+flankR
    composite are the SAME node for both arms; only HD1 differs. Mirrors
    AG5-shape with shared var content.

    flankL — HD1a — link — HD2_shared — link — HD2_flankR_shared
           \\ HD1b /
    """
    net = TestNetwork(name="closed_bubble_var_spine_shared",
                       expected="closed_bubble")
    net.add_node("flankL", label="flankL")
    net.add_node("HD1a", label="HD1+flankL", vars={"HD1"})
    net.add_node("HD1b", label="HD1+flankL", vars={"HD1"})
    net.add_node("link_1")
    net.add_node("HD2_shared", label="HD2", vars={"HD2"})
    net.add_node("link_2")
    net.add_node("HD2_R_shared", label="HD2+flankR", vars={"HD2"})
    for e in [("flankL", "HD1a"), ("flankL", "HD1b"),
              ("HD1a", "link_1"), ("HD1b", "link_1"),
              ("link_1", "HD2_shared"),
              ("HD2_shared", "link_2"),
              ("link_2", "HD2_R_shared")]:
        net.add_edge(*e)
    return net


def single_one_composite_node() -> TestNetwork:
    """A single composite var+flankL+flankR node with pure-flank chains on
    each side. The graph genuinely has one allele; picker forces 2-allele
    pair so picks-heuristic would mislabel. Mirrors AJB36-shape.

    flankL — flankL — composite(HD1+HD2+flankL+flankR) — flankR — flankR
    """
    net = TestNetwork(name="single_one_composite_node", expected="single")
    net.add_node("flankL_1", label="flankL")
    net.add_node("flankL_2", label="flankL")
    net.add_node("composite", label="HD1+HD2+flankL+flankR",
                  vars={"HD1", "HD2"})
    net.add_node("flankR_1", label="flankR")
    net.add_node("flankR_2", label="flankR")
    for e in [("flankL_1", "flankL_2"), ("flankL_2", "composite"),
              ("composite", "flankR_1"), ("flankR_1", "flankR_2")]:
        net.add_edge(*e)
    return net


def separate() -> TestNetwork:
    """Two completely disjoint flank-HD-flank chains."""
    net = TestNetwork(name="separate", expected="separate")
    for s in "ab":
        net.add_node(f"flankL{s}", label="flankL")
        net.add_node(f"flankR{s}", label="flankR")
        net.add_node(f"HD1{s}", label="HD1", vars={"HD1"})
        net.add_node(f"HD2{s}", label="HD2", vars={"HD2"})
        net.add_edge(f"flankL{s}", f"HD1{s}")
        net.add_edge(f"HD1{s}", f"HD2{s}")
        net.add_edge(f"HD2{s}", f"flankR{s}")
    return net


# All clean cases, grouped by ground truth
CLEAN_CASES = [
    closed_bubble,
    closed_bubble_with_connector_joints,
    closed_bubble_ag17_shape,
    closed_bubble_long_flank_chain,
    closed_bubble_long_unlabeled_to_flank,
    closed_bubble_var_spine_shared,
    open_bubble_case1,
    open_bubble_case2,
    complex_case1, complex_case2, complex_case3, complex_case4, complex_case5,
    separate,
    single_one_composite_node,
]
