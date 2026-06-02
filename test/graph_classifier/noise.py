"""Graph perturbations that preserve topology class:
  - noisy:   add random dangling chains hanging off existing nodes (outward only)
  - linker:  insert 1-5 unlabeled connector nodes inside existing edges
  - reality: combine adjacent edges (collapse interior node) + subdivide some
  - degvar:  drop one var-gene tag (loses one HD copy)

Each transformation returns a deep copy; the original stays untouched.
"""
from __future__ import annotations
import random
from .network import TestNetwork
from .classifier import MAX_LINKER_PADDING


def _gen_id(rng: random.Random, used: set[str], prefix: str) -> str:
    while True:
        nid = f"{prefix}_{rng.randrange(10**6, 10**7)}"
        if nid not in used:
            return nid


def add_noise(net: TestNetwork, max_chains: int = 5, max_len: int = 5,
              rng: random.Random | None = None) -> TestNetwork:
    """Add up to `max_chains` dangling chains of length 1..max_len hanging off random
    existing nodes. The chains extend OUTWARD only (no back-links into clean network).
    Chain nodes carry no label and no var-gene tag — pure noise."""
    rng = rng or random.Random(0)
    out = net.deep_copy(new_name=net.name + "+noise")
    original_nodes = list(out.nodes)
    n_chains = rng.randint(1, max_chains)
    for _ in range(n_chains):
        anchor = rng.choice(original_nodes)
        length = rng.randint(1, max_len)
        prev = anchor
        for _ in range(length):
            nid = _gen_id(rng, out.nodes, "noise")
            out.add_edge(prev, nid)
            prev = nid
    return out


def add_linkers(net: TestNetwork, max_subdivs: int = 5,
                max_chain: int = MAX_LINKER_PADDING,
                rng: random.Random | None = None) -> TestNetwork:
    """Subdivide some clean-network edges with chains of 1..max_chain anonymous linker
    nodes. (a — b) becomes (a — link1 — ... — linkN — b). Linker nodes are unlabeled
    and carry no var-gene tag. Topology class is preserved.

    Default `max_chain` is bound to `MAX_LINKER_PADDING` so the test data never
    produces a flank-to-var linker chain longer than the classifier's validity
    horizon — keeps the two coupled."""
    rng = rng or random.Random(1)
    out = net.deep_copy(new_name=net.name + "+linker")
    edges_to_subdivide = list(out.edges)
    rng.shuffle(edges_to_subdivide)
    n = min(rng.randint(1, max_subdivs), len(edges_to_subdivide))
    for e in edges_to_subdivide[:n]:
        a, b = tuple(e)
        chain_len = rng.randint(1, max_chain)
        out.edges.discard(e)
        prev = a
        for _ in range(chain_len):
            nid = _gen_id(rng, out.nodes, "link")
            out.add_edge(prev, nid)
            prev = nid
        out.add_edge(prev, b)
    return out


def apply_reality(net: TestNetwork, rng: random.Random | None = None) -> TestNetwork:
    """Random structural noise that mimics real-graph quirks:
    - collapse some interior degree-2 unlabeled nodes (merge a-x-b into a-b)
    - subdivide some edges (a-b becomes a-x-y-b)
    Class-preserving by construction."""
    rng = rng or random.Random(2)
    out = net.deep_copy(new_name=net.name + "+reality")
    adj = out.adj()
    # Find degree-2 unlabeled non-var nodes (candidates to collapse)
    collapsibles = [n for n in out.nodes
                    if len(adj.get(n, ())) == 2
                    and not out.labels.get(n)
                    and not out.var_per.get(n)]
    rng.shuffle(collapsibles)
    n_collapse = rng.randint(0, min(3, len(collapsibles)))
    for n in collapsibles[:n_collapse]:
        nbrs = list(adj.get(n, ()))
        if len(nbrs) != 2 or n not in out.nodes:
            continue
        a, b = nbrs
        out.edges = {e for e in out.edges if n not in e}
        out.nodes.discard(n)
        out.labels.pop(n, None)
        out.var_per.pop(n, None)
        out.add_edge(a, b)
        adj = out.adj()
    # Now subdivide a random edge
    edges = list(out.edges)
    rng.shuffle(edges)
    for e in edges[:rng.randint(0, 2)]:
        a, b = tuple(e)
        out.edges.discard(e)
        x = _gen_id(rng, out.nodes, "real")
        y = _gen_id(rng, out.nodes, "real")
        out.add_edge(a, x)
        out.add_edge(x, y)
        out.add_edge(y, b)
    return out


def apply_reality_full(net: TestNetwork, rng: random.Random | None = None,
                        passes: int = 1) -> TestNetwork:
    """Random mixed perturbation — combines all class-preserving transformations
    into one pass. Per pass, randomly applies SOME of:
      - add_noise (dangling chains)
      - add_linkers (subdivide edges)
      - apply_reality (collapse degree-2 unlabeled + subdivide)
      - lump_var_flank (var + adjacent flank → composite)
      - lump_var_var (two adjacent vars → composite)
      - partial_lump_flank (insert flank+gene composite between flank and var)
      - fragment_var / fragment_flank (split into adjacent same-tag nodes)
      - extra_flank_L / extra_flank_R (add an isolated extra flank node)

    Class-preserving except in degenerate cases where structure-collapsing
    perturbations (e.g. lump_var_var on a complex case) simplify the topology.
    Test runner should treat the verdict as 'log only' for thoroughness."""
    from . import lumping
    rng = rng or random.Random(20)
    out = net
    OPS = [
        ("noise",          lambda n, r: add_noise(n, rng=r)),
        ("linkers",        lambda n, r: add_linkers(n, rng=r)),
        ("reality",        lambda n, r: apply_reality(n, rng=r)),
        ("lump_var_flank", lambda n, r: lumping.lump_var_flank(n, rng=r)),
        ("lump_var_var",   lambda n, r: lumping.lump_var_var(n, rng=r)),
        ("partial_lump",   lambda n, r: lumping.partial_lump_flank(n, rng=r)),
        ("fragment_var",   lambda n, r: lumping.fragment_var(n, rng=r)),
        ("fragment_flank", lambda n, r: lumping.fragment_flank(n, rng=r)),
        ("extra_L",        lambda n, r: lumping.extra_flank(n, side="L", rng=r)),
        ("extra_R",        lambda n, r: lumping.extra_flank(n, side="R", rng=r)),
    ]
    applied = []
    for _ in range(passes):
        k = rng.randint(2, max(2, len(OPS) // 2))
        chosen = rng.sample(OPS, k)
        for tag, fn in chosen:
            try:
                out = fn(out, rng)
                applied.append(tag)
            except Exception:
                continue
    out.name = net.name + "+REALITY[" + ",".join(applied) + "]"
    return out


def drop_var_copy(net: TestNetwork, gene: str | None = None,
                  rng: random.Random | None = None) -> TestNetwork:
    """Remove the var-gene tag from one randomly chosen copy of one gene. The node
    survives as an unlabeled non-var node (the "degenerated" copy). Topology may
    shift class because fewer var-bearing arms exist now — caller decides the
    expected class for the result."""
    rng = rng or random.Random(3)
    out = net.deep_copy(new_name=net.name + "+degvar", new_expected="unknown")
    # Group var nodes by gene name
    by_gene: dict[str, list[str]] = {}
    for n, genes in out.var_per.items():
        for g in genes:
            by_gene.setdefault(g, []).append(n)
    if gene is None:
        candidates = [g for g, copies in by_gene.items() if len(copies) >= 2]
        if not candidates: return out
        gene = rng.choice(candidates)
    if gene not in by_gene or len(by_gene[gene]) < 1: return out
    drop_node = rng.choice(by_gene[gene])
    # Remove the specific gene from var_per[drop_node]; if empty, drop the entry.
    s = out.var_per.get(drop_node, set())
    s.discard(gene)
    if s: out.var_per[drop_node] = s
    else: out.var_per.pop(drop_node, None)
    # Also strip the gene name from the label
    lab = out.labels.get(drop_node, "")
    toks = [t for t in lab.split("+") if t and t != gene]
    if toks: out.labels[drop_node] = "+".join(toks)
    else: out.labels.pop(drop_node, None)
    return out
