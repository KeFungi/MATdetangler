"""Additional perturbations: things that real GFAs do.

  - lump_var_flank:   merge an HD node and its adjacent flank into ONE composite-tagged node
                      (e.g. label "HD1+flankL", var_per={"HD1"})
  - lump_var_var:     merge two adjacent HD nodes into one composite (e.g. "HD1+HD2")
  - fragment_var:     split a single HD node into two adjacent nodes that BOTH carry that gene
  - fragment_flank:   same idea for flankL / flankR
  - extra_flank:      add a paralog flank in a disconnected subgraph (see function doc)
  - cut_edge:         delete a random edge (can shift class — caller decides new expected)

All return a deep copy. Class is preserved unless explicitly noted.

Allele-identity rule
--------------------
All var-gene nodes in test cases carry an allele letter suffix in their ID
(HD1a, HD2a → allele "a"; HD1b, HD2b → allele "b"; etc.). Perturbations that
COMBINE two var-bearing nodes into one composite (lump_var_var, and any
cascade of lumps via lump_var_flank → adjacent-var merge) are restricted to
SAME-ALLELE pairs. A chimeric composite like "HD1a+HD1b" doesn't occur in
real assemblies — the two alleles are distinct genomic entities — and
prohibiting it removes a class of artificial collapses in the metatest.
"""
from __future__ import annotations
import random, re
from .network import TestNetwork


_ALLELE_RX = re.compile(r"HD\d+([a-z])")


def _node_allele(node_id: str) -> str | None:
    """Extract the allele letter from a node ID. Returns the letter if all HD
    occurrences share the same allele; None if there are none (e.g. flanks) or
    if multiple distinct alleles appear (chimera — already-merged across
    alleles, which we DON'T want to merge further)."""
    found = _ALLELE_RX.findall(node_id)
    if not found:
        return None
    s = set(found)
    return found[0] if len(s) == 1 else None


def _same_allele(a: str, b: str) -> bool:
    """True iff both node IDs encode a known allele AND they match."""
    aa, ab = _node_allele(a), _node_allele(b)
    return aa is not None and ab is not None and aa == ab


def _merge_nodes(net: TestNetwork, a: str, b: str, merged_id: str | None = None) -> str:
    """Replace nodes a and b with a single composite node carrying the union of labels
    and var-gene tags. The edge (a, b) is removed; other edges to a and b are redirected
    to the new merged node."""
    if a not in net.nodes or b not in net.nodes:
        raise ValueError(f"merge: missing {a} or {b}")
    new = merged_id or f"{a}__{b}"
    labels = "+".join(t for t in (net.labels.get(a, ""), net.labels.get(b, "")) if t)
    vars_ = (net.var_per.get(a, set()) | net.var_per.get(b, set()))
    new_edges = set()
    for e in net.edges:
        if {a, b} == set(e): continue                                  # drop a-b edge itself
        rewritten = frozenset(new if x in (a, b) else x for x in e)
        if len(rewritten) == 2: new_edges.add(rewritten)               # skip self-loops
    net.edges = new_edges
    net.nodes.discard(a); net.nodes.discard(b)
    net.labels.pop(a, None); net.labels.pop(b, None)
    net.var_per.pop(a, None); net.var_per.pop(b, None)
    net.add_node(new, label=labels, vars=vars_ if vars_ else None)
    return new


def _is_pure_var(out: TestNetwork, n: str) -> bool:
    """True iff n carries a var gene AND no flank tag."""
    if not out.var_per.get(n): return False
    return not any(t.startswith("flank")
                   for t in out.labels.get(n, "").split("+"))


def _allele_pure_var_count(out: TestNetwork, allele: str) -> int:
    """Number of pure-var (no flank tag) nodes belonging to a given allele."""
    return sum(1 for n in out.var_per
               if _node_allele(n) == allele and _is_pure_var(out, n))


def lump_var_flank(net: TestNetwork, rng: random.Random | None = None) -> TestNetwork:
    """Merge one var node with an adjacent flank into a composite (HD1+flankL etc.).

    SEPARATION rule: the var node being merged must NOT be the last pure-var
    (no-flank-tag) node of its allele. Each allele has to keep at least one
    purely-var node so the two alleles stay distinguishable in the var-only
    subgraph after the merge."""
    rng = rng or random.Random(10)
    out = net.deep_copy(new_name=net.name + "+lumpVF")
    adj = out.adj()
    candidates = []
    for v in list(out.var_per):
        v_allele = _node_allele(v)
        if v_allele is None: continue
        # Refuse if v is the last pure-var of its allele
        if _allele_pure_var_count(out, v_allele) <= 1: continue
        for f in adj.get(v, set()):
            if any(t.startswith("flank")
                   for t in out.labels.get(f, "").split("+")):
                candidates.append((v, f))
    if not candidates: return out
    v, f = rng.choice(candidates)
    _merge_nodes(out, v, f)
    return out


def lump_var_var(net: TestNetwork, rng: random.Random | None = None) -> TestNetwork:
    """Merge two adjacent var nodes into one composite (HD1+HD2). The two genes get
    co-located on one segment — common when both genes fit inside one SPAdes contig.

    No same-allele constraint: ANY two adjacent var nodes can be merged,
    including cross-allele pairs (HD1a+HD1b). Cross-allele composites are
    structurally meaningful — they don't normally occur in SPAdes output,
    but stressing the classifier with them is informative."""
    rng = rng or random.Random(11)
    out = net.deep_copy(new_name=net.name + "+lumpVV")
    adj = out.adj()
    pairs = [(a, b) for a in list(out.var_per) for b in adj.get(a, set())
             if b in out.var_per and a < b]
    if not pairs: return out
    a, b = rng.choice(pairs)
    _merge_nodes(out, a, b)
    return out


def fragment_var(net: TestNetwork, rng: random.Random | None = None) -> TestNetwork:
    """Split a random var node V into V_p1 — V_p2; both carry V's gene tag.
    Edges to V are split: incoming neighbors of V get redistributed so the chain
    stays connected (one neighbor goes to V_p1, the rest to V_p2)."""
    rng = rng or random.Random(12)
    out = net.deep_copy(new_name=net.name + "+fragV")
    candidates = list(out.var_per)
    if not candidates: return out
    v = rng.choice(candidates)
    label = out.labels.get(v, "")
    vars_ = set(out.var_per.get(v, set()))
    adj = out.adj()
    nbrs = list(adj.get(v, set()))
    if len(nbrs) < 2:
        return out
    p1, p2 = f"{v}_p1", f"{v}_p2"
    # Pick a random partition of v's neighbors into two non-empty halves
    rng.shuffle(nbrs)
    split = max(1, rng.randint(1, len(nbrs) - 1))
    nbrs_p1, nbrs_p2 = nbrs[:split], nbrs[split:]
    out.nodes.discard(v)
    out.labels.pop(v, None); out.var_per.pop(v, None)
    out.edges = {e for e in out.edges if v not in e}
    out.add_node(p1, label=label, vars=vars_)
    out.add_node(p2, label=label, vars=vars_)
    out.add_edge(p1, p2)
    for n in nbrs_p1: out.add_edge(p1, n)
    for n in nbrs_p2: out.add_edge(p2, n)
    return out


def fragment_flank(net: TestNetwork, rng: random.Random | None = None) -> TestNetwork:
    """Same idea for a flank node: split flankL into flankL_p1 — flankL_p2."""
    rng = rng or random.Random(13)
    out = net.deep_copy(new_name=net.name + "+fragF")
    candidates = [n for n in out.nodes
                  if any(t.startswith("flank") for t in out.labels.get(n, "").split("+"))]
    if not candidates: return out
    f = rng.choice(candidates)
    label = out.labels.get(f, "")
    adj = out.adj()
    nbrs = list(adj.get(f, set()))
    if len(nbrs) < 2:
        return out
    p1, p2 = f"{f}_p1", f"{f}_p2"
    rng.shuffle(nbrs)
    split = max(1, rng.randint(1, len(nbrs) - 1))
    out.nodes.discard(f)
    out.labels.pop(f, None)
    out.edges = {e for e in out.edges if f not in e}
    out.add_node(p1, label=label)
    out.add_node(p2, label=label)
    out.add_edge(p1, p2)
    for n in nbrs[:split]: out.add_edge(p1, n)
    for n in nbrs[split:]: out.add_edge(p2, n)
    return out


def extra_flank(net: TestNetwork, side: str = "L",
                rng: random.Random | None = None,
                chain_len: int = 2) -> TestNetwork:
    """Add a paralog flank-tagged node in a DISCONNECTED subgraph (no path to
    the main locus). Optionally extends a short unlabeled noise chain off it.

    Why disconnected: an `extra_flank` simulates a paralog flank somewhere
    outside the locus that the BFS happened to pull in. The classifier's
    validity rule (a flank region is valid only if some flank in it reaches
    a var within MAX_LINKER_PADDING unlabeled intermediates) correctly
    excludes this — there's no var reachable from the disconnected subgraph
    at all. Attaching the extra via a chain shorter than the validity
    horizon would let it sneak through as a "real" flank, which violates
    the spec that this shouldn't shift class."""
    rng = rng or random.Random(14)
    out = net.deep_copy(new_name=net.name + f"+extraFlank{side}")
    base_id = f"extra_flank{side}_{rng.randrange(10**6, 10**7)}"
    out.add_node(base_id, label=f"flank{side}")
    prev = base_id
    for _ in range(chain_len):
        nid = f"extra_noise_{rng.randrange(10**6, 10**7)}"
        out.add_edge(prev, nid)
        prev = nid
    return out


def partial_lump_flank(net: TestNetwork, rng: random.Random | None = None) -> TestNetwork:
    """Insert a composite (flank + adjacent-gene) node between a flank F and its
    var neighbor V. After: F — composite(flankX+geneY) — V. The pure flank F and
    pure var V both still exist; the composite is a third node carrying BOTH tags
    AND the gene. This mirrors the realistic SPAdes case where the gene/flank
    boundary lands inside one segment — you end up with flankL AND flankL+HD1."""
    rng = rng or random.Random(16)
    out = net.deep_copy(new_name=net.name + "+partialLumpVF")
    adj = out.adj()
    candidates = []
    for f in list(out.nodes):
        toks = [t for t in out.labels.get(f, "").split("+") if t]
        if not toks or not any(t.startswith("flank") for t in toks): continue
        if out.var_per.get(f): continue                        # already composite — skip
        for v in adj.get(f, set()):
            if out.var_per.get(v):
                candidates.append((f, v))
    if not candidates: return out
    f, v = rng.choice(candidates)
    f_label = out.labels.get(f, "")
    v_genes = out.var_per.get(v, set())
    # New composite node carries flank label AND one gene tag from v
    pick_gene = sorted(v_genes)[0]
    composite_id = f"{f}__plus__{pick_gene}"
    composite_label = f_label + "+" + pick_gene
    out.add_node(composite_id, label=composite_label, vars={pick_gene})
    # Rewire: remove the direct F-V edge, add F-composite and composite-V edges
    out.edges.discard(frozenset((f, v)))
    out.add_edge(f, composite_id)
    out.add_edge(composite_id, v)
    return out


def cut_edge(net: TestNetwork, rng: random.Random | None = None) -> TestNetwork:
    """Delete a random edge. Caller takes responsibility for the new expected class —
    this is a destructive perturbation. Returns with expected='unknown'."""
    rng = rng or random.Random(15)
    out = net.deep_copy(new_name=net.name + "+cut", new_expected="unknown")
    if not out.edges: return out
    e = rng.choice(list(out.edges))
    out.edges.discard(e)
    return out
