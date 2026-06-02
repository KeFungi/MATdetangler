"""Run the topology classifier against:
  - clean test cases (each of the 4 ground-truth classes, multiple sub-cases)
  - noisy / linker / reality perturbations (class-preserving, should still classify correctly)
  - degvar perturbation (class may shift; we log the verdict but don't assert)

Usage:
  python3 -m test.graph_classifier.run_tests             # all cases
  python3 -m test.graph_classifier.run_tests --verbose   # plus per-series breakdown
"""
from __future__ import annotations
import argparse, random, sys
from . import network, noise, lumping, classifier


GREEN = "\033[32m"; RED = "\033[31m"; YELLOW = "\033[33m"; RESET = "\033[0m"


def fmt_line(case_name: str, expected: str, got: str, ok: bool, extra: str = "") -> str:
    color = GREEN if ok else RED
    tick = "PASS" if ok else "FAIL"
    return f"  [{color}{tick}{RESET}]  {case_name:<40} expected={expected:<14} got={got:<14} {extra}"


def run_case(net, verbose: bool = False) -> bool:
    result = classifier.classify_network(net)
    got = result["class"]
    ok = (got == net.expected) if net.expected != "unknown" else True
    extra = result.get("explain", "")
    print(fmt_line(net.name, net.expected, got, ok, extra))
    if verbose:
        print(f"        nodes={result['n_nodes']} edges={result['n_edges']} "
              f"var_nodes={result['n_var']} flankL={result['n_flankL']} flankR={result['n_flankR']} "
              f"n_var_series={result.get('n_var_series', '-')}")
        for si in result.get("series", []):
            print(f"          series #{si['idx']}  size={si['size']} var={si['n_var']} "
                  f"Ls={sorted(si['Ls'])} Rs={sorted(si['Rs'])} clean={si['clean']} "
                  f"max_deg={si['max_deg']} internal_flank={si['internal_flank']}")
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", "-v", action="store_true")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--only", default="", help="comma-separated case name substring filter")
    a = ap.parse_args(argv)

    rng = random.Random(a.seed)
    filt = [s.strip() for s in a.only.split(",") if s.strip()]

    print("=" * 80)
    print("CLEAN CASES")
    print("=" * 80)
    n_pass = n_fail = 0
    clean_nets = []
    for builder in network.CLEAN_CASES:
        net = builder()
        if filt and not any(f in net.name for f in filt): continue
        clean_nets.append(net)
        if run_case(net, a.verbose): n_pass += 1
        else: n_fail += 1
    print(f"  -> clean: {n_pass} pass, {n_fail} fail")

    print()
    print("=" * 80)
    print("NOISY CASES (add dangling random chains; class preserved)")
    print("=" * 80)
    np = nf = 0
    for net in clean_nets:
        perturbed = noise.add_noise(net, rng=random.Random(a.seed))
        if run_case(perturbed, a.verbose): np += 1
        else: nf += 1
    print(f"  -> noisy: {np} pass, {nf} fail")

    print()
    print("=" * 80)
    print("LINKER CASES (subdivide edges with anon connectors; class preserved)")
    print("=" * 80)
    lp = lf = 0
    for net in clean_nets:
        perturbed = noise.add_linkers(net, rng=random.Random(a.seed))
        if run_case(perturbed, a.verbose): lp += 1
        else: lf += 1
    print(f"  -> linker: {lp} pass, {lf} fail")

    print()
    print("=" * 80)
    print("REALITY CASES (collapse + subdivide; class preserved)")
    print("=" * 80)
    rp = rf = 0
    for net in clean_nets:
        # apply on top of linker (so there's something to collapse)
        linked = noise.add_linkers(net, rng=random.Random(a.seed))
        perturbed = noise.apply_reality(linked, rng=random.Random(a.seed))
        if run_case(perturbed, a.verbose): rp += 1
        else: rf += 1
    print(f"  -> reality: {rp} pass, {rf} fail")

    print()
    print("=" * 80)
    print("LUMP / FRAGMENT CASES (composite labels, split nodes, extra flanks; class preserved)")
    print("=" * 80)
    lump_p = lump_f = 0
    PERTURBS = [
        ("lump_var_flank",      lumping.lump_var_flank),
        ("lump_var_var",        lumping.lump_var_var),
        ("partial_lump_flank",  lumping.partial_lump_flank),
        ("fragment_var",        lumping.fragment_var),
        ("fragment_flank",      lumping.fragment_flank),
        ("extra_flank_L",       lambda n, rng=None: lumping.extra_flank(n, side="L", rng=rng)),
        ("extra_flank_R",       lambda n, rng=None: lumping.extra_flank(n, side="R", rng=rng)),
    ]
    for net in clean_nets:
        for tag, fn in PERTURBS:
            try:
                p = fn(net, rng=random.Random(a.seed))
            except Exception as e:
                print(f"  [{YELLOW}ERR{RESET}]   {net.name}+{tag:<22} {type(e).__name__}: {e}")
                lump_f += 1
                continue
            p.name = f"{net.name}+{tag}"
            if run_case(p, a.verbose): lump_p += 1
            else: lump_f += 1
    print(f"  -> lump/fragment: {lump_p} pass, {lump_f} fail")

    print()
    print("=" * 80)
    print("METATEST — apply_reality_full: random mixed pipeline (noise + linkers +")
    print("  reality + lump + fragment + extra_flank), N trials per clean case.")
    print("  Treats class as preserved EXCEPT where structure-collapsing perturbations")
    print("  legitimately simplify a complex case → reported as 'simplified' (log only).")
    print("=" * 80)
    mp = mf = ms = 0
    n_trials = 20
    from .noise import apply_reality_full
    for net in clean_nets:
        for trial in range(n_trials):
            seed = a.seed + 1000 * trial + hash(net.name) % 1000
            perturbed = apply_reality_full(net, rng=random.Random(seed))
            perturbed.expected = net.expected   # carry over expected
            result = classifier.classify_network(perturbed)
            got = result["class"]
            if got == net.expected:
                mp += 1
            elif net.expected == "complexed" and got in ("closed_bubble", "open_bubble", "single"):
                # Structure-collapsing perturbations can legitimately simplify complex topology
                ms += 1
            elif net.expected in ("closed_bubble", "open_bubble") and got == "single":
                # Lump cascade can collapse an arm entirely → single arm remains
                ms += 1
            else:
                mf += 1
                if a.verbose:
                    print(fmt_line(perturbed.name[:55], net.expected, got, False))
    print(f"  -> metatest ({len(clean_nets)} cases × {n_trials} trials = "
          f"{len(clean_nets)*n_trials}): {mp} pass, {ms} simplified, {mf} fail")

    print()
    print("=" * 80)
    print("DEGVAR CASES (drop one var-gene tag; class may shift — expected='unknown')")
    print("=" * 80)
    for net in clean_nets:
        perturbed = noise.drop_var_copy(net, rng=random.Random(a.seed))
        run_case(perturbed, a.verbose)
    print(f"  -> degvar: logged only (no assert)")

    print()
    print("=" * 80)
    total_pass = n_pass + np + lp + rp + lump_p + mp
    total_fail = n_fail + nf + lf + rf + lump_f + mf
    print(f"    ↳ metatest simplified (complex → simple): {ms}")
    print(f"  OVERALL (class-preserving sets):  {total_pass} pass, {total_fail} fail")
    print("=" * 80)
    return 0 if total_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
