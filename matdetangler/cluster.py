"""Mating-type clustering across many samples (`MATdetangler cluster`).

Two-step pipeline (port of cluster_align.py + cluster_cut.py from the Pcub project):

  1) align — extract each picked allele's "locus core" (span from leftmost-to-rightmost variable-gene
     tblastn hit +/- 50 bp; conserved flanks excluded), run ONE MAFFT alignment on all cores, and save
     the pairwise similarity matrix (sim = identical aligned columns / shorter-core length).
  2) cut   — single-linkage at THRESH on the saved similarity matrix; name clusters 1, 2, ... by size;
     emit allele_classification.tsv + allele_distance_matrix.tsv.

The align step is the expensive part (one MAFFT call); the cut step is cheap and re-runnable with
any threshold without re-aligning.
"""
from __future__ import annotations
import os, sys, subprocess, argparse, tempfile, collections, re
from .input_process import read_fasta
from . import blast_utils as bu

MAFFT = os.environ.get("MAFFT", "mafft")

def _all_picked_alleles(results_dir: str) -> dict[str, str]:
    """Collect every per-sample primary_alleles.fasta under results/<sample>/ and pool them."""
    out = {}
    for s in sorted(os.listdir(results_dir)):
        fa = os.path.join(results_dir, s, "primary_alleles.fasta")
        if not os.path.exists(fa): continue
        out.update(read_fasta(fa))
    return out

def _core_spans(alleles: dict[str, str], variable_proteins_fa: str,
                tblastn_pid: float = 30.0, tblastn_aa: int = 50) -> dict[str, str]:
    """For each allele, extract the locus core (variable-gene span +/- 50bp)."""
    if not alleles: return {}
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "a.fa")
        with open(sf, "w") as o:
            for n, s in alleles.items(): o.write(f">{n}\n{s}\n")
        db = bu.fasta_to_db(sf, t, name="db")
        rows = bu.tblastn_hits(variable_proteins_fa, db, min_pid=tblastn_pid,
                                min_aa=tblastn_aa, extra_outfmt="sstart send")
    pts = collections.defaultdict(list)
    for r in rows:
        try: pts[r[0]] += [int(r[4]), int(r[5])]
        except (IndexError, ValueError): pass
    out = {}
    for n, ps in pts.items():
        seq = alleles[n]
        lo, hi = max(0, min(ps) - 50), min(len(seq), max(ps) + 50)
        out[n] = seq[lo:hi]
    return out

def align(results_dir: str, queries_dir: str, out_dir: str,
          tblastn_pid: float = 30.0, tblastn_aa: int = 50, threads: int = 4) -> dict:
    """Run align step. Returns dict with paths to alignment + pairs + meta files."""
    os.makedirs(out_dir, exist_ok=True)
    proteins = os.path.join(queries_dir, "variable_proteins.fasta")
    alleles = _all_picked_alleles(results_dir)
    if not alleles:
        raise SystemExit(f"no primary_alleles.fasta found under {results_dir}/*/")
    cores = _core_spans(alleles, proteins, tblastn_pid, tblastn_aa)
    if not cores:
        raise SystemExit("no allele cores extracted (no variable-gene tblastn hits)")
    cores_fa = os.path.join(out_dir, "cores.fasta")
    with open(cores_fa, "w") as o:
        for n, s in cores.items():
            o.write(f">{n}\n{s}\n")
    # one MAFFT alignment
    aln_fa = os.path.join(out_dir, "cores.aln.fasta")
    with open(aln_fa, "w") as o:
        subprocess.run([MAFFT, "--auto", "--thread", str(threads), "--quiet", cores_fa],
                       stdout=o, check=True)
    # parse MSA + compute pairwise sim
    msa = read_fasta(aln_fa)
    names = list(msa)
    lens = {n: len(cores[n]) for n in names}
    pairs_tsv = os.path.join(out_dir, "cores_pairs.tsv")
    meta_tsv = os.path.join(out_dir, "cores_meta.tsv")
    with open(meta_tsv, "w") as o:
        o.write("allele\tcore_len\n")
        for n in names: o.write(f"{n}\t{lens[n]}\n")
    with open(pairs_tsv, "w") as o:
        o.write("a\tb\tsim\talnid\n")
        for i in range(len(names)):
            sa = msa[names[i]]
            for j in range(i + 1, len(names)):
                sb = msa[names[j]]
                m = al = 0
                for x, y in zip(sa, sb):
                    if x != "-" and y != "-":
                        al += 1
                        if x.upper() == y.upper(): m += 1
                sim = m / min(lens[names[i]], lens[names[j]]) if min(lens[names[i]], lens[names[j]]) else 0.0
                alnid = m / al if al else 0.0
                o.write(f"{names[i]}\t{names[j]}\t{sim:.4f}\t{alnid:.4f}\n")
    print(f"[cluster.align] {len(names)} cores, {len(names)*(len(names)-1)//2} pairs -> {out_dir}/")
    return {"cores_fa": cores_fa, "aln_fa": aln_fa, "pairs_tsv": pairs_tsv, "meta_tsv": meta_tsv}

def _degHD_of(allele_name: str, deghd_set: set[str]) -> bool:
    """Allele is degHD if its sample/<sample>/picks.tsv labeled it so."""
    return allele_name in deghd_set

def _collect_degHD_set(results_dir: str) -> set[str]:
    """Read each picks.tsv; collect names with type=element-like (degHD-labeled)."""
    deg = set()
    for s in sorted(os.listdir(results_dir)):
        pt = os.path.join(results_dir, s, "picks.tsv")
        if not os.path.exists(pt): continue
        with open(pt) as fh:
            hdr = next(fh, "").rstrip("\n").split("\t")
            try:
                i_type = hdr.index("type"); i_name = hdr.index("allele")
            except ValueError: continue
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) > max(i_type, i_name) and f[i_type] == "element-like":
                    # picks.tsv has "allele1" / "allele2"; allele_name in pairs.tsv is "Pcub_<sample>_allele<i>"
                    deg.add(f"Pcub_{f[0]}_{f[i_name]}")
    return deg

def cut(align_dir: str, results_dir: str, out_dir: str, thresh: float = 0.90,
        deghd_label: str = "degHD") -> dict:
    """Read pairs/meta from align_dir; cluster real alleles single-linkage at thresh."""
    os.makedirs(out_dir, exist_ok=True)
    pairs_tsv = os.path.join(align_dir, "cores_pairs.tsv")
    meta_tsv  = os.path.join(align_dir, "cores_meta.tsv")
    # nodes
    cores = []
    with open(meta_tsv) as fh:
        next(fh, None)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            cores.append(f[0])
    # degHD labels
    deghd = _collect_degHD_set(results_dir)
    # pairs
    sim = {}; alnid = {}
    with open(pairs_tsv) as fh:
        next(fh, None)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            key = tuple(sorted([f[0], f[1]]))
            sim[key] = float(f[2]); alnid[key] = float(f[3])
    # single-linkage on REAL alleles
    real = [c for c in cores if c not in deghd]
    parent = {c: c for c in real}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for (a, b), v in sim.items():
        if v >= thresh and a not in deghd and b not in deghd:
            parent[find(a)] = find(b)
    clusters = collections.defaultdict(list)
    for c in real: clusters[find(c)].append(c)
    order = sorted(clusters.values(), key=lambda m: (-len(m), sorted(m)[0]))
    mt = {}
    for i, members in enumerate(order, 1):
        for c in members: mt[c] = str(i)
    for c in cores:
        if c in deghd: mt[c] = deghd_label
    # write classification + distance matrix
    cls = os.path.join(out_dir, "allele_classification.tsv")
    dm  = os.path.join(out_dir, "allele_distance_matrix.tsv")
    def sample_of(allele):
        m = re.match(r"Pcub_(.+)_allele\d+", allele)
        return m.group(1) if m else "-"
    with open(cls, "w") as o:
        o.write("allele\tsample\tcluster\n")
        for c in sorted(cores):
            o.write(f"{c}\t{sample_of(c)}\t{mt.get(c, 'NA')}\n")
    # group rows numerically by cluster
    def mtkey(a):
        m = mt.get(a, "zz")
        return (int(m) if m.isdigit() else 9999, m, a)
    allo = sorted(cores, key=mtkey)
    with open(dm, "w") as o:
        o.write("allele\t" + "\t".join(allo) + "\n")
        for a in allo:
            row = [a]
            for b in allo:
                d = 0.0 if a == b else 1.0 - sim.get(tuple(sorted([a, b])), 0.0)
                row.append(f"{d:.3f}")
            o.write("\t".join(row) + "\n")
    n_mt = len(order); n_singletons = sum(1 for m in order if len(m) == 1)
    print(f"[cluster.cut] THRESH={thresh}: {n_mt} real clusters ({n_singletons} singletons) + {deghd_label}")
    print(f"  wrote {cls}\n  wrote {dm}")
    return {"classification_tsv": cls, "distance_matrix_tsv": dm, "n_mating_types": n_mt}

def _cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="action", required=True)
    a = sub.add_parser("align"); a.add_argument("--results-dir", required=True)
    a.add_argument("--queries-dir", required=True); a.add_argument("--out-dir", required=True)
    a.add_argument("--threads", type=int, default=4)
    c = sub.add_parser("cut"); c.add_argument("--align-dir", required=True)
    c.add_argument("--results-dir", required=True); c.add_argument("--out-dir", required=True)
    c.add_argument("--thresh", type=float, default=0.90)
    ns = p.parse_args(argv)
    if ns.action == "align":
        align(ns.results_dir, ns.queries_dir, ns.out_dir, threads=ns.threads)
    else:
        cut(ns.align_dir, ns.results_dir, ns.out_dir, thresh=ns.thresh)

if __name__ == "__main__":
    _cli()
