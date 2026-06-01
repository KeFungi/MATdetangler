"""Per-allele coverage that masks repeat / non-locus content.

Why: under competitive `-k 1` mapping, an allele that carries a high-copy repeat (e.g. a
MITE inside a stitched allele) gobbles all the genome-wide repeat reads, inflating its
mean depth by 2-3x while the other allele looks normal. That looks like a dikaryon imbalance
but is an artifact of the reference itself.

Fix: report depth restricted to the variable-gene span (the "HD-core") and, optionally, to
the non-repeat portion of the allele. This matches the Pcub `finalqc_one.py` HD-core depth
that drove the original deliverables.

Outputs (overwrite coverage.tsv):
   #allele  len  mapped_reads  whole_breadth_pct  whole_meandepth  core_bp  core_meandepth
"""
from __future__ import annotations
import os, sys, argparse, subprocess, collections

def _samtools_view_h(bam: str) -> list[tuple[str, int]]:
    out = []
    o = subprocess.run(["samtools", "view", "-H", bam], stdout=subprocess.PIPE, text=True).stdout
    for ln in o.splitlines():
        if not ln.startswith("@SQ"): continue
        f = dict(x.split(":", 1) for x in ln.split("\t")[1:] if ":" in x)
        out.append((f["SN"], int(f["LN"])))
    return out

def _depth_per_pos(bam: str, name: str) -> dict[int, int]:
    o = subprocess.run(["samtools", "depth", "-a", "-r", name, bam],
                       stdout=subprocess.PIPE, text=True).stdout
    d = {}
    for ln in o.splitlines():
        p = ln.split("\t")
        if len(p) >= 3: d[int(p[1])] = int(p[2])
    return d

def _mapped_reads(bam: str, name: str) -> int:
    o = subprocess.run(["samtools", "view", "-c", "-F", "4", bam, name],
                       stdout=subprocess.PIPE, text=True).stdout.strip()
    try: return int(o)
    except ValueError: return 0

def _hits_intervals(query_fa: str, alleles_fa: str, kind: str,
                    min_len: int = 50, min_id: float = 25.0) -> dict[str, list[tuple[int, int]]]:
    """kind: 'tblastn' for protein->NT, 'blastn' for NT->NT."""
    if not query_fa or not os.path.exists(query_fa): return {}
    cmd = [kind, "-query", query_fa, "-subject", alleles_fa,
           "-outfmt", "6 sseqid sstart send length pident"]
    if kind == "blastn": cmd += ["-dust", "no"]
    o = subprocess.run(cmd, stdout=subprocess.PIPE, text=True).stdout
    iv: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    for ln in o.splitlines():
        f = ln.split("\t")
        if len(f) < 5: continue
        L = int(f[3]); pid = float(f[4])
        if L < min_len or pid < min_id: continue
        a, b = int(f[1]), int(f[2])
        iv[f[0]].append((min(a, b), max(a, b)))
    return iv

def _to_set(iv: list[tuple[int, int]]) -> set[int]:
    s: set[int] = set()
    for lo, hi in iv: s.update(range(lo, hi + 1))
    return s

def write_coverage(bam: str, alleles_fa: str, queries_dir: str,
                   repeats_fa: str | None, out_tsv: str) -> None:
    variable_prot = os.path.join(queries_dir, "variable_proteins.fasta")
    hd_iv = _hits_intervals(variable_prot, alleles_fa, "tblastn", min_len=50, min_id=25.0)
    repeat_iv = _hits_intervals(repeats_fa, alleles_fa, "blastn", min_len=50, min_id=80.0) if repeats_fa else {}
    with open(out_tsv, "w") as o:
        o.write("#allele\tlen\tmapped_reads\twhole_breadth_pct\twhole_meandepth"
                "\tcore_bp\tcore_meandepth\trepeat_bp\n")
        for name, L in _samtools_view_h(bam):
            dep = _depth_per_pos(bam, name)
            covered = sum(1 for v in dep.values() if v > 0)
            whole_breadth = 100.0 * covered / max(1, L)
            whole_mean = sum(dep.values()) / max(1, L)
            mapped = _mapped_reads(bam, name)
            hc = _to_set(hd_iv.get(name, []))
            rp = _to_set(repeat_iv.get(name, []))
            core_pos = sorted(p for p in hc if p not in rp)
            core_dep = [dep.get(p, 0) for p in core_pos]
            core_mean = sum(core_dep) / max(1, len(core_dep)) if core_dep else 0.0
            o.write(f"{name}\t{L}\t{mapped}\t{whole_breadth:.4f}\t{whole_mean:.3f}"
                    f"\t{len(core_pos)}\t{core_mean:.3f}\t{len(rp)}\n")
    print(f"[coverage_core] wrote {out_tsv}")

def _cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--bam", required=True)
    p.add_argument("--alleles", required=True)
    p.add_argument("--queries-dir", required=True)
    p.add_argument("--repeats", default=None)
    p.add_argument("--out-tsv", required=True)
    a = p.parse_args(argv)
    write_coverage(a.bam, a.alleles, a.queries_dir, a.repeats, a.out_tsv)

if __name__ == "__main__":
    _cli()
