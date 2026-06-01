"""Locate the variable genes inside the locus reference (via tblastn), and derive flanks
automatically as the locus sequence outside the HD-hit envelope.

Replaces the previous BED-driven design which silently produced broken proteins on
intron-containing genes. The user now supplies:

    --locus-ref   locus.fasta      single-record DNA spanning HD region + both flanks
    --proteins    HDs.fasta        curated HD/variable proteins (one per gene)

The pipeline:

    1. tblastn(proteins -> locus_ref): per protein, take the best hit and record its
       nucleotide range on locus_ref.
    2. HD envelope = [min(hit.start), max(hit.end)] across all variable genes.
    3. flankL = locus[0 : envelope.start]                       (everything 5' of the leftmost HD)
       flankR = locus[envelope.end : end_of_locus]              (everything 3' of the rightmost HD)
    4. Sanitize protein IDs so blast (which requires unique sseqids) is happy.
    5. Write variable_proteins.fasta + flankL.fasta + flankR.fasta + variable_nt.fasta
       (variable_nt is the NT span of each gene's hit on the locus, useful for the optional
        blastn fallback).

Outputs (in <outdir>/):
    variable_proteins.fasta   ID-sanitized copy of --proteins
    variable_nt.fasta         NT spans of the gene hits (per gene, strand-corrected)
    flankL.fasta              auto-derived 5' flank NT
    flankR.fasta              auto-derived 3' flank NT
"""
from __future__ import annotations
import os, re, statistics, subprocess, tempfile, argparse

# ---------- small helpers used by other modules (kept stable) ----------


# SPAdes-output regex. Header is "NODE_<id>_length_<L>_cov_<F>" where F is
# ALWAYS "\d+\.\d{6}" across every SPAdes run inspected in this project (83k+
# contigs, six SPAdes versions / sample dirs sampled). The `(?:\.\d+)?` part
# is a safety fallback for hypothetical integer-only cov; full decimal
# precision is preserved by float() on the captured string.
_SPADES_CONTIG_HEADER_RE = re.compile(r"length_(\d+)_cov_(\d+(?:\.\d+)?)")


def estimate_genome_cov_from_contigs(contigs_fa: str, min_len: int = 5000) -> float:
    """Median of the `cov_` field from SPAdes contigs.fasta headers, restricted
    to contigs >= min_len bp. Units: bp-equivalent (matches GFA `DP:f:` segment
    depths used downstream for repeat detection — `cov_repeat_factor * genome_cov`
    becomes apples-to-apples).

    Median (not mean) chosen because it's robust to the very long tail of
    high-cov collapsed-repeat contigs (e.g. cov ~13e6 on 46 bp tips in real data).

    Returns 0.0 if no contig is >= min_len (e.g. tiny/poor assembly). The
    pipeline treats genome_cov=0 as "disable coverage-based repeat detection",
    which degrades gracefully.
    """
    covs: list[float] = []
    with open(contigs_fa) as f:
        for ln in f:
            if not ln.startswith(">"): continue
            m = _SPADES_CONTIG_HEADER_RE.search(ln)
            if not m: continue
            if int(m.group(1)) < min_len: continue
            covs.append(float(m.group(2)))
    if not covs: return 0.0
    return statistics.median(covs)

def rc(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

def read_fasta(path: str) -> dict[str, str]:
    d, name = {}, None
    for ln in open(path):
        if not ln.strip(): continue
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            d[name] = []
        elif name is not None:
            d[name].append(ln.strip())
    return {k: "".join(v) for k, v in d.items()}

# ---------- protein ID sanitization ----------

_BAD = re.compile(r"[^A-Za-z0-9_.-]+")

def _sanitize_protein_ids(proteins_fa: str, out_fa: str) -> list[str]:
    """Write a sanitized copy of proteins_fa to out_fa, ensuring unique IDs.
    Returns the list of new (sanitized) IDs in input order. The original record DESCRIPTION
    (after the first space) is preserved verbatim — we only rewrite the ID portion.

    Strategy: take the first token of the header, replace illegal chars with '_', and append
    a 1-based suffix when IDs collide. So two records both headed `>NC_062999_cut - A-alpha`
    and `>NC_062999_cut - beta1-1` become `>NC_062999_cut_1 - A-alpha` and
    `>NC_062999_cut_2 - beta1-1`.
    """
    out_ids: list[str] = []
    base_counts: dict[str, int] = {}
    headers, bodies = [], []
    cur_body: list[str] = []
    with open(proteins_fa) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln: continue
            if ln.startswith(">"):
                if headers: bodies.append(cur_body); cur_body = []
                headers.append(ln[1:])
            else:
                cur_body.append(ln)
        if headers: bodies.append(cur_body)
    # plan IDs
    plan = []
    for h in headers:
        first = h.split()[0] if h.split() else "p"
        first = _BAD.sub("_", first).strip("_") or "p"
        n = base_counts.get(first, 0) + 1
        base_counts[first] = n
        plan.append((first, n))
    # apply suffix when the same base appears more than once
    final_ids = []
    seen_final: set[str] = set()
    base_total = {k: v for k, v in base_counts.items()}
    counters = {k: 0 for k in base_total}
    for base, _ in plan:
        if base_total[base] > 1:
            counters[base] += 1
            new = f"{base}_{counters[base]}"
        else:
            new = base
        # in case sanitization still collided across different originals, append another suffix
        while new in seen_final:
            new = f"{new}_x"
        final_ids.append(new); seen_final.add(new)
    with open(out_fa, "w") as out:
        for new_id, h, body in zip(final_ids, headers, bodies):
            # preserve the original description (everything after the first whitespace token)
            parts = h.split(maxsplit=1)
            desc = (" " + parts[1]) if len(parts) > 1 else ""
            out.write(f">{new_id}{desc}\n")
            for ln in body: out.write(ln + "\n")
    return final_ids

# ---------- tblastn → locus, derive HD coords ----------

def tblastn_locate(proteins_fa: str, locus_fa: str, threads: int = 4,
                   min_pid: float = 30.0, min_aa: int = 50) -> tuple[dict[str, dict], list[int]]:
    """Run tblastn(proteins → locus) and resolve paralog leakage.

    Returns:
        per_gene_map: {protein_id: {chrom, start, end, strand, pid, aln_aa,
                                     n_hsps_accepted}}
        intergenic_ranges: list of dicts describing the intergenic span between each
                            pair of adjacent genes (left-to-right), with fields
                            {left_gene, left_end, right_gene, right_start,
                             intergenic_len}. We do NOT split intergenic space into
                            either gene's "territory" — it's genuine non-gene sequence
                            and is left in the envelope as-is.

    Algorithm — confidence-ordered greedy with cross-protein non-overlap:
        1. tblastn all proteins → locus; keep HSPs that pass (pid, aln_aa) cutoffs.
        2. Rank HSPs by confidence: score = pid x aln_aa, descending.
        3. Walk in order. Accept each HSP unless it overlaps an already-accepted
           HSP from a DIFFERENT protein by > 0.5 x shorter-HSP length. Same-protein
           HSPs never conflict, so divergent N- or C-terminus HSPs extend the
           gene's span without inflating aln_aa.
        4. Per gene, span = (min(start), max(end)) over accepted same-protein HSPs.
           Annotation fields (pid, aln_aa, strand) come from the BEST accepted HSP.
        5. Inter-gene boundary = midpoint between adjacent accepted spans.

    Why confidence order. HD1 and HD2 are paralogs and each protein's query produces
    weak partial HSPs in the sister's region. Walking in confidence order means the
    high-quality real hits land their claims first; the low-quality leakage HSPs
    arrive later and get rejected because their range already belongs to a different
    protein. Earlier versions that took the naive union of all HSPs per protein
    inflated HD1's reported aln_aa from 643 aa (true protein length) to 720 aa.
    """
    import collections
    with tempfile.TemporaryDirectory() as t:
        db = os.path.join(t, "locus_db")
        subprocess.run(["makeblastdb", "-in", locus_fa, "-dbtype", "nucl", "-out", db],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        out = subprocess.run(
            ["tblastn", "-query", proteins_fa, "-db", db,
             "-evalue", "1e-5", "-num_threads", str(threads),
             "-outfmt", "6 qseqid sseqid pident length sstart send"],
            stdout=subprocess.PIPE, text=True).stdout
    by_qid: dict[str, list[dict]] = collections.defaultdict(list)
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 6: continue
        qid, sid, pid, ln_aa, ss, se = f[0], f[1], float(f[2]), int(f[3]), int(f[4]), int(f[5])
        if pid < min_pid or ln_aa < min_aa: continue
        strand = "+" if se >= ss else "-"
        s, e = sorted([ss, se])
        by_qid[qid].append({"chrom": sid, "start": s, "end": e,
                             "strand": strand, "pid": pid, "aln_aa": ln_aa,
                             "score": pid * ln_aa})
    # Algorithm: walk HSPs in CONFIDENCE order (score = pid * aln_aa, descending). Each
    # HSP makes a territorial claim for its protein. A lower-confidence HSP is rejected
    # if it overlaps an already-claimed HSP from a DIFFERENT protein by more than
    # overlap_frac of the shorter HSP's length. Same-protein HSPs never conflict, so
    # divergent N/C-terminus HSPs from the SAME protein still extend the span. The
    # non-overlapping boundary between adjacent genes is then derived from the accepted
    # claims (midpoint between the rightmost accepted HSP of gene A and the leftmost
    # accepted HSP of gene B).
    overlap_frac = 0.5
    all_hsps_ranked: list[tuple[str, dict]] = []
    for qid, hsps in by_qid.items():
        for h in hsps:
            all_hsps_ranked.append((qid, h))
    all_hsps_ranked.sort(key=lambda kv: -kv[1]["score"])
    accepted: list[tuple[str, dict]] = []
    for qid, h in all_hsps_ranked:
        bad = False
        for k_qid, k in accepted:
            if k_qid == qid or k["chrom"] != h["chrom"]: continue
            ov = max(0, min(h["end"], k["end"]) - max(h["start"], k["start"]) + 1)
            shorter = min(h["end"] - h["start"] + 1, k["end"] - k["start"] + 1)
            if shorter > 0 and ov / shorter > overlap_frac:
                bad = True; break
        if not bad: accepted.append((qid, h))
    accepted_by_qid: dict[str, list[dict]] = collections.defaultdict(list)
    for qid, h in accepted: accepted_by_qid[qid].append(h)
    # Per gene span: best HSP for annotation (chrom, strand, pid, aln_aa); span coords
    # extended by min(start) / max(end) over ALL accepted same-protein HSPs.
    raw_per_qid: dict[str, dict] = {}
    for qid, hsps in accepted_by_qid.items():
        best = max(hsps, key=lambda h: h["score"])
        chrom = best["chrom"]
        same_chr = [h for h in hsps if h["chrom"] == chrom]
        s = min(h["start"] for h in same_chr)
        e = max(h["end"]   for h in same_chr)
        raw_per_qid[qid] = {"chrom": chrom, "start": s, "end": e,
                             "strand": best["strand"], "pid": best["pid"],
                             "aln_aa": best["aln_aa"],
                             "n_hsps_accepted": len(same_chr)}
    # Intergenic regions between adjacent accepted gene spans (sorted left-to-right).
    # Each entry is (left_gene_qid, left_end, right_gene_qid, right_start). When the
    # two adjacent spans abut or overlap, the entry's length is 0 or negative — we
    # still report it (length <= 0 means "no intergenic seq between these genes").
    sorted_genes = sorted(raw_per_qid.items(), key=lambda kv: kv[1]["start"])
    intergenic: list[dict] = []
    for i in range(len(sorted_genes) - 1):
        a_qid, a = sorted_genes[i]
        b_qid, b = sorted_genes[i + 1]
        intergenic.append({
            "left_gene":    a_qid, "left_end":     a["end"],
            "right_gene":   b_qid, "right_start":  b["start"],
            "intergenic_len": max(0, b["start"] - a["end"] - 1),
        })
    return raw_per_qid, intergenic

def _write_record(path: str, name: str, seq: str) -> None:
    with open(path, "w") as o:
        o.write(f">{name}\n")
        for i in range(0, len(seq), 80): o.write(seq[i:i+80] + "\n")

# ---------- main ----------

def _files_equal(a: str, b: str) -> bool:
    """Byte-for-byte file comparison."""
    if not (os.path.exists(a) and os.path.exists(b)): return False
    if os.path.getsize(a) != os.path.getsize(b): return False
    with open(a, "rb") as fa, open(b, "rb") as fb:
        while True:
            ba, bb = fa.read(65536), fb.read(65536)
            if ba != bb: return False
            if not ba: return True


def write_queries(locus_fa: str, proteins_fa: str, outdir: str,
                  threads: int = 4, max_flank_len: int = 2000,
                  envelope_padding: int = 500,
                  min_flank_len: int = 100,
                  cache_dir: str | None = None) -> dict:
    """Run tblastn(proteins → locus), derive flanks TRIMMED to `max_flank_len` each side,
    derive the path-enumerator bp budget, and write the manifest.

    Cache behavior — OPT-IN, only when `cache_dir` is given (matches the wrapper's
    `--caching DIR` / blast_cache pattern). When cache_dir is set:
        cache_dir/input_process/
            locus.fasta      — exact bytes of the locus input from the producing run
            proteins.fasta   — exact bytes of the proteins input from the producing run
            manifest.json    — the manifest the producing run wrote
            variable_proteins.fasta, variable_nt.fasta, flankL.fasta, flankR.fasta
    On the NEXT call with the same `cache_dir`, if the incoming locus.fasta + proteins.fasta
    byte-match the saved copies AND the same params (max_flank_len/envelope_padding/
    min_flank_len) are in the saved manifest, we copy the saved outputs to `outdir` and
    return the saved manifest — skipping tblastn + envelope + flank work.

    When `cache_dir` is None (default), every call runs fresh.

    `max_flank_len` (default 2000 bp per side) caps how much locus we keep flanking the HD
    envelope. It also drives:
        derived_max_locus_len = envelope_size + len(trimmed_flankL) + len(trimmed_flankR)
    which downstream modules (graph_path_search) read off the manifest as their bp budget
    when enumerating allele paths in the GFA. Smaller flank ⇒ tighter search.
    """
    import json, shutil
    os.makedirs(outdir, exist_ok=True)
    cache_entry = os.path.join(cache_dir, "input_process") if cache_dir else None
    if cache_entry:
        saved_locus     = os.path.join(cache_entry, "locus.fasta")
        saved_proteins  = os.path.join(cache_entry, "proteins.fasta")
        saved_manifest  = os.path.join(cache_entry, "manifest.json")
        if (os.path.exists(saved_manifest)
            and _files_equal(locus_fa, saved_locus)
            and _files_equal(proteins_fa, saved_proteins)):
            try:
                m = json.load(open(saved_manifest))
            except (json.JSONDecodeError, OSError):
                m = None
            if (m
                and m.get("max_flank_len_arg")    == max_flank_len
                and m.get("envelope_padding_arg") == envelope_padding
                and m.get("min_flank_len_arg",     min_flank_len) == min_flank_len):
                # Hit: copy the cached outputs into outdir, return the cached manifest.
                for fname in ("variable_proteins.fasta", "variable_nt.fasta",
                              "flankL.fasta", "flankR.fasta", "manifest.json"):
                    src = os.path.join(cache_entry, fname)
                    if os.path.exists(src):
                        shutil.copyfile(src, os.path.join(outdir, fname))
                print(f"[input_process] cache HIT  ({cache_entry}/) -> copied to {outdir}/, skipping tblastn")
                return m
            print(f"[input_process] cache STALE (param drift) at {cache_entry}/, recomputing")
        elif os.path.exists(saved_locus) or os.path.exists(saved_proteins):
            print(f"[input_process] cache STALE (input bytes differ) at {cache_entry}/, recomputing")
    # 1. sanitize protein IDs (handles duplicate sseqid case like NC_062999_cut x2)
    proteins_out = os.path.join(outdir, "variable_proteins.fasta")
    new_ids = _sanitize_protein_ids(proteins_fa, proteins_out)
    if not new_ids:
        raise SystemExit(f"--proteins: no records found in {proteins_fa}")
    # 2. tblastn(proteins → locus); returns per-gene span (extended via accepted same-
    #    protein HSPs after confidence-ordered cross-protein non-overlap resolution)
    #    AND the intergenic ranges between adjacent genes (descriptive only — we don't
    #    collapse intergenic seq into either gene's territory).
    hits, intergenic = tblastn_locate(proteins_out, locus_fa, threads=threads)
    missing = [i for i in new_ids if i not in hits]
    if missing:
        raise SystemExit(
            f"tblastn: these proteins didn't hit the locus reference at >=30% pid / >=50 aa: "
            f"{missing}. Either the locus is wrong, or the proteins are too divergent."
        )
    # 3. load locus, compute HD envelope
    locus_recs = read_fasta(locus_fa)
    if not locus_recs:
        raise SystemExit(f"--locus-ref: no records in {locus_fa}")
    if len(locus_recs) > 1:
        print(f"  warn: --locus-ref has {len(locus_recs)} records; using the first ({next(iter(locus_recs))})")
    chrom = next(iter(locus_recs))
    locus_seq = locus_recs[chrom]
    starts = [hits[i]["start"] for i in new_ids if hits[i]["chrom"] == chrom]
    ends   = [hits[i]["end"]   for i in new_ids if hits[i]["chrom"] == chrom]
    if not starts:
        raise SystemExit(f"none of the protein hits landed on the locus chrom {chrom}")
    # apply envelope_padding to compensate for tblastn missing the divergent N/C ends of a
    # protein at the pid cutoff (e.g., HD1's variable N-terminus when only the conserved
    # homeodomain hits at >= 30% pid)
    env_start = max(1, min(starts) - envelope_padding)              # 1-based, inclusive
    env_end   = min(len(locus_seq), max(ends) + envelope_padding)   # 1-based, inclusive
    envelope_size = env_end - env_start + 1
    # 4. TRIM flanks to max_flank_len each side, starting right at the padded envelope
    flankL_start = max(0, env_start - 1 - max_flank_len)
    flankL_end   = env_start - 1                                    # 0-based exclusive
    flankR_start = env_end                                          # 0-based inclusive
    flankR_end   = min(len(locus_seq), env_end + max_flank_len)
    flankL_nt = locus_seq[flankL_start:flankL_end]
    flankR_nt = locus_seq[flankR_start:flankR_end]
    if len(flankL_nt) < min_flank_len:
        print(f"  warn: trimmed flankL is only {len(flankL_nt)} bp (< {min_flank_len}); "
              f"the locus reference may not extend far enough 5' of the variable region")
    if len(flankR_nt) < min_flank_len:
        print(f"  warn: trimmed flankR is only {len(flankR_nt)} bp; locus may not extend far enough 3'")
    _write_record(os.path.join(outdir, "flankL.fasta"), "flankL", flankL_nt)
    _write_record(os.path.join(outdir, "flankR.fasta"), "flankR", flankR_nt)
    # 5. variable_nt.fasta — per-gene NT span on the locus, strand-corrected to the protein
    nt_path = os.path.join(outdir, "variable_nt.fasta")
    with open(nt_path, "w") as nt_out:
        for gid in new_ids:
            h = hits[gid]
            seg = locus_seq[h["start"] - 1 : h["end"]]
            if h["strand"] == "-": seg = rc(seg)
            nt_out.write(f">{gid}\n")
            for i in range(0, len(seg), 80): nt_out.write(seg[i:i+80] + "\n")
    # 6. derive the path-enumerator bp budget from envelope + trimmed flanks
    derived_max_locus_len = envelope_size + len(flankL_nt) + len(flankR_nt)
    manifest = {
        "chrom": chrom, "chrom_len": len(locus_seq),
        "max_flank_len_arg":     max_flank_len,
        "envelope_padding_arg":  envelope_padding,
        "min_flank_len_arg":     min_flank_len,
        "envelope_start":        env_start,
        "envelope_end":          env_end,
        "envelope_size":         envelope_size,
        "derived_max_locus_len": derived_max_locus_len,
        "variable_genes": [{"name": gid, "start": hits[gid]["start"], "end": hits[gid]["end"],
                            "strand": hits[gid]["strand"], "pid": hits[gid]["pid"],
                            "aln_aa": hits[gid]["aln_aa"],
                            "n_hsps_accepted": hits[gid].get("n_hsps_accepted", 1)}
                           for gid in new_ids],
        "intergenic_ranges": intergenic,
        "flankL": {"name": "flankL", "len": len(flankL_nt),
                   "locus_start": flankL_start + 1, "locus_end": flankL_end},
        "flankR": {"name": "flankR", "len": len(flankR_nt),
                   "locus_start": flankR_start + 1, "locus_end": flankR_end},
        "variable_proteins_fasta": proteins_out,
        "variable_nt_fasta":        nt_path,
        "flankL_fasta":             os.path.join(outdir, "flankL.fasta"),
        "flankR_fasta":             os.path.join(outdir, "flankR.fasta"),
    }
    with open(os.path.join(outdir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    # Cache write — only when cache_dir was supplied. Save inputs + outputs to the cache
    # entry so the next call can detect "inputs unchanged" and skip.
    if cache_entry:
        os.makedirs(cache_entry, exist_ok=True)
        shutil.copyfile(locus_fa,    saved_locus)
        shutil.copyfile(proteins_fa, saved_proteins)
        for fname in ("variable_proteins.fasta", "variable_nt.fasta",
                      "flankL.fasta", "flankR.fasta", "manifest.json"):
            src = os.path.join(outdir, fname)
            if os.path.exists(src): shutil.copyfile(src, os.path.join(cache_entry, fname))
        print(f"[input_process] cache WRITE -> {cache_entry}/")
    return manifest

# ---------- CLI ----------

def _cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--locus-ref", required=True,
                   help="single-record DNA reference: HD region + both flanks")
    p.add_argument("--proteins", required=True,
                   help="curated protein fasta (HD1, HD2, ...). Real proteins, not "
                        "BED-translated DNA. Records with duplicate IDs get a numeric suffix.")
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--max-flank-len", type=int, default=2000,
                   help="trim each flank to at most this many bp (default 2000). Also drives "
                        "the path-enumerator bp budget: envelope_size + len(flankL) + len(flankR)")
    p.add_argument("--envelope-padding", type=int, default=500,
                   help="bp padding added to each side of the tblastn HD envelope (default 500). "
                        "Compensates for tblastn missing the divergent N/C ends of a protein "
                        "below the pid cutoff (e.g. HD1's variable N-terminus).")
    p.add_argument("--min-flank-len", type=int, default=100,
                   help="warn if a trimmed flank is shorter than this (default 100 bp)")
    p.add_argument("--cache-dir", default=None,
                   help="if set, enable input_process caching: save locus.fasta + "
                        "proteins.fasta + outputs under <cache_dir>/input_process/ on the "
                        "first call, and skip on subsequent calls when the inputs and "
                        "params match byte-for-byte. Default off (always runs fresh).")
    a = p.parse_args(argv)
    m = write_queries(a.locus_ref, a.proteins, a.outdir,
                      threads=a.threads, max_flank_len=a.max_flank_len,
                      envelope_padding=a.envelope_padding,
                      min_flank_len=a.min_flank_len,
                      cache_dir=a.cache_dir)
    print(f"locus chrom: {m['chrom']} ({m['chrom_len']} bp)")
    print(f"variable genes ({len(m['variable_genes'])}):")
    for g in m["variable_genes"]:
        print(f"  {g['name']}: locus {g['start']}-{g['end']} strand={g['strand']} "
              f"tblastn pid={g['pid']:.1f} aln={g['aln_aa']}aa")
    print(f"HD envelope: {m['envelope_start']} .. {m['envelope_end']}  ({m['envelope_size']} bp)")
    print(f"flankL: {m['flankL']['len']} bp (locus {m['flankL']['locus_start']} .. {m['flankL']['locus_end']})")
    print(f"flankR: {m['flankR']['len']} bp (locus {m['flankR']['locus_start']} .. {m['flankR']['locus_end']})")
    print(f"derived max_locus_len = {m['derived_max_locus_len']} bp  (envelope + 2 trimmed flanks; "
          f"INFORMATIONAL — the wrapper's --max-locus-len CLI flag is authoritative for downstream stages)")
    print(f"wrote queries + manifest.json to {a.outdir}/")

if __name__ == "__main__":
    _cli()
