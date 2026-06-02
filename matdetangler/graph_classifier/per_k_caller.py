"""Per-k allele caller — one GFA + one label-hits TSV in, one allele set out.

Contract: a single k (Suilu k33 / k45 / k53; Pcub k55 / k77). Cross-k
integration is a downstream selection step that picks the best per-k
verdict + allele set.

Algorithm
=========

Phase 1 (main BFS — var-seeded):

    while nhop <= max_nhop:
        1. neighborhood = N-hop BFS from var segments in the GFA.
        2. Apply coverage filter [lo*genome_cov, hi*genome_cov] on segments.
        3. Classify the resulting bubble.
           If closed_bubble AND arms are divergent: break.
        4. Fallback: remove the coverage filter, classify again.
           If closed_bubble AND arms are divergent: break.
        5. Else: nhop += 1, loop.

Phase 2 (flank fallback — var + flank seeded):

    Same loop, but seed BFS from var nodes AND flank nodes.

Emission rule (uniform dedup across all topologies)
===================================================

Every topology produces a list of candidate sequences (one per arm /
component / simple-path). The candidate list is then collapsed by
pairwise divergence — sequences that are NOT divergent from an
already-kept representative drop out, and the longest representative
wins. The verdict **always preserves the classifier's origin**; only
the emitted record naming reflects the post-dedup count `n`:

    n == 1  → allele1                                    (verdict = origin)
    n == 2  → allele1 + allele2                          (verdict = origin)
    n >= 3  → chimera1, chimera2, … (one per dedup seq)  (verdict = origin)

Examples:
    single        (1 candidate)                          → allele1            verdict=single
    closed_bubble (2 arms, dedup→1)                      → allele1            verdict=closed_bubble
    closed_bubble (2 arms, dedup→2)                      → allele1+allele2    verdict=closed_bubble
    open_bubble   (2 arms, dedup→2)                      → allele1+allele2    verdict=open_bubble
    separate      (3 components, dedup→2)                → allele1+allele2    verdict=separate
    separate      (3 components, dedup→3)                → chimera1..3        verdict=separate
    complexed     (5 paths, dedup→1)                     → allele1            verdict=complexed
    complexed     (5 paths, dedup→2)                     → allele1+allele2    verdict=complexed
    complexed     (5 paths, dedup→4)                     → chimera1..4        verdict=complexed

Divergence
==========

Two reconstructed arm sequences are "divergent" if pairwise identity is
below `1 - divergence_threshold`. The default threshold is 0.05 (i.e.
arms must differ by ≥5% to count as distinct alleles).
"""
from __future__ import annotations
from collections import deque
from .seg_processor import directional_split
from .bubble_classifier import classify
from .bubble_bfs import build_adj
from .arm_sequences import parse_gfa_sequences, reconstruct_arms, reverse_complement
from .labeler import read_seg_label_hits, read_segment_lengths_from_hits


# -----------------------------------------------------------------------------
# GFA-level helpers
# -----------------------------------------------------------------------------

def parse_gfa(gfa_path: str
              ) -> tuple[dict[str, str], dict[str, float],
                          set[frozenset], dict[tuple, tuple[str, str]]]:
    """Parse a GFA. Returns (sequences, depths, edges, edge_endpoints)."""
    sequences: dict[str, str] = {}
    depths: dict[str, float] = {}
    edges: set[frozenset] = set()
    endpoints: dict[tuple, tuple[str, str]] = {}
    with open(gfa_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f: continue
            if f[0] == "S":
                if len(f) >= 3:
                    sequences[f[1]] = f[2]
                    for tag in f[3:]:
                        if tag.startswith("DP:f:"):
                            try: depths[f[1]] = float(tag[5:])
                            except ValueError: pass
                            break
            elif f[0] == "L" and len(f) >= 5:
                a, oa, b, ob = f[1], f[2], f[3], f[4]
                end_a = "R" if oa == "+" else "L"
                end_b = "L" if ob == "+" else "R"
                edges.add(frozenset((a, b)))
                key = tuple(sorted([a, b]))
                if a == key[0]:
                    endpoints[key] = (end_a, end_b)
                else:
                    endpoints[key] = (end_b, end_a)
    return sequences, depths, edges, endpoints


def parse_contigs_paths(path: str) -> dict[str, list[str]]:
    """Parse SPAdes contigs.paths. Returns {contig_name: [seg_id, ...]}.

    Format (modern SPAdes):
        NODE_1_length_5000_cov_45.6
        12+,3-,45+;
        78+,9-;
        NODE_1_length_5000_cov_45.6'                   ← RC variant
        45-,3+,12-
        ...

    Each contig has 1+ subpaths (semicolon-separated) of `seg_id±` tokens.
    We strip orientation marks and return the union of all visited seg IDs
    (orientation-insensitive). RC entries (trailing `'`) are merged into
    their parent contig.
    """
    walks: dict[str, list[str]] = {}
    cur: str | None = None
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line: continue
            if cur is None:
                cur = line.rstrip("'")                          # strip RC marker
                walks.setdefault(cur, [])
            else:
                for chunk in line.replace(";", ",").split(","):
                    seg = chunk.rstrip("+-").strip()
                    if seg:
                        walks[cur].append(seg)
                cur = None
    return walks


def contig_seeds_for_locus(
        contig_blast_tsvs: list[str],
        contigs_paths_file: str,
        min_pid: float = 80.0,
        min_alnlen: int = 100,
) -> set[str]:
    """Identify GFA segments anchored by any contig that BLAST-hits the locus.

    Reads outfmt-6 contig-BLAST TSVs (qseqid sseqid pident length ...),
    keeps subjects (sseqid = contig names) with pident ≥ min_pid AND
    length ≥ min_alnlen. Maps those contigs → segments via contigs.paths.
    Returns the union of segments across all locus-bearing contigs.
    """
    import os
    if not os.path.exists(contigs_paths_file):
        return set()
    walks = parse_contigs_paths(contigs_paths_file)
    locus_contigs: set[str] = set()
    for path in contig_blast_tsvs:
        if not os.path.exists(path): continue
        with open(path) as fh:
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) < 4: continue
                try:
                    pid = float(f[2]); ln_ = int(f[3])
                except ValueError: continue
                if pid >= min_pid and ln_ >= min_alnlen:
                    locus_contigs.add(f[1])
    segs: set[str] = set()
    for c in locus_contigs:
        # SPAdes contigs.paths may name contigs without prefix; also
        # strip the trailing _cov_X.X if header has it. Try both forms.
        if c in walks:
            segs.update(walks[c])
            continue
        # Heuristic: tolerate suffix differences
        cn = c.rstrip("'")
        for k, v in walks.items():
            if k == cn or k.startswith(cn + "_") or cn.startswith(k + "_"):
                segs.update(v)
                break
    return segs


def bfs_expand_segments(seeds: set[str], adj_und: dict[str, set[str]],
                         n_hops: int) -> set[str]:
    """N-hop BFS in the GFA undirected adjacency, returns the reachable set."""
    if n_hops <= 0: return set(seeds)
    seen = set(seeds)
    frontier = set(seeds)
    for _ in range(n_hops):
        nxt = set()
        for n in frontier:
            for m in adj_und.get(n, ()):
                if m not in seen:
                    nxt.add(m)
        if not nxt: break
        seen |= nxt
        frontier = nxt
    return seen


def filter_segs_by_depth(segs: set[str], depths: dict[str, float],
                          genome_cov: float,
                          lo_mult: float = 0.25,
                          hi_mult: float = 2.0) -> set[str]:
    """Keep only segments whose depth is within [lo_mult * genome_cov,
    hi_mult * genome_cov]. Segments without a depth value are kept (we don't
    want to throw them away just because we couldn't look up depth)."""
    lo, hi = lo_mult * genome_cov, hi_mult * genome_cov
    return {s for s in segs if s not in depths or lo <= depths[s] <= hi}


def restrict_subgraph(nodes: set[str], edges: set[frozenset],
                       endpoints: dict[tuple, tuple[str, str]]
                       ) -> tuple[set[frozenset], dict[tuple, tuple[str, str]]]:
    """Drop edges whose endpoints aren't both in `nodes`."""
    kept = {e for e in edges if all(x in nodes for x in tuple(e))}
    kept_ep = {k: v for k, v in endpoints.items()
                if k[0] in nodes and k[1] in nodes}
    return kept, kept_ep


# -----------------------------------------------------------------------------
# Divergence
# -----------------------------------------------------------------------------

def _identity(a: str, b: str) -> float:
    """Pairwise nucleotide identity via edlib infix (HW) edit distance,
    strand-agnostic. Query = shorter seq, target = longer (end gaps free).
    We try the target both forward AND reverse-complemented and take the
    closer alignment — graph walks through the same bubble can be emitted
    in either orientation depending on which anchor BFS started from, so
    two alleles that are reverse complements of each other must collapse
    in dedup. Returns 1.0 for identical (or fully contained), lower for
    divergent. Handles indels and end-fragmentation properly."""
    import edlib
    if not a or not b: return 0.0
    q, t = (a, b) if len(a) <= len(b) else (b, a)
    r_fwd = edlib.align(q, t, task="distance", mode="HW")
    r_rc  = edlib.align(q, reverse_complement(t), task="distance", mode="HW")
    best_ed = min(r_fwd["editDistance"], r_rc["editDistance"])
    return 1.0 - best_ed / max(len(q), 1)


def is_divergent(seq_a: str, seq_b: str, threshold: float = 0.05) -> bool:
    """True if the two sequences differ by at least `threshold` (edit-distance
    identity below 1 − threshold). Default 0.05 → pairs with ≥ 5% edit distance
    count as divergent alleles, < 5% collapse as duplicates."""
    if not seq_a or not seq_b: return False
    return (1.0 - _identity(seq_a, seq_b)) >= threshold


# -----------------------------------------------------------------------------
# Orchestrator
# -----------------------------------------------------------------------------

def _try_one_pass(seeds: set[str], all_edges: set[frozenset],
                   endpoints: dict[tuple, tuple[str, str]],
                   adj_und: dict[str, set[str]], n_hops: int,
                   seg_labels: dict[str, list], seg_length: dict[str, int],
                   depths: dict[str, float], gfa_sequences: dict[str, str],
                   genome_cov: float | None,
                   divergence_threshold: float,
                   apply_cov_filter: bool,
                   lo_mult: float, hi_mult: float) -> dict:
    """One pass: BFS-expand seeds N hops, optionally cov-filter, classify,
    return verdict dict with arm sequences attached if it classified."""
    nhood = bfs_expand_segments(seeds, adj_und, n_hops)
    if apply_cov_filter and genome_cov:
        nhood = filter_segs_by_depth(nhood, depths, genome_cov, lo_mult, hi_mult)
    pass_edges, pass_endpoints = restrict_subgraph(nhood, all_edges, endpoints)
    pass_seg_length = {s: seg_length.get(s, len(gfa_sequences.get(s, "")))
                       for s in nhood}
    pass_seg_labels = {s: seg_labels.get(s, []) for s in nhood}
    nodes, edges_pp, labels, var_per, provenance = directional_split(
        pass_seg_labels, pass_seg_length, pass_edges, pass_endpoints,
    )
    res = classify(nodes, edges_pp, labels, var_per)
    res["_provenance"] = provenance
    res["_nhood"] = nhood
    res["_var_per_node"] = var_per
    res["_label_per_node"] = labels
    # Adjacency over post-P1 nodes — needed by _emit_result to look up the
    # exterior neighbors of each candidate path's endpoints (for the
    # innermost-flank record).
    _adj_pp: dict[str, set[str]] = {}
    for e in edges_pp:
        t = tuple(e)
        if len(t) == 1: continue
        a, b = t
        _adj_pp.setdefault(a, set()).add(b)
        _adj_pp.setdefault(b, set()).add(a)
    res["_adj"] = _adj_pp
    return res


def find_alleles(
        seg_label_hits_tsv: str,
        gfa_path: str,
        genome_cov: float | None = None,
        init_nhop: int = 3,
        max_nhop: int = 10,
        divergence_threshold: float = 0.05,
        lo_mult: float = 0.25,
        hi_mult: float = 2.0,
        k: int | str | None = None,
        var_proteins_ref: str | None = None,
        expected_var_tags: set[str] | None = None,
        locus_padding: int = 1500,
        contig_seeds: set[str] | None = None,
        queries_dir: str | None = None,
) -> dict:
    """Run the full orchestrator. Returns a dict with the final classification
    and any emitted allele/chimera sequences:

      {
        "verdict": "closed_bubble" | "open_bubble" | "single" | "separate" | "complexed",
        "divergent": bool,
        "n_hops_used": int,
        "phase": "main" | "flank_fallback",
        "alleles": list[(name, sequence)],
        "info": <last classify info>,
      }
    """
    gfa_seqs, depths, all_edges, endpoints = parse_gfa(gfa_path)
    seg_labels = read_seg_label_hits(seg_label_hits_tsv)
    seg_length = read_segment_lengths_from_hits(seg_label_hits_tsv)
    for s, seq in gfa_seqs.items():
        seg_length.setdefault(s, len(seq))

    adj_und: dict[str, set[str]] = {}
    for e in all_edges:
        t = tuple(e)
        if len(t) == 1: continue                            # self-loop in GFA
        a, b = t
        adj_und.setdefault(a, set()).add(b)
        adj_und.setdefault(b, set()).add(a)

    # Seeds for the two phases. Optional contig_seeds = GFA segments walked
    # by SPAdes contigs that BLAST-hit the locus (HD/flank); used to anchor
    # the BFS in samples where the labeler's segment-level hits split the
    # locus across disconnected GFA components.
    contig_seeds = set(contig_seeds or ())
    var_seeds = {s for s, hits in seg_labels.items()
                 if any(h.kind == "var" for h in hits)} | contig_seeds
    flank_seeds = {s for s, hits in seg_labels.items()
                   if any(h.kind == "flank" for h in hits)} | contig_seeds

    def _accept(res: dict) -> bool:
        if res["class"] != "closed_bubble": return False
        arms = res.get("closed_arms", [])
        if len(arms) < 2: return False
        seqs = [_arm_sequence_for(a, res["_provenance"], gfa_seqs) for a in arms[:2]]
        if is_divergent(seqs[0], seqs[1], divergence_threshold):
            res["_emit_arms"] = arms[:2]
            res["_emit_seqs"] = seqs
            res["divergent"] = True
            return True
        return False

    def _log(phase: str, nhop: int, use_cov: bool, res: dict) -> None:
        n_arms = res.get("n_arms", 0)
        n_var  = res.get("n_var", 0)
        nhood  = len(res.get("_nhood", ()))
        print(f"  [{phase} nhop={nhop} cov={'on' if use_cov else 'off':>3}] "
              f"|nhood|={nhood:<6} var={n_var:<3} cls={res['class']:<14} arms={n_arms}",
              flush=True)

    last_res = None
    fallback_cov_on = None                     # cov-on @ max_nhop, for the no-acceptance path

    def _phase_loop(phase: str, seeds: set[str]) -> dict | None:
        nonlocal last_res, fallback_cov_on
        for nhop in range(init_nhop, max_nhop + 1):
            for use_cov in (True, False):      # cov-on first → if no acceptance, we replay cov-on at max below
                res = _try_one_pass(
                    seeds, all_edges, endpoints, adj_und, nhop,
                    seg_labels, seg_length, depths, gfa_seqs,
                    genome_cov, divergence_threshold,
                    apply_cov_filter=use_cov,
                    lo_mult=lo_mult, hi_mult=hi_mult,
                )
                res["_phase"]       = phase
                res["_n_hops"]      = nhop
                res["_cov_filter"]  = use_cov
                _log(phase, nhop, use_cov, res)
                last_res = res
                if use_cov and nhop == max_nhop:
                    fallback_cov_on = res      # remember the cov-filtered @ max for emission fallback
                if _accept(res):
                    return res
        return None

    def _finalize(res: dict, divergence_threshold: float) -> dict:
        # _emit_result now does the locus-trim + completeness check inline.
        out = _emit_result(
            res, gfa_seqs, depths=depths,
            divergence_threshold=divergence_threshold,
            var_proteins_ref=var_proteins_ref,
            locus_padding=locus_padding,
            expected_var_tags=expected_var_tags,
            genome_cov=genome_cov, lo_mult=lo_mult, hi_mult=hi_mult,
            queries_dir=queries_dir,
        )
        out["k"] = k
        out["genome_cov"] = genome_cov
        out["bubble_type"] = out["verdict"]                            # alias for TSV clarity
        # Segments + segment-label list, from the result's emitted candidates.
        prov = res.get("_provenance", {})
        emitted_subnodes: set[str] = set()
        for p in res.get("closed_arms", []):    emitted_subnodes.update(p)
        for p in res.get("dangling_arms", []):  emitted_subnodes.update(p)
        for c in res.get("var_components") or []: emitted_subnodes.update(c)
        emitted_segs: set[str] = {prov[n][0] for n in emitted_subnodes if n in prov}
        labels_by_seg: dict[str, str] = {}
        for s in emitted_segs:
            tags = sorted({h.tag for h in seg_labels.get(s, [])})
            if tags: labels_by_seg[s] = "+".join(tags)
        out["component_list"]   = [n for n, _ in out["alleles"]]
        out["segments"]         = sorted(emitted_segs)
        out["segments_labeled"] = [(s, labels_by_seg.get(s, "")) for s in sorted(emitted_segs)]
        out["basepair"]         = sum(len(s) for _, s in out["alleles"])
        return out

    accepted = _phase_loop("main", var_seeds)
    if accepted is not None:
        return _finalize(accepted, divergence_threshold)

    accepted = _phase_loop("flank_fallback", var_seeds | flank_seeds)
    if accepted is not None:
        return _finalize(accepted, divergence_threshold)

    # Loop exhausted — fallback chain:
    #   1) cov-on @ max_nhop (preferred — selective)
    #   2) if cov-on yielded 0 candidates, fall back to cov-off @ max_nhop (last_res)
    def _has_candidates(r: dict) -> bool:
        return (bool(r.get("closed_arms")) or bool(r.get("dangling_arms"))
                or bool(r.get("var_components")))

    emit_from = fallback_cov_on or last_res
    if emit_from is fallback_cov_on and not _has_candidates(fallback_cov_on):
        emit_from = last_res
        print(f"  [exhausted] cov-on @ max_nhop emitted 0 candidates; falling back to cov-off",
              flush=True)
    print(f"  [exhausted] emitting from phase={emit_from.get('_phase')} "
          f"nhop={emit_from.get('_n_hops')} cov={emit_from.get('_cov_filter')}",
          flush=True)
    return _finalize(emit_from, divergence_threshold)


def _arm_sequence_for(arm_path: list[str], provenance: dict,
                       gfa_seqs: dict[str, str]) -> str:
    """Build a sequence from one arm path, without joint-skipping.
    Used for the orchestrator's _accept divergence check (full arm content)."""
    parts = []
    for n in arm_path:
        if n not in provenance: continue
        seg, start, end, strand = provenance[n]
        sub = gfa_seqs.get(seg, "")[start:end]
        if strand == "-":
            sub = reverse_complement(sub)
        parts.append(sub)
    return "".join(parts)


def _arm_sequence_for_skip(arm_path: list[str], provenance: dict,
                            gfa_seqs: dict[str, str],
                            skip_nodes: set[str]) -> str:
    """Build a sequence from one arm path, omitting nodes in `skip_nodes`
    (joint/fork nodes shared with other candidate paths). Used by emission
    so that each candidate's sequence reflects only its distinguishing
    (non-shared) content."""
    parts = []
    for n in arm_path:
        if n in skip_nodes: continue
        if n not in provenance: continue
        seg, start, end, strand = provenance[n]
        sub = gfa_seqs.get(seg, "")[start:end]
        if strand == "-":
            sub = reverse_complement(sub)
        parts.append(sub)
    return "".join(parts)


def _dedup_sequences(seqs: list[str],
                      divergence_threshold: float = 0.05) -> list[str]:
    """Collapse non-divergent sequences. Walk in length-desc order so the
    longest representative of each equivalence class is kept."""
    nonempty = [s for s in seqs if s]
    nonempty.sort(key=len, reverse=True)
    kept: list[str] = []
    for s in nonempty:
        if any(not is_divergent(s, k, threshold=divergence_threshold) for k in kept):
            continue
        kept.append(s)
    return kept


def _tblastn_trim_each(
        named_seqs: list[tuple[str, str]],
        var_proteins_ref: str,
        padding: int,
        min_pid: float = 30.0,
        min_aa: int = 50,
) -> tuple[list[tuple[str, str]], dict[str, set[str]]]:
    """Run a single tblastn(var proteins → all candidates concat-as-multi-fasta).
    For each candidate, trim to [min hit start − padding, max hit end + padding].
    Drop candidates with no hits.

    Returns (kept_named_seqs, found_tags_by_name).
    """
    import subprocess, tempfile, os
    if not named_seqs: return [], {}
    with tempfile.TemporaryDirectory() as td:
        sfa = os.path.join(td, "all_cands.fa")
        with open(sfa, "w") as fh:
            for name, seq in named_seqs:
                fh.write(f">{name}\n{seq}\n")
        out = subprocess.run(
            ["tblastn", "-query", var_proteins_ref, "-subject", sfa,
             "-evalue", "1e-5", "-outfmt",
             "6 qseqid sseqid pident length sstart send"],
            capture_output=True, text=True).stdout
        spans_by: dict[str, list[tuple[int, int]]] = {}
        tags_by:  dict[str, set[str]]               = {}
        for ln in out.splitlines():
            f = ln.split("\t")
            if len(f) >= 6 and float(f[2]) >= min_pid and int(f[3]) >= min_aa:
                qid, sid = f[0], f[1]
                ss, se = int(f[4]), int(f[5])
                spans_by.setdefault(sid, []).append((min(ss, se), max(ss, se)))
                tags_by .setdefault(sid, set()).add(qid)
        kept: list[tuple[str, str]] = []
        for name, seq in named_seqs:
            spans = spans_by.get(name)
            if not spans: continue
            lo = max(0, min(s for s, _ in spans) - padding)
            hi = min(len(seq), max(e for _, e in spans) + padding)
            kept.append((name, seq[lo:hi]))
    return kept, tags_by


def _dedup_named(named_seqs: list[tuple[str, str]],
                  divergence_threshold: float = 0.05
                  ) -> list[tuple[str, str]]:
    """Length-desc walk, drop seqs whose k-mer Jaccard distance < threshold from
    any already-kept seq. Like _dedup_sequences but preserves the (name, seq) pairing."""
    nonempty = [(n, s) for n, s in named_seqs if s]
    nonempty.sort(key=lambda ns: len(ns[1]), reverse=True)
    kept: list[tuple[str, str]] = []
    for n, s in nonempty:
        if any(not is_divergent(s, ks, threshold=divergence_threshold)
               for _, ks in kept):
            continue
        kept.append((n, s))
    return kept


def _dedup_named_ranked(named_seqs: list[tuple[str, str]],
                         divergence_threshold: float,
                         key_fn) -> list[tuple[str, str]]:
    """Same as _dedup_named, but the walk order is determined by `key_fn`
    (smaller value = higher priority) instead of length-desc. This lets us
    rank by completeness FIRST so the more-complete representative of each
    equivalence class survives, not just the longest one. The divergence
    check itself (5% edit-distance, RC-aware) is unchanged."""
    nonempty = [(n, s) for n, s in named_seqs if s]
    nonempty.sort(key=lambda ns: key_fn(ns[0], ns[1]))
    kept: list[tuple[str, str]] = []
    for n, s in nonempty:
        if any(not is_divergent(s, ks, threshold=divergence_threshold)
               for _, ks in kept):
            continue
        kept.append((n, s))
    return kept


def _blastn_flank_presence(named_seqs: list[tuple[str, str]],
                            queries_dir: str | None
                            ) -> dict[str, tuple[bool, bool]]:
    """Run a single blastn pass with queries/flankL.fasta + flankR.fasta
    against the candidate sequences (one DB built per call). Returns
    {name: (has_flankL, has_flankR)} where each bool is True iff at least
    one hit at pid >= 85% and aln length >= 100 bp survives.

    Returns {name: (False, False)} for all candidates when queries_dir is
    missing, the flank FASTAs are absent, or BLAST exits non-zero — so
    completeness ranking degrades gracefully (every candidate looks
    "incomplete at flanks" and the next tiers do the ranking)."""
    import subprocess, tempfile, os
    presence = {n: (False, False) for n, _ in named_seqs}
    if not queries_dir or not named_seqs: return presence
    flankL_q = os.path.join(queries_dir, "flankL.fasta")
    flankR_q = os.path.join(queries_dir, "flankR.fasta")
    if not (os.path.exists(flankL_q) and os.path.exists(flankR_q)):
        return presence
    with tempfile.TemporaryDirectory() as td:
        sfa = os.path.join(td, "cands.fa")
        with open(sfa, "w") as fh:
            for n, s in named_seqs:
                fh.write(f">{n}\n{s}\n")
        db = os.path.join(td, "candb")
        r = subprocess.run(
            ["makeblastdb", "-in", sfa, "-dbtype", "nucl", "-out", db],
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0: return presence
        for q_path, idx in [(flankL_q, 0), (flankR_q, 1)]:
            r = subprocess.run(
                ["blastn", "-query", q_path, "-db", db, "-outfmt", "6",
                 "-evalue", "1e-10", "-dust", "no"],
                capture_output=True, text=True)
            if r.returncode != 0: continue
            for ln in r.stdout.splitlines():
                f = ln.split("\t")
                if len(f) < 12: continue
                sname = f[1]
                try:
                    pid = float(f[2]); aln = int(f[3])
                except ValueError: continue
                if pid >= 85.0 and aln >= 100 and sname in presence:
                    cur = list(presence[sname]); cur[idx] = True
                    presence[sname] = tuple(cur)
    return presence


def _path_mean_cov(path: list[str], provenance: dict,
                    depths: dict[str, float]) -> float:
    """Length-weighted mean DP:f: depth across a path's nodes."""
    total_bp = 0
    total_w  = 0.0
    for n in path:
        if n not in provenance: continue
        seg, s, e, _ = provenance[n]
        L = e - s
        d = depths.get(seg)
        if d is None or L <= 0: continue
        total_bp += L
        total_w  += L * d
    return (total_w / total_bp) if total_bp else 0.0


def _emit_result(res: dict, gfa_seqs: dict[str, str],
                  depths: dict[str, float] | None = None,
                  divergence_threshold: float = 0.05,
                  var_proteins_ref: str | None = None,
                  locus_padding: int = 1500,
                  expected_var_tags: set[str] | None = None,
                  genome_cov: float | None = None,
                  lo_mult: float = 0.25,
                  hi_mult: float = 2.0,
                  queries_dir: str | None = None,
                  force: bool = False) -> dict:
    """Pipe every var-bearing candidate through trim → dedup → emit.

    1. Build raw candidate sequences (full path content, no node-level trim).
    2. If `var_proteins_ref`: tblastn(var_proteins → each candidate), trim each
       to [min_hit − padding, max_hit + padding]. Candidates with no var hits
       are dropped here. This is the locus-trim step, inline with dedup.
    3. Dedup the trimmed sequences (k-mer Jaccard at `divergence_threshold`).
    4. Emit by post-dedup count `n`:
         n == 1 → allele1
         n == 2 → allele1 + allele2
         n >= 3 → chimera1, chimera2, …

    Candidate sources (the classifier populates whichever apply):
      - closed_arms + dangling_arms   (single, open/closed_bubble, complexed)
      - var_components                (separate — disjoint var subgraphs)
    """
    cls = res["class"]
    prov = res.get("_provenance", {})
    var_per_node = res.get("_var_per_node", {})
    var_nodes = {n for n, v in var_per_node.items() if v}
    # Label tokens per post-P1 node (set by directional_split). Carries
    # the per-node tag string like "HD1+HD2" or "flankL"; empty for Uvars.
    label_per_node = res.get("_label_per_node", {})
    adj_pp = res.get("_adj", {})

    def _is_flankL(n: str) -> bool:
        return "flankL" in (label_per_node.get(n, "")).split("+")

    def _is_flankR(n: str) -> bool:
        return "flankR" in (label_per_node.get(n, "")).split("+")

    def _innermost_flank_bounds(p: list[str]) -> tuple[str | None, str | None]:
        """Find innermost flank-labeled nodes bracketing the candidate path.
        Path nodes are bubble-internal (var + unlabeled); flank-labeled
        nodes sit OUTSIDE the bubble, adjacent to its boundary. So we look
        at each endpoint's exterior neighbors:
            L_bound = a flankL-labeled neighbor of path[0]   (or None)
            R_bound = a flankR-labeled neighbor of path[-1]  (or None)
        Emission stays var-trimmed; this is a metadata record only."""
        if not p: return (None, None)
        L_bound = R_bound = None
        for nb in adj_pp.get(p[0], ()):
            if _is_flankL(nb): L_bound = nb; break
        for nb in adj_pp.get(p[-1], ()):
            if _is_flankR(nb): R_bound = nb; break
        # Symmetric check (in case the path's "first" is actually flankR-side
        # due to enumeration order from a different anchor): try swapping.
        if L_bound is None and R_bound is None:
            for nb in adj_pp.get(p[-1], ()):
                if _is_flankL(nb): L_bound = nb; break
            for nb in adj_pp.get(p[0], ()):
                if _is_flankR(nb): R_bound = nb; break
        return (L_bound, R_bound)

    def _seq_full(p): return _arm_sequence_for(p, prov, gfa_seqs)

    def _trim_path_to_var(p: list[str]) -> list[str]:
        """Node-level trim: ordered path → sub-path [first var, last var]."""
        if not p: return p
        var_idx = [i for i, n in enumerate(p) if n in var_nodes]
        if not var_idx: return p
        return p[var_idx[0]:var_idx[-1] + 1]

    def _trim_component_to_var(c: list[str]) -> list[str]:
        """Node-level trim: unordered component → only var nodes."""
        return [n for n in c if n in var_nodes]

    var_per_node = res.get("_var_per_node", {})

    # The component cov-band check is GATED by the same `cov_filter` state
    # the BFS iteration used (res["_cov_filter"]). When the BFS ran with
    # cov-filter on, we apply it to var_components too; when the BFS ran
    # with cov-filter off (the lenient fallback), we don't.
    cov_filter_active = bool(res.get("_cov_filter", False))

    def _component_passes_filter(comp_nodes: list[str]) -> bool:
        """Filter for raw `var_components` emission (rarely used now that
        we recurse classify per network).
        (a) ≥ 2 DISTINCT var tags — always applies (the orphan-fragment guard).
        (b) mean depth in cov-filter band — applied only when the BFS this
            iteration also ran with the cov filter on (consistent state)."""
        tags: set[str] = set()
        for n in comp_nodes:
            t = var_per_node.get(n)
            if t: tags |= set(t)
        if len(tags) < 2:
            return False
        if cov_filter_active and genome_cov and depths:
            lo, hi = lo_mult * genome_cov, hi_mult * genome_cov
            total_bp = 0; weighted = 0.0
            for n in comp_nodes:
                if n not in prov: continue
                seg, s, e, _ = prov[n]
                L = e - s
                d = depths.get(seg)
                if d is None or L <= 0: continue
                total_bp += L; weighted += L * d
            mean_dp = (weighted / total_bp) if total_bp else 0.0
            if mean_dp and not (lo <= mean_dp <= hi):
                return False
        return True

    # ---- Per-network pool builder -------------------------------------------
    #
    # For a `separate` verdict the classifier recurses into each disjoint
    # network and returns its sub-results in res["sub_results"]. We then run
    # the full trim → dedup → emit pipeline INDEPENDENTLY per network and
    # tag each network's emitted alleles with an "n{i}_" prefix, so the
    # output for an N-network separate sample is the concatenation of N
    # per-network allele sets. Dedup does NOT cross networks (alleles from
    # disjoint graph regions stay as distinct records even if their
    # sequences happen to be similar — they represent different loci).
    #
    # For non-separate verdicts there's a single "pool" with prefix = "".

    def _process_pool(closed_arms, dangling_arms, var_components, prefix):
        """Build raw candidates → locus-trim → dedup → emit names.
        Returns dict with keys: alleles, allele_cov, allele_segments,
        extend_bounds, n_raw, n_dedup, found_tags_surviving."""
        local_raw: list[tuple[str, str]] = []
        local_path: dict[str, list[str]] = {}
        for j, p in enumerate(closed_arms):
            tp = _trim_path_to_var(list(p))
            nm = f"{prefix}cl{j}"
            local_raw.append((nm, _seq_full(tp))); local_path[nm] = tp
        for j, p in enumerate(dangling_arms):
            tp = _trim_path_to_var(list(p))
            nm = f"{prefix}da{j}"
            local_raw.append((nm, _seq_full(tp))); local_path[nm] = tp
        for j, c in enumerate(var_components or []):
            tc = _trim_component_to_var(list(c))
            if not _component_passes_filter(tc): continue
            nm = f"{prefix}vc{j}"
            local_raw.append((nm, _arm_sequence_for(tc, prov, gfa_seqs)))
            local_path[nm] = tc

        n_raw_pool = len(local_raw)

        # Locus trim — drops no-HD candidates
        local_found: dict[str, set[str]] = {}
        if var_proteins_ref:
            local_raw, local_found = _tblastn_trim_each(
                local_raw, var_proteins_ref, locus_padding)

        depths_d = depths or {}
        _cov  = lambda rn: _path_mean_cov(local_path.get(rn, []), prov, depths_d)
        _bnds = lambda rn: _innermost_flank_bounds(local_path.get(rn, []))

        # Completeness ranking for dedup: fresh blastn over the trimmed
        # candidate set against queries/flankL.fasta + queries/flankR.fasta,
        # combined with the tblastn var-tag info from _tblastn_trim_each.
        # Per candidate, compute:
        #   complete_var   = expected_var_tags ⊆ found_var_tags(this candidate)
        #   complete_locus = complete_var AND has_flankL AND has_flankR
        # Then re-sort candidates by (complete_locus DESC, complete_var DESC,
        # length DESC, diploid_cov_distance ASC) before walking the divergence
        # check — the more-complete representative survives each equivalence
        # class, instead of the longest one.
        flank_presence = _blastn_flank_presence(local_raw, queries_dir) if local_raw else {}

        def _rank_key(name: str, seq: str) -> tuple:
            fv = local_found.get(name, set())
            cv = bool(expected_var_tags) and (expected_var_tags <= fv)
            hL, hR = flank_presence.get(name, (False, False))
            cl = cv and hL and hR
            bp = len(seq)
            cov = _cov(name)
            diploid_dist = abs(cov / genome_cov - 0.5) if genome_cov else 1.0
            # sort ascending: smaller is better → negate for DESC tiers
            return (not cl, not cv, -bp, diploid_dist)

        # Dedup within this pool only (no cross-network collapse), now with
        # completeness-first ranking. is_divergent threshold unchanged.
        local_dedup = _dedup_named_ranked(local_raw, divergence_threshold, _rank_key)
        nd = len(local_dedup)

        # seg_processor.directional_split now emits clean position-indexed
        # sub-node IDs ("{parent}#1", "{parent}#2", "{parent}#N"). The
        # provenance dict carries the coords for sequence materialization;
        # the IDs themselves stay coord-free for clean display in bubble.*.
        def _segs(rn):
            out, seen = [], set()
            for sub in local_path.get(rn, ()):
                if sub not in prov: continue
                if sub in seen: continue
                seen.add(sub); out.append(sub)
            return out

        # Per-pool label scheme: n=1→allele1, n=2→allele1/2, n>=3→chimeraN
        if nd <= 1:
            if local_dedup:
                rn, rs = local_dedup[0]
                alleles_p     = [(f"{prefix}allele1", rs)]
                cov_p         = [_cov(rn)]
                bounds_p      = [_bnds(rn)]
                segs_p        = [_segs(rn)]
            else:
                alleles_p = []; cov_p = []; bounds_p = []; segs_p = []
        elif nd == 2:
            alleles_p = [(f"{prefix}allele1", local_dedup[0][1]),
                          (f"{prefix}allele2", local_dedup[1][1])]
            cov_p     = [_cov(local_dedup[0][0]), _cov(local_dedup[1][0])]
            bounds_p  = [_bnds(local_dedup[0][0]), _bnds(local_dedup[1][0])]
            segs_p    = [_segs(local_dedup[0][0]), _segs(local_dedup[1][0])]
        else:
            alleles_p = [(f"{prefix}chimera{k+1}", s)
                         for k, (_, s) in enumerate(local_dedup)]
            cov_p    = [_cov(rn) for rn, _ in local_dedup]
            bounds_p = [_bnds(rn) for rn, _ in local_dedup]
            segs_p   = [_segs(rn) for rn, _ in local_dedup]

        surviving = set()
        for rn, _ in local_dedup:
            surviving |= local_found.get(rn, set())

        return {
            "alleles": alleles_p, "allele_cov": cov_p,
            "extend_bounds": bounds_p, "allele_segments": segs_p,
            "n_raw": n_raw_pool, "n_dedup": nd,
            "found_tags_surviving": surviving,
        }

    # Build pools — one per network for `separate`, one global pool otherwise.
    sub_results = res.get("sub_results")
    if sub_results:
        pools = []
        for i, sub in enumerate(sub_results):
            pools.append(_process_pool(
                sub.get("closed_arms", []) or [],
                sub.get("dangling_arms", []) or [],
                sub.get("var_components", []) or [],
                prefix=f"n{i+1}_",
            ))
    else:
        pools = [_process_pool(
            res.get("closed_arms", []) or [],
            res.get("dangling_arms", []) or [],
            res.get("var_components", []) or [],
            prefix="",
        )]

    alleles         = [a for p in pools for a in p["alleles"]]
    allele_cov      = [c for p in pools for c in p["allele_cov"]]
    extend_bounds   = [b for p in pools for b in p["extend_bounds"]]
    allele_segments = [s for p in pools for s in p["allele_segments"]]
    n_raw           = sum(p["n_raw"]   for p in pools)
    n               = sum(p["n_dedup"] for p in pools)
    surviving_tags: set[str] = set()
    for p in pools: surviving_tags |= p["found_tags_surviving"]

    depths = depths or {}

    if expected_var_tags is not None:
        complete_var = expected_var_tags <= surviving_tags
        locus_coverage = len(expected_var_tags & surviving_tags) / max(1, len(expected_var_tags))
    else:
        complete_var = None
        locus_coverage = None

    # Decide the final verdict from the actual network sourcing of emissions:
    #   * alleles span ≥ 2 distinct networks    → "separate"  (HARD OVERRIDE,
    #                                              wins over any sub-verdict like
    #                                              open_bubble / closed_bubble /
    #                                              single — the sample really is
    #                                              split across disconnected
    #                                              graph regions)
    #   * alleles all from ONE network          → promote that network's
    #                                              sub-verdict (so a separate
    #                                              sample where only one network
    #                                              had real alleles is reported
    #                                              as the network's own shape:
    #                                              closed_bubble / open_bubble /
    #                                              complexed / single)
    #   * no networks (no sub_results, normal)  → use the classifier's top-level
    #                                              verdict as-is
    final_verdict = cls
    if sub_results:
        emit_nets = set()
        for name, _ in alleles:
            if name.startswith("n") and "_" in name:
                emit_nets.add(name.split("_", 1)[0])
        if len(emit_nets) >= 2:
            final_verdict = "separate"
        elif len(emit_nets) == 1:
            net_idx = int(next(iter(emit_nets))[1:]) - 1
            if 0 <= net_idx < len(sub_results):
                final_verdict = sub_results[net_idx].get("class", cls)
        # else: no network emitted anything — keep classifier's "separate"

    # Materialize each unique sub-node ID referenced by any emitted allele.
    # The IDs themselves are coord-free ("{parent}#N"); coords live in `prov`.
    # Returned as a small {sub_id: sub_sequence} dict so run_per_k can
    # write a per-k FASTA that graph_paths reads when assembling
    # bubble.txt / bubble.gfa / bubble.png.
    sub_seqs: dict[str, str] = {}
    referenced_ids: set[str] = set()
    for seg_list in allele_segments:
        referenced_ids.update(seg_list)
    for sid in referenced_ids:
        if sid not in prov: continue
        parent, s, e, strand = prov[sid]
        if parent not in gfa_seqs: continue
        sub = gfa_seqs[parent][s:e]
        if not sub: continue
        if strand == "-": sub = reverse_complement(sub)
        sub_seqs[sid] = sub

    return {
        "verdict": final_verdict,
        "topology": cls,                                                # raw classifier top-level (always "separate" when sub_results present)
        "complete_var": complete_var,
        "complete_locus": complete_var,                                # same notion: did we recover all var tags
        "locus_coverage": locus_coverage,
        "found_var_tags": sorted(surviving_tags),
        "n_candidates": n_raw,                                         # before dedup (pre-trim count)
        "n_after_dedup": n,
        "divergent": n >= 2,                                           # at least 2 distinct allele sequences
        "n_hops_used": res.get("_n_hops"),
        "phase": res.get("_phase"),
        "cov_filter_used": res.get("_cov_filter"),
        "alleles": alleles,
        "allele_cov": allele_cov,                                       # length-weighted DP:f: mean per allele
        "extend_bounds": extend_bounds,                                 # per allele (innermost_flankL, innermost_flankR)
        "allele_segments": allele_segments,                             # per allele: list of GFA seg IDs walked (path order, deduped)
        "subnode_seqs": sub_seqs,                                       # {sub_id: materialized sub-region sequence} for split parents
        "info": {k: v for k, v in res.items() if not k.startswith("_")
                 and k not in ("closed_arms", "dangling_arms", "var_components")},
    }


def write_fasta(alleles: list[tuple[str, str]], out_path: str,
                 sample: str = "") -> int:
    """Write the emitted alleles/chimeras to a FASTA file."""
    with open(out_path, "w") as fh:
        for name, seq in alleles:
            fh.write(f">{sample + '_' if sample else ''}{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    return len(alleles)
