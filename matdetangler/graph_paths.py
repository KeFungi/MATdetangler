"""Locate the picked alleles in the SPAdes GFA and emit annotated path views.

Outputs (per sample):
  bubble.txt   ASCII drawing of the closed bubble (or single arm for haploid)
  bubble.gfa   sub-GFA: just the bubble nodes + their L-links (loads in Bandage)
  bubble.dot   Graphviz DOT (loads in Cytoscape / Gephi / graphviz / networkx)
  bubble.tsv   edge list: node_a  node_b  allele_arm  labels   (any network library)
  bubble.png   matplotlib figure: each allele on its own row, flank-bearing nodes linked by dashed lines

Algorithm:
  1. label every GFA node by content (variable genes via tblastn, flankL/flankR via blastn, degHD
     via blastn against --known-degHD if given, repeat via blastn against --repeats if given)
  2. for each picked allele a in {a1, a2}:
       seed_nodes  = nodes whose sequence matches a at >=99% over >=1 kb (blastn)
       L_anchor    = pure-flankL node nearest a seed (BFS)
       R_anchor    = pure-flankR node nearest a seed
       path(a)     = shortest BFS path from L_anchor to R_anchor through the seeds, preferring
                     variable-gene-bearing nodes and depth within +/- cov_dev_frac of expected
  3. emit the four files
"""
from __future__ import annotations
import os, sys, subprocess, argparse, tempfile, collections, re, math
from .input_process import read_fasta
from . import blast_utils as bu

def _gfa_segments(gfa: str) -> dict[str, str]:
    d = {}
    for ln in open(gfa):
        if ln[0] == "S":
            f = ln.split("\t"); d[f[1]] = f[2]
    return d

def _gfa_adj(gfa: str) -> dict[str, set[str]]:
    a = collections.defaultdict(set)
    for ln in open(gfa):
        if ln[0] == "L":
            f = ln.split("\t"); a[f[1]].add(f[3]); a[f[3]].add(f[1])
    return a

def _gfa_depth(gfa: str) -> dict[str, float]:
    """Parse depth from KC:i:/dp:f:/RC:i: tags if present (best-effort)."""
    d = {}
    for ln in open(gfa):
        if ln[0] != "S": continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 3: continue
        name, seq = f[1], f[2]
        depth = None
        for tag in f[3:]:
            if tag.startswith("dp:f:"):
                depth = float(tag[5:]); break
            if tag.startswith("KC:i:"):
                kc = int(tag[5:]); depth = kc / max(1, len(seq)); break
            if tag.startswith("RC:i:"):
                rc = int(tag[5:]); depth = rc / max(1, len(seq)); break
        if depth is not None: d[name] = depth
    return d

def label_nodes(seqs: dict[str, str], queries_dir: str,
                known_degHD: str | None, repeats: str | None,
                tblastn_pid: float = 30.0, tblastn_aa: int = 50,
                flank_pid: float = 85.0, flank_minlen: int = 100,
                degHD_pid: float = 95.0, degHD_minlen: int = 1000) -> dict[str, str]:
    """Return {node_id: label_string} where label is a "+"-joined union of detected features."""
    if not seqs: return {}
    proteins = os.path.join(queries_dir, "variable_proteins.fasta")
    flankL = os.path.join(queries_dir, "flankL.fasta")
    flankR = os.path.join(queries_dir, "flankR.fasta")
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "n.fa")
        with open(sf, "w") as o:
            for i, s in seqs.items(): o.write(f">{i}\n{s}\n")
        db = bu.fasta_to_db(sf, t, name="db")
        feats: dict[str, list[str]] = collections.defaultdict(list)
        for r in bu.tblastn_hits(proteins, db, min_pid=tblastn_pid, min_aa=tblastn_aa):
            feats[r[0]].append(r[1])
        for nm in bu.hits_ids(bu.blastn_hits(flankL, db, min_pid=flank_pid, min_len=flank_minlen)):
            feats[nm].append("flankL")
        for nm in bu.hits_ids(bu.blastn_hits(flankR, db, min_pid=flank_pid, min_len=flank_minlen)):
            feats[nm].append("flankR")
        if known_degHD and os.path.exists(known_degHD):
            for nm in bu.hits_ids(bu.blastn_hits(known_degHD, db, min_pid=degHD_pid, min_len=degHD_minlen)):
                feats[nm].append("degHD")
        if repeats and os.path.exists(repeats):
            for nm in bu.hits_ids(bu.blastn_hits(repeats, db, min_pid=80, min_len=50)):
                feats[nm].append("repeat")
    out = {}
    for i in seqs:
        fs = list(dict.fromkeys(feats.get(i, [])))   # dedup, preserve order
        # collapse degHD-marked nodes to label "degHD"
        if "degHD" in fs: out[i] = "degHD"
        else:             out[i] = "+".join(fs)
    return out

def _bfs(adj: dict[str, set[str]], src: str, dst: set[str], blocked: set[str]) -> list[str] | None:
    if src in dst: return [src]
    prev = {src: None}; q = collections.deque([src])
    while q:
        x = q.popleft()
        for y in adj[x]:
            if y in prev or (y in blocked and y not in dst): continue
            prev[y] = x
            if y in dst:
                p = [y]
                while prev[p[-1]] is not None: p.append(prev[p[-1]])
                return list(reversed(p))
            q.append(y)
    return None

def _seed_nodes_for(allele_seq: str, gfa_seqs: dict[str, str],
                    min_pid: float = 99.0, min_len: int = 1000) -> list[tuple[str, int]]:
    """Which GFA nodes match this allele closely? Returns [(node_id, hit_length)] sorted by length desc."""
    if not allele_seq or not gfa_seqs: return []
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "n.fa")
        with open(sf, "w") as o:
            for i, s in gfa_seqs.items(): o.write(f">{i}\n{s}\n")
        db = bu.fasta_to_db(sf, t, name="db")
        qf = os.path.join(t, "q.fa")
        open(qf, "w").write(">a\n" + allele_seq + "\n")
        hits = bu.blastn_hits(qf, db, min_pid=min_pid, min_len=min_len)
    # blastn_hits row format: (sseqid, length, pident, ...) — varies; here we treat col[3] as length
    best: dict[str, int] = {}
    for r in hits:
        nid = r[0]
        try:    L = int(r[3])
        except (IndexError, ValueError): L = min_len
        if L > best.get(nid, 0): best[nid] = L
    return sorted(best.items(), key=lambda kv: -kv[1])

def _trace_path_by_alignment(allele_seq: str, gfa_seqs: dict[str, str],
                              min_pid: float = 90.0, min_len: int = 200) -> list[str]:
    """Trace the path of an allele through a GFA by ordering blastn hits along the allele.

    Process:
      1. blastn allele -> GFA segments with `qstart qend sstart send pident length sseqid`.
      2. Sort hits by qstart (position along the allele).
      3. Greedy non-overlapping cover: take hits in order, skip any whose qstart < prev qend - 50
         (50 bp overlap tolerance for k-mer-like joins).
      4. The chosen segments in order are the allele's path.

    Two distinct alleles produce distinct paths because the blastn hit set itself differs
    wherever the allele sequences differ. Shared content (e.g. flanks) produces shared
    segments; unique content (e.g. MITE in a stitched allele) produces unique segments.
    Lower thresholds than `_seed_nodes_for` (`min_pid=90`, `min_len=200`) so short, divergent
    interior segments still register."""
    if not allele_seq or not gfa_seqs: return []
    with tempfile.TemporaryDirectory() as t:
        sf = os.path.join(t, "g.fa")
        with open(sf, "w") as o:
            for i, s in gfa_seqs.items(): o.write(f">{i}\n{s}\n")
        db = bu.fasta_to_db(sf, t, name="db")
        qf = os.path.join(t, "a.fa"); open(qf, "w").write(">a\n" + allele_seq + "\n")
        # custom outfmt with positional columns
        proc = subprocess.run(
            ["blastn", "-query", qf, "-db", db, "-task", "blastn",
             "-outfmt", "6 qstart qend sseqid sstart send pident length",
             "-evalue", "1e-10", "-dust", "no", "-perc_identity", str(min_pid)],
            stdout=subprocess.PIPE, text=True)
    rows = []
    for ln in proc.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 7: continue
        try: rows.append((int(f[0]), int(f[1]), f[2], int(f[3]), int(f[4]), float(f[5]), int(f[6])))
        except ValueError: continue
    # filter by min_len
    rows = [r for r in rows if r[6] >= min_len]
    if not rows: return []
    # for each query position, keep only the best (longest) hit
    rows.sort(key=lambda r: (r[0], -r[6]))
    path: list[str] = []
    cur_end = -10**9
    for qs, qe, sid, ss, se, pid, ln in rows:
        # allow up to 50 bp overlap between consecutive hits (k-mer overlap)
        if qs < cur_end - 50: continue
        if path and path[-1] == sid: continue  # collapse consecutive same-segment hits
        path.append(sid)
        cur_end = qe
    return path

def _astr(path: list[str], labels: dict[str, str]) -> str:
    return " <-> ".join(f"{n} ({labels[n]})" if labels.get(n) else n for n in path)

def _merge_tokens(arm1: list[str], arm2: list[str], labels: dict[str, str]
                  ) -> list[tuple[str | None, str | None, str]]:
    """LCS-align two walks by node-label content. Returns ordered tokens (id_arm1, id_arm2, label),
    with id_armN = None when only the other arm has that position. Empty label is wildcard so
    unlabeled hub nodes don't block alignment."""
    if not arm1 and not arm2: return []
    if not arm2: return [(n, None, labels.get(n, "")) for n in arm1]
    if not arm1: return [(None, n, labels.get(n, "")) for n in arm2]
    la = [labels.get(n, "") for n in arm1]
    lb = [labels.get(n, "") for n in arm2]
    def matchable(x, y): return x == y or x == "" or y == ""
    n_, m_ = len(la), len(lb)
    dp = [[0] * (m_ + 1) for _ in range(n_ + 1)]
    for i in range(n_):
        for j in range(m_):
            dp[i + 1][j + 1] = (dp[i][j] + 1) if matchable(la[i], lb[j]) else max(dp[i + 1][j], dp[i][j + 1])
    i, j = n_, m_
    tokens: list[tuple[str | None, str | None, str]] = []
    while i > 0 and j > 0:
        if matchable(la[i - 1], lb[j - 1]) and dp[i][j] == dp[i - 1][j - 1] + 1:
            tokens.append((arm1[i - 1], arm2[j - 1], la[i - 1] or lb[j - 1])); i -= 1; j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            tokens.append((arm1[i - 1], None, la[i - 1])); i -= 1
        else:
            tokens.append((None, arm2[j - 1], lb[j - 1])); j -= 1
    while i > 0: tokens.append((arm1[i - 1], None, la[i - 1])); i -= 1
    while j > 0: tokens.append((None, arm2[j - 1], lb[j - 1])); j -= 1
    return list(reversed(tokens))

def _tokens_to_str(tokens: list[tuple[str | None, str | None, str]]) -> str:
    out = []
    for a, b, lab in tokens:
        tok = f"{a or '-'}/{b or '-'}"
        out.append(f"{tok} ({lab})" if lab else tok)
    return " <-> ".join(out)

def _merge_walks(arm1: list[str], arm2: list[str], labels: dict[str, str]) -> str:
    return _tokens_to_str(_merge_tokens(arm1, arm2, labels))

def _dedup_loops(path: list[str]) -> list[str]:
    """If any node appears more than once, collapse the loop:
        [..., X, ..., X, ...]  ->  [..., X, ...]   (keep prefix up to first X, then suffix after last X)
    Repeats until the path is loop-free. Preserves order and endpoints."""
    while True:
        seen: dict[str, int] = {}
        cut = None
        for i, n in enumerate(path):
            if n in seen:
                cut = (seen[n], i); break
            seen[n] = i
        if cut is None: return path
        first, last = cut
        path = path[: first + 1] + path[last + 1 :]

def _pure_flank(side: str, labels: dict[str, str]) -> list[str]:
    """Prefer nodes labeled ONLY with the flank (not co-labeled with variable genes / degHD)."""
    pure = [n for n, l in labels.items()
            if l == side or (side in l and "degHD" not in l and not any(g in l for g in ("HD1", "HD2")))]
    any_ = [n for n, l in labels.items() if side in l]
    return pure or any_

def emit_ascii(arm1: list[str], arm2: list[str] | None, labels: dict[str, str], out_path: str) -> None:
    """Write the two allele walks as separate labeled lines (mirrors the two-row layout in
    bubble.png). Same content as summary.tsv's allele1_path / allele2_path columns."""
    with open(out_path, "w") as o:
        if not arm1:
            o.write("(no allele walks found in the GFA)\n"); return
        o.write("allele1: " + _astr(arm1, labels) + "\n")
        o.write("allele2: " + _astr(arm2 or [], labels) + "\n")

def emit_sub_gfa(arm1: list[str], arm2: list[str] | None, gfa: str, out_path: str) -> None:
    """Extract S-lines for the bubble nodes + L-lines among them."""
    keep = set(arm1) | (set(arm2) if arm2 else set())
    with open(gfa) as fh, open(out_path, "w") as o:
        for ln in fh:
            if ln[0] == "S":
                f = ln.split("\t")
                if f[1] in keep: o.write(ln)
            elif ln[0] == "L":
                f = ln.split("\t")
                if f[1] in keep and f[3] in keep: o.write(ln)

def emit_dot_tsv(arm1: list[str], arm2: list[str] | None, labels: dict[str, str],
                 out_dot: str, out_tsv: str) -> None:
    nodes = list(dict.fromkeys((arm1 or []) + (arm2 or [])))
    with open(out_dot, "w") as o, open(out_tsv, "w") as t:
        o.write("graph bubble {\n")
        o.write('  graph [rankdir=LR, layout=neato];\n')
        for n in nodes:
            lab = (labels.get(n, "") or "").replace('"', "")
            color = ("lightblue" if "flank" in lab and "degHD" not in lab
                     else "lightcoral" if lab == "degHD"
                     else "lightyellow")
            o.write(f'  "{n}" [label="{n}\\n{lab}", style=filled, fillcolor={color}];\n')
        t.write("node_a\tnode_b\tarm\tlabel_a\tlabel_b\n")
        def add_edges(path, arm_name):
            for a, b in zip(path[:-1], path[1:]):
                o.write(f'  "{a}" -- "{b}" [color={"royalblue" if arm_name=="arm1" else "darkorange"}];\n')
                t.write(f"{a}\t{b}\t{arm_name}\t{labels.get(a,'')}\t{labels.get(b,'')}\n")
        if arm1: add_edges(arm1, "arm1")
        if arm2: add_edges(arm2, "arm2")
        o.write("}\n")

def emit_png_paired(arm1: list[str], arm2: list[str] | None, labels: dict[str, str],
                    out_path: str, title: str | None = None) -> None:
    """Two horizontal walks (arm1 top, arm2 bottom), each retaining its own node IDs.
    Flank-bearing nodes (labels containing 'flank') across arms are linked by a thin
    dashed gray line — a visual homology hint without forcing the two arms onto a single
    backbone. If the two walks happen to converge (the picker chose two distinct alleles
    but the GFA path-finder traced them to the same nodes), they'll just be drawn as two
    identical-looking rows; that's the honest representation."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [png] matplotlib not available; skipping bubble.png", file=sys.stderr); return
    if not arm1 and not arm2: return
    # Allow singleton mode: one arm (the picked single allele) rendered centered.
    if not arm2: arm2 = []
    if not arm1: arm1, arm2 = arm2, []
    n1, n2 = len(arm1), len(arm2)
    singleton = (n2 == 0)
    fig_h = 3.0 if singleton else 4.5
    fig, ax = plt.subplots(figsize=(max(9, 2 + 1.8 * max(n1, n2)), fig_h))
    def positions(path, y):
        n = len(path)
        if n < 2: return [(0.5, y)]
        return [(i / (n - 1), y) for i in range(n)]
    pos1 = positions(arm1, 0.0 if singleton else +0.45)
    pos2 = positions(arm2, -0.45) if not singleton else []
    def face_for(lab):
        # Variable-gene (HD) takes precedence over flank; pure-number / empty
        # labels are white. A node tagged "HD1+flankL" is colored as HD-bearing
        # because the HD content is what diverges between alleles.
        tokens = [t for t in lab.split("+") if t]
        has_flank = any(t.startswith("flank") for t in tokens)
        has_var   = any(not t.startswith("flank") for t in tokens)
        return ("#fff3b0" if has_var      # variable gene → yellow
                else "#cfe8ff" if has_flank  # flank only   → blue
                else "#ffffff")              # pure number  → white
    # Cross-arm homology lines:
    #   1. For EVERY GFA segment ID that appears in both arms (any label),
    #      draw a dashed line linking the two occurrences. Same ID = same
    #      segment in the GFA = definite homology.
    #   2. As a flank-end fallback: if NO flank-only node (label is flankL or
    #      flankR with NO variable-gene tag) is shared between the arms for
    #      that flank type, draw ONE extra line connecting the outermost
    #      flank-only node of each arm. Outermost: leftmost flankL (arm
    #      entry), rightmost flankR (arm exit).
    # The pure-pairwise variant (linking every flank-bearing node to every
    # other in the opposite arm) turned bubbles with many flank-tagged
    # segments into a hairball; this rule preserves homology evidence
    # without spurious cross-links.
    def _is_flank_only(lab: str, want: str) -> bool:
        toks = [t for t in lab.split("+") if t]
        return bool(toks) and all(t.startswith("flank") for t in toks) and want in toks
    if not singleton:
        ids1: dict[str, int] = {}
        for i, n in enumerate(arm1): ids1.setdefault(n, i)
        ids2: dict[str, int] = {}
        for j, n in enumerate(arm2): ids2.setdefault(n, j)
        shared_ids = set(ids1) & set(ids2)
        for node in shared_ids:
            i, j = ids1[node], ids2[node]
            ax.plot([pos1[i][0], pos2[j][0]], [pos1[i][1], pos2[j][1]],
                    color="#9aa0a6", lw=0.9, ls="--", zorder=0)
        # outermost flank-only fallback for whichever flank has no shared-ID anchor
        OUTER = {"flankL": "first", "flankR": "last"}
        for want in ("flankL", "flankR"):
            idx1 = [i for i, n in enumerate(arm1) if _is_flank_only(labels.get(n, ""), want)]
            idx2 = [j for j, n in enumerate(arm2) if _is_flank_only(labels.get(n, ""), want)]
            if not idx1 or not idx2: continue
            if any(arm1[i] in shared_ids for i in idx1): continue
            if any(arm2[j] in shared_ids for j in idx2): continue
            pick = (lambda xs: xs[0]) if OUTER[want] == "first" else (lambda xs: xs[-1])
            i = pick(idx1); j = pick(idx2)
            ax.plot([pos1[i][0], pos2[j][0]], [pos1[i][1], pos2[j][1]],
                    color="#9aa0a6", lw=0.9, ls="--", zorder=0)
    # arm 1 (blue) — backbone + boxes
    for (xa, ya), (xb, yb) in zip(pos1[:-1], pos1[1:]):
        ax.plot([xa, xb], [ya, yb], color="#3b6db8", lw=1.8, zorder=1)
    for (x, y), n in zip(pos1, arm1):
        la = labels.get(n, "")
        ax.annotate(f"{n}\n{la}" if la else n, (x, y), ha="center", va="center", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.30", fc=face_for(la), ec="0.4"), zorder=2)
    # arm 2 (orange) — skipped in singleton mode
    if not singleton:
        for (xa, ya), (xb, yb) in zip(pos2[:-1], pos2[1:]):
            ax.plot([xa, xb], [ya, yb], color="#d97a3a", lw=1.8, zorder=1)
        for (x, y), n in zip(pos2, arm2):
            la = labels.get(n, "")
            ax.annotate(f"{n}\n{la}" if la else n, (x, y), ha="center", va="center", fontsize=8,
                        bbox=dict(boxstyle="round,pad=0.30", fc=face_for(la), ec="0.4"), zorder=2)
    ax.set_xlim(-0.05, 1.05)
    if singleton: ax.set_ylim(-0.5, 0.5)
    else: ax.set_ylim(-0.9, 0.9)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ("top", "right", "bottom", "left"): ax.spines[sp].set_visible(False)
    if title: ax.set_title(title, fontsize=10)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)

def _trace_allele_in_gfa(seq: str, gfa: str, queries_dir: str,
                         known_degHD: str | None, repeats: str | None,
                         recorded_path: list[str] | None = None
                         ) -> tuple[list[str], dict[str, str], dict[str, str]]:
    """Return (path, labels, seqs) for one allele in one GFA.

    If `recorded_path` is given (i.e. segment_alleles wrote the walk into picks.tsv), use it
    verbatim — that's the authoritative path from the candidate-generation step. Otherwise
    fall back to blastn-ordering against the GFA segments (re-derives the path from sequence,
    used only when picks.tsv was written by a legacy contig-based picker)."""
    seqs = _gfa_segments(gfa)
    labels = label_nodes(seqs, queries_dir, known_degHD, repeats)
    if recorded_path is not None:
        # filter to segments that actually exist in this GFA (defensive — different k = different segs)
        path = [n for n in recorded_path if n in seqs]
        return path, labels, seqs
    return _trace_path_by_alignment(seq, seqs), labels, seqs

def run(sample: str, gfa: str, primary_alleles_fa: str, queries_dir: str, outdir: str,
        known_degHD: str | None = None, repeats: str | None = None,
        per_allele_gfas: dict[str, str] | None = None,
        recorded_paths: dict[str, list[str]] | None = None) -> dict:
    """Trace each picked allele's walk through a GFA and emit ASCII/DOT/PNG views.

    `per_allele_gfas`: optional mapping allele_name -> GFA path. When given (and the allele's
    GFA exists), each allele is traced in its own native k's GFA — node IDs and label
    placements differ between the two walks, which is precisely how the user sees that the
    alleles came from different assemblies.

    `recorded_paths`: optional mapping allele_name -> list of segment IDs (the path written by
    `segment_alleles.py` into `picks.tsv:segments`). When given, used verbatim — no
    re-derivation needed. Strands (`+`/`-`) in segment IDs are stripped before lookup.
    """
    os.makedirs(outdir, exist_ok=True)
    alleles = read_fasta(primary_alleles_fa)
    if not alleles:
        print(f"[graph_paths] {sample}: no picked alleles; nothing to draw")
        return {}
    # Two-stage optimization for the labeling cost:
    # (1) PER-GFA CACHE: multiple alleles from the same k share the same GFA;
    #     parse it once and reuse.
    # (2) RECORDED-PATH SUBSET: if recorded_paths covers every allele in this
    #     GFA, we know upfront which segments will be drawn — only THOSE need
    #     content labels. Drops labeling work from O(|GFA|=100k-300k segs) to
    #     O(|union of paths|=~tens of segs). When any allele lacks a recorded
    #     path, we fall back to labeling the full GFA (alignment-based trace
    #     needs to know labels for every candidate segment).
    needed_segs_per_gfa: dict[str, set[str] | None] = {}
    for nm in alleles:
        g = (per_allele_gfas or {}).get(nm) or gfa
        rec = (recorded_paths or {}).get(nm)
        if rec is None:
            needed_segs_per_gfa[g] = None   # mark "must label all"
        else:
            if needed_segs_per_gfa.get(g) is None and g in needed_segs_per_gfa:
                continue  # already marked all-labeling for this GFA
            needed_segs_per_gfa.setdefault(g, set()).update(rec)
    gfa_cache: dict[str, tuple[dict[str, str], dict[str, str]]] = {}
    def _seqs_labels_for(g: str):
        if g in gfa_cache: return gfa_cache[g]
        ss = _gfa_segments(g)
        needed = needed_segs_per_gfa.get(g)
        if needed is None:
            ll = label_nodes(ss, queries_dir, known_degHD, repeats)
        else:
            subset = {sid: ss[sid] for sid in needed if sid in ss}
            ll = label_nodes(subset, queries_dir, known_degHD, repeats)
        gfa_cache[g] = (ss, ll)
        return ss, ll
    # per-allele trace
    arm_data = []  # list of (allele_name, path, labels, seqs, gfa_used)
    for nm, seq in alleles.items():
        g = (per_allele_gfas or {}).get(nm) or gfa
        rec = (recorded_paths or {}).get(nm)
        seqs, labels = _seqs_labels_for(g)
        if rec is not None:
            path = [n for n in rec if n in seqs]
        else:
            path = _trace_path_by_alignment(seq, seqs)
        arm_data.append((nm, path, labels, seqs, g))
    # merge labels + sequences across arms so the renderers don't need to know about
    # multiple graphs (node IDs are k-prefixed when graphs differ, to avoid collisions)
    same_gfa = len(set(a[4] for a in arm_data)) <= 1
    merged_labels: dict[str, str] = {}
    merged_seqs: dict[str, str] = {}
    arm_paths: list[list[str]] = []
    for nm, path, labels, seqs, g in arm_data:
        if same_gfa:
            tag = ""
        else:
            # tag node IDs with the GFA's parent dir basename so allele1 and allele2 walks have
            # distinct IDs even if they're drawn together
            tag = "K" + os.path.basename(os.path.dirname(g)).lstrip("kK") + ":"
        tagged_path = [f"{tag}{n}" for n in path]
        for n, tagn in zip(path, tagged_path):
            if tagn not in merged_labels: merged_labels[tagn] = labels.get(n, "")
            if tagn not in merged_seqs:   merged_seqs[tagn]   = seqs.get(n, "")
        arm_paths.append(tagged_path)
    arm1 = arm_paths[0] if arm_paths else []
    arm2 = arm_paths[1] if len(arm_paths) > 1 else None
    emit_ascii(arm1, arm2, merged_labels, os.path.join(outdir, "bubble.txt"))
    # sub-GFA: write per-GFA S-lines and L-lines with the same K-tag prefix used everywhere else.
    # Cross-graph links don't exist, so L-lines are only within each source GFA.
    used_gfas = list(dict.fromkeys(a[4] for a in arm_data))
    def _tag_for(g):
        return "" if same_gfa else "K" + os.path.basename(os.path.dirname(g)).lstrip("kK") + ":"
    with open(os.path.join(outdir, "bubble.gfa"), "w") as o:
        for g in used_gfas:
            tag = _tag_for(g)
            tagged_ids_for_this_g = {n for a in arm_paths for n in a if (not tag) or n.startswith(tag)}
            raw_ids = {n[len(tag):] if tag else n for n in tagged_ids_for_this_g}
            for ln in open(g):
                if ln[0] == "S":
                    f = ln.split("\t")
                    if f[1] in raw_ids:
                        o.write("\t".join([f[0], f"{tag}{f[1]}"] + f[2:]))
                elif ln[0] == "L":
                    f = ln.split("\t")
                    if f[1] in raw_ids and f[3] in raw_ids:
                        f[1] = f"{tag}{f[1]}"; f[3] = f"{tag}{f[3]}"
                        o.write("\t".join(f))
    emit_dot_tsv(arm1, arm2, merged_labels,
                  os.path.join(outdir, "bubble.dot"),
                  os.path.join(outdir, "bubble.tsv"))
    if same_gfa:
        _png_title = f"{sample}: locus walks (allele1=blue, allele2=orange)"
    else:
        ks_label = ' vs '.join('K' + os.path.basename(os.path.dirname(g)).lstrip('kK') for g in used_gfas)
        _png_title = f"{sample}: per-allele walks ({ks_label})"
    emit_png_paired(arm1, arm2, merged_labels,
                    os.path.join(outdir, "bubble.png"), title=_png_title)
    arm1_str = _astr(arm1, merged_labels)
    arm2_str = _astr(arm2 or [], merged_labels)
    import json
    with open(os.path.join(outdir, "_graph.json"), "w") as o:
        json.dump({"arm1_str": arm1_str or "-", "arm2_str": arm2_str or "-"}, o)
    # for the caller's convenience, expose the merged structures under the legacy names
    labels = merged_labels
    print(f"[graph_paths] {sample}: arm1={len(arm1)} nodes, arm2={len(arm2) if arm2 else 0} nodes -> {outdir}/bubble.*")
    return {"arm1": arm1, "arm2": arm2, "labels": labels,
            "arm1_str": arm1_str, "arm2_str": arm2_str}

def _paths_from_cand_ann(cand_ann_tsv: str, spades_dir: str, primary_alleles_fa: str
                          ) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Same shape as _per_allele_gfas_and_paths but reads bubble_alleles.ann.tsv
    (the step-3 emission). 10-col schema: name, len, k, segments, var_genes,
    flankL, flankR, cov, is_degHD, has_repeat. We use cols 1 (name), 3 (k),
    4 (segments). Used when picks.tsv doesn't exist (--skip-pick mode).
    """
    out_gfas: dict[str, str] = {}
    out_paths: dict[str, list[str]] = {}
    if not (cand_ann_tsv and os.path.exists(cand_ann_tsv)):
        return out_gfas, out_paths
    from .paths import spades_k_paths
    # read names actually present in the primary fasta
    names = set()
    for ln in open(primary_alleles_fa):
        if ln.startswith(">"): names.add(ln[1:].split()[0])
    for ln in open(cand_ann_tsv):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4: continue
        nm, _len, k, segs = f[0], f[1], f[2], f[3]
        if nm not in names: continue
        if spades_dir:
            _, g = spades_k_paths(spades_dir, k)
            if g: out_gfas[nm] = g
        if segs and segs != "-":
            out_paths[nm] = [s.rstrip("+-") for s in segs.split(",")]
    return out_gfas, out_paths


def _per_allele_gfas_and_paths(picks_tsv: str, spades_dir: str, primary_alleles_fa: str
                                ) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Read picks.tsv and resolve each allele's (a) source-k GFA and (b) recorded segment path.

    Returns (gfas, paths) both keyed by the FASTA name of the allele.
    `paths[allele_name]` is the list of segment IDs (strands stripped) from the picks `segments`
    column, written by `segment_alleles.py`. If the column is missing or "-", that allele's
    entry is absent from `paths` and the caller falls back to alignment.
    """
    if not (picks_tsv and os.path.exists(picks_tsv)): return {}, {}
    from .paths import spades_k_paths
    short_to_k: dict[str, str] = {}
    short_to_segs: dict[str, list[str]] = {}
    with open(picks_tsv) as fh:
        hdr = next(fh, "").rstrip("\n").split("\t")
        try:
            i_allele = hdr.index("allele"); i_k = hdr.index("k")
        except ValueError:
            return {}, {}
        i_segs = hdr.index("segments") if "segments" in hdr else -1
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) <= max(i_allele, i_k): continue
            short_to_k[f[i_allele]] = f[i_k]
            if i_segs >= 0 and len(f) > i_segs and f[i_segs] not in ("", "-"):
                # strand-stripped segment IDs (e.g. "975421+" -> "975421")
                short_to_segs[f[i_allele]] = [s.rstrip("+-") for s in f[i_segs].split(",")]
    fasta_names = []
    for ln in open(primary_alleles_fa):
        if ln.startswith(">"): fasta_names.append(ln[1:].split()[0])
    gfas: dict[str, str] = {}
    paths: dict[str, list[str]] = {}
    for nm in fasta_names:
        for short, k in short_to_k.items():
            if nm.endswith("_" + short) or nm == short or nm.endswith(short):
                _, gfa = spades_k_paths(spades_dir, k)
                if gfa: gfas[nm] = gfa
                if short in short_to_segs: paths[nm] = short_to_segs[short]
                break
    return gfas, paths

# legacy alias for any external caller
def _per_allele_gfas_from_picks(picks_tsv: str, spades_dir: str,
                                primary_alleles_fa: str) -> dict[str, str]:
    return _per_allele_gfas_and_paths(picks_tsv, spades_dir, primary_alleles_fa)[0]

def _cli(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--sample", required=True)
    p.add_argument("--gfa", required=True, help="fallback GFA used when --picks-tsv/--spades-dir aren't supplied or an allele's source-k GFA can't be resolved")
    p.add_argument("--primary-alleles", required=True)
    p.add_argument("--queries-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--known-degHD", default=None)
    p.add_argument("--repeats", default=None)
    p.add_argument("--picks-tsv", default=None, help="if given with --spades-dir, each allele is traced in its native-k GFA")
    p.add_argument("--cand-ann", default=None,
                   help="bubble_alleles.ann.tsv from step 3; fallback source of per-allele "
                        "(source-k, segment-walk) when picks.tsv is unavailable (--skip-pick mode). "
                        "Lets step 7 use the recorded path verbatim — no whole-GFA labeling needed.")
    p.add_argument("--spades-dir", default=None)
    a = p.parse_args(argv)
    if a.picks_tsv and a.spades_dir and os.path.exists(a.picks_tsv):
        pag, rec = _per_allele_gfas_and_paths(a.picks_tsv, a.spades_dir, a.primary_alleles)
    elif a.cand_ann and a.spades_dir:
        pag, rec = _paths_from_cand_ann(a.cand_ann, a.spades_dir, a.primary_alleles)
    else:
        pag, rec = ({}, {})
    run(a.sample, a.gfa, a.primary_alleles, a.queries_dir, a.outdir,
        a.known_degHD, a.repeats, per_allele_gfas=pag, recorded_paths=rec)

if __name__ == "__main__":
    _cli()
