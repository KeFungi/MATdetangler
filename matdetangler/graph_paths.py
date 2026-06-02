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


# Sub-node ID convention emitted by per_k_caller._segs when P1 split a GFA
# segment into multiple sub-nodes: "<parent>#<N>" (1-based index, no coords).
# Lets bubble.txt / bubble.gfa / bubble.png distinguish HD1-bearing vs
# HD2-bearing sub-regions of the same parent segment instead of collapsing
# them to a single ambiguous node label. Bare IDs (no "#") are un-split
# parents. The actual materialized sub-sequence for each sub-node ID is
# stored by run_per_k in <sample_dir>/<k>/subnode_seqs.fasta, which the
# resolver reads to populate the seqs dict when needed.
def _parent_of(sid: str) -> str:
    return sid.split("#", 1)[0] if "#" in sid else sid

def _load_subnode_seqs(k_dir: str) -> dict[str, str]:
    """Read <k_dir>/subnode_seqs.fasta if present. Returns {} otherwise."""
    p = os.path.join(k_dir, "subnode_seqs.fasta")
    if not os.path.exists(p): return {}
    out: dict[str, str] = {}
    cur = None; buf: list[str] = []
    with open(p) as fh:
        for ln in fh:
            ln = ln.rstrip()
            if ln.startswith(">"):
                if cur is not None: out[cur] = "".join(buf)
                cur = ln[1:].split()[0]; buf = []
            else:
                buf.append(ln)
        if cur is not None: out[cur] = "".join(buf)
    return out

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

def emit_sub_gfa(arm1: list[str], arm2: list[str] | None, gfa: str, out_path: str,
                 sub_seqs: dict[str, str] | None = None) -> None:
    """Extract S-lines for the bubble nodes + L-lines among them.

    If `sub_seqs` is given, S-lines are written directly from it (one per
    sub-node ID, supports decorated `parent#start-end±` IDs). L-lines are
    synthesized from arm adjacency since the source GFA only has parent-level
    edges. Otherwise we fall back to copying S-lines verbatim from the GFA
    (legacy bare-ID case) and copying parent-to-parent L-lines.
    """
    keep = set(arm1) | (set(arm2) if arm2 else set())
    with open(out_path, "w") as o:
        if sub_seqs:
            # S-lines: write each retained node with its sub-region sequence
            for nm in keep:
                seq = sub_seqs.get(nm, "")
                if not seq: continue
                o.write(f"S\t{nm}\t{seq}\n")
            # L-lines: synthesized from arm-walk adjacency. Each arm
            # contributes (n_i, n_{i+1}) pairs; strand is unknown post-trim
            # so we mark + by convention. Duplicates collapsed.
            seen_links: set[tuple[str, str]] = set()
            for arm in (arm1, arm2 or []):
                for a, b in zip(arm[:-1], arm[1:]):
                    if (a, b) in seen_links or (b, a) in seen_links: continue
                    seen_links.add((a, b))
                    o.write(f"L\t{a}\t+\t{b}\t+\t0M\n")
            return
        # Fallback: copy from source GFA (bare-ID only)
        with open(gfa) as fh:
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

def emit_png_paired(arms: list[list[str]], arm_names: list[str], labels: dict[str, str],
                    out_path: str, title: str | None = None) -> None:
    """Render N walks as N stacked horizontal rows (top to bottom). Each row
    keeps its own node IDs. For N=1 the single row is centered (singleton
    mode). For N=2 the original blue/orange paired layout is preserved.
    For N>=3 a cycling palette is used and cross-row homology lines are
    drawn between every PAIR of consecutive rows that share a node ID.
    A title (if given) sits above the figure — bash wrapper uses this to
    surface the bubble verdict + allele count.

    Backward compatibility shim: legacy callers passed (arm1, arm2, labels, out, title).
    The wrapper at module bottom translates that to the new signature.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [png] matplotlib not available; skipping bubble.png", file=sys.stderr); return
    arms = [a for a in arms if a]
    if not arms: return
    # Cap visible rows at 4 — above that the plot stops being a useful
    # picture (overlapping boxes, illegible labels). Extras are dropped from
    # the PNG only; bubble.txt / picks.tsv still carry the full set.
    truncated = max(0, len(arms) - 4)
    arm_names = list(arm_names) + [""] * max(0, len(arms) - len(arm_names))
    if truncated:
        arms = arms[:4]; arm_names = arm_names[:4]
    N = len(arms)
    max_n = max(len(a) for a in arms)
    # 1 inch per row + header padding, no compression at higher N
    fig_h = 2.0 + 1.0 * N
    fig_w = max(9, 2 + 1.8 * max_n)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    # Row layout: 1.0 data-unit spacing per row, centered around y=0.
    # Top row = (N-1)/2, bottom row = -(N-1)/2. Plenty of headroom so the
    # outermost rows don't clip on annotations.
    ys = [(N - 1) / 2.0 - i for i in range(N)]
    def positions(path, y):
        n = len(path)
        if n < 2: return [(0.5, y)]
        return [(i / (n - 1), y) for i in range(n)]
    row_pos = [positions(arm, ys[i]) for i, arm in enumerate(arms)]
    palette = ["#3b6db8", "#d97a3a", "#5e9c64", "#a463b5", "#c7503f",
               "#1f8a8a", "#8a6f2e", "#5b5f96"]
    row_colors = [palette[i % len(palette)] for i in range(N)]

    def face_for(lab):
        tokens = [t for t in lab.split("+") if t]
        has_flank = any(t.startswith("flank") for t in tokens)
        has_var   = any(not t.startswith("flank") for t in tokens)
        return ("#fff3b0" if has_var
                else "#cfe8ff" if has_flank
                else "#ffffff")

    def _is_flank_only(lab: str, want: str) -> bool:
        toks = [t for t in lab.split("+") if t]
        return bool(toks) and all(t.startswith("flank") for t in toks) and want in toks

    # Cross-row homology lines — drawn between EVERY consecutive pair of
    # rows. Same logic as the old 2-arm version, just applied repeatedly.
    if N >= 2:
        OUTER = {"flankL": "first", "flankR": "last"}
        for k in range(N - 1):
            a, b = arms[k], arms[k + 1]
            pa, pb = row_pos[k], row_pos[k + 1]
            ids_a: dict[str, int] = {}
            for i, n in enumerate(a): ids_a.setdefault(n, i)
            ids_b: dict[str, int] = {}
            for j, n in enumerate(b): ids_b.setdefault(n, j)
            shared = set(ids_a) & set(ids_b)
            for node in shared:
                i, j = ids_a[node], ids_b[node]
                ax.plot([pa[i][0], pb[j][0]], [pa[i][1], pb[j][1]],
                        color="#9aa0a6", lw=0.9, ls="--", zorder=0)
            for want in ("flankL", "flankR"):
                ia = [i for i, n in enumerate(a) if _is_flank_only(labels.get(n, ""), want)]
                ib = [j for j, n in enumerate(b) if _is_flank_only(labels.get(n, ""), want)]
                if not ia or not ib: continue
                if any(a[i] in shared for i in ia): continue
                if any(b[j] in shared for j in ib): continue
                pick = (lambda xs: xs[0]) if OUTER[want] == "first" else (lambda xs: xs[-1])
                i = pick(ia); j = pick(ib)
                ax.plot([pa[i][0], pb[j][0]], [pa[i][1], pb[j][1]],
                        color="#9aa0a6", lw=0.9, ls="--", zorder=0)

    # Draw each row: backbone + boxes + row label
    for i, (arm, name, color, pos) in enumerate(zip(arms, arm_names, row_colors, row_pos)):
        for (xa, ya), (xb, yb) in zip(pos[:-1], pos[1:]):
            ax.plot([xa, xb], [ya, yb], color=color, lw=1.8, zorder=1)
        for (x, y), n in zip(pos, arm):
            la = labels.get(n, "")
            ax.annotate(f"{n}\n{la}" if la else n, (x, y), ha="center", va="center", fontsize=8,
                        bbox=dict(boxstyle="round,pad=0.30", fc=face_for(la), ec="0.4"), zorder=2)
        # Row label (e.g. "chimera3") above its segment row, NOT on the
        # same y as the segment boxes — keeps the label out of the box
        # area when arms are dense and avoids overlap with cross-row
        # homology lines. Centered along the row so it's visible regardless
        # of how many boxes are on this row.
        ax.text(-0.05, ys[i] + 0.32, name, ha="left", va="bottom", fontsize=9,
                color=color, weight="bold")

    ax.set_xlim(-0.08, 1.05)
    pad_top, pad_bot = 0.9, 0.6   # extra headroom for the per-row labels
    ax.set_ylim(min(ys) - pad_bot, max(ys) + pad_top)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ("top", "right", "bottom", "left"): ax.spines[sp].set_visible(False)
    full_title = title or ""
    if truncated:
        # Footnote-style suffix when extras were dropped
        full_title = (full_title + f"   [showing 4 of {N + truncated}]").strip()
    if full_title: ax.set_title(full_title, fontsize=11, weight="bold")
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
        recorded_paths: dict[str, list[str]] | None = None,
        subnode_seqs_per_gfa: dict[str, dict[str, str]] | None = None,
        bubble_type: str | None = None) -> dict:
    """Trace each picked allele's walk through a GFA and emit ASCII/DOT/PNG views.

    `per_allele_gfas`: optional mapping allele_name -> GFA path. When given (and the allele's
    GFA exists), each allele is traced in its own native k's GFA — node IDs and label
    placements differ between the two walks, which is precisely how the user sees that the
    alleles came from different assemblies.

    `recorded_paths`: optional mapping allele_name -> list of segment IDs (the path written by
    `segment_alleles.py` into `picks.tsv:segments`). When given, used verbatim — no
    re-derivation needed. Strands (`+`/`-`) in segment IDs are stripped before lookup.

    `subnode_seqs_per_gfa`: optional mapping GFA path -> {sub_id: sub_sequence}. Provides
    materialized sub-region sequences for "{parent}#N" decorated IDs. When given, these
    sub-IDs are added to the per-GFA seqs dict so labels reflect per-sub-region content.

    `bubble_type`: optional verdict label (closed_bubble, open_bubble, complexed, ...) for
    the PNG title. When given together with N alleles, the figure title shows
    "{sample}: {bubble_type} (N alleles)".
    """
    os.makedirs(outdir, exist_ok=True)
    alleles = read_fasta(primary_alleles_fa)
    if not alleles:
        print(f"[graph_paths] {sample}: no picked alleles; nothing to draw")
        return {}
    sub_seqs_per_gfa = subnode_seqs_per_gfa or {}
    # Two-stage optimization for the labeling cost:
    # (1) PER-GFA CACHE: multiple alleles from the same k share the same GFA;
    #     parse it once and reuse.
    # (2) RECORDED-PATH SUBSET: if recorded_paths covers every allele in this
    #     GFA, we know upfront which segments will be drawn — only THOSE need
    #     content labels. Drops labeling work from O(|GFA|=100k-300k segs) to
    #     O(|union of paths|=~tens of segs). When any allele lacks a recorded
    #     path, we fall back to labeling the full GFA (alignment-based trace
    #     needs to know labels for every candidate segment).
    # Per GFA: collect the set of recorded sub-node IDs needed. Decorated
    # sub-node IDs ("{parent}#N") have their materialized sequence in the
    # per-k subnode_seqs.fasta; bare IDs come from the GFA directly. We
    # augment the seqs dict so label_nodes() sees the actual sub-region
    # content (so an HD1-only sub-node labels as "HD1" rather than the
    # parent's joint "HD1+HD2"), and so downstream lookups for these IDs
    # succeed.
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
        parent_seqs = _gfa_segments(g)
        # Sub-node seqs live next to the GFA's per-k output dir, NOT next
        # to the GFA file itself. Per-allele resolver passes (per_allele_gfas)
        # the SPAdes graph path, but the materialized sub-node FASTA was
        # written by run_per_k to <results>/<sample>/<k>/subnode_seqs.fasta.
        # We can't infer that from the spades path alone, so the resolver
        # below attaches it via a side dict (sub_seqs_per_gfa).
        sub_seqs = sub_seqs_per_gfa.get(g, {})
        needed = needed_segs_per_gfa.get(g)
        if needed is None:
            ss = dict(parent_seqs); ss.update(sub_seqs)
            ll = label_nodes(ss, queries_dir, known_degHD, repeats)
        else:
            ss = {}
            for sid in needed:
                if "#" in sid:
                    if sid in sub_seqs: ss[sid] = sub_seqs[sid]
                elif sid in parent_seqs:
                    ss[sid] = parent_seqs[sid]
            ll = label_nodes(ss, queries_dir, known_degHD, repeats)
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
    # Orient each row so its var-tag fingerprint matches the first row's
    # (the leftmost var-bearing node on every row carries the same gene tag
    # — matching the File-S1-style canonical convention). A row whose
    # var-tag order is the EXACT REVERSE of row 0's gets walked backwards.
    # Rows with disjoint or ambiguous fingerprints are left alone.
    def _vartag_seq(arm: list[str]) -> list[str]:
        out: list[str] = []
        for n in arm:
            lab = merged_labels.get(n, "")
            toks = [t for t in lab.split("+") if t and not t.startswith("flank")]
            if toks and (not out or out[-1] != toks[0]):
                out.append(toks[0])
        return out
    if arm_paths:
        canon = _vartag_seq(arm_paths[0])
        if canon:
            for i in range(1, len(arm_paths)):
                vs = _vartag_seq(arm_paths[i])
                if vs and vs == list(reversed(canon)):
                    arm_paths[i] = list(reversed(arm_paths[i]))
    arm1 = arm_paths[0] if arm_paths else []
    arm2 = arm_paths[1] if len(arm_paths) > 1 else None
    emit_ascii(arm1, arm2, merged_labels, os.path.join(outdir, "bubble.txt"))
    # sub-GFA: write per-GFA S-lines and L-lines with the same K-tag prefix used everywhere else.
    # Cross-graph links don't exist, so L-lines are only within each source GFA.
    used_gfas = list(dict.fromkeys(a[4] for a in arm_data))
    def _tag_for(g):
        return "" if same_gfa else "K" + os.path.basename(os.path.dirname(g)).lstrip("kK") + ":"
    # If any allele's path contains decorated sub-node IDs (anything with
    # "#" — emitted by per_k_caller when P1 split a parent into multiple
    # sub-nodes), the source GFA can't supply those S-lines verbatim; emit
    # them from the merged sub-seqs dict instead, and synthesize L-lines
    # from arm-walk adjacency (parent-to-parent L-lines from the GFA no
    # longer have the right node names anyway).
    any_subnode = any("#" in n for a in arm_paths for n in a)
    with open(os.path.join(outdir, "bubble.gfa"), "w") as o:
        if any_subnode:
            for n in dict.fromkeys(n for a in arm_paths for n in a):
                seq = merged_seqs.get(n, "")
                if seq: o.write(f"S\t{n}\t{seq}\n")
            seen_links: set[tuple[str, str]] = set()
            for arm in arm_paths:
                for a, b in zip(arm[:-1], arm[1:]):
                    if (a, b) in seen_links or (b, a) in seen_links: continue
                    seen_links.add((a, b))
                    o.write(f"L\t{a}\t+\t{b}\t+\t0M\n")
        else:
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
    # PNG title: "{sample}: {bubble_type} (N alleles)" on top — surfaces
    # the classifier verdict + the allele count right in the figure.
    n_drawn = len(arm_paths)
    allele_names_list = list(alleles.keys())
    if bubble_type:
        _png_title = f"{sample}: {bubble_type} ({n_drawn} allele{'s' if n_drawn != 1 else ''})"
    elif same_gfa:
        _png_title = f"{sample}: {n_drawn} allele walk{'s' if n_drawn != 1 else ''}"
    else:
        ks_label = ' vs '.join('K' + os.path.basename(os.path.dirname(g)).lstrip('kK') for g in used_gfas)
        _png_title = f"{sample}: per-allele walks ({ks_label})"
    emit_png_paired(arm_paths, allele_names_list, merged_labels,
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


def _picks_summary_bubble_type(picks_tsv: str) -> str | None:
    """Read bubble_type from the sibling picks_summary.tsv (single row).
    Returns None when the file is absent or unreadable."""
    if not picks_tsv: return None
    summary = picks_tsv.replace("picks.tsv", "picks_summary.tsv")
    if not os.path.exists(summary): return None
    with open(summary) as fh:
        hdr = next(fh, "").rstrip("\n").split("\t")
        if "bubble_type" not in hdr: return None
        i = hdr.index("bubble_type")
        row = next(fh, "").rstrip("\n").split("\t")
        if len(row) <= i: return None
        v = row[i].strip()
        return v if v and v not in ("-", "NA", "NO_RESULT") else None


def _subnode_seqs_per_gfa(per_allele_gfas: dict[str, str]) -> dict[str, dict[str, str]]:
    """For each unique GFA path in per_allele_gfas, load <gfa_dir>/subnode_seqs.fasta
    if present. Returns {gfa_path: {sub_id: sub_seq}}. The directory holding
    the GFA is also where run_per_k writes the subnode FASTA, so we look there."""
    out: dict[str, dict[str, str]] = {}
    seen_gfas = set(per_allele_gfas.values())
    # Each GFA path is .../spades/<sample>/<k>/assembly_graph_after_simplification.gfa
    # but subnode_seqs.fasta is in the RESULTS sample/k dir, not next to the GFA.
    # The caller is expected to populate this via _subnode_seqs_per_gfa_from_picks
    # — this default helper just covers the (rare) case where the FASTA sits
    # next to the GFA itself.
    for g in seen_gfas:
        k_dir = os.path.dirname(g)
        seqs = _load_subnode_seqs(k_dir)
        if seqs: out[g] = seqs
    return out


def _subnode_seqs_per_gfa_from_picks(picks_tsv: str, per_allele_gfas: dict[str, str]
                                      ) -> dict[str, dict[str, str]]:
    """Resolve subnode_seqs.fasta from the RESULTS dir: each allele's source
    k is in picks.tsv, and run_per_k wrote subnode_seqs.fasta to
    <picks_dir>/<k>/subnode_seqs.fasta. We aggregate across alleles, keyed
    by the per-allele GFA path so the resolver can look up directly."""
    if not picks_tsv or not os.path.exists(picks_tsv): return {}
    results_sample_dir = os.path.dirname(picks_tsv)
    # parse picks.tsv: allele -> k
    allele_to_k: dict[str, str] = {}
    with open(picks_tsv) as fh:
        hdr = next(fh, "").rstrip("\n").split("\t")
        try:
            i_a = hdr.index("allele"); i_k = hdr.index("k")
        except ValueError:
            return {}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) > max(i_a, i_k):
                allele_to_k[f[i_a]] = f[i_k]
    # for each GFA in per_allele_gfas, derive the k and load that k's
    # subnode_seqs.fasta from the results dir
    out: dict[str, dict[str, str]] = {}
    for nm, g in per_allele_gfas.items():
        # nm in primary fasta is like "<SAMPLE>_k53_allele1"; the picks
        # entries use short names. We need to find an entry that matches.
        k = None
        for short, kk in allele_to_k.items():
            if nm.endswith("_" + short) or nm == short or nm.endswith(short):
                k = kk; break
        if not k: continue
        k_dir = os.path.join(results_sample_dir, k)
        seqs = _load_subnode_seqs(k_dir)
        if seqs: out[g] = seqs
    return out


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
                # Strip trailing strand marks ONLY on bare parent IDs
                # (e.g. "975421+" -> "975421"). Decorated sub-node IDs of
                # the form "parent#start-end[+|-]" keep their full form —
                # the trailing strand is part of the sub-region encoding.
                def _strip_bare(s: str) -> str:
                    return s.rstrip("+-") if "#" not in s else s
                short_to_segs[f[i_allele]] = [_strip_bare(s) for s in f[i_segs].split(",")]
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
        sub_per_gfa = _subnode_seqs_per_gfa_from_picks(a.picks_tsv, pag)
        btype = _picks_summary_bubble_type(a.picks_tsv)
    elif a.cand_ann and a.spades_dir:
        pag, rec = _paths_from_cand_ann(a.cand_ann, a.spades_dir, a.primary_alleles)
        sub_per_gfa = _subnode_seqs_per_gfa(pag)
        btype = None
    else:
        pag, rec, sub_per_gfa, btype = ({}, {}, {}, None)
    run(a.sample, a.gfa, a.primary_alleles, a.queries_dir, a.outdir,
        a.known_degHD, a.repeats, per_allele_gfas=pag, recorded_paths=rec,
        subnode_seqs_per_gfa=sub_per_gfa, bubble_type=btype)

if __name__ == "__main__":
    _cli()
