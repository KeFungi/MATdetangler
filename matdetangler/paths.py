"""Resolve per-k contigs + GFA paths under a SPAdes output directory.

Recommended layout (per-k SPAdes runs):

  spades_dir/k21/{contigs.fasta, assembly_graph_after_simplification.gfa}
  spades_dir/k33/...
  spades_dir/k55/...

Produced by looping `spades.py -k <k> -o spades_dir/k<k>/` once per k.
This is the only layout that exposes a usable GFA for every k in the sweep.

Also accepted:

  (a) `spades_dir/K21/...` (uppercase) — only relevant if a tool other than the SPAdes
      multi-k pipeline produced the per-k subdirs; SPAdes' own `-k 21,33,55` mode strips
      K21/K33 to internal binary state.

  (b) Top-level final-k fallback: if no per-k subdir is found for the requested k but the
      spades_dir itself contains an `assembly_graph_after_simplification.gfa` + `contigs.fasta`,
      we return those (with `top_level=True`) — this is what `spades.py -k 21,33,55 -o spades_dir/`
      leaves behind, and it represents the final k only.

Returns (contigs_fasta, gfa) for the requested k, or (None, None) if nothing found.
The caller can detect the fallback by comparing the returned paths against `spades_dir`.
"""
from __future__ import annotations
import os

CONTIG_NAMES = ("contigs.fasta", "final_contigs.fasta", "before_rr.fasta")
GFA_NAME = "assembly_graph_after_simplification.gfa"

def spades_k_paths(spades_dir: str, k: str) -> tuple[str | None, str | None]:
    knum = k.lstrip("kK")
    for d in (os.path.join(spades_dir, f"k{knum}"),
              os.path.join(spades_dir, f"K{knum}")):
        if not os.path.isdir(d):
            continue
        gfa = os.path.join(d, GFA_NAME)
        if not os.path.exists(gfa):
            continue
        for cf in CONTIG_NAMES:
            p = os.path.join(d, cf)
            if os.path.exists(p):
                return p, gfa
    # final-k-only fallback: top-level GFA + contigs (what `spades.py -k 21,33,55` leaves)
    top_gfa = os.path.join(spades_dir, GFA_NAME)
    if os.path.exists(top_gfa):
        for cf in CONTIG_NAMES:
            p = os.path.join(spades_dir, cf)
            if os.path.exists(p):
                return p, top_gfa
    return None, None
