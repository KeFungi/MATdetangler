"""Reconstruct allele DNA sequences from classifier arm paths.

Pipeline:
  1. classify() returns each arm as a path of post-P1 node IDs.
  2. seg_processor.directional_split() supplies provenance:
       provenance[node] = (original_seg_id, start, end, strand)
     For unsplit segments, the range is the whole segment.
     For sub-segments of a split, ranges partition the original segment
     at midpoints between consecutive BLAST hits.
  3. The GFA gives the original segment's stored-strand DNA.

This module joins all that into per-arm sequences:
  - Walk each arm's path
  - For each node, pull the sub-region from the original segment
  - Apply RC if strand is "-"
  - Concatenate
  - Optionally skip joint/fork nodes (nodes shared between arms)
"""
from __future__ import annotations


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def parse_gfa_sequences(gfa_path: str) -> dict[str, str]:
    """Return {seg_id: stored-strand DNA}. Parses S-lines from a GFA file."""
    out: dict[str, str] = {}
    with open(gfa_path) as fh:
        for line in fh:
            if not line.startswith("S\t"): continue
            parts = line.rstrip("\n").split("\t", 3)
            if len(parts) >= 3:
                out[parts[1]] = parts[2]
    return out


def find_joint_nodes(arms: list[list[str]]) -> set[str]:
    """A joint/fork node = a node that appears in 2+ arms.

    For closed bubbles where arms share L-anchor / R-anchor (the divergence
    and convergence points), those anchors are joint nodes by this rule.
    """
    count: dict[str, int] = {}
    for arm in arms:
        for n in arm:
            count[n] = count.get(n, 0) + 1
    return {n for n, c in count.items() if c >= 2}


def reconstruct_arm_sequence(
        arm_path: list[str],
        provenance: dict[str, tuple[str, int, int, str]],
        gfa_sequences: dict[str, str],
        skip_nodes: set[str] | None = None,
) -> str:
    """Build the DNA for one arm by walking its path.

    Each node contributes the sub-region of its original segment, RC'd if
    the sub-region's stored strand is "-". Nodes in `skip_nodes` (e.g. joints)
    are omitted from the concatenation.
    """
    skip = skip_nodes or set()
    parts: list[str] = []
    for n in arm_path:
        if n in skip: continue
        if n not in provenance: continue
        seg, start, end, strand = provenance[n]
        seq = gfa_sequences.get(seg, "")
        if not seq: continue
        sub = seq[start:end]
        if strand == "-":
            sub = reverse_complement(sub)
        parts.append(sub)
    return "".join(parts)


def reconstruct_arms(
        arms: list[list[str]],
        provenance: dict[str, tuple[str, int, int, str]],
        gfa_sequences: dict[str, str],
        skip_joints: bool = True,
) -> list[dict]:
    """Top-level: reconstruct every arm's sequence.

    Returns a list of dicts, one per arm:
        {
            "path":      [node_ids],
            "trimmed":   [node_ids with joints removed],
            "sequence":  "ACGT...",
            "length":    int,
            "n_joints_skipped": int,
        }
    """
    joints = find_joint_nodes(arms) if skip_joints else set()
    out = []
    for arm in arms:
        trimmed = [n for n in arm if n not in joints]
        seq = reconstruct_arm_sequence(arm, provenance, gfa_sequences, joints)
        out.append({
            "path": list(arm),
            "trimmed": trimmed,
            "sequence": seq,
            "length": len(seq),
            "n_joints_skipped": len(arm) - len(trimmed),
        })
    return out
