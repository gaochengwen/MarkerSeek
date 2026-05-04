"""Plotly-ready visualisation helpers for feature detail pages."""

from __future__ import annotations

from collections import defaultdict
from html import escape
from itertools import combinations

import networkx as nx
import numpy as np

from .diagnostics import _p_distance, assign_haplotypes


CANONICAL_OR_GAP = frozenset("ACGT-")


def compute_haplotype_network(aligned_block: dict[str, str], species_map: dict[str, str]) -> dict:
    """Return a minimum-spanning haplotype network for an aligned feature block."""

    if not aligned_block:
        return {"nodes": [], "edges": []}

    sample_to_haplotype = assign_haplotypes(aligned_block)
    haplotype_samples: dict[str, list[str]] = defaultdict(list)
    haplotype_sequences: dict[str, str] = {}
    for sample_name, haplotype_id in sample_to_haplotype.items():
        haplotype_samples[haplotype_id].append(sample_name)
        haplotype_sequences.setdefault(haplotype_id, _normalise_haplotype(aligned_block[sample_name]))

    graph = nx.Graph()
    for haplotype_id in sorted(haplotype_samples):
        graph.add_node(haplotype_id)
    for left, right in combinations(sorted(haplotype_sequences), 2):
        graph.add_edge(
            left,
            right,
            weight=_hamming_distance(haplotype_sequences[left], haplotype_sequences[right]),
        )

    mst = nx.minimum_spanning_tree(graph, weight="weight") if graph.number_of_nodes() else graph
    positions = nx.spring_layout(mst, seed=42, iterations=200) if mst.number_of_nodes() else {}

    nodes = []
    for haplotype_id in sorted(haplotype_samples):
        samples = haplotype_samples[haplotype_id]
        x, y = positions.get(haplotype_id, (0.0, 0.0))
        nodes.append(
            {
                "id": haplotype_id,
                "frequency": len(samples),
                "species": sorted({species_map.get(sample_name, "Unknown") or "Unknown" for sample_name in samples}),
                "x": float(x),
                "y": float(y),
            }
        )

    edges = [
        {"source": str(left), "target": str(right), "weight": int(data.get("weight", 0))}
        for left, right, data in sorted(mst.edges(data=True), key=lambda item: (str(item[0]), str(item[1])))
    ]
    return {"nodes": nodes, "edges": edges}


def compute_species_pca(aligned_block: dict[str, str], species_map: dict[str, str]) -> dict:
    """Project samples onto PC1/PC2 from a centred pairwise p-distance matrix."""

    samples = list(aligned_block)
    if not samples:
        return {"samples": [], "explained_variance": [0.0, 0.0]}

    distance_matrix = np.zeros((len(samples), len(samples)), dtype=float)
    for left_index, right_index in combinations(range(len(samples)), 2):
        distance, _valid_sites = _p_distance(aligned_block[samples[left_index]], aligned_block[samples[right_index]])
        distance_value = 0.0 if distance is None else float(distance)
        distance_matrix[left_index, right_index] = distance_value
        distance_matrix[right_index, left_index] = distance_value

    centred = (
        distance_matrix
        - np.mean(distance_matrix, axis=0)
        - np.mean(distance_matrix, axis=1)[:, None]
        + float(np.mean(distance_matrix))
    )
    u, singular_values, _vt = np.linalg.svd(centred)

    pc1_values = u[:, 0] * singular_values[0] if singular_values.size >= 1 else np.zeros(len(samples))
    pc2_values = u[:, 1] * singular_values[1] if singular_values.size >= 2 else np.zeros(len(samples))
    total_variance = float(np.sum(singular_values))
    explained = [
        float(singular_values[0] / total_variance) if total_variance and singular_values.size >= 1 else 0.0,
        float(singular_values[1] / total_variance) if total_variance and singular_values.size >= 2 else 0.0,
    ]

    return {
        "samples": [
            {
                "sample_name": sample_name,
                "species": species_map.get(sample_name, "Unknown") or "Unknown",
                "pc1": float(pc1_values[index]),
                "pc2": float(pc2_values[index]),
            }
            for index, sample_name in enumerate(samples)
        ],
        "explained_variance": explained,
    }


def render_alignment_svg_block(aligned_block: dict[str, str], *, columns_per_row: int = 60) -> str:
    """Render a server-side HTML alignment viewer fragment."""

    if columns_per_row <= 0:
        raise ValueError("columns_per_row must be a positive integer.")
    if not aligned_block:
        return '<div class="alignment-viewer"></div>'

    max_length = max(len(sequence) for sequence in aligned_block.values())
    parts = ['<div class="alignment-viewer">']
    for row_start in range(0, max_length, columns_per_row):
        row_end = min(row_start + columns_per_row, max_length)
        parts.append(_render_position_ruler(row_start, row_end, columns_per_row))
        for sample_name, sequence in aligned_block.items():
            segment = sequence[row_start:row_end]
            sample_label = escape(sample_name[:24].ljust(24))
            rendered_segment = "".join(_render_base(base) for base in segment)
            parts.append(
                f'<pre class="alignment-row" data-sequence="{escape(segment.upper(), quote=True)}">'
                f"{sample_label} {rendered_segment}</pre>"
            )
    parts.append("</div>")
    return "\n".join(parts)


def _normalise_haplotype(sequence: str) -> str:
    return "".join(base if base in CANONICAL_OR_GAP else "N" for base in sequence.upper())


def _hamming_distance(left: str, right: str) -> int:
    return sum(1 for left_base, right_base in zip(left, right, strict=False) if left_base != right_base)


def _render_base(base: str) -> str:
    base = base.upper()
    css_class = f"b-{base}" if base in CANONICAL_OR_GAP else "b-N"
    if base == "-":
        css_class = "b-gap"
    label = base if base else " "
    return f'<span class="{css_class}">{escape(label)}</span>'


def _render_position_ruler(row_start: int, row_end: int, columns_per_row: int) -> str:
    ruler = [" "] * columns_per_row
    for position in range(row_start + 1, row_end + 1):
        if position == row_start + 1 or (position - 1) % columns_per_row == 0:
            label = str(position)
            offset = position - row_start - 1
            for index, character in enumerate(label):
                if offset + index < len(ruler):
                    ruler[offset + index] = character
    return f'<pre class="alignment-ruler">{" " * 25}{"".join(ruler[: row_end - row_start])}</pre>'
