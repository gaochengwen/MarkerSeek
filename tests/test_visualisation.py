from __future__ import annotations

import re

from markerseek.visualisation import (
    compute_haplotype_network,
    compute_species_pca,
    render_alignment_svg_block,
)


def test_haplotype_network_minimum_spanning_tree_has_n_minus_one_edges() -> None:
    aligned_block = {
        "sample_a": "A" * 40,
        "sample_b": "C" + ("A" * 39),
        "sample_c": "CC" + ("A" * 38),
        "sample_d": "CCC" + ("A" * 37),
    }
    species_map = {sample_name: "species" for sample_name in aligned_block}

    network = compute_haplotype_network(aligned_block, species_map)

    assert len(network["nodes"]) == 4
    assert len(network["edges"]) == 3
    assert sum(edge["weight"] for edge in network["edges"]) == 3


def test_species_pca_separates_two_clusters() -> None:
    aligned_block = {
        "alpha_1": "A" * 100,
        "alpha_2": "A" * 100,
        "alpha_3": "A" * 100,
        "beta_1": ("C" * 5) + ("A" * 95),
        "beta_2": ("C" * 5) + ("A" * 95),
        "beta_3": ("C" * 5) + ("A" * 95),
    }
    species_map = {sample_name: sample_name.split("_")[0] for sample_name in aligned_block}

    pca = compute_species_pca(aligned_block, species_map)
    alpha_pc1 = {sample["pc1"] > 0 for sample in pca["samples"] if sample["species"] == "alpha"}
    beta_pc1 = {sample["pc1"] > 0 for sample in pca["samples"] if sample["species"] == "beta"}

    assert len(alpha_pc1) == 1
    assert len(beta_pc1) == 1
    assert alpha_pc1 != beta_pc1
    assert 0.0 <= pca["explained_variance"][0] <= 1.0


def test_render_alignment_svg_wraps_at_100_columns() -> None:
    sequence = ("A" * 100) + ("C" * 100) + ("GTACGTACGT")
    html = render_alignment_svg_block({"sample": sequence})

    rows = re.findall(r'<pre class="alignment-row" data-sequence="([^"]+)"', html)

    assert len(rows) == 3
    assert rows[0] == "A" * 100
    assert rows[1] == "C" * 100
    assert rows[2] == "GTACGTACGT"
    assert html.count('class="alignment-ruler"') == 3
    assert "101" in html
    assert '<span class="b-A">A</span>' in html
