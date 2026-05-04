from __future__ import annotations

from markerseek.diagnostics import alignment_reliability_block, assign_haplotypes, divergence_stats
from markerseek.scoring import compute_score


def test_haplotype_assignment_groups_identical_blocks() -> None:
    haplotypes = assign_haplotypes({"s1": "ACGT", "s2": "ACGT", "s3": "ACGA", "s4": "TCGT"})

    assert len(set(haplotypes.values())) == 3
    assert haplotypes["s1"] == haplotypes["s2"]


def test_alignment_reliability_penalises_gappy_columns() -> None:
    pure_block = {f"s{index}": "A" * 50 for index in range(4)}
    gappy_block = {
        "s1": "-" * 30 + "A" * 20,
        "s2": "-" * 30 + "A" * 20,
        "s3": "A" * 50,
        "s4": "A" * 50,
    }

    assert alignment_reliability_block(pure_block) >= 0.95
    assert alignment_reliability_block(gappy_block) < 0.5


def test_markerseek_score_extremes() -> None:
    best_metrics = {
        "pi": 0.05,
        "variable_site_density": 0.10,
        "indel_density": 0.05,
        "flanking_conservation_min": 1.0,
        "missing_ambig_ratio": 0.0,
        "alignment_reliability": 1.0,
        "species_resolution": 1.0,
        "barcoding_gap": 0.05,
        "nearest_neighbor_discrimination": 1.0,
        "length_suitability": 1.0,
    }

    assert compute_score(best_metrics) >= 95
    assert compute_score({}) <= 5


def test_divergence_barcoding_gap_positive_for_well_separated_species() -> None:
    aligned_block = {
        "sp1_a": "A" * 100,
        "sp1_b": "A" * 100,
        "sp2_a": "T" * 5 + "A" * 95,
        "sp2_b": "T" * 5 + "A" * 95,
    }
    species_map = {
        "sp1_a": "species_1",
        "sp1_b": "species_1",
        "sp2_a": "species_2",
        "sp2_b": "species_2",
    }

    stats = divergence_stats(aligned_block, species_map)

    assert stats["barcoding_gap"] is not None
    assert stats["barcoding_gap"] > 0
    assert stats["nearest_neighbor_discrimination"] == 1.0
