from __future__ import annotations

from markerseek.diagnostics import divergence_stats


def test_single_sample_per_species_returns_none_for_intra_dependent_metrics() -> None:
    aligned_block = {
        "s1": "A" * 100,
        "s2": "C" * 100,
        "s3": "G" * 100,
        "s4": "T" * 100,
    }
    species_map = {
        "s1": "species_1",
        "s2": "species_2",
        "s3": "species_3",
        "s4": "species_4",
    }

    stats = divergence_stats(aligned_block, species_map)

    assert stats["intraspecific"] is None
    assert stats["nearest_neighbor_discrimination"] is None
    assert stats["barcoding_gap"] is None
    assert stats["misclassification_risk"] is None
