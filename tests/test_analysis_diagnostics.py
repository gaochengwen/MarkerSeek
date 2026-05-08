from __future__ import annotations

import pytest

from markerseek.analysis import (
    MarkerSeekError,
    infer_regions,
    infer_regions_or_unpartitioned,
    promote_label_window_candidates,
    promote_top_score_candidates,
    region_for_position,
)
from markerseek.diagnostics import alignment_reliability_block, assign_haplotypes, divergence_stats
from markerseek.models import AnnotatedInterval, FeatureResult, RegionSegment, WindowResult
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


def test_promote_top_score_candidates_keeps_hotspots_and_adds_top_five() -> None:
    features = [
        FeatureResult(
            feature_id=f"F{i}",
            feature_type="gene",
            parent_gene=f"gene{i}",
            label_name=f"gene{i}",
            start=i * 10 + 1,
            end=i * 10 + 9,
            strand=1,
            length_bp=9,
            region="LSC",
            pi=0.01,
            valid_sites=9,
            spans_origin=False,
            is_hotspot=(i == 0),
            hotspot_rank=(1 if i == 0 else None),
            markerseek_score=float(i),
        )
        for i in range(8)
    ]

    promote_top_score_candidates(features, extra_count=5)

    selected = {feature.feature_id for feature in features if feature.is_hotspot}
    assert selected == {"F0", "F3", "F4", "F5", "F6", "F7"}


def test_promote_label_window_candidates_uses_pi_plot_labels_before_score_topups() -> None:
    features = [
        FeatureResult(
            feature_id="ycf1",
            feature_type="gene",
            parent_gene="ycf1",
            label_name="ycf1",
            start=101,
            end=160,
            strand=1,
            length_bp=60,
            region="SSC",
            pi=0.01,
            valid_sites=60,
            spans_origin=False,
            markerseek_score=10.0,
        ),
        FeatureResult(
            feature_id="ndhF",
            feature_type="gene",
            parent_gene="ndhF",
            label_name="ndhF",
            start=301,
            end=360,
            strand=1,
            length_bp=60,
            region="SSC",
            pi=0.02,
            valid_sites=60,
            spans_origin=False,
            markerseek_score=20.0,
        ),
        FeatureResult(
            feature_id="top_score",
            feature_type="igs",
            parent_gene="a|b",
            label_name="a-b",
            start=501,
            end=560,
            strand=0,
            length_bp=60,
            region="LSC",
            pi=0.03,
            valid_sites=60,
            spans_origin=False,
            markerseek_score=99.0,
        ),
    ]
    windows = [
        WindowResult("W1", 105, 154, 130, 0.04, 50, "SSC", "ycf1", is_hotspot=True),
        WindowResult("W2", 305, 354, 330, 0.05, 50, "SSC", "ndhF", is_hotspot=True),
    ]

    promote_label_window_candidates(features, windows, genome_length=1000, label_mode="all")
    promote_top_score_candidates(features, extra_count=1)

    assert [(feature.feature_id, feature.hotspot_rank) for feature in features if feature.is_hotspot] == [
        ("ycf1", 1),
        ("ndhF", 2),
        ("top_score", 3),
    ]


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


def test_ir_detection_failure_falls_back_to_unpartitioned_genome() -> None:
    sequence = "ACGT" * 250

    with pytest.raises(MarkerSeekError, match="Automatic IR detection failed"):
        infer_regions(sequence, [])

    regions = infer_regions_or_unpartitioned(sequence, [])

    assert regions == [RegionSegment("Genome", 0, len(sequence))]
    assert region_for_position(123, regions, len(sequence)) == "Genome"


def test_ir_detection_fallback_only_catches_missing_ir_candidates() -> None:
    def interval(gene: str, start: int, end: int) -> AnnotatedInterval:
        return AnnotatedInterval(
            feature_id=f"{gene}_{start}",
            label=gene,
            feature_type="gene",
            parent_gene=gene,
            start=start,
            end=end,
        )

    sequence = "".join("ACGT"[(index * index + index * 7) % 4] for index in range(6000))
    features = []
    for index, gene in enumerate(("gene1", "gene2", "gene3", "gene4")):
        features.append(interval(gene, 100 + index * 200, 150 + index * 200))
        features.append(interval(gene, 3000 + index * 200, 3050 + index * 200))

    with pytest.raises(MarkerSeekError, match="implausibly short"):
        infer_regions_or_unpartitioned(sequence, features)
