from __future__ import annotations

from markerseek.models import AnalysisResult, FeatureResult, RegionSegment, WindowResult
from markerseek.output import write_analysis_outputs, write_marker_features_tsv, write_primers_tsv


def test_candidate_marker_features_tsv_has_25_columns(tmp_path) -> None:
    result = AnalysisResult(
        reference_name="reference",
        genome_length=100,
        sample_count=2,
        regions=[],
        position_pi=[],
        windows=[],
        features=[
            FeatureResult(
                feature_id="IGS_001",
                feature_type="igs",
                parent_gene="geneA|geneB",
                label_name="geneA-geneB",
                start=1,
                end=100,
                strand=0,
                length_bp=100,
                region="LSC",
                pi=0.01,
                valid_sites=100,
                spans_origin=False,
            )
        ],
    )
    path = tmp_path / "candidate_marker_features.tsv"

    write_marker_features_tsv(path, result.features)

    header = path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert header == [
        "feature_id",
        "feature_type",
        "parent_gene",
        "label_name",
        "start",
        "end",
        "strand",
        "length_bp",
        "region",
        "pi",
        "variable_sites",
        "indel_sites",
        "conserved_left_bp",
        "conserved_right_bp",
        "primer_available",
        "species_resolution",
        "unique_haplotype_count",
        "species_specific_haplotype_ratio",
        "interspecific_divergence",
        "intraspecific_divergence",
        "nearest_neighbor_discrimination",
        "barcoding_gap",
        "misclassification_risk",
        "alignment_reliability",
        "markerseek_score",
    ]
    assert len(header) == 25


def test_candidate_marker_features_filename(tmp_path) -> None:
    result = AnalysisResult(
        reference_name="reference",
        genome_length=100,
        sample_count=2,
        regions=[RegionSegment("LSC", 0, 100)],
        position_pi=[0.0] * 100,
        windows=[
            WindowResult(
                window_id="W0001",
                start=1,
                end=100,
                midpoint=50,
                pi=0.01,
                valid_sites=100,
                region="LSC",
                label_name="geneA-geneB",
            )
        ],
        features=[
            FeatureResult(
                feature_id="IGS_001",
                feature_type="igs",
                parent_gene="geneA|geneB",
                label_name="geneA-geneB",
                start=1,
                end=100,
                strand=0,
                length_bp=100,
                region="LSC",
                pi=0.01,
                valid_sites=100,
                spans_origin=False,
            )
        ],
    )

    write_analysis_outputs(
        result,
        tmp_path,
        hotspot_mode="top-percent",
        hotspot_value=3.0,
        label_mode="peak-only",
        label_max=None,
        label_min_distance_bp=0,
        similarity_window=200,
        similarity_step=60,
        similarity_floor=0.5,
        include_similarity_plot=False,
        primer_design=False,
    )

    assert (tmp_path / "candidate_marker_features.tsv").exists()
    assert not (tmp_path / "Marker_features.tsv").exists()


def test_primers_tsv_has_expected_columns_and_writes_when_empty(tmp_path) -> None:
    result = AnalysisResult(
        reference_name="reference",
        genome_length=100,
        sample_count=2,
        regions=[],
        position_pi=[],
        windows=[],
        features=[],
    )
    path = tmp_path / "primers.tsv"

    write_primers_tsv(path, result.primers)

    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    assert len(header) == 19
    assert header[:2] == ["primer_id", "label_name"]
    assert "amplicon_size" in header
    assert "amplicon_min_len" not in header
    assert "amplicon_variable_sites" not in header
    assert lines == ["\t".join(header)]
