from __future__ import annotations

from markerseek.models import AnalysisResult, FeatureResult
from markerseek.output import write_marker_features_tsv


def test_marker_features_tsv_has_25_columns(tmp_path) -> None:
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
    path = tmp_path / "Marker_features.tsv"

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
