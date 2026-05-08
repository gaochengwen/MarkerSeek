from __future__ import annotations

from pathlib import Path

import pytest

from markerseek import plotting
from markerseek.models import AnalysisResult, FeatureResult, RegionSegment, WindowResult


def _plot_result(regions: list[RegionSegment]) -> AnalysisResult:
    sequence = "ACGT" * 20
    return AnalysisResult(
        reference_name="reference",
        genome_length=len(sequence),
        sample_count=2,
        regions=regions,
        position_pi=[0.0] * len(sequence),
        windows=[
            WindowResult("W0001", 1, 40, 20, 0.01, 40, regions[0].name, "geneA", is_hotspot=True),
            WindowResult("W0002", 41, 80, 60, 0.02, 40, regions[0].name, "geneA", is_hotspot=False),
        ],
        features=[
            FeatureResult(
                feature_id="geneA",
                feature_type="gene",
                parent_gene="geneA",
                label_name="geneA",
                start=5,
                end=45,
                strand=1,
                length_bp=41,
                region=regions[0].name,
                pi=0.01,
                valid_sites=41,
                spans_origin=False,
            )
        ],
        aligned_sequences={
            "reference": sequence,
            "sample": sequence[:20] + "T" + sequence[21:],
        },
        sample_order=["reference", "sample"],
    )


def test_unpartitioned_genome_region_is_not_an_informative_plot_partition() -> None:
    genome_result = _plot_result([RegionSegment("Genome", 0, 80)])
    partitioned_result = _plot_result([RegionSegment("LSC", 0, 80)])

    assert not plotting._has_informative_region_partition(genome_result)
    assert plotting._has_informative_region_partition(partitioned_result)


def test_no_ir_fallback_plots_skip_genome_region_bar(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    result = _plot_result([RegionSegment("Genome", 0, 80)])

    def fail_region_bar(*_args, **_kwargs) -> None:
        raise AssertionError("The no-IR Genome fallback should not draw a region bar.")

    monkeypatch.setattr(plotting, "PLOT_DPI", 80)
    monkeypatch.setattr(plotting, "_plot_region_bar", fail_region_bar)

    plotting.plot_pi_figure(result, tmp_path)
    plotting.plot_similarity_figure(result, tmp_path)

    assert (tmp_path / "pi_plot.png").exists()
    assert (tmp_path / "similarity_plot.png").exists()
