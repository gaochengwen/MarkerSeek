from __future__ import annotations

from markerseek.analysis import compute_position_pi, infer_regions
from markerseek.models import AnalysisResult, AnnotatedInterval, WindowResult
from markerseek.plotting import select_label_windows


def test_compute_position_pi_ignores_missing_symbols():
    pi_values, valid_mask = compute_position_pi(
        [
            "ACGT",
            "AGGT",
            "A-GT",
            "NNGT",
        ]
    )

    assert valid_mask == [1, 1, 1, 1]
    assert pi_values[0] == 0.0
    assert pi_values[1] == 1.0
    assert pi_values[2] == 0.0
    assert pi_values[3] == 0.0


def test_compute_position_pi_matches_pairwise_diversity():
    pi_values, valid_mask = compute_position_pi(
        [
            "A",
            "A",
            "G",
            "G",
            "G",
        ]
    )

    assert valid_mask == [1]
    assert pi_values[0] == 0.6


def test_infer_regions_expands_from_marker_pairs():
    irb = "ATGCGTAC" * 160
    ira = reverse_complement(irb)
    sequence = ("A" * 400) + irb + ("C" * 320) + ira + ("G" * 300)

    irb_start = 400
    irb_end = irb_start + len(irb)
    ira_start = irb_end + 320
    ira_end = ira_start + len(ira)

    markers = [
        AnnotatedInterval("rrn16S_a", "rrn16S", "rRNA", "rrn16S", irb_start + 60, irb_start + 200),
        AnnotatedInterval("rrn16S_b", "rrn16S", "rRNA", "rrn16S", ira_start + 60, ira_start + 200),
        AnnotatedInterval("rpl2_a", "rpl2", "gene", "rpl2", irb_start + 260, irb_start + 420),
        AnnotatedInterval("rpl2_b", "rpl2", "gene", "rpl2", ira_start + 260, ira_start + 420),
        AnnotatedInterval("ndhB_a", "ndhB", "gene", "ndhB", irb_start + 520, irb_start + 700),
        AnnotatedInterval("ndhB_b", "ndhB", "gene", "ndhB", ira_start + 520, ira_start + 700),
        AnnotatedInterval("ycf2_a", "ycf2", "gene", "ycf2", irb_start + 860, irb_start + 1040),
        AnnotatedInterval("ycf2_b", "ycf2", "gene", "ycf2", ira_start + 860, ira_start + 1040),
    ]

    regions = infer_regions(sequence, markers)
    region_map = {region.name: region for region in regions}

    assert region_map["IRb"].start == irb_start
    assert region_map["IRb"].end == irb_end
    assert region_map["IRa"].start == ira_start
    assert region_map["IRa"].end == ira_end
    assert region_map["SSC"].length(len(sequence)) == 320
    assert region_map["LSC"].length(len(sequence)) == 700


def test_select_label_windows_keeps_cluster_peaks_and_spacing():
    result = AnalysisResult(
        reference_name="ref",
        genome_length=2000,
        sample_count=3,
        regions=[],
        position_pi=[],
        windows=[
            WindowResult("W0001", 0, 200, 100, 0.020, 200, "LSC", "peak_a_left", is_hotspot=True, hotspot_rank=4),
            WindowResult("W0002", 100, 300, 200, 0.050, 200, "LSC", "peak_a", is_hotspot=True, hotspot_rank=2),
            WindowResult("W0003", 200, 400, 300, 0.030, 200, "LSC", "peak_a_right", is_hotspot=True, hotspot_rank=3),
            WindowResult("W0004", 420, 620, 520, 0.025, 200, "LSC", "peak_b_left", is_hotspot=True, hotspot_rank=5),
            WindowResult("W0005", 460, 660, 560, 0.040, 200, "LSC", "peak_b", is_hotspot=True, hotspot_rank=1),
            WindowResult("W0006", 1180, 1380, 1280, 0.060, 200, "SSC", "peak_c", is_hotspot=True, hotspot_rank=6),
            WindowResult("W0007", 1500, 1700, 1600, 0.010, 200, "IRa", "background", is_hotspot=False),
        ],
        features=[],
    )

    labels = select_label_windows(
        result,
        label_mode="peak-only",
        label_max=10,
        label_min_distance_bp=400,
    )

    assert [window.label_name for window in labels] == ["peak_a", "peak_c"]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
