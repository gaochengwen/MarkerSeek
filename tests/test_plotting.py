from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from markerseek.models import WindowResult
from markerseek.plotting import (
    _annotate_peak_windows,
    _bbox_overlap_area,
    _peak_text_bbox,
    _pi_y_axis_limits,
    compute_similarity_tracks,
)


def test_peak_annotations_avoid_overlapping_text_boxes():
    fig, ax = plt.subplots(figsize=(4.0, 2.5), dpi=180)
    ax.set_xlim(1, 2000)
    ax.set_ylim(0, 0.08)
    windows = [
        WindowResult("W0001", 100, 700, 700, 0.070, 600, "LSC", "peak_alpha", True, 1),
        WindowResult("W0002", 300, 900, 760, 0.068, 600, "LSC", "peak_beta", True, 2),
        WindowResult("W0003", 500, 1100, 820, 0.066, 600, "LSC", "peak_gamma", True, 3),
        WindowResult("W0004", 700, 1300, 880, 0.064, 600, "LSC", "peak_delta", True, 4),
    ]
    ax.plot([window.midpoint for window in windows], [window.pi for window in windows])

    annotations = _annotate_peak_windows(ax, windows)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [
        _peak_text_bbox(ax, window, *annotation.get_position(), renderer).expanded(1.02, 1.08)
        for window, annotation in zip(windows, annotations)
    ]
    for left_index, left in enumerate(boxes):
        for right in boxes[left_index + 1 :]:
            assert _bbox_overlap_area(left, right) < 1.0
    plt.close(fig)


def test_compute_similarity_tracks_handles_gaps_and_mismatches():
    aligned = {
        "ref": "ACGTACGTAC",
        "match": "ACGTACGTAC",
        "with_mismatch": "ACGTAGGTAC",
        "with_gap": "ACGT-CGTAC",
    }

    tracks = compute_similarity_tracks(aligned, "ref", window=4, step=2)

    assert "ref" not in tracks
    assert set(tracks) == {"match", "with_mismatch", "with_gap"}

    match_positions, match_identities = tracks["match"]
    assert all(value == 1.0 for value in match_identities)
    assert match_positions[0] == 3

    _, mismatch_identities = tracks["with_mismatch"]
    assert max(mismatch_identities) == 1.0
    assert min(mismatch_identities) == 0.75

    _, gap_identities = tracks["with_gap"]
    for value in gap_identities:
        assert value is not None
        assert value == 1.0


def test_compute_similarity_tracks_skips_windows_without_valid_pairs():
    aligned = {
        "ref": "----ACGT",
        "sample": "----ACGT",
    }

    _, identities = compute_similarity_tracks(aligned, "ref", window=4, step=2)["sample"]
    assert identities[0] is None
    assert identities[-1] == 1.0


def test_pi_y_axis_limits_scale_with_low_diversity_values():
    hotspot = WindowResult("W0001", 1, 600, 300, 0.012, 600, "LSC", "peak", True, 1)

    _, y_top = _pi_y_axis_limits(
        [0.0, 0.006, 0.012],
        [hotspot],
        hotspot_mode="top-percent",
        hotspot_value=3.0,
    )

    assert y_top < 0.016


def test_pi_y_axis_limits_include_threshold_when_it_is_higher_than_data():
    threshold = 0.02

    _, y_top = _pi_y_axis_limits(
        [0.0, 0.006, 0.012],
        [],
        hotspot_mode="threshold",
        hotspot_value=threshold,
    )

    assert y_top > threshold
