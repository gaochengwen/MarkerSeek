"""MarkerSeek score normalization and weighted aggregation."""

from __future__ import annotations

DEFAULT_WEIGHTS: dict[str, float] = {
    "pi": 0.18,
    "variable_site_density": 0.10,
    "indel_density": 0.05,
    "flanking_conservation_min": 0.12,
    "missing_ambig_ratio": 0.05,
    "alignment_reliability": 0.10,
    "species_resolution": 0.15,
    "barcoding_gap": 0.10,
    "nearest_neighbor_discrimination": 0.10,
    "length_suitability": 0.05,
}

METRIC_NORMALIZATION: dict[str, tuple[str, float, float]] = {
    "pi": ("higher_better", 0.0, 0.05),
    "variable_site_density": ("higher_better", 0.0, 0.10),
    "indel_density": ("higher_better", 0.0, 0.05),
    "flanking_conservation_min": ("higher_better", 0.0, 1.0),
    "missing_ambig_ratio": ("lower_better", 0.0, 0.20),
    "alignment_reliability": ("higher_better", 0.0, 1.0),
    "species_resolution": ("higher_better", 0.0, 1.0),
    "barcoding_gap": ("higher_better", -0.02, 0.05),
    "nearest_neighbor_discrimination": ("higher_better", 0.0, 1.0),
    "length_suitability": ("higher_better", 0.0, 1.0),
}

assert abs(sum(DEFAULT_WEIGHTS.values()) - 1.0) < 1e-6


def normalize_metric(value: float | None, direction: str, lo: float, hi: float) -> float:
    """Normalize a metric to 0..1, optionally inverting lower-is-better values."""

    if value is None:
        return 0.0
    if direction not in {"higher_better", "lower_better"}:
        raise ValueError(f"Unsupported metric direction: {direction}")
    if hi <= lo:
        raise ValueError("Metric upper bound must be greater than lower bound.")
    clipped = min(max(float(value), lo), hi)
    normalized = (clipped - lo) / (hi - lo)
    if direction == "lower_better":
        return 1.0 - normalized
    return normalized


def length_suitability(length_bp: int) -> float:
    """Score feature length using the MarkerSeek trapezoidal suitability curve."""

    length = max(0, int(length_bp))
    if length < 300:
        return 0.0
    if length < 600:
        return (length - 300) / 300
    if length <= 1500:
        return 1.0
    if length < 3000:
        return 1.0 - ((length - 1500) / 1500) * 0.5
    if length < 5000:
        return 0.5 - ((length - 3000) / 2000) * 0.5
    return 0.0


def compute_score(metrics: dict[str, float | None], weights: dict[str, float] = DEFAULT_WEIGHTS) -> float:
    """Compute a 0..100 weighted MarkerSeek score rounded to one decimal place."""

    score = 0.0
    for metric_name, weight in weights.items():
        direction, lo, hi = METRIC_NORMALIZATION[metric_name]
        score += weight * normalize_metric(metrics.get(metric_name), direction, lo, hi)
    return round(min(max(score * 100.0, 0.0), 100.0), 1)
