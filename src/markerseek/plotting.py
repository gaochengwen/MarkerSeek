"""Plotting utilities for MarkerSeek."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch, Polygon, Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib.transforms import Bbox, IdentityTransform

from .models import AnalysisResult, AnnotatedInterval
from .models import WindowResult

CANONICAL_BASES = {"A", "C", "G", "T"}
SIMILARITY_FILL_COLOR = "#2f3236"
SIMILARITY_LINE_COLOR = "#101216"
SIMILARITY_BG_ALPHA = 0.18

MM_PER_INCH = 25.4
FIGURE_WIDTH_MM = 183
FIGURE_HEIGHT_MM = 89
PLOT_DPI = 600
REGION_COLORS = {
    "LSC": "#f1f4f7",
    "IRb": "#eef4eb",
    "SSC": "#faeadb",
    "IRa": "#f3e9f2",
}
PI_LINE_COLOR = "#0d4596"
HOTSPOT_COLOR = "#252525"
ANNOTATION_LINE_COLOR = "#3a3a3a"
GRID_COLOR = "#d6dde3"
FEATURE_COLORS = {
    "Photosynthesis": "#5b8cc0",
    "ATP synthase": "#55b4c3",
    "Ribosomal proteins": "#f0a05f",
    "tRNA": "#6f9a60",
    "rRNA": "#aac987",
    "RNA polymerase": "#da77a3",
    "NADH dehydrogenase": "#8872af",
    "Other / hypothetical": "#a2a2a2",
}
LABEL_MODE_CHOICES = {"peak-only", "all", "none"}


def _label_for_region(region_name: str) -> str:
    return "IR" if region_name.startswith("IR") else region_name


def _preferred_font_family() -> str:
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    for candidate in ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"):
        if candidate in available_fonts:
            return candidate
    return "sans-serif"


def _window_span(window: WindowResult, genome_length: int) -> int:
    if window.end > window.start:
        return window.end - window.start
    return (genome_length - window.start) + window.end


def _group_hotspot_clusters(hotspots: list[WindowResult], genome_length: int) -> list[list[WindowResult]]:
    if not hotspots:
        return []

    sorted_hotspots = sorted(hotspots, key=lambda window: (window.midpoint, window.start, window.end))
    clusters: list[list[WindowResult]] = [[sorted_hotspots[0]]]
    previous = sorted_hotspots[0]
    for window in sorted_hotspots[1:]:
        overlap_gap = max(1, _window_span(previous, genome_length) // 2)
        if window.start <= previous.end or (window.midpoint - previous.midpoint) <= overlap_gap:
            clusters[-1].append(window)
        else:
            clusters.append([window])
        previous = window
    return clusters


def _select_peak_windows(
    hotspots: list[WindowResult],
    genome_length: int,
) -> list[WindowResult]:
    clusters = _group_hotspot_clusters(hotspots, genome_length)
    return [
        max(
            cluster,
            key=lambda window: (
                window.pi if window.pi is not None else float("-inf"),
                -(window.hotspot_rank or 10**9),
                -window.valid_sites,
                -window.midpoint,
            ),
        )
        for cluster in clusters
    ]


def select_label_windows(
    result: AnalysisResult,
    *,
    hotspot_mode: str = "top-percent",
    label_mode: str,
    label_max: int | None,
    label_min_distance_bp: int,
) -> list[WindowResult]:
    hotspots = [window for window in result.windows if window.is_hotspot and window.pi is not None]
    if not hotspots or label_mode == "none":
        return []
    if label_mode == "all":
        candidates = sorted(hotspots, key=lambda window: (window.midpoint, window.start, window.end))
    elif label_mode == "peak-only":
        candidates = _select_peak_windows(hotspots, result.genome_length)
    else:
        raise ValueError(f"Unsupported label_mode: {label_mode}")

    ranked_candidates = sorted(
        candidates,
        key=lambda window: (
            -(window.pi if window.pi is not None else float("-inf")),
            window.hotspot_rank or 10**9,
            window.midpoint,
        ),
    )
    chosen: list[WindowResult] = []
    chosen_labels: set[str] = set()
    for window in ranked_candidates:
        if window.label_name and window.label_name in chosen_labels:
            continue
        nearby_selected = [
            selected
            for selected in chosen
            if abs(window.midpoint - selected.midpoint) < label_min_distance_bp
        ]
        if nearby_selected:
            continue
        chosen.append(window)
        if window.label_name:
            chosen_labels.add(window.label_name)
        if label_max is not None and len(chosen) >= label_max:
            break
    return sorted(chosen, key=lambda window: window.midpoint)


def _hotspot_threshold_label(hotspot_mode: str, hotspot_value: float) -> str | None:
    if hotspot_mode != "threshold":
        return None
    return f"Threshold Pi = {hotspot_value:.3f}"


def _format_bp_tick(value: float, _) -> str:
    if value <= 1:
        return "0"
    if value >= 1000:
        return f"{int(round(value / 1000))}k"
    return f"{int(value)}"


def _position_ticks(genome_length: int) -> list[int]:
    ticks = [1]
    ticks.extend(range(20000, genome_length, 20000))
    if genome_length > 1:
        ticks.append(genome_length)
    return ticks


def _shade_regions(ax, result: AnalysisResult) -> None:
    for region in sorted(result.regions, key=lambda item: item.start):
        for start, end in region.spans(result.genome_length):
            ax.axvspan(
                start + 1,
                end,
                color=REGION_COLORS.get(region.name, "#f2f2f2"),
                alpha=0.64,
                zorder=0,
                linewidth=0,
            )


def _region_spans_one_based(region, genome_length: int) -> list[tuple[int, int]]:
    return [(start + 1, end) for start, end in region.spans(genome_length)]


def _plot_region_bar(ax, result: AnalysisResult) -> None:
    ax.set_xlim(1, result.genome_length)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.tick_params(axis="x", bottom=False, labelbottom=False, top=False, labeltop=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for region in sorted(result.regions, key=lambda item: item.start):
        spans = _region_spans_one_based(region, result.genome_length)
        for start, end in spans:
            ax.add_patch(
                Rectangle(
                    (start, 0.12),
                    max(end - start + 1, 1),
                    0.76,
                    facecolor=REGION_COLORS.get(region.name, "#f2f2f2"),
                    edgecolor="#b9c2ca",
                    linewidth=0.35,
                )
            )
        label_start, label_end = max(spans, key=lambda span: span[1] - span[0])
        if label_end - label_start > result.genome_length * 0.018:
            ax.text(
                (label_start + label_end) / 2,
                0.50,
                region.name,
                ha="center",
                va="center",
                fontsize=5.7,
                color="#20272d",
                fontstyle="italic",
            )


def _region_boundary_positions(result: AnalysisResult) -> list[int]:
    boundaries = {1}
    for region in result.regions:
        for start, end in region.spans(result.genome_length):
            if start > 0:
                boundaries.add(start + 1)
            if 0 < end < result.genome_length:
                boundaries.add(end + 1)
    return sorted(position for position in boundaries if 1 < position < result.genome_length)


def _draw_region_boundaries(axes, result: AnalysisResult) -> None:
    for axis in axes:
        for position in _region_boundary_positions(result):
            axis.axvline(
                position,
                color="#9aa1a8",
                linewidth=0.42,
                linestyle=(0, (3, 3)),
                zorder=1,
                alpha=0.62,
            )


def _feature_category(feature) -> str:
    gene_name = (feature.parent_gene or feature.label_name or feature.feature_id).lower()
    feature_type = feature.feature_type.lower()
    if feature_type == "trna" or gene_name.startswith("trn"):
        return "tRNA"
    if feature_type == "rrna" or gene_name.startswith("rrn"):
        return "rRNA"
    if gene_name.startswith(("psa", "psb", "pet", "rbc")):
        return "Photosynthesis"
    if gene_name.startswith("atp"):
        return "ATP synthase"
    if gene_name.startswith(("rpl", "rps")):
        return "Ribosomal proteins"
    if gene_name.startswith("ndh"):
        return "NADH dehydrogenase"
    if gene_name.startswith("rpo"):
        return "RNA polymerase"
    return "Other / hypothetical"


def _feature_spans_one_based(feature, genome_length: int) -> list[tuple[int, int]]:
    if not feature.spans_origin:
        return [(feature.start, feature.end)]
    return [(feature.start, genome_length), (1, feature.end)]


def _add_feature_arrow(ax, start: int, end: int, y: float, strand: int, color: str, genome_length: int) -> None:
    if end <= start:
        return

    min_width = max(70.0, genome_length * 0.00045)
    if end - start < min_width:
        if strand < 0:
            start = max(1, int(end - min_width))
        else:
            end = min(genome_length, int(start + min_width))

    height = 0.070
    span = max(end - start, 1)
    head = min(max(span * 0.30, genome_length * 0.00045), genome_length * 0.0019)
    if strand < 0:
        points = [
            (end, y - height / 2),
            (start + head, y - height / 2),
            (start, y),
            (start + head, y + height / 2),
            (end, y + height / 2),
        ]
    elif strand > 0:
        points = [
            (start, y - height / 2),
            (end - head, y - height / 2),
            (end, y),
            (end - head, y + height / 2),
            (start, y + height / 2),
        ]
    else:
        points = [
            (start, y - height / 2),
            (end, y - height / 2),
            (end, y + height / 2),
            (start, y + height / 2),
        ]
    ax.add_patch(Polygon(points, closed=True, facecolor=color, edgecolor="none", alpha=0.86))


def _feature_label_allowed(feature, span_start: int, span_end: int, genome_length: int) -> bool:
    return False


def _plot_feature_tracks(ax, result: AnalysisResult, *, show_xaxis: bool = True) -> None:
    ax.set_xlim(1, result.genome_length)
    ax.set_ylim(-0.21, 0.21)
    ax.set_ylabel("Genes", fontsize=6.2, labelpad=10)
    ax.set_yticks([-0.105, 0.105], labels=["(-)", "(+)"])
    ax.tick_params(axis="y", length=0, labelsize=5.4, colors="#202020")
    if show_xaxis:
        ax.tick_params(axis="x", bottom=True, labelbottom=True, width=0.45, length=2.2, pad=1.5, colors="#202020")
    else:
        ax.tick_params(axis="x", bottom=False, labelbottom=False, top=False, labeltop=False)
    ax.axhline(0, color="#303030", linewidth=0.34, zorder=0)
    for side, spine in ax.spines.items():
        spine.set_visible(side == "bottom" and show_xaxis)
        spine.set_linewidth(0.45)
        spine.set_color("#202020")
    if show_xaxis:
        ax.set_xticks(_position_ticks(result.genome_length))
        ax.xaxis.set_major_formatter(FuncFormatter(_format_bp_tick))
        ax.set_xlabel("Nucleotide position (bp)", fontsize=6.5, labelpad=3)
    else:
        ax.set_xticks([])
        ax.set_xlabel("")

    drawable_features = [
        feature
        for feature in result.features
        if feature.feature_type in {"gene", "tRNA", "rRNA"} and feature.length_bp > 0
    ]
    for feature in sorted(drawable_features, key=lambda item: (item.start, item.end, item.label_name)):
        category = _feature_category(feature)
        color = FEATURE_COLORS[category]
        y = 0.105 if feature.strand >= 0 else -0.105
        for start, end in _feature_spans_one_based(feature, result.genome_length):
            _add_feature_arrow(ax, start, end, y, feature.strand, color, result.genome_length)
            if _feature_label_allowed(feature, start, end, result.genome_length):
                ax.text(
                    (start + end) / 2,
                    y,
                    feature.label_name,
                    ha="center",
                    va="center",
                    fontsize=5.0,
                    color="#111111",
                    fontweight="normal",
                    bbox={"boxstyle": "round,pad=0.08", "facecolor": "white", "edgecolor": "none", "alpha": 0.70},
                    clip_on=True,
                )


def _plot_function_legend(ax) -> None:
    ax.axis("off")
    handles = [Patch(facecolor=color, edgecolor="none", label=label) for label, color in FEATURE_COLORS.items()]
    ax.legend(
        handles=handles,
        loc="center",
        ncol=8,
        frameon=False,
        fontsize=5.3,
        handlelength=0.82,
        handleheight=0.82,
        handletextpad=0.46,
        columnspacing=1.18,
        borderaxespad=0,
    )


def _draw_y_axis_top_tick(ax) -> None:
    x_min, x_max = ax.get_xlim()
    y_top = ax.get_ylim()[1]
    tick_len = (x_max - x_min) * 0.006
    ax.plot(
        [x_min, x_min - tick_len],
        [y_top, y_top],
        color="#202020",
        linewidth=0.45,
        clip_on=False,
        solid_capstyle="butt",
    )


def _pi_y_axis_limits(
    y_values: list[float],
    labeled_hotspots: list[WindowResult],
    *,
    hotspot_mode: str,
    hotspot_value: float,
) -> tuple[float, float]:
    if not y_values:
        return 0.0, 0.1

    y_max = max(y_values)
    reference_top = y_max
    if hotspot_mode == "threshold" and hotspot_value > reference_top:
        reference_top = hotspot_value

    margin_fraction = 0.18 if labeled_hotspots else 0.08
    top_margin = max(reference_top * margin_fraction, 0.001)
    return 0.0, reference_top + top_margin


def _annotation_candidate_offsets() -> list[tuple[int, int]]:
    return [
        (0, 5),
        (0, 7),
        (0, 9),
        (-5, 6),
        (5, 6),
        (0, 12),
        (-8, 9),
        (8, 9),
        (-12, 8),
        (12, 8),
        (0, 16),
        (-14, 12),
        (14, 12),
        (-20, 10),
        (20, 10),
        (-10, 18),
        (10, 18),
        (0, 22),
        (-26, 14),
        (26, 14),
        (-18, 24),
        (18, 24),
        (0, 30),
        (-34, 18),
        (34, 18),
        (-26, 30),
        (26, 30),
        (-44, 22),
        (44, 22),
        (-36, 36),
        (36, 36),
        (0, 42),
        (-56, 28),
        (56, 28),
        (-46, 44),
        (46, 44),
        (0, 52),
        (-68, 32),
        (68, 32),
    ]


def _bbox_overlap_area(left: Bbox, right: Bbox) -> float:
    overlap_width = min(left.x1, right.x1) - max(left.x0, right.x0)
    overlap_height = min(left.y1, right.y1) - max(left.y0, right.y0)
    if overlap_width <= 0 or overlap_height <= 0:
        return 0.0
    return overlap_width * overlap_height


def _bbox_outside_penalty(bbox: Bbox, container: Bbox) -> float:
    left = max(container.x0 - bbox.x0, 0)
    right = max(bbox.x1 - container.x1, 0)
    bottom = max(container.y0 - bbox.y0, 0)
    top = max(bbox.y1 - container.y1, 0)
    return (left * left) + (right * right) + (bottom * bottom) + (top * top)


def _segments_cross(
    seg_a: tuple[tuple[float, float], tuple[float, float]],
    seg_b: tuple[tuple[float, float], tuple[float, float]],
) -> bool:
    (ax_, ay_), (bx_, by_) = seg_a
    (cx_, cy_), (dx_, dy_) = seg_b

    def _orient(px: float, py: float, qx: float, qy: float, rx: float, ry: float) -> float:
        return (qx - px) * (ry - py) - (qy - py) * (rx - px)

    o1 = _orient(ax_, ay_, bx_, by_, cx_, cy_)
    o2 = _orient(ax_, ay_, bx_, by_, dx_, dy_)
    o3 = _orient(cx_, cy_, dx_, dy_, ax_, ay_)
    o4 = _orient(cx_, cy_, dx_, dy_, bx_, by_)
    return ((o1 > 0) != (o2 > 0)) and ((o3 > 0) != (o4 > 0))


def _peak_annotation_alignment(x_offset: int) -> str:
    if x_offset < -2:
        return "right"
    if x_offset > 2:
        return "left"
    return "center"


def _peak_text_bbox(ax, window: WindowResult, x_offset: int, y_offset: int, renderer) -> Bbox:
    anchor_x, anchor_y = ax.transData.transform((window.midpoint, window.pi))
    points_to_pixels = ax.figure.dpi / 72.0
    text = ax.text(
        anchor_x + (x_offset * points_to_pixels),
        anchor_y + (y_offset * points_to_pixels),
        window.label_name or window.window_id,
        transform=IdentityTransform(),
        ha=_peak_annotation_alignment(x_offset),
        va="bottom",
        fontsize=5.0,
        fontstyle="italic",
        color="#111111",
        alpha=0,
    )
    bbox = text.get_window_extent(renderer=renderer)
    text.remove()
    return bbox


def _connector_needed(x_offset: int, y_offset: int) -> bool:
    return abs(x_offset) > 4 or y_offset > 12


def _make_peak_annotation(ax, window: WindowResult, x_offset: int, y_offset: int):
    arrowprops = None
    if _connector_needed(x_offset, y_offset):
        arrowprops = {
            "arrowstyle": "-",
            "color": ANNOTATION_LINE_COLOR,
            "linewidth": 0.42,
            "linestyle": (0, (2.2, 2.2)),
            "shrinkA": 0,
            "shrinkB": 2,
        }
    return ax.annotate(
        window.label_name or window.window_id,
        xy=(window.midpoint, window.pi),
        xytext=(x_offset, y_offset),
        textcoords="offset points",
        ha=_peak_annotation_alignment(x_offset),
        va="bottom",
        fontsize=5.0,
        fontstyle="italic",
        color="#111111",
        arrowprops=arrowprops,
        annotation_clip=True,
        zorder=5,
    )


def _annotate_peak_windows(ax, windows: list[WindowResult], reserved_bboxes: list[Bbox] | None = None) -> list:
    if not windows:
        return []

    figure = ax.figure
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    axes_bbox = ax.get_window_extent(renderer=renderer)
    safe_bbox = Bbox.from_extents(axes_bbox.x0 + 4, axes_bbox.y0 + 4, axes_bbox.x1 - 4, axes_bbox.y1 - 4)
    occupied = list(reserved_bboxes or [])
    placed_lines: list[tuple[tuple[float, float], tuple[float, float]]] = []
    annotations_by_x = []
    points_to_pixels = ax.figure.dpi / 72.0

    ranked_windows = sorted(
        [window for window in windows if window.pi is not None],
        key=lambda window: (-(window.pi or 0.0), window.midpoint),
    )
    for window in ranked_windows:
        anchor_pixel = ax.transData.transform((window.midpoint, window.pi))
        anchor_xy = (float(anchor_pixel[0]), float(anchor_pixel[1]))
        best_offset = _annotation_candidate_offsets()[0]
        best_score = float("inf")
        best_line: tuple[tuple[float, float], tuple[float, float]] | None = None
        best_connector_drawn = False
        for x_offset, y_offset in _annotation_candidate_offsets():
            bbox = _peak_text_bbox(ax, window, x_offset, y_offset, renderer).expanded(1.05, 1.18)
            overlap = sum(_bbox_overlap_area(bbox, placed) for placed in occupied)
            outside = _bbox_outside_penalty(bbox, safe_bbox)
            text_xy = (
                anchor_xy[0] + x_offset * points_to_pixels,
                anchor_xy[1] + y_offset * points_to_pixels,
            )
            candidate_segment = (anchor_xy, text_xy)
            connector_drawn = _connector_needed(x_offset, y_offset)
            crossings = (
                sum(1 for placed in placed_lines if _segments_cross(candidate_segment, placed))
                if connector_drawn
                else 0
            )
            line_length = (x_offset * x_offset) + (y_offset * y_offset)
            score = (
                (outside * 100000.0)
                + (overlap * 3000.0)
                + (crossings * 60000.0)
                + (line_length * 0.08)
            )
            if score < best_score:
                best_score = score
                best_offset = (x_offset, y_offset)
                best_line = candidate_segment
                best_connector_drawn = connector_drawn

        annotation = _make_peak_annotation(ax, window, best_offset[0], best_offset[1])
        occupied.append(_peak_text_bbox(ax, window, best_offset[0], best_offset[1], renderer).expanded(1.05, 1.18))
        if best_connector_drawn and best_line is not None:
            placed_lines.append(best_line)
        annotations_by_x.append((window.midpoint, annotation))

    return [annotation for _, annotation in sorted(annotations_by_x, key=lambda item: item[0])]


def plot_pi_figure(
    result: AnalysisResult,
    outdir: Path,
    *,
    hotspot_mode: str = "top-percent",
    hotspot_value: float = 5.0,
    label_mode: str = "peak-only",
    label_max: int | None = None,
    label_min_distance_bp: int = 0,
) -> tuple[Path, Path]:
    """Create publication-style Pi plots in PDF and PNG formats."""

    if label_mode not in LABEL_MODE_CHOICES:
        raise ValueError(f"Unsupported label_mode: {label_mode}")

    plt.rcParams.update(
        {
            "font.family": _preferred_font_family(),
            "font.size": 6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.labelsize": 6.5,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.8,
            "axes.linewidth": 0.5,
        }
    )

    width_in = FIGURE_WIDTH_MM / MM_PER_INCH
    height_in = FIGURE_HEIGHT_MM / MM_PER_INCH
    fig = plt.figure(figsize=(width_in, height_in), constrained_layout=False)
    fig.patch.set_facecolor("white")
    grid = fig.add_gridspec(5, 1, height_ratios=[49.2, 4.8, 9.0, 3.0, 10.0], hspace=0.025)
    ax = fig.add_subplot(grid[0])
    region_ax = fig.add_subplot(grid[1], sharex=ax)
    feature_ax = fig.add_subplot(grid[2], sharex=ax)
    spacer_ax = fig.add_subplot(grid[3])
    legend_ax = fig.add_subplot(grid[4])
    spacer_ax.axis("off")

    _plot_region_bar(region_ax, result)
    _plot_feature_tracks(feature_ax, result)
    _plot_function_legend(legend_ax)

    x_values = [window.midpoint for window in result.windows if window.pi is not None]
    y_values = [window.pi for window in result.windows if window.pi is not None]
    ax.plot(
        x_values,
        y_values,
        color=PI_LINE_COLOR,
        linewidth=0.50,
        solid_joinstyle="round",
        solid_capstyle="round",
        zorder=2,
    )

    labeled_hotspots = select_label_windows(
        result,
        hotspot_mode=hotspot_mode,
        label_mode=label_mode,
        label_max=label_max,
        label_min_distance_bp=label_min_distance_bp,
    )
    if labeled_hotspots:
        ax.scatter(
            [window.midpoint for window in labeled_hotspots],
            [window.pi for window in labeled_hotspots],
            color="#ffffff",
            s=7,
            linewidths=0.42,
            edgecolors=HOTSPOT_COLOR,
            zorder=4,
        )

    ax.set_xlim(1, result.genome_length)
    y_bottom, y_top = _pi_y_axis_limits(
        y_values,
        labeled_hotspots,
        hotspot_mode=hotspot_mode,
        hotspot_value=hotspot_value,
    )
    ax.set_ylim(y_bottom, y_top)

    threshold_label = _hotspot_threshold_label(hotspot_mode, hotspot_value)
    threshold_artist = None
    if threshold_label is not None:
        ax.axhline(hotspot_value, color="#d65f5f", linewidth=0.55, linestyle=(0, (4, 3)), zorder=1)
        threshold_artist = ax.text(
            0.995,
            0.94,
            threshold_label,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=5.2,
            color="#b13a3a",
        )

    reserved_bboxes: list[Bbox] = []
    if threshold_artist is not None:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        reserved_bboxes.append(threshold_artist.get_window_extent(renderer=renderer).expanded(1.05, 1.18))
    _annotate_peak_windows(ax, labeled_hotspots, reserved_bboxes=reserved_bboxes)

    ax.set_xlabel("")
    ax.set_ylabel("Nucleotide diversity", labelpad=6)
    ax.set_facecolor("white")
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.38, linestyle=(0, (3, 3)))
    for side, spine in ax.spines.items():
        spine.set_visible(side in {"left", "bottom"})
        spine.set_linewidth(0.55 if side == "left" else 0.45)
        spine.set_color("#202020")
    ax.tick_params(width=0.45, length=2.3, color="#202020")
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.margins(x=0)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.3f}"))
    _draw_y_axis_top_tick(ax)

    fig.subplots_adjust(left=0.082, right=0.986, top=0.961, bottom=0.101)

    outdir.mkdir(parents=True, exist_ok=True)
    pdf_path = outdir / "pi_plot.pdf"
    png_path = outdir / "pi_plot.png"
    fig.savefig(pdf_path, dpi=PLOT_DPI)
    fig.savefig(png_path, dpi=PLOT_DPI)
    plt.close(fig)
    return pdf_path, png_path


def compute_similarity_tracks(
    aligned_sequences: dict[str, str],
    reference_name: str,
    *,
    window: int = 100,
    step: int = 20,
) -> dict[str, tuple[list[int], list[float | None]]]:
    """Sliding-window pairwise identity vs the reference (mVISTA-style).

    For each non-reference sample, walks the reference-projected alignment in
    fixed-size windows and reports per-window pairwise identity::

        identity = (# positions where ref == sample, both A/C/G/T)
                 / (# positions where both ref and sample are A/C/G/T)

    Gaps and IUPAC ambiguity codes are excluded from both numerator and
    denominator, so a window full of gaps yields ``None`` rather than 0%.

    Returns
    -------
    dict[sample_name, (positions, identities)]
        ``positions`` are 1-based midpoints in reference coordinates;
        ``identities`` are floats in [0.0, 1.0] or ``None`` for windows with
        no comparable sites.
    """

    if reference_name not in aligned_sequences:
        raise KeyError(f"Reference {reference_name!r} not in alignment")
    if window <= 0 or step <= 0:
        raise ValueError("Similarity window and step must be positive integers.")

    reference = aligned_sequences[reference_name]
    length = len(reference)
    tracks: dict[str, tuple[list[int], list[float | None]]] = {}

    for name, sequence in aligned_sequences.items():
        if name == reference_name:
            continue
        if len(sequence) != length:
            raise ValueError("Aligned sequences must share the reference projected length.")

        prefix_match = [0] * (length + 1)
        prefix_valid = [0] * (length + 1)
        for i, (ref_base, sample_base) in enumerate(zip(reference, sequence), start=1):
            valid = ref_base in CANONICAL_BASES and sample_base in CANONICAL_BASES
            prefix_valid[i] = prefix_valid[i - 1] + (1 if valid else 0)
            prefix_match[i] = prefix_match[i - 1] + (1 if (valid and ref_base == sample_base) else 0)

        positions: list[int] = []
        identities: list[float | None] = []
        if length <= window:
            valid_total = prefix_valid[length]
            match_total = prefix_match[length]
            positions.append(length // 2 + 1)
            identities.append(match_total / valid_total if valid_total else None)
        else:
            starts = list(range(0, length - window + 1, step))
            if starts[-1] != length - window:
                starts.append(length - window)
            for start in starts:
                end = start + window
                v = prefix_valid[end] - prefix_valid[start]
                m = prefix_match[end] - prefix_match[start]
                positions.append(start + window // 2 + 1)
                identities.append(m / v if v else None)
        tracks[name] = (positions, identities)
    return tracks


def _plot_xaxis_strip(ax, genome_length: int) -> None:
    ax.set_xlim(1, genome_length)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    for side, spine in ax.spines.items():
        spine.set_visible(side == "top")
        spine.set_linewidth(0.45)
        spine.set_color("#202020")
    ticks = _position_ticks(genome_length)
    ax.set_xticks(ticks)
    ax.set_xticklabels([])
    ax.tick_params(
        axis="x",
        top=True,
        labeltop=False,
        bottom=False,
        labelbottom=False,
        direction="in",
        width=0.45,
        length=2.4,
        colors="#202020",
    )
    for index, tick_position in enumerate(ticks):
        if index == 0:
            ha = "left"
        elif index == len(ticks) - 1:
            ha = "right"
        else:
            ha = "center"
        ax.text(
            tick_position,
            0.74,
            _format_bp_tick(tick_position, None),
            ha=ha,
            va="top",
            fontsize=5.4,
            color="#202020",
        )
    ax.text(
        0.5,
        0.04,
        "Nucleotide position (bp)",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.0,
        color="#202020",
    )


def _feature_category_from_interval(interval: AnnotatedInterval) -> str:
    gene_name = (interval.parent_gene or interval.label or interval.feature_id).lower()
    feature_type = interval.feature_type.lower()
    if feature_type == "trna" or gene_name.startswith("trn"):
        return "tRNA"
    if feature_type == "rrna" or gene_name.startswith("rrn"):
        return "rRNA"
    if gene_name.startswith(("psa", "psb", "pet", "rbc")):
        return "Photosynthesis"
    if gene_name.startswith("atp"):
        return "ATP synthase"
    if gene_name.startswith(("rpl", "rps")):
        return "Ribosomal proteins"
    if gene_name.startswith("ndh"):
        return "NADH dehydrogenase"
    if gene_name.startswith("rpo"):
        return "RNA polymerase"
    return "Other / hypothetical"


def _shade_exons(ax, exon_intervals: list[AnnotatedInterval], genome_length: int) -> None:
    for interval in exon_intervals:
        category = _feature_category_from_interval(interval)
        color = FEATURE_COLORS.get(category, FEATURE_COLORS["Other / hypothetical"])
        for span_start, span_end in interval.spans(genome_length):
            ax.axvspan(
                span_start + 1,
                span_end,
                color=color,
                alpha=SIMILARITY_BG_ALPHA,
                linewidth=0,
                zorder=0,
            )


def _plot_similarity_row(
    ax,
    sample_name: str,
    positions: list[int],
    identities: list[float | None],
    exon_intervals: list[AnnotatedInterval],
    genome_length: int,
    *,
    floor: float = 0.5,
    is_last: bool = False,
) -> None:
    ax.set_xlim(1, genome_length)
    ax.set_ylim(floor, 1.0)
    ax.set_yticks([])
    ax.text(
        1.006,
        1.0,
        "100%",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=3.8,
        color="#404040",
        clip_on=False,
    )
    ax.text(
        1.006,
        0.0,
        f"{int(floor * 100)}%",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=3.8,
        color="#404040",
        clip_on=False,
    )
    ax.tick_params(
        axis="y",
        which="both",
        left=False,
        right=False,
        labelleft=False,
        labelright=False,
    )
    ax.tick_params(axis="x", bottom=False, labelbottom=False, top=False, labeltop=False)
    for side, spine in ax.spines.items():
        spine.set_visible(side in {"left", "right", "bottom" if is_last else "top"})
        spine.set_linewidth(0.32)
        spine.set_color("#9aa0a6")
    ax.set_facecolor("white")

    _shade_exons(ax, exon_intervals, genome_length)

    if positions:
        clipped = [
            (pos, max(floor, min(1.0, identity if identity is not None else floor)))
            for pos, identity in zip(positions, identities)
        ]
        xs = [point[0] for point in clipped]
        ys = [point[1] for point in clipped]
        ax.plot(xs, ys, color=SIMILARITY_LINE_COLOR, linewidth=0.18, zorder=3)

    pretty_name = sample_name.rsplit("_", 1)[0].replace("_", " ")
    ax.set_ylabel(
        pretty_name,
        rotation=0,
        ha="right",
        va="center",
        fontsize=5.6,
        fontstyle="italic",
        labelpad=8,
    )


def plot_similarity_figure(
    result: AnalysisResult,
    outdir: Path,
    *,
    similarity_window: int = 100,
    similarity_step: int = 20,
    similarity_floor: float = 0.5,
) -> tuple[Path, Path]:
    """Create an mVISTA-style similarity figure (PDF + PNG).

    Layout (top to bottom):
      gene-annotation track  →  N similarity rows (one per non-reference
      sample, ``floor``–100% identity)  →  shared x-axis tick labels  →
      LSC/IRb/SSC/IRa region bar  →  feature-category legend.

    Per-sample background tinting follows the same FEATURE_COLORS palette as
    the legend, so exon spans (CDS / tRNA / rRNA) appear as faint coloured
    stripes in each sample row. This figure is independent from the existing
    Pi figure (``pi_plot.*``) and is written to ``similarity_plot.{pdf,png}``.
    """

    aligned = result.aligned_sequences
    if not aligned:
        raise ValueError(
            "AnalysisResult.aligned_sequences is empty; rerun run_analysis() to populate it."
        )

    if result.sample_order:
        ordered = [name for name in result.sample_order if name != result.reference_name]
    else:
        ordered = sorted(name for name in aligned if name != result.reference_name)
    sample_order = [name for name in ordered if name in aligned]
    if not sample_order:
        raise ValueError("Similarity figure requires at least one non-reference sample.")

    tracks = compute_similarity_tracks(
        aligned,
        result.reference_name,
        window=similarity_window,
        step=similarity_step,
    )

    plt.rcParams.update(
        {
            "font.family": _preferred_font_family(),
            "font.size": 6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.labelsize": 6.0,
            "xtick.labelsize": 5.6,
            "ytick.labelsize": 5.0,
            "axes.linewidth": 0.5,
        }
    )

    n_samples = len(sample_order)
    gene_track_mm = 6.5
    spacer_top_mm = 0.4
    per_sample_mm = 8.05 if n_samples <= 8 else (6.21 if n_samples <= 16 else 5.06)
    xaxis_strip_mm = 7.6
    spacer_xaxis_mm = 0.8
    region_bar_mm = 4.8
    spacer_region_mm = 1.8
    legend_mm = 7.5
    bottom_padding_mm = 1.5

    width_in = FIGURE_WIDTH_MM / MM_PER_INCH
    height_mm = (
        gene_track_mm
        + spacer_top_mm
        + per_sample_mm * n_samples
        + xaxis_strip_mm
        + spacer_xaxis_mm
        + region_bar_mm
        + spacer_region_mm
        + legend_mm
        + bottom_padding_mm
    )
    height_in = height_mm / MM_PER_INCH

    fig = plt.figure(figsize=(width_in, height_in), constrained_layout=False)
    fig.patch.set_facecolor("white")

    height_ratios = [
        gene_track_mm,
        spacer_top_mm,
        *([per_sample_mm] * n_samples),
        xaxis_strip_mm,
        spacer_xaxis_mm,
        region_bar_mm,
        spacer_region_mm,
        legend_mm,
    ]
    grid = fig.add_gridspec(len(height_ratios), 1, height_ratios=height_ratios, hspace=0.0)

    cursor = 0
    feature_ax = fig.add_subplot(grid[cursor])
    cursor += 1
    fig.add_subplot(grid[cursor]).axis("off")
    cursor += 1
    sample_axes = []
    for _ in range(n_samples):
        sample_axes.append(fig.add_subplot(grid[cursor]))
        cursor += 1
    xtick_ax = fig.add_subplot(grid[cursor])
    cursor += 1
    fig.add_subplot(grid[cursor]).axis("off")
    cursor += 1
    region_ax = fig.add_subplot(grid[cursor])
    cursor += 1
    fig.add_subplot(grid[cursor]).axis("off")
    cursor += 1
    legend_ax = fig.add_subplot(grid[cursor])

    _plot_feature_tracks(feature_ax, result, show_xaxis=False)
    for index, (ax, sample) in enumerate(zip(sample_axes, sample_order)):
        positions, identities = tracks[sample]
        _plot_similarity_row(
            ax,
            sample,
            positions,
            identities,
            result.exon_intervals,
            result.genome_length,
            floor=similarity_floor,
            is_last=(index == n_samples - 1),
        )
    _draw_region_boundaries(sample_axes, result)
    _plot_xaxis_strip(xtick_ax, result.genome_length)
    _plot_region_bar(region_ax, result)
    _plot_function_legend(legend_ax)

    bottom_fraction = bottom_padding_mm / height_mm
    fig.subplots_adjust(left=0.155, right=0.955, top=0.985, bottom=bottom_fraction)

    outdir.mkdir(parents=True, exist_ok=True)
    pdf_path = outdir / "similarity_plot.pdf"
    png_path = outdir / "similarity_plot.png"
    fig.savefig(pdf_path, dpi=PLOT_DPI)
    fig.savefig(png_path, dpi=PLOT_DPI)
    plt.close(fig)
    return pdf_path, png_path
