"""Output helpers shared by the CLI and web application."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from .models import AnalysisResult
from .plotting import plot_pi_figure, plot_similarity_figure

RESULT_FILENAMES = (
    "pi_windows.tsv",
    "pi_features.tsv",
    "pi_plot.pdf",
    "pi_plot.png",
    "similarity_plot.pdf",
    "similarity_plot.png",
)


def format_float(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.6f}"


def write_windows_tsv(path: Path, windows) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "window_id",
                "start",
                "end",
                "midpoint",
                "pi",
                "valid_sites",
                "region",
                "label_name",
                "is_hotspot",
            ]
        )
        for row in windows:
            writer.writerow(
                [
                    row.window_id,
                    row.start,
                    row.end,
                    row.midpoint,
                    format_float(row.pi),
                    row.valid_sites,
                    row.region,
                    row.label_name,
                    "yes" if row.is_hotspot else "no",
                ]
            )


def write_features_tsv(path: Path, features) -> None:
    parent_genes_with_parts = {
        row.parent_gene
        for row in features
        if row.feature_id.startswith(f"{row.parent_gene}_part")
    }
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
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
            ]
        )
        for row in features:
            if (
                row.feature_id == row.parent_gene
                and row.parent_gene in parent_genes_with_parts
            ):
                continue
            writer.writerow(
                [
                    row.feature_id,
                    row.feature_type,
                    row.parent_gene,
                    row.label_name,
                    row.start,
                    row.end,
                    row.strand,
                    row.length_bp,
                    row.region,
                    format_float(row.pi),
                ]
            )


def write_analysis_outputs(
    result: AnalysisResult,
    outdir: str | Path,
    *,
    hotspot_mode: str,
    hotspot_value: float,
    label_mode: str,
    label_max: int | None,
    label_min_distance_bp: int,
    similarity_window: int,
    similarity_step: int,
    similarity_floor: float,
    include_similarity_plot: bool = True,
) -> None:
    output_dir = Path(outdir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_windows_tsv(output_dir / "pi_windows.tsv", result.windows)
    write_features_tsv(output_dir / "pi_features.tsv", result.features)
    plot_pi_figure(
        result,
        output_dir,
        hotspot_mode=hotspot_mode,
        hotspot_value=hotspot_value,
        label_mode=label_mode,
        label_max=label_max,
        label_min_distance_bp=label_min_distance_bp,
    )
    if include_similarity_plot:
        plot_similarity_figure(
            result,
            output_dir,
            similarity_window=similarity_window,
            similarity_step=similarity_step,
            similarity_floor=similarity_floor,
        )


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def create_results_zip(outdir: str | Path, archive_path: str | Path) -> Path:
    output_dir = Path(outdir).resolve()
    archive = Path(archive_path).resolve()
    archive.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(archive, "w", compression=ZIP_DEFLATED) as handle:
        for name in RESULT_FILENAMES:
            path = output_dir / name
            if path.exists() and path.is_file():
                handle.write(path, arcname=name)
    return archive
