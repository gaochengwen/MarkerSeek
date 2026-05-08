"""Output helpers shared by the CLI and web application."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from .diagnostics import CANONICAL_BASES
from .models import AnalysisResult, FeatureResult, PrimerResult
from .plotting import plot_pi_figure, plot_primer_summary, plot_similarity_figure
from .visualisation import compute_haplotype_network, compute_species_pca

RESULT_FILENAMES = (
    "pi_windows.tsv",
    "candidate_marker_features.tsv",
    "haplotype_assignments.tsv",
    "sample_metadata.tsv",
    "pi_plot.pdf",
    "pi_plot.png",
    "similarity_plot.pdf",
    "similarity_plot.png",
    "primers.tsv",
    "primer_summary.png",
)

MARKER_FEATURE_COLUMNS = [
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

PRIMER_COLUMNS = [
    "primer_id",
    "label_name",
    "rank",
    "fwd_seq",
    "rev_seq",
    "fwd_len",
    "rev_len",
    "fwd_gc",
    "rev_gc",
    "fwd_tm",
    "rev_tm",
    "fwd_self_any_th",
    "rev_self_any_th",
    "primer3_penalty",
    "target_start",
    "target_end",
    "amplicon_len",
    "cross_species_success_rate",
    "primer_score",
]


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


def marker_feature_rows(features):
    parent_genes_with_parts = {
        row.parent_gene
        for row in features
        if row.feature_id.startswith(f"{row.parent_gene}_part")
    }
    for row in features:
        if (
            row.feature_id == row.parent_gene
            and row.parent_gene in parent_genes_with_parts
        ):
            continue
        yield row


def write_marker_features_tsv(path: Path, features) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(MARKER_FEATURE_COLUMNS)
        for row in marker_feature_rows(features):
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
                    str(int(row.variable_sites)),
                    str(int(row.indel_sites)),
                    str(int(row.conserved_left_bp)),
                    str(int(row.conserved_right_bp)),
                    row.primer_available,
                    format_float(row.species_resolution),
                    str(int(row.unique_haplotype_count)),
                    format_float(row.species_specific_haplotype_ratio),
                    format_float(row.interspecific_divergence),
                    format_float(row.intraspecific_divergence),
                    format_float(row.nearest_neighbor_discrimination),
                    format_float(row.barcoding_gap),
                    format_float(row.misclassification_risk),
                    format_float(row.alignment_reliability),
                    format_float(row.markerseek_score),
                ]
            )


def write_haplotype_assignments_tsv(path: Path, features, sample_order: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["feature_id", "sample_name", "haplotype_id"])
        for feature in features:
            for index, sample_name in enumerate(sample_order):
                haplotype_id = feature.haplotypes[index] if index < len(feature.haplotypes) else "NA"
                writer.writerow([feature.feature_id, sample_name, haplotype_id])


def write_sample_metadata_tsv(path: Path, sample_metadata) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_name", "species", "source_path"])
        for metadata in sample_metadata:
            writer.writerow([metadata.sample_name, metadata.species, metadata.source_path])


def write_primers_tsv(path: Path, primers) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(PRIMER_COLUMNS)
        for row in primers:
            writer.writerow(
                [
                    row.primer_id,
                    row.label_name,
                    row.rank,
                    row.fwd_seq,
                    row.rev_seq,
                    row.fwd_len,
                    row.rev_len,
                    f"{row.fwd_gc:.1f}",
                    f"{row.rev_gc:.1f}",
                    f"{row.fwd_tm:.1f}",
                    f"{row.rev_tm:.1f}",
                    f"{row.fwd_self_any_th:.1f}",
                    f"{row.rev_self_any_th:.1f}",
                    f"{row.primer3_penalty:.6f}",
                    row.target_start,
                    row.target_end,
                    int(round(row.amplicon_mean_len)),
                    f"{row.cross_species_success_rate:.6f}",
                    f"{row.primer_score:.1f}",
                ]
            )


def write_feature_payloads(outdir: str | Path, result: AnalysisResult) -> None:
    """Persist internal JSON payloads for hotspot feature detail pages."""

    output_dir = Path(outdir).resolve()
    payload_dir = output_dir / "feature_payload"
    payload_dir.mkdir(parents=True, exist_ok=True)
    for stale_payload in payload_dir.glob("*.json"):
        stale_payload.unlink()

    species_map = {metadata.sample_name: metadata.species for metadata in result.sample_metadata}
    for feature in marker_feature_rows(result.features):
        if not feature.is_hotspot:
            continue

        spans = feature.spans(result.genome_length)
        aligned_block = _slice_aligned_block(result.aligned_sequences, spans)
        positions = _feature_positions(spans)
        snp_positions, indel_positions = _site_positions(aligned_block, positions)
        haplotypes = {
            sample_name: feature.haplotypes[index] if index < len(feature.haplotypes) else "NA"
            for index, sample_name in enumerate(result.sample_order)
        }
        primers = [_primer_payload(primer) for primer in result.primers if primer.label_name == feature.label_name]
        payload = {
            "feature": _feature_payload(feature),
            "pi_curve": {
                "positions": positions,
                "pi": [
                    0.0 if result.position_pi[position - 1] is None else float(result.position_pi[position - 1])
                    for position in positions
                ],
            },
            "snp_positions": snp_positions,
            "indel_positions": indel_positions,
            "aligned_block": aligned_block,
            "haplotypes": haplotypes,
            "species_map": species_map,
            "haplotype_network": compute_haplotype_network(aligned_block, species_map),
            "species_pca": compute_species_pca(aligned_block, species_map),
            "primers": primers,
        }
        write_json(payload_dir / f"{feature_id_safe(feature.feature_id)}.json", payload)


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
    primer_design: bool = False,
    mafft_bin: str = "mafft",
) -> None:
    output_dir = Path(outdir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_windows_tsv(output_dir / "pi_windows.tsv", result.windows)
    write_marker_features_tsv(output_dir / "candidate_marker_features.tsv", result.features)
    write_haplotype_assignments_tsv(output_dir / "haplotype_assignments.tsv", result.features, result.sample_order)
    write_sample_metadata_tsv(output_dir / "sample_metadata.tsv", result.sample_metadata)
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
    if primer_design:
        write_primers_tsv(output_dir / "primers.tsv", result.primers)
        plot_primer_summary(result, output_dir)
    write_feature_payloads(output_dir, result)


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def feature_id_safe(feature_id: str) -> str:
    return feature_id.replace("/", "_").replace(".", "_").replace(":", "_")


def _slice_aligned_block(aligned_sequences: dict[str, str], spans: list[tuple[int, int]]) -> dict[str, str]:
    return {
        sample_name: "".join(sequence[start:end] for start, end in spans)
        for sample_name, sequence in aligned_sequences.items()
    }


def _feature_positions(spans: list[tuple[int, int]]) -> list[int]:
    return [position + 1 for start, end in spans for position in range(start, end)]


def _site_positions(aligned_block: dict[str, str], positions: list[int]) -> tuple[list[int], list[int]]:
    snp_positions: list[int] = []
    indel_positions: list[int] = []
    sequences = [sequence.upper() for sequence in aligned_block.values()]
    for position, column in zip(positions, zip(*(sequence for sequence in sequences), strict=False), strict=False):
        counts = Counter(column)
        canonical = {base for base in counts if base in CANONICAL_BASES}
        if len(canonical) >= 2:
            snp_positions.append(position)
        if counts["-"] and canonical:
            indel_positions.append(position)
    return snp_positions, indel_positions


def _feature_payload(feature: FeatureResult) -> dict:
    return {
        "feature_id": feature.feature_id,
        "feature_type": feature.feature_type,
        "parent_gene": feature.parent_gene,
        "label_name": feature.label_name,
        "start": feature.start,
        "end": feature.end,
        "strand": feature.strand,
        "length_bp": feature.length_bp,
        "region": feature.region,
        "pi": feature.pi,
        "variable_sites": feature.variable_sites,
        "indel_sites": feature.indel_sites,
        "conserved_left_bp": feature.conserved_left_bp,
        "conserved_right_bp": feature.conserved_right_bp,
        "primer_available": feature.primer_available,
        "species_resolution": feature.species_resolution,
        "unique_haplotype_count": feature.unique_haplotype_count,
        "species_specific_haplotype_ratio": feature.species_specific_haplotype_ratio,
        "interspecific_divergence": feature.interspecific_divergence,
        "intraspecific_divergence": feature.intraspecific_divergence,
        "nearest_neighbor_discrimination": feature.nearest_neighbor_discrimination,
        "barcoding_gap": feature.barcoding_gap,
        "misclassification_risk": feature.misclassification_risk,
        "alignment_reliability": feature.alignment_reliability,
        "markerseek_score": feature.markerseek_score,
    }


def _primer_payload(primer: PrimerResult) -> dict:
    return {
        "primer_id": primer.primer_id,
        "label_name": primer.label_name,
        "rank": primer.rank,
        "fwd_seq": primer.fwd_seq,
        "rev_seq": primer.rev_seq,
        "fwd_len": primer.fwd_len,
        "rev_len": primer.rev_len,
        "fwd_gc": primer.fwd_gc,
        "rev_gc": primer.rev_gc,
        "fwd_tm": primer.fwd_tm,
        "rev_tm": primer.rev_tm,
        "fwd_self_any_th": primer.fwd_self_any_th,
        "rev_self_any_th": primer.rev_self_any_th,
        "primer3_penalty": primer.primer3_penalty,
        "target_start": primer.target_start,
        "target_end": primer.target_end,
        "amplicon_min_len": primer.amplicon_min_len,
        "amplicon_max_len": primer.amplicon_max_len,
        "amplicon_mean_len": primer.amplicon_mean_len,
        "amplicon_len": int(round(primer.amplicon_mean_len)),
        "cross_species_success_rate": primer.cross_species_success_rate,
        "amplicon_variable_sites": primer.amplicon_variable_sites,
        "amplicon_indel_sites": primer.amplicon_indel_sites,
        "primer_score": primer.primer_score,
    }


RESULTS_ZIP_EXCLUDED_SUFFIXES = (".png",)


def create_results_zip(outdir: str | Path, archive_path: str | Path) -> Path:
    output_dir = Path(outdir).resolve()
    archive = Path(archive_path).resolve()
    archive.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(archive, "w", compression=ZIP_DEFLATED) as handle:
        for name in RESULT_FILENAMES:
            if name.endswith(RESULTS_ZIP_EXCLUDED_SUFFIXES):
                continue
            path = output_dir / name
            if path.exists() and path.is_file():
                handle.write(path, arcname=name)
    return archive
