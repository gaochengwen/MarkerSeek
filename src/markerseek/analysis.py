"""Core analysis pipeline for MarkerSeek."""

from __future__ import annotations

import math
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from io import StringIO
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqRecord import SeqRecord

from .diagnostics import (
    alignment_reliability_block,
    assign_haplotypes,
    conserved_column_mask,
    count_variable_and_indel_sites,
    divergence_stats,
    flanking_conservation_from_mask,
    species_resolution,
    species_specific_haplotype_ratio,
)
from .models import (
    AnalysisResult,
    AnnotatedInterval,
    FeatureResult,
    PrimerResult,
    RegionSegment,
    SampleMetadata,
    WindowResult,
)
from .primer import (
    AmpliconHit,
    FlankCandidate,
    PrimerPair,
    design_primers_for_region,
    design_primers_for_template,
    primer_binds_in_sample,
    find_conserved_flanks,
    in_silico_pcr,
    score_primer_pair,
)
from .plotting import select_label_windows
from .scoring import DEFAULT_WEIGHTS, compute_score, length_suitability
from .species import build_sample_metadata

SUPPORTED_SUFFIXES = {".gb", ".gbk", ".genbank"}
CANONICAL_BASES = {"A", "C", "G", "T"}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
KNOWN_IR_MARKERS = (
    "rrn16S",
    "rrn23S",
    "rrn4.5S",
    "rrn5S",
    "trnA-UGC",
    "trnI-CAU",
    "trnI-GAU",
    "trnN-GUU",
    "rpl2",
    "rpl23",
    "rps7",
    "rps12",
    "ndhB",
    "ycf2",
    "ycf15",
)
FEATURE_PRIORITY = {"gene": 1, "tRNA": 2, "rRNA": 3}
PRIMER_SUBTARGET_LENGTHS = (1000, 800, 600, 400, 250, 150, 120)
PRIMER_TARGETS_PER_LENGTH = 3
PRIMER_MAX_TARGET_CANDIDATES = 14
PRIMER_RELAXED_FLANK_IDENTITY = 0.75
PRIMER_RELAXED_FLANK_SCAN_BP = 1200
PRIMER_TEMPLATE_PADDING_FACTORS = (0.5, 1.0)
AUTOMATIC_IR_DETECTION_FAILED = (
    "Automatic IR detection failed. Provide --lsc/--irb/--ssc/--ira to override the region boundaries."
)
UNPARTITIONED_REGION_NAME = "Genome"


@dataclass(frozen=True)
class PrimerTargetCandidate:
    start: int
    end: int
    pi: float
    full_feature: bool = False


class MarkerSeekError(RuntimeError):
    """Raised when MarkerSeek cannot complete an analysis."""


def reverse_complement_text(sequence: str) -> str:
    """Return the reverse complement of an uppercase DNA string."""

    return sequence.translate(COMPLEMENT)[::-1]


def run_analysis(
    inputs: Iterable[str | Path],
    *,
    outdir: str | Path | None = None,
    reference: str | Path | None = None,
    window: int = 600,
    step: int = 200,
    hotspot_mode: str = "top-percent",
    hotspot_value: float = 3.0,
    mafft_bin: str = "mafft",
    mafft_threads: int | None = None,
    manual_regions: dict[str, tuple[int, int]] | None = None,
    label_mode: str = "peak-only",
    label_max: int | None = None,
    label_min_distance_bp: int = 0,
    primer_design: bool = False,
    primer_settings: dict | None = None,
    in_silico_settings: dict | None = None,
) -> AnalysisResult:
    """Run the full GenBank-to-Pi workflow and return structured results."""

    if window <= 0:
        raise MarkerSeekError("--window must be a positive integer.")
    if step <= 0:
        raise MarkerSeekError("--step must be a positive integer.")
    if mafft_threads is not None and mafft_threads <= 0:
        raise MarkerSeekError("--mafft-threads must be a positive integer.")

    genome_paths = discover_genbank_files(inputs)
    if len(genome_paths) < 2:
        raise MarkerSeekError("MarkerSeek requires at least two GenBank inputs.")

    genome_records = [load_genbank(path, index=index) for index, path in enumerate(genome_paths, start=1)]
    reference_record = choose_reference_record(genome_records, reference)
    reference_sequence = str(reference_record.record.seq).upper()
    genome_length = len(reference_sequence)

    atomic_features = extract_atomic_features(reference_record.record)
    exon_intervals = extract_exon_intervals(reference_record.record)
    regions = (
        parse_manual_regions(manual_regions, genome_length)
        if manual_regions
        else infer_regions_or_unpartitioned(reference_sequence, atomic_features)
    )
    atomic_features = assign_regions_to_features(atomic_features, regions, genome_length)
    feature_catalog = build_feature_catalog(atomic_features, genome_length)
    feature_catalog = assign_regions_to_features(feature_catalog, regions, genome_length)

    aligned_records = run_mafft_alignment(genome_records, mafft_bin, mafft_threads=mafft_threads)
    aligned_sequences = project_alignment_to_reference(aligned_records, reference_record.sample_name)
    if not aligned_sequences:
        raise MarkerSeekError("MAFFT produced no aligned sequences.")
    if len(next(iter(aligned_sequences.values()))) != genome_length:
        raise MarkerSeekError("Projected alignment length does not match the reference genome length.")

    position_pi, valid_mask = compute_position_pi(aligned_sequences.values())
    prefix_pi, prefix_valid = build_prefix_arrays(position_pi, valid_mask)

    windows = build_window_results(
        genome_length=genome_length,
        window=window,
        step=step,
        prefix_pi=prefix_pi,
        prefix_valid=prefix_valid,
        regions=regions,
        labels=feature_catalog,
    )
    features = build_feature_results(
        feature_catalog,
        genome_length=genome_length,
        prefix_pi=prefix_pi,
        prefix_valid=prefix_valid,
    )

    apply_hotspots(windows, hotspot_mode, hotspot_value)

    sample_order = [reference_record.sample_name] + [
        record.sample_name
        for record in genome_records
        if record.sample_name != reference_record.sample_name
    ]

    result = AnalysisResult(
        reference_name=reference_record.sample_name,
        genome_length=genome_length,
        sample_count=len(aligned_sequences),
        regions=regions,
        position_pi=position_pi,
        windows=windows,
        features=features,
        aligned_sequences=aligned_sequences,
        exon_intervals=exon_intervals,
        sample_order=sample_order,
    )
    sample_metadata = build_sample_metadata(genome_records, result.sample_order)
    result.sample_metadata = sample_metadata
    result.score_weights = dict(DEFAULT_WEIGHTS)
    conserved_mask = enrich_features_with_diagnostics(result, result.aligned_sequences, sample_metadata)
    apply_marker_scores(result)
    promote_label_window_candidates(
        result.features,
        result.windows,
        genome_length=result.genome_length,
        label_mode=label_mode,
        label_max=label_max,
        label_min_distance_bp=label_min_distance_bp,
    )
    promote_top_score_candidates(result.features, extra_count=10)
    if primer_design:
        effective_in_silico_settings = dict(in_silico_settings or {})
        effective_in_silico_settings.setdefault("mafft_bin", mafft_bin)
        design_primers_for_hotspots(
            result,
            genome_records,
            primer_settings or {},
            effective_in_silico_settings,
            conserved_mask=conserved_mask,
        )
    return result


def enrich_features_with_diagnostics(
    result: AnalysisResult,
    aligned_sequences: dict[str, str],
    sample_metadata: list[SampleMetadata],
) -> list[bool]:
    """Fill feature-level SNP, haplotype, divergence, flank, and score diagnostics."""

    species_map = {metadata.sample_name: metadata.species for metadata in sample_metadata}
    conserved_mask = conserved_column_mask(aligned_sequences)
    for feature in result.features:
        aligned_block = slice_aligned_block(aligned_sequences, feature.spans(result.genome_length))
        variable_sites, indel_sites = count_variable_and_indel_sites(aligned_block)
        haplotypes = assign_haplotypes(aligned_block)
        divergences = divergence_stats(aligned_block, species_map)
        conserved_left_bp, conserved_right_bp = flanking_conservation_from_mask(
            conserved_mask,
            region_start=feature.start - 1,
            region_end=feature.end,
            genome_length=result.genome_length,
        )

        feature.variable_sites = variable_sites
        feature.indel_sites = indel_sites
        feature.conserved_left_bp = conserved_left_bp
        feature.conserved_right_bp = conserved_right_bp
        feature.species_resolution = species_resolution(haplotypes, species_map)
        feature.unique_haplotype_count = len(set(haplotypes.values()))
        feature.species_specific_haplotype_ratio = species_specific_haplotype_ratio(haplotypes, species_map)
        feature.interspecific_divergence = divergences["interspecific"]
        feature.intraspecific_divergence = divergences["intraspecific"]
        feature.nearest_neighbor_discrimination = divergences["nearest_neighbor_discrimination"]
        feature.barcoding_gap = divergences["barcoding_gap"]
        feature.misclassification_risk = divergences["misclassification_risk"]
        feature.alignment_reliability = alignment_reliability_block(aligned_block)
        feature.missing_ambig_ratio = missing_ambig_ratio(aligned_block)
        feature.haplotypes = [haplotypes.get(sample_name, "NA") for sample_name in result.sample_order]
    return conserved_mask


REPLICATE_DEPENDENT_SCORE_METRICS = ("barcoding_gap", "nearest_neighbor_discrimination")


def _feature_score_metrics(feature: FeatureResult) -> dict[str, float | None]:
    """Build the input dict for ``compute_score`` from a feature's stored fields."""

    length_bp = max(feature.length_bp, 1)
    return {
        "pi": feature.pi,
        "variable_site_density": feature.variable_sites / length_bp,
        "indel_density": feature.indel_sites / length_bp,
        "flanking_conservation_min": min(feature.conserved_left_bp, feature.conserved_right_bp) / 200.0,
        "missing_ambig_ratio": feature.missing_ambig_ratio,
        "alignment_reliability": feature.alignment_reliability,
        "species_resolution": feature.species_resolution,
        "barcoding_gap": feature.barcoding_gap,
        "nearest_neighbor_discrimination": feature.nearest_neighbor_discrimination,
        "length_suitability": length_suitability(feature.length_bp),
    }


def effective_score_weights(
    features: list[FeatureResult], base_weights: dict[str, float]
) -> dict[str, float]:
    """Drop replicate-dependent metrics that are NA across the **entire dataset**
    (i.e., one individual per species) and proportionally redistribute their
    weight to the surviving metrics so ``sum(weights) == 1.0``. Without this, a
    one-individual-per-species plastome collection has a structural ~20-point
    score ceiling (``barcoding_gap`` + ``nearest_neighbor_discrimination``
    contribute 20% of the weight but are computable only with replicates)."""

    drop = [
        name
        for name in REPLICATE_DEPENDENT_SCORE_METRICS
        if name in base_weights
        and all(getattr(feature, name) is None for feature in features)
    ]
    if not drop:
        return dict(base_weights)
    keep = {k: v for k, v in base_weights.items() if k not in drop}
    total = sum(keep.values())
    if total <= 0:
        return keep
    return {k: v / total for k, v in keep.items()}


def apply_marker_scores(result: AnalysisResult) -> None:
    """Compute every feature's MarkerSeek score using the dataset's effective
    weight set, and persist that weight set on ``result.score_weights``."""

    weights = effective_score_weights(result.features, DEFAULT_WEIGHTS)
    result.score_weights = weights
    for feature in result.features:
        feature.markerseek_score = compute_score(_feature_score_metrics(feature), weights)


def feature_window_overlap(feature: FeatureResult, window: WindowResult, genome_length: int) -> int:
    overlap = 0
    window_spans = (
        [(window.start - 1, window.end)]
        if window.end >= window.start
        else [(window.start - 1, genome_length), (0, window.end)]
    )
    for feature_start, feature_end in feature.spans(genome_length):
        for window_start, window_end in window_spans:
            left = max(feature_start, window_start)
            right = min(feature_end, window_end)
            if right > left:
                overlap += right - left
    return overlap


def feature_midpoint(feature: FeatureResult, genome_length: int) -> int:
    return ((feature.start - 1) + (feature.length_bp // 2)) % genome_length + 1


def feature_for_labeled_window(
    features: list[FeatureResult],
    window: WindowResult,
    genome_length: int,
    selected_feature_ids: set[str],
) -> FeatureResult | None:
    exact_label_matches = [feature for feature in features if feature.label_name == window.label_name]
    candidates = exact_label_matches or [
        feature
        for feature in features
        if feature_window_overlap(feature, window, genome_length) > 0
    ]
    if not candidates:
        return None

    return min(
        candidates,
        key=lambda feature: (
            feature.feature_id in selected_feature_ids,
            -feature_window_overlap(feature, window, genome_length),
            abs(feature_midpoint(feature, genome_length) - window.midpoint),
            feature.start,
            feature.end,
            feature.feature_id,
        ),
    )


def promote_label_window_candidates(
    features: list[FeatureResult],
    windows: list[WindowResult],
    *,
    genome_length: int,
    label_mode: str = "peak-only",
    label_max: int | None = None,
    label_min_distance_bp: int = 0,
) -> None:
    """Promote the feature regions that correspond to labels drawn on the Pi plot."""

    if not windows or label_mode == "none":
        return
    label_result = AnalysisResult(
        reference_name="",
        genome_length=genome_length,
        sample_count=0,
        regions=[],
        position_pi=[],
        windows=windows,
        features=features,
    )
    labeled_windows = select_label_windows(
        label_result,
        label_mode=label_mode,
        label_max=label_max,
        label_min_distance_bp=label_min_distance_bp,
    )
    selected_feature_ids: set[str] = set()
    rank = 1
    for window in labeled_windows:
        feature = feature_for_labeled_window(features, window, genome_length, selected_feature_ids)
        if feature is None or feature.feature_id in selected_feature_ids:
            continue
        feature.is_hotspot = True
        feature.hotspot_rank = rank
        selected_feature_ids.add(feature.feature_id)
        rank += 1


def promote_top_score_candidates(features: list[FeatureResult], *, extra_count: int = 10) -> None:
    """Keep Pi hotspot candidates and add the highest-scoring remaining features."""

    if extra_count <= 0:
        return
    existing_ranks = [feature.hotspot_rank for feature in features if feature.hotspot_rank is not None]
    next_rank = (max(existing_ranks) + 1) if existing_ranks else 1
    remaining = [
        feature
        for feature in features
        if not feature.is_hotspot and feature.markerseek_score is not None
    ]
    remaining.sort(
        key=lambda feature: (
            -float(feature.markerseek_score or 0.0),
            feature.start,
            feature.end,
            feature.label_name,
        )
    )
    for feature in remaining[:extra_count]:
        feature.is_hotspot = True
        feature.hotspot_rank = next_rank
        next_rank += 1


def design_primers_for_hotspots(
    result: AnalysisResult,
    raw_records: list["GenomeInput"],
    primer_settings: dict | None,
    in_silico_settings: dict | None,
    *,
    conserved_mask: list[bool] | None = None,
) -> None:
    """Design primers for every hotspot feature and attach scored results."""

    result.primers.clear()
    result.primer_amplicon_samples.clear()
    for feature in result.features:
        feature.primer_available = "no"

    hotspots = [feature for feature in result.features if feature.is_hotspot]
    if not hotspots:
        return

    reference_record = next(
        (record for record in raw_records if record.sample_name == result.reference_name),
        None,
    )
    if reference_record is None:
        raise MarkerSeekError("Reference record is unavailable for primer design.")

    reference_seq = str(reference_record.record.seq).upper().replace("-", "").replace(".", "")
    reference_name = result.reference_name
    sample_seqs = {
        record.sample_name: str(record.record.seq).upper().replace("-", "").replace(".", "")
        for record in raw_records
    }
    total_samples = len(sample_seqs)
    if total_samples == 0:
        return
    other_sample_names = [name for name in sample_seqs if name != reference_name]

    effective_primer_settings = dict(primer_settings or {})
    effective_in_silico_settings = dict(in_silico_settings or {})
    max_pairs = int(effective_primer_settings.get("PRIMER_NUM_RETURN", 5))
    min_amplicon = int(effective_in_silico_settings.get("min_amplicon", 80))
    max_amplicon = int(effective_in_silico_settings.get("max_amplicon", 3000))
    primer_min_len = int(effective_primer_settings.get("PRIMER_MIN_SIZE", 18))
    strict_primer_settings = _screening_primer_settings(effective_primer_settings, max_pairs, relaxed=False)
    relaxed_primer_settings = _screening_primer_settings(effective_primer_settings, max_pairs, relaxed=True)
    max_mismatch = int(effective_in_silico_settings.get("max_mismatch", 1))
    anchor_3prime = int(effective_in_silico_settings.get("anchor_3prime", 5))
    flank_max_len = max(primer_min_len, int(effective_in_silico_settings.get("flank_max_len", 120)))
    mafft_bin = str(effective_in_silico_settings.get("mafft_bin", "mafft"))
    mask = conserved_mask if conserved_mask is not None else conserved_column_mask(result.aligned_sequences)

    for feature in hotspots:
        if feature.spans_origin:
            continue
        target_candidates = _primer_target_candidates(
            feature,
            result.position_pi,
            max_amplicon=max_amplicon,
            primer_min_len=primer_min_len,
        )
        scored_pairs: list[tuple[PrimerPair, dict[str, str]]] = []
        seen_primer_sequences: set[tuple[str, str]] = set()
        for target in target_candidates:
            primer_pairs = _primer_pairs_for_target(
                reference_seq,
                result.aligned_sequences,
                mask,
                feature,
                target,
                genome_length=result.genome_length,
                primer_min_len=primer_min_len,
                flank_max_len=flank_max_len,
                max_amplicon=max_amplicon,
                max_pairs=max_pairs,
                strict_primer_settings=strict_primer_settings,
                relaxed_primer_settings=relaxed_primer_settings,
            )
            for pair in primer_pairs:
                sequence_key = (pair.fwd_seq, pair.rev_seq)
                if sequence_key in seen_primer_sequences:
                    continue
                seen_primer_sequences.add(sequence_key)
                # Full in-silico PCR only on the reference — that is what
                # determines amplicon length and whether the primer pair is
                # usable at all. Skip if it does not amplify cleanly.
                ref_hits = in_silico_pcr(
                    pair.fwd_seq,
                    pair.rev_seq,
                    {reference_name: reference_seq},
                    min_amplicon=min_amplicon,
                    max_amplicon=max_amplicon,
                    max_mismatch=max_mismatch,
                    anchor_3prime=anchor_3prime,
                )
                ref_hit = ref_hits.get(reference_name)
                if ref_hit is None or ref_hit.multiple_hits:
                    score_primer_pair(pair, ref_hits, total_samples)
                    continue
                # Cross-species universality: confirm both primers bind in
                # every non-reference sample with a compatible amplicon
                # window. We collect the per-sample amplicon length from the
                # binding check (cheap) and reuse the genome-wide MAFFT
                # alignment for variable/indel counting (no per-pair MAFFT).
                successful_samples = {reference_name}
                sample_lengths: dict[str, int] = {reference_name: ref_hit.length}
                synthetic_hits: dict[str, AmpliconHit | None] = {reference_name: ref_hit}
                for sample_name in other_sample_names:
                    binding = primer_binds_in_sample(
                        pair.fwd_seq,
                        pair.rev_seq,
                        sample_seqs[sample_name],
                        min_amplicon=min_amplicon,
                        max_amplicon=max_amplicon,
                        max_mismatch=max_mismatch,
                        anchor_3prime=anchor_3prime,
                    )
                    if binding is None:
                        synthetic_hits[sample_name] = None
                        continue
                    successful_samples.add(sample_name)
                    sample_lengths[sample_name] = binding[2]
                    synthetic_hits[sample_name] = AmpliconHit(
                        sample_name=sample_name,
                        amplicon_seq="",
                        fwd_pos=binding[0],
                        rev_pos=binding[1],
                        length=binding[2],
                        fwd_mismatches=0,
                        rev_mismatches=0,
                        multiple_hits=False,
                    )

                # Slice the genome-wide projected alignment using reference
                # coordinates. ``project_alignment_to_reference`` made each
                # sample column-equivalent to the reference's gap-free index
                # space, so [fwd_pos, rev_pos + reverse_len) maps straight to
                # the homologous amplicon region in every sample.
                ref_start = ref_hit.fwd_pos
                ref_end = ref_hit.rev_pos + len(pair.rev_seq)
                aligned_block = {
                    name: result.aligned_sequences[name][ref_start:ref_end]
                    for name in successful_samples
                    if name in result.aligned_sequences
                }
                if len(aligned_block) >= 2:
                    var_sites, indel_sites = count_variable_and_indel_sites(aligned_block)
                else:
                    var_sites = indel_sites = 0

                lengths = list(sample_lengths.values())
                pair.amplicon_min_len = min(lengths)
                pair.amplicon_max_len = max(lengths)
                pair.amplicon_mean_len = sum(lengths) / len(lengths)
                pair.amplicon_variable_sites = var_sites
                pair.amplicon_indel_sites = indel_sites
                score_primer_pair(pair, synthetic_hits, total_samples)
                if pair.primer_score > 0:
                    scored_pairs.append((pair, successful_samples))
                    if len(scored_pairs) >= max_pairs:
                        break
            if len(scored_pairs) >= max_pairs:
                break

        if not scored_pairs:
            continue

        scored_pairs.sort(
            key=lambda item: (
                -item[0].primer_score,
                -item[0].cross_species_success_rate,
                item[0].primer3_penalty,
                item[0].fwd_seq,
                item[0].rev_seq,
            )
        )
        feature.primer_available = "yes"
        for rank, (pair, samples_with_amplicon) in enumerate(scored_pairs[:max_pairs], start=1):
            pair.rank = rank
            primer_result = _primer_result(feature.label_name, pair)
            result.primers.append(primer_result)
            result.primer_amplicon_samples[primer_result.primer_id] = set(samples_with_amplicon)


def _primer_pairs_for_target(
    reference_seq: str,
    aligned_sequences: dict[str, str],
    conserved_mask: list[bool],
    feature: FeatureResult,
    target: PrimerTargetCandidate,
    *,
    primer_min_len: int,
    flank_max_len: int,
    genome_length: int,
    max_amplicon: int,
    max_pairs: int,
    strict_primer_settings: dict,
    relaxed_primer_settings: dict,
) -> list[PrimerPair]:
    """Generate primer pairs using strict flanks, broad templates, then relaxed settings."""

    pairs: list[PrimerPair] = []
    seen_contexts: set[tuple[str, int, int, int, int]] = set()

    def add_region_pairs(
        *,
        settings_name: str,
        settings: dict,
        min_len: int,
        identity_threshold: float,
        max_scan: int,
    ) -> None:
        left_flank, right_flank = find_conserved_flanks(
            aligned_sequences,
            conserved_mask,
            target.start,
            target.end,
            max_scan=max_scan,
            min_len=min_len,
            max_len=max(flank_max_len, min_len),
            identity_threshold=identity_threshold,
        )
        if left_flank is None or right_flank is None:
            return
        context = (settings_name, left_flank.start, right_flank.end, target.start, target.end)
        if context in seen_contexts:
            return
        seen_contexts.add(context)
        pairs.extend(
            design_primers_for_region(
                reference_seq,
                target.start,
                target.end,
                left_flank,
                right_flank,
                max_pairs=max_pairs,
                primer_settings=settings,
            )
        )

    def add_template_pairs(*, settings_name: str, settings: dict) -> None:
        for template_start, template_end in _template_bounds_for_target(
            feature,
            target,
            genome_length=genome_length,
            max_amplicon=max_amplicon,
        ):
            context = (settings_name, template_start, template_end, target.start, target.end)
            if context in seen_contexts:
                continue
            seen_contexts.add(context)
            pairs.extend(
                design_primers_for_template(
                    reference_seq,
                    template_start,
                    template_end,
                    target.start,
                    target.end,
                    max_pairs=max_pairs,
                    primer_settings=settings,
                )
            )

    add_region_pairs(
        settings_name="strict",
        settings=strict_primer_settings,
        min_len=primer_min_len,
        identity_threshold=0.85,
        max_scan=500,
    )
    add_template_pairs(settings_name="strict", settings=strict_primer_settings)
    add_region_pairs(
        settings_name="strict_relaxed_flanks",
        settings=strict_primer_settings,
        min_len=primer_min_len,
        identity_threshold=PRIMER_RELAXED_FLANK_IDENTITY,
        max_scan=PRIMER_RELAXED_FLANK_SCAN_BP,
    )
    relaxed_min_len = int(relaxed_primer_settings.get("PRIMER_MIN_SIZE", primer_min_len))
    add_region_pairs(
        settings_name="relaxed",
        settings=relaxed_primer_settings,
        min_len=relaxed_min_len,
        identity_threshold=PRIMER_RELAXED_FLANK_IDENTITY,
        max_scan=PRIMER_RELAXED_FLANK_SCAN_BP,
    )
    add_template_pairs(settings_name="relaxed", settings=relaxed_primer_settings)
    return pairs


def _template_bounds_for_target(
    feature: FeatureResult,
    target: PrimerTargetCandidate,
    *,
    genome_length: int,
    max_amplicon: int,
) -> list[tuple[int, int]]:
    feature_start = feature.start - 1
    feature_end = feature.end
    broad_padding = max(max_amplicon // 2, 120)
    if target.full_feature:
        bounds = [(max(0, feature_start - broad_padding), min(genome_length, feature_end + broad_padding))]
    elif feature.feature_type in {"gene", "tRNA", "rRNA"}:
        bounds = [
            (feature_start, feature_end),
            (
                max(0, feature_start - max(120, max_amplicon // 4)),
                min(genome_length, feature_end + max(120, max_amplicon // 4)),
            ),
        ]
    else:
        bounds = [
            (feature_start, feature_end),
            (max(0, feature_start - broad_padding), min(genome_length, feature_end + broad_padding)),
        ]

    templates: list[tuple[int, int]] = []
    seen_templates: set[tuple[int, int]] = set()
    for bound_start, bound_end in bounds:
        if bound_end <= bound_start:
            continue
        for factor in PRIMER_TEMPLATE_PADDING_FACTORS:
            padding = max(80, int(max_amplicon * factor / 2))
            template_start = max(bound_start, target.start - padding)
            template_end = min(bound_end, target.end + padding)
            if template_start > target.start or template_end < target.end:
                continue
            key = (template_start, template_end)
            if key in seen_templates:
                continue
            seen_templates.add(key)
            templates.append(key)
    return templates


def _screening_primer_settings(base_settings: dict, max_pairs: int, *, relaxed: bool) -> dict:
    settings = dict(base_settings)
    requested = int(settings.get("PRIMER_NUM_RETURN", max_pairs) or max_pairs)
    settings["PRIMER_NUM_RETURN"] = max(requested, max_pairs * (4 if relaxed else 3), 24 if relaxed else 15)
    if not relaxed:
        return settings

    settings["PRIMER_MIN_TM"] = min(float(settings.get("PRIMER_MIN_TM", 52.0)), 50.0)
    settings["PRIMER_MAX_TM"] = max(float(settings.get("PRIMER_MAX_TM", 70.0)), 72.0)
    settings["PRIMER_MIN_SIZE"] = min(int(settings.get("PRIMER_MIN_SIZE", 18)), 16)
    settings["PRIMER_MAX_SIZE"] = max(int(settings.get("PRIMER_MAX_SIZE", 27)), 30)
    settings["PRIMER_MIN_GC"] = min(float(settings.get("PRIMER_MIN_GC", 25.0)), 20.0)
    settings["PRIMER_MAX_GC"] = max(float(settings.get("PRIMER_MAX_GC", 70.0)), 75.0)
    return settings


def _primer_target_candidates(
    feature: FeatureResult,
    position_pi: list[float | None],
    *,
    max_amplicon: int,
    primer_min_len: int,
) -> list[PrimerTargetCandidate]:
    """Return full or high-Pi subtargets that can fit inside the amplicon bounds."""

    if feature.spans_origin:
        return []
    feature_start = feature.start - 1
    feature_end = feature.end
    feature_len = feature_end - feature_start
    if feature_len <= 0:
        return []

    core_len_cap = max_amplicon - (2 * max(primer_min_len, 1))
    if core_len_cap <= 0:
        core_len_cap = max(1, min(max_amplicon, feature_len))
    max_target_len = max(1, min(feature_len, core_len_cap))
    max_subtarget_len = max(1, min(max_target_len, 1000))
    lengths = _primer_candidate_lengths(
        feature_len,
        max_target_len=max_subtarget_len,
    )
    pi_prefix, valid_prefix = _feature_pi_prefix(position_pi, feature_start, feature_len)

    full_candidate: PrimerTargetCandidate | None = None
    if feature_len <= max_target_len:
        full_candidate = PrimerTargetCandidate(
            start=feature_start,
            end=feature_end,
            pi=_interval_pi(pi_prefix, valid_prefix, 0, feature_len),
            full_feature=True,
        )

    subtargets: list[PrimerTargetCandidate] = []
    for target_len in lengths:
        if target_len == feature_len:
            continue
        for offset, pi_value in _top_pi_offsets(
            pi_prefix,
            valid_prefix,
            feature_len=feature_len,
            target_len=target_len,
            limit=PRIMER_TARGETS_PER_LENGTH,
        ):
            subtargets.append(
                PrimerTargetCandidate(
                    start=feature_start + offset,
                    end=feature_start + offset + target_len,
                    pi=pi_value,
                )
            )

    deduplicated = _deduplicate_primer_targets(subtargets)
    deduplicated.sort(key=lambda item: (-round(item.pi, 12), -(item.end - item.start), item.start, item.end))
    if full_candidate is not None:
        deduplicated = [
            target
            for target in deduplicated
            if (target.start, target.end) != (full_candidate.start, full_candidate.end)
        ]
        return [full_candidate, *deduplicated[: PRIMER_MAX_TARGET_CANDIDATES - 1]]
    return deduplicated[:PRIMER_MAX_TARGET_CANDIDATES]


def _primer_candidate_lengths(feature_len: int, *, max_target_len: int) -> list[int]:
    lengths = {max_target_len}
    lengths.add(min(max_target_len, 120))
    for length in PRIMER_SUBTARGET_LENGTHS:
        if length <= max_target_len:
            lengths.add(length)
    return sorted((length for length in lengths if 0 < length <= feature_len), reverse=True)


def _feature_pi_prefix(
    position_pi: list[float | None],
    feature_start: int,
    feature_len: int,
) -> tuple[list[float], list[int]]:
    pi_prefix = [0.0]
    valid_prefix = [0]
    for offset in range(feature_len):
        position = feature_start + offset
        pi_value = position_pi[position] if 0 <= position < len(position_pi) else None
        if pi_value is None or math.isnan(pi_value):
            pi_prefix.append(pi_prefix[-1])
            valid_prefix.append(valid_prefix[-1])
        else:
            pi_prefix.append(pi_prefix[-1] + float(pi_value))
            valid_prefix.append(valid_prefix[-1] + 1)
    return pi_prefix, valid_prefix


def _interval_pi(pi_prefix: list[float], valid_prefix: list[int], start: int, end: int) -> float:
    valid_sites = valid_prefix[end] - valid_prefix[start]
    if valid_sites <= 0:
        return 0.0
    return (pi_prefix[end] - pi_prefix[start]) / valid_sites


def _top_pi_offsets(
    pi_prefix: list[float],
    valid_prefix: list[int],
    *,
    feature_len: int,
    target_len: int,
    limit: int,
) -> list[tuple[int, float]]:
    if target_len <= 0 or target_len > feature_len or limit <= 0:
        return []
    ranked = [
        (_interval_pi(pi_prefix, valid_prefix, offset, offset + target_len), offset)
        for offset in range(0, feature_len - target_len + 1)
    ]
    ranked.sort(key=lambda item: (-round(item[0], 12), item[1]))

    selected: list[tuple[int, float]] = []
    min_spacing = max(1, target_len // 2)
    for pi_value, offset in ranked:
        if any(abs(offset - existing_offset) < min_spacing for existing_offset, _ in selected):
            continue
        selected.append((offset, pi_value))
        if len(selected) >= limit:
            break
    return selected


def _deduplicate_primer_targets(targets: list[PrimerTargetCandidate]) -> list[PrimerTargetCandidate]:
    best_by_interval: dict[tuple[int, int], PrimerTargetCandidate] = {}
    for target in targets:
        key = (target.start, target.end)
        existing = best_by_interval.get(key)
        if existing is None or target.pi > existing.pi:
            best_by_interval[key] = target
    return list(best_by_interval.values())


def _primer_result(label_name: str, pair: PrimerPair) -> PrimerResult:
    return PrimerResult(
        label_name=label_name,
        rank=pair.rank,
        fwd_seq=pair.fwd_seq,
        rev_seq=pair.rev_seq,
        fwd_len=pair.fwd_len,
        rev_len=pair.rev_len,
        fwd_gc=round(pair.fwd_gc, 1),
        rev_gc=round(pair.rev_gc, 1),
        fwd_tm=round(pair.fwd_tm, 1),
        rev_tm=round(pair.rev_tm, 1),
        fwd_self_any_th=round(pair.fwd_self_any_th, 1),
        rev_self_any_th=round(pair.rev_self_any_th, 1),
        primer3_penalty=pair.primer3_penalty,
        target_start=pair.target_start + 1,
        target_end=pair.target_end,
        amplicon_min_len=pair.amplicon_min_len,
        amplicon_max_len=pair.amplicon_max_len,
        amplicon_mean_len=round(pair.amplicon_mean_len, 1),
        cross_species_success_rate=round(pair.cross_species_success_rate, 6),
        amplicon_variable_sites=pair.amplicon_variable_sites,
        amplicon_indel_sites=pair.amplicon_indel_sites,
        primer_score=round(pair.primer_score, 1),
    )


def slice_aligned_block(aligned_sequences: dict[str, str], spans: list[tuple[int, int]]) -> dict[str, str]:
    """Slice and concatenate reference-coordinate spans from each aligned sequence."""

    return {
        sample_name: "".join(sequence[start:end] for start, end in spans)
        for sample_name, sequence in aligned_sequences.items()
    }


def missing_ambig_ratio(aligned_block: dict[str, str]) -> float | None:
    """Return the fraction of feature cells that are gaps or non-canonical bases."""

    total = 0
    missing_or_ambiguous = 0
    for sequence in aligned_block.values():
        for base in sequence.upper():
            total += 1
            if base not in CANONICAL_BASES:
                missing_or_ambiguous += 1
    if total == 0:
        return None
    return missing_or_ambiguous / total


def discover_genbank_files(inputs: Iterable[str | Path]) -> list[Path]:
    """Expand paths and directories into a deterministic GenBank file list."""

    resolved: list[Path] = []
    seen: set[Path] = set()
    for raw_path in inputs:
        path = Path(raw_path).expanduser().resolve()
        if not path.exists():
            raise MarkerSeekError(f"Input path does not exist: {path}")
        if path.is_dir():
            matches = sorted(
                candidate
                for candidate in path.rglob("*")
                if candidate.is_file() and candidate.suffix.lower() in SUPPORTED_SUFFIXES
            )
            resolved.extend(matches)
        elif path.suffix.lower() in SUPPORTED_SUFFIXES:
            resolved.append(path)
        else:
            raise MarkerSeekError(f"Unsupported input file type: {path}")

    unique_paths: list[Path] = []
    for path in resolved:
        if path not in seen:
            unique_paths.append(path)
            seen.add(path)
    return unique_paths


class GenomeInput:
    """Loaded GenBank input with a stable sample identifier."""

    def __init__(self, path: Path, sample_name: str, record: SeqRecord) -> None:
        self.path = path
        self.sample_name = sample_name
        self.record = record


def load_genbank(path: Path, index: int) -> GenomeInput:
    records = list(SeqIO.parse(str(path), "genbank"))
    if not records:
        raise MarkerSeekError(f"No GenBank records were found in {path}")
    record = records[0]
    sample_name = sanitize_sample_name(path.stem, index)
    return GenomeInput(path=path, sample_name=sample_name, record=record)


def sanitize_sample_name(raw_name: str, index: int) -> str:
    cleaned = "".join(char if char.isalnum() or char in {"_", "-"} else "_" for char in raw_name)
    cleaned = cleaned.strip("_") or f"sample_{index}"
    return f"{cleaned}_{index}"


def choose_reference_record(
    genome_records: list[GenomeInput],
    reference: str | Path | None,
) -> GenomeInput:
    if reference is None:
        return genome_records[0]
    target = Path(reference).expanduser().resolve()
    for record in genome_records:
        if record.path == target:
            return record
    raise MarkerSeekError(f"Reference file is not part of the provided inputs: {target}")


def extract_atomic_features(record: SeqRecord) -> list[AnnotatedInterval]:
    """Flatten GenBank gene features into non-overlapping atomic intervals."""

    deduplicated: dict[tuple[str, int, int, int], AnnotatedInterval] = {}
    for feature in record.features:
        if feature.type not in FEATURE_PRIORITY:
            continue
        gene_name = (
            feature.qualifiers.get("gene", [None])[0]
            or feature.qualifiers.get("locus_tag", [None])[0]
            or feature.qualifiers.get("product", [feature.type])[0]
        )
        location_parts = feature.location.parts if isinstance(feature.location, CompoundLocation) else [feature.location]
        for index, part in enumerate(location_parts, start=1):
            start = int(part.start)
            end = int(part.end)
            if start == end:
                continue
            strand = int(part.strand or 0)
            feature_id = gene_name if len(location_parts) == 1 else f"{gene_name}_part{index}"
            candidate = AnnotatedInterval(
                feature_id=feature_id,
                label=gene_name,
                feature_type=feature.type,
                parent_gene=gene_name,
                start=start,
                end=end,
                strand=strand,
            )
            key = (gene_name, start, end, strand)
            existing = deduplicated.get(key)
            if existing is None or FEATURE_PRIORITY[candidate.feature_type] > FEATURE_PRIORITY[existing.feature_type]:
                deduplicated[key] = candidate
    return sorted(deduplicated.values(), key=lambda item: (item.start, item.end, item.label))


def extract_exon_intervals(record: SeqRecord) -> list[AnnotatedInterval]:
    """Flatten CDS/tRNA/rRNA exon parts into atomic intervals for shading.

    Unlike ``extract_atomic_features`` (which uses gene/tRNA/rRNA features and
    therefore reports a single contiguous span for spliced genes), this
    function reads CDS features so that multi-exon genes contribute one
    interval per exon (CompoundLocation parts). tRNA and rRNA features are
    also included so that ribosomal and tRNA bodies show up as exonic
    background where no overlapping CDS exists.
    """

    exon_types = {"CDS", "tRNA", "rRNA"}
    priority = {"CDS": 1, "tRNA": 2, "rRNA": 3}
    deduplicated: dict[tuple[str, int, int, int], AnnotatedInterval] = {}
    for feature in record.features:
        if feature.type not in exon_types:
            continue
        gene_name = (
            feature.qualifiers.get("gene", [None])[0]
            or feature.qualifiers.get("locus_tag", [None])[0]
            or feature.qualifiers.get("product", [feature.type])[0]
        )
        location_parts = feature.location.parts if isinstance(feature.location, CompoundLocation) else [feature.location]
        for index, part in enumerate(location_parts, start=1):
            start = int(part.start)
            end = int(part.end)
            if start == end:
                continue
            strand = int(part.strand or 0)
            feature_id = gene_name if len(location_parts) == 1 else f"{gene_name}_exon{index}"
            candidate = AnnotatedInterval(
                feature_id=feature_id,
                label=gene_name,
                feature_type=feature.type,
                parent_gene=gene_name,
                start=start,
                end=end,
                strand=strand,
            )
            key = (gene_name, start, end, strand)
            existing = deduplicated.get(key)
            if existing is None or priority[candidate.feature_type] < priority[existing.feature_type]:
                deduplicated[key] = candidate
    return sorted(deduplicated.values(), key=lambda item: (item.start, item.end, item.label))


def infer_regions(sequence: str, features: list[AnnotatedInterval]) -> list[RegionSegment]:
    """Infer LSC, SSC, and the two IR copies from reference annotations."""

    copies_by_gene: dict[str, list[AnnotatedInterval]] = defaultdict(list)
    for feature in features:
        copies_by_gene[feature.parent_gene].append(feature)

    candidate_pairs: list[tuple[AnnotatedInterval, AnnotatedInterval]] = []
    preferred_pairs: list[tuple[AnnotatedInterval, AnnotatedInterval]] = []
    for gene_name, copies in copies_by_gene.items():
        ordered = sorted({(item.start, item.end): item for item in copies}.values(), key=lambda item: item.start)
        if len(ordered) != 2:
            continue
        pair = (ordered[0], ordered[1])
        candidate_pairs.append(pair)
        if gene_name in KNOWN_IR_MARKERS:
            preferred_pairs.append(pair)

    pair_pool = preferred_pairs or candidate_pairs
    if len(pair_pool) < 4:
        raise MarkerSeekError(AUTOMATIC_IR_DETECTION_FAILED)

    candidate_boundaries: list[tuple[int, int, int, int]] = []
    coarse_candidates = realign_ir_cores(
        sequence,
        min(pair[0].start for pair in pair_pool),
        max(pair[0].end for pair in pair_pool),
        min(pair[1].start for pair in pair_pool),
        max(pair[1].end for pair in pair_pool),
    )
    candidate_boundaries.extend(
        collect_valid_ir_candidates(
            sequence,
            coarse_candidates,
        )
    )
    for irb_feature, ira_feature in pair_pool:
        candidate_boundaries.extend(
            collect_valid_ir_candidates(
                sequence,
                [
                    (irb_feature.start, irb_feature.end, ira_feature.start, ira_feature.end),
                ],
            )
        )

    if not candidate_boundaries:
        raise MarkerSeekError(
            "The inferred IR regions are implausibly short. Provide manual boundaries with --lsc/--irb/--ssc/--ira."
        )

    # Expand each duplicated marker pair independently, then keep the consensus
    # boundary set supported by the most markers.
    irb_start, irb_end, ira_start, ira_end = max(
        Counter(candidate_boundaries).items(),
        key=lambda item: (
            item[1],
            min(item[0][1] - item[0][0], item[0][3] - item[0][2]),
            -(item[0][0] + item[0][2]),
        ),
    )[0]

    genome_length = len(sequence)
    between_ir = RegionSegment("segment_between_ir", irb_end, ira_start)
    wrap_segment = RegionSegment("segment_wrap", ira_end, irb_start)
    if between_ir.length(genome_length) <= 0 or wrap_segment.length(genome_length) <= 0:
        raise MarkerSeekError(
            "Failed to assign LSC and SSC from the inferred IR boundaries. Provide manual boundaries instead."
        )

    if between_ir.length(genome_length) <= wrap_segment.length(genome_length):
        ssc = RegionSegment("SSC", between_ir.start, between_ir.end)
        lsc = RegionSegment("LSC", wrap_segment.start, wrap_segment.end)
    else:
        lsc = RegionSegment("LSC", between_ir.start, between_ir.end)
        ssc = RegionSegment("SSC", wrap_segment.start, wrap_segment.end)

    return [
        lsc,
        RegionSegment("IRb", irb_start, irb_end),
        ssc,
        RegionSegment("IRa", ira_start, ira_end),
    ]


def infer_regions_or_unpartitioned(sequence: str, features: list[AnnotatedInterval]) -> list[RegionSegment]:
    """Infer plastome regions, falling back only for references with no detectable IR."""

    try:
        return infer_regions(sequence, features)
    except MarkerSeekError as exc:
        if str(exc) != AUTOMATIC_IR_DETECTION_FAILED:
            raise
        return [RegionSegment(UNPARTITIONED_REGION_NAME, 0, len(sequence))]


def collect_valid_ir_candidates(
    sequence: str,
    seeds: list[tuple[int, int, int, int]],
) -> list[tuple[int, int, int, int]]:
    """Expand candidate IR seeds and keep only plausible interval pairs."""

    valid: list[tuple[int, int, int, int]] = []
    genome_length = len(sequence)
    for irb_start, irb_end, ira_start, ira_end in seeds:
        for candidate in (
            expand_ir_boundaries(sequence, irb_start, irb_end, ira_start, ira_end),
            expand_ir_boundaries_parallel(sequence, irb_start, irb_end, ira_start, ira_end),
        ):
            expanded_irb_start, expanded_irb_end, expanded_ira_start, expanded_ira_end = candidate
            if expanded_irb_end <= expanded_irb_start or expanded_ira_end <= expanded_ira_start:
                continue
            if min(expanded_irb_end - expanded_irb_start, expanded_ira_end - expanded_ira_start) < 1000:
                continue
            single_copy_lengths = (
                RegionSegment("segment_between_ir", expanded_irb_end, expanded_ira_start).length(genome_length),
                RegionSegment("segment_wrap", expanded_ira_end, expanded_irb_start).length(genome_length),
            )
            if min(single_copy_lengths) < 50:
                continue
            valid.append(candidate)
    return valid


def realign_ir_cores(
    sequence: str,
    irb_start: int,
    irb_end: int,
    ira_start: int,
    ira_end: int,
) -> list[tuple[int, int, int, int]]:
    """Realign rough IR cores before boundary expansion."""

    coarse_span = max(irb_end - irb_start, ira_end - ira_start)
    left_window_start = max(0, irb_start - coarse_span)
    left_window_end = min(len(sequence), irb_end + coarse_span)
    right_window_start = max(0, ira_start - coarse_span)
    right_window_end = min(len(sequence), ira_end + coarse_span)

    left_window = sequence[left_window_start:left_window_end].upper()
    right_window_rc = reverse_complement_text(sequence[right_window_start:right_window_end].upper())
    match = longest_common_substring_positions(left_window, right_window_rc)
    if match is None:
        return [(irb_start, irb_end, ira_start, ira_end)]

    left_match_start, right_match_start, match_size = match
    aligned_irb_start = left_window_start + left_match_start
    aligned_irb_end = aligned_irb_start + match_size
    aligned_ira_end = right_window_end - right_match_start
    aligned_ira_start = aligned_ira_end - match_size
    return [(aligned_irb_start, aligned_irb_end, aligned_ira_start, aligned_ira_end)]


def longest_common_substring_positions(left: str, right: str) -> tuple[int, int, int] | None:
    """Return the start offsets and size of the longest common substring."""

    if not left or not right:
        return None

    base = 911382323
    mask = (1 << 64) - 1

    def build_rolling_hash(text: str) -> tuple[list[int], list[int]]:
        prefix = [0]
        powers = [1]
        for char in text.encode("ascii"):
            prefix.append(((prefix[-1] * base) + char) & mask)
            powers.append((powers[-1] * base) & mask)
        return prefix, powers

    def get_hash(prefix: list[int], powers: list[int], start: int, end: int) -> int:
        return (prefix[end] - (prefix[start] * powers[end - start])) & mask

    left_prefix, left_powers = build_rolling_hash(left)
    right_prefix, right_powers = build_rolling_hash(right)

    def find_match(size: int) -> tuple[int, int] | None:
        if size <= 0:
            return 0, 0
        seen: dict[int, list[int]] = {}
        for left_start in range(len(left) - size + 1):
            seen.setdefault(get_hash(left_prefix, left_powers, left_start, left_start + size), []).append(left_start)
        for right_start in range(len(right) - size + 1):
            candidates = seen.get(get_hash(right_prefix, right_powers, right_start, right_start + size))
            if not candidates:
                continue
            right_slice = right[right_start : right_start + size]
            for left_start in candidates:
                if left[left_start : left_start + size] == right_slice:
                    return left_start, right_start
        return None

    best: tuple[int, int, int] | None = None
    low = 1
    high = min(len(left), len(right))
    while low <= high:
        mid = (low + high) // 2
        match = find_match(mid)
        if match is None:
            high = mid - 1
            continue
        best = (match[0], match[1], mid)
        low = mid + 1
    return best


def expand_ir_boundaries(
    sequence: str,
    irb_start: int,
    irb_end: int,
    ira_start: int,
    ira_end: int,
) -> tuple[int, int, int, int]:
    """Expand IR cores while sequence identity supports the mirrored boundaries."""

    uppercase = sequence.upper()
    while irb_start > 0 and ira_end < len(uppercase):
        left_base = uppercase[irb_start - 1]
        right_base = uppercase[ira_end]
        if left_base.translate(COMPLEMENT) != right_base:
            break
        irb_start -= 1
        ira_end += 1

    while irb_end < ira_start and ira_start > 0:
        left_base = uppercase[irb_end]
        right_base = uppercase[ira_start - 1]
        if right_base.translate(COMPLEMENT) != left_base:
            break
        irb_end += 1
        ira_start -= 1

    return irb_start, irb_end, ira_start, ira_end


def expand_ir_boundaries_parallel(
    sequence: str,
    irb_start: int,
    irb_end: int,
    ira_start: int,
    ira_end: int,
) -> tuple[int, int, int, int]:
    """Expand IR cores when duplicated annotations track the same-side offsets."""

    uppercase = sequence.upper()
    while irb_start > 0 and ira_start > 0:
        left_base = uppercase[irb_start - 1]
        right_base = uppercase[ira_start - 1]
        if right_base.translate(COMPLEMENT) != left_base:
            break
        irb_start -= 1
        ira_start -= 1

    while irb_end < len(uppercase) and ira_end < len(uppercase):
        left_base = uppercase[irb_end]
        right_base = uppercase[ira_end]
        if right_base.translate(COMPLEMENT) != left_base:
            break
        irb_end += 1
        ira_end += 1

    return irb_start, irb_end, ira_start, ira_end


def parse_manual_regions(
    manual_regions: dict[str, tuple[int, int]],
    genome_length: int,
) -> list[RegionSegment]:
    """Validate manual region overrides passed in from the CLI."""

    expected_keys = ("LSC", "IRb", "SSC", "IRa")
    if set(manual_regions) != set(expected_keys):
        raise MarkerSeekError("Manual regions must include LSC, IRb, SSC, and IRa.")

    segments = [RegionSegment(name=key, start=value[0], end=value[1]) for key, value in manual_regions.items()]
    total = sum(segment.length(genome_length) for segment in segments)
    if total != genome_length:
        raise MarkerSeekError("Manual regions must partition the full circular genome without gaps.")

    occupancy = [0] * genome_length
    for segment in segments:
        for start, end in segment.spans(genome_length):
            for position in range(start, end):
                occupancy[position] += 1
    if any(value != 1 for value in occupancy):
        raise MarkerSeekError("Manual regions overlap or leave uncovered positions.")
    return segments


def assign_regions_to_features(
    features: list[AnnotatedInterval],
    regions: list[RegionSegment],
    genome_length: int,
) -> list[AnnotatedInterval]:
    assigned: list[AnnotatedInterval] = []
    for feature in features:
        midpoint = feature.midpoint(genome_length)
        region_name = region_for_position(midpoint, regions, genome_length)
        assigned.append(replace(feature, region=region_name))
    return assigned


def region_for_position(position: int, regions: list[RegionSegment], genome_length: int) -> str:
    for region in regions:
        if region.contains(position, genome_length):
            return region.name
    raise MarkerSeekError(f"No region covers position {position + 1}.")


def build_feature_catalog(features: list[AnnotatedInterval], genome_length: int) -> list[AnnotatedInterval]:
    """Add intergenic spacers to the atomic gene catalogue."""

    catalog = list(features)
    ordered = sorted(features, key=lambda item: (item.start, item.end, item.feature_id))
    for index, current in enumerate(ordered):
        next_item = ordered[(index + 1) % len(ordered)]
        if index < len(ordered) - 1:
            if next_item.start <= current.end:
                continue
            start, end = current.end, next_item.start
        else:
            gap_length = (genome_length - current.end) + next_item.start
            if gap_length <= 0:
                continue
            start, end = current.end, next_item.start

        spacer = AnnotatedInterval(
            feature_id=f"IGS_{index + 1:03d}",
            label=f"{current.label}-{next_item.label}",
            feature_type="igs",
            parent_gene=f"{current.label}|{next_item.label}",
            start=start,
            end=end,
        )
        if spacer.length(genome_length) > 0:
            catalog.append(spacer)
    return sorted(catalog, key=lambda item: (item.start, item.end, item.feature_id))


def run_mafft_alignment(
    genome_records: list[GenomeInput],
    mafft_bin: str,
    *,
    mafft_threads: int | None = None,
) -> list[SeqRecord]:
    """Align all GenBank sequences with MAFFT."""

    if shutil.which(mafft_bin) is None and not Path(mafft_bin).exists():
        raise MarkerSeekError(
            f"MAFFT executable was not found: {mafft_bin}. Install MAFFT or pass --mafft-bin explicitly."
        )

    with tempfile.TemporaryDirectory(prefix="markerseek_") as tmpdir:
        fasta_path = Path(tmpdir) / "inputs.fasta"
        seq_records: list[SeqRecord] = []
        for genome_record in genome_records:
            copy = genome_record.record[:]
            copy.id = genome_record.sample_name
            copy.name = genome_record.sample_name
            copy.description = genome_record.path.name
            seq_records.append(copy)
        SeqIO.write(seq_records, fasta_path, "fasta")

        mafft_path = Path(mafft_bin)
        # Force FFT-NS-2 (progressive only) so the runtime stays roughly O(N log N)
        # for hundreds of plastomes; the default --auto switches to FFT-NS-i for
        # large N which is O(N^2) in the disttbfast stage and runs out of memory
        # on small servers.
        if mafft_path.exists() and mafft_path.suffix == ".py":
            command = [sys.executable, str(mafft_path), "--retree", "2", "--maxiterate", "0", str(fasta_path)]
        else:
            command = [mafft_bin]
            if mafft_threads is not None:
                command.extend(["--thread", str(mafft_threads)])
            command.extend(["--retree", "2", "--maxiterate", "0", str(fasta_path)])
        try:
            completed = subprocess.run(command, check=True, capture_output=True, text=True)
        except FileNotFoundError as exc:
            raise MarkerSeekError(f"MAFFT executable was not found: {mafft_bin}") from exc
        except subprocess.CalledProcessError as exc:
            stderr = exc.stderr.strip() or exc.stdout.strip()
            raise MarkerSeekError(f"MAFFT failed: {stderr}") from exc
        aligned = list(SeqIO.parse(StringIO(completed.stdout), "fasta"))
    if len(aligned) != len(genome_records):
        raise MarkerSeekError("MAFFT did not return one alignment entry per input sequence.")
    lengths = {len(record.seq) for record in aligned}
    if len(lengths) != 1:
        raise MarkerSeekError("MAFFT produced inconsistent alignment lengths.")
    return aligned


def project_alignment_to_reference(
    aligned_records: list[SeqRecord],
    reference_name: str,
) -> dict[str, str]:
    """Drop alignment columns where the reference has a gap."""

    aligned_map = {record.id: str(record.seq).upper() for record in aligned_records}
    if reference_name not in aligned_map:
        raise MarkerSeekError(f"Reference sequence {reference_name} is missing from the MAFFT alignment.")
    reference_alignment = aligned_map[reference_name]
    keep_indices = [index for index, base in enumerate(reference_alignment) if base not in {"-", "."}]
    projected = {}
    for name, sequence in aligned_map.items():
        projected[name] = "".join(sequence[index] for index in keep_indices)
    return projected


def compute_position_pi(projected_sequences: Iterable[str]) -> tuple[list[float | None], list[int]]:
    """Per-site nucleotide diversity (Nei & Li 1979) on reference-projected sequences.

    For each alignment column, count canonical bases A/C/G/T only — gaps ("-",
    ".") and IUPAC ambiguity codes (N, Y, R, W, S, K, M, B, D, H, V) are
    treated as missing data and excluded from the per-site sample size n.
    Sites with n < 2 cannot define a pair and are returned as ``None`` with
    valid_mask = 0; sites with n >= 2 are valid.

    For a valid site with allele counts c_i (Σ c_i = n) the per-site
    diversity is computed as the fraction of pairwise comparisons that
    differ::

        π_l = (C(n,2) − Σ C(c_i,2)) / C(n,2)
            = (n² − Σ c_i²) / (n(n−1))

    which is algebraically identical to Nei's (n/(n-1))·(1 − Σ p_i²) with
    p_i = c_i / n. Each site uses its own n, so missing data at one sample
    does not invalidate the whole column.

    Inputs must already be projected to reference coordinates (see
    ``project_alignment_to_reference``); columns where the reference carries
    a gap are dropped beforehand, so insertions present only in non-reference
    samples do not contribute to π.

    Returns
    -------
    position_pi : list[float | None]
        Per-site π for valid sites, ``None`` for sites with n < 2.
    valid_mask : list[int]
        1 for valid sites, 0 otherwise. Same length as the projected
        alignment.
    """

    sequences = [sequence.upper() for sequence in projected_sequences]
    if not sequences:
        return [], []

    sequence_length = len(sequences[0])
    for sequence in sequences:
        if len(sequence) != sequence_length:
            raise MarkerSeekError("Projected sequences do not share the same length.")

    position_pi: list[float | None] = []
    valid_mask: list[int] = []
    for column in zip(*sequences):
        counts: dict[str, int] = defaultdict(int)
        for base in column:
            if base in CANONICAL_BASES:
                counts[base] += 1
        sample_count = sum(counts.values())
        if sample_count < 2:
            position_pi.append(None)
            valid_mask.append(0)
            continue
        total_pairs = sample_count * (sample_count - 1) / 2
        same_pairs = sum(count * (count - 1) / 2 for count in counts.values())
        diff_pairs = total_pairs - same_pairs
        position_pi.append(diff_pairs / total_pairs)
        valid_mask.append(1)
    return position_pi, valid_mask


def build_prefix_arrays(position_pi: list[float | None], valid_mask: list[int]) -> tuple[list[float], list[int]]:
    prefix_pi = [0.0]
    prefix_valid = [0]
    for pi_value, valid in zip(position_pi, valid_mask):
        prefix_pi.append(prefix_pi[-1] + (pi_value or 0.0))
        prefix_valid.append(prefix_valid[-1] + valid)
    return prefix_pi, prefix_valid


def interval_summary(
    start: int,
    end: int,
    genome_length: int,
    prefix_pi: list[float],
    prefix_valid: list[int],
) -> tuple[float | None, int]:
    """Mean per-site π over an interval, as the standard window/feature estimator.

    Aggregates the per-site π values produced by ``compute_position_pi``
    across the half-open interval [start, end) on a circular genome
    (wrap-around handled when end <= start). The interval estimator is the
    unweighted arithmetic mean of valid-site π_l::

        π_interval = Σ_{l ∈ valid sites} π_l / (number of valid sites)

    Sites with n < 2 (missing or non-canonical bases) do not contribute to
    either the numerator or the denominator. This matches the convention
    used by DnaSP and most population-genetic tools. Note that each valid
    site is weighted equally regardless of how many samples were callable
    at that site, so a window where most sites had small n is comparable in
    weight to one where most sites had full n; this is a deliberate
    standard choice, not a bug.

    Returns ``(None, 0)`` if the interval contains no valid sites.
    """

    total_pi = 0.0
    total_valid = 0
    spans = [(start, end)] if end > start else [(start, genome_length), (0, end)]
    for span_start, span_end in spans:
        total_pi += prefix_pi[span_end] - prefix_pi[span_start]
        total_valid += prefix_valid[span_end] - prefix_valid[span_start]
    if total_valid == 0:
        return None, 0
    return total_pi / total_valid, total_valid


def window_starts(genome_length: int, window: int, step: int) -> list[int]:
    if genome_length <= window:
        return [0]
    starts = list(range(0, genome_length - window + 1, step))
    final_start = genome_length - window
    if starts[-1] != final_start:
        starts.append(final_start)
    return starts


def build_window_results(
    *,
    genome_length: int,
    window: int,
    step: int,
    prefix_pi: list[float],
    prefix_valid: list[int],
    regions: list[RegionSegment],
    labels: list[AnnotatedInterval],
) -> list[WindowResult]:
    results: list[WindowResult] = []
    for index, start in enumerate(window_starts(genome_length, window, step), start=1):
        end = min(start + window, genome_length)
        pi_value, valid_sites = interval_summary(start, end, genome_length, prefix_pi, prefix_valid)
        midpoint = start + ((end - start) // 2)
        label_name = choose_window_label(start, end, labels, genome_length)
        results.append(
            WindowResult(
                window_id=f"W{index:04d}",
                start=start + 1,
                end=end,
                midpoint=midpoint + 1,
                pi=pi_value,
                valid_sites=valid_sites,
                region=region_for_position(midpoint % genome_length, regions, genome_length),
                label_name=label_name,
            )
        )
    return results


def choose_window_label(
    window_start: int,
    window_end: int,
    labels: list[AnnotatedInterval],
    genome_length: int,
) -> str:
    best_label = ""
    best_overlap = -1
    for feature in labels:
        overlap = overlap_length(window_start, window_end, feature, genome_length)
        if overlap > best_overlap:
            best_overlap = overlap
            best_label = feature.label
    return best_label


def overlap_length(window_start: int, window_end: int, feature: AnnotatedInterval, genome_length: int) -> int:
    overlap = 0
    for span_start, span_end in feature.spans(genome_length):
        left = max(window_start, span_start)
        right = min(window_end, span_end)
        if right > left:
            overlap += right - left
    return overlap


def build_feature_results(
    catalog: list[AnnotatedInterval],
    *,
    genome_length: int,
    prefix_pi: list[float],
    prefix_valid: list[int],
) -> list[FeatureResult]:
    results: list[FeatureResult] = []
    for feature in catalog:
        pi_value, valid_sites = interval_summary(feature.start, feature.end, genome_length, prefix_pi, prefix_valid)
        results.append(
            FeatureResult(
                feature_id=feature.feature_id,
                feature_type=feature.feature_type,
                parent_gene=feature.parent_gene,
                label_name=feature.label,
                start=feature.start + 1,
                end=feature.end if feature.end > 0 else genome_length,
                strand=feature.strand,
                length_bp=feature.length(genome_length),
                region=feature.region,
                pi=pi_value,
                valid_sites=valid_sites,
                spans_origin=feature.wraps_origin,
            )
        )
    return results


def apply_hotspots(rows: list[WindowResult] | list[FeatureResult], mode: str, value: float) -> None:
    valid_rows = [row for row in rows if row.pi is not None]
    ordered = sorted(valid_rows, key=lambda item: (-float(item.pi), item.start, item.end))

    selected: set[int] = set()
    if mode == "threshold":
        selected = {id(row) for row in valid_rows if float(row.pi) >= value}
    elif mode == "top-n":
        top_n = max(1, int(value))
        selected = {id(row) for row in ordered[:top_n]}
    elif mode == "top-percent":
        count = max(1, math.ceil(len(ordered) * (value / 100.0)))
        selected = {id(row) for row in ordered[:count]}
    else:
        raise MarkerSeekError(f"Unsupported hotspot mode: {mode}")

    rank = 1
    for row in ordered:
        if id(row) in selected:
            row.is_hotspot = True
            row.hotspot_rank = rank
            rank += 1
