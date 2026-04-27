"""Core analysis pipeline for MarkerSeek."""

from __future__ import annotations

import math
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import replace
from io import StringIO
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
from Bio.SeqRecord import SeqRecord

from .models import AnalysisResult, AnnotatedInterval, FeatureResult, RegionSegment, WindowResult

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
    manual_regions: dict[str, tuple[int, int]] | None = None,
) -> AnalysisResult:
    """Run the full GenBank-to-Pi workflow and return structured results."""

    if window <= 0:
        raise MarkerSeekError("--window must be a positive integer.")
    if step <= 0:
        raise MarkerSeekError("--step must be a positive integer.")

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
        else infer_regions(reference_sequence, atomic_features)
    )
    atomic_features = assign_regions_to_features(atomic_features, regions, genome_length)
    feature_catalog = build_feature_catalog(atomic_features, genome_length)
    feature_catalog = assign_regions_to_features(feature_catalog, regions, genome_length)

    aligned_records = run_mafft_alignment(genome_records, mafft_bin)
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
    apply_hotspots(features, hotspot_mode, hotspot_value)

    sample_order = [reference_record.sample_name] + [
        record.sample_name
        for record in genome_records
        if record.sample_name != reference_record.sample_name
    ]

    return AnalysisResult(
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
        raise MarkerSeekError(
            "Automatic IR detection failed. Provide --lsc/--irb/--ssc/--ira to override the region boundaries."
        )

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


def run_mafft_alignment(genome_records: list[GenomeInput], mafft_bin: str) -> list[SeqRecord]:
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
        if mafft_path.exists() and mafft_path.suffix == ".py":
            command = [sys.executable, str(mafft_path), "--auto", str(fasta_path)]
        else:
            command = [mafft_bin, "--auto", str(fasta_path)]
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
