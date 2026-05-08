"""Feature-level diagnostic metrics for MarkerSeek."""

from __future__ import annotations

import hashlib
import math
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np

CANONICAL_BASES = frozenset("ACGT")


def _sequence_length(aligned_block: dict[str, str]) -> int:
    lengths = {len(sequence) for sequence in aligned_block.values()}
    if not lengths:
        return 0
    if len(lengths) != 1:
        raise ValueError("Aligned block sequences must have equal length.")
    return lengths.pop()


def count_variable_and_indel_sites(aligned_block: dict[str, str]) -> tuple[int, int]:
    """Count variable canonical-base columns and mixed gap/base indel sites."""

    _sequence_length(aligned_block)
    variable_sites = 0
    indel_sites = 0
    sequences = [sequence.upper() for sequence in aligned_block.values()]
    for column in zip(*(sequence for sequence in sequences), strict=False):
        canonical = {base for base in column if base in CANONICAL_BASES}
        if len(canonical) >= 2:
            variable_sites += 1
        if "-" in column and canonical:
            indel_sites += 1
    return variable_sites, indel_sites


def assign_haplotypes(aligned_block: dict[str, str]) -> dict[str, str]:
    """Group normalized sequence strings and assign H001-style haplotype ids."""

    groups: dict[str, list[str]] = defaultdict(list)
    first_seen: dict[str, int] = {}
    for index, (sample_name, sequence) in enumerate(aligned_block.items()):
        normalized = "".join(base if base in "ACGT-" else "N" for base in sequence.upper())
        digest = hashlib.blake2b(normalized.encode("ascii"), digest_size=16).hexdigest()
        groups[digest].append(sample_name)
        first_seen.setdefault(digest, index)

    ordered_digests = sorted(groups, key=lambda digest: (-len(groups[digest]), first_seen[digest]))
    labels = {digest: f"H{index:03d}" for index, digest in enumerate(ordered_digests, start=1)}
    return {sample_name: labels[digest] for digest in ordered_digests for sample_name in groups[digest]}


def species_resolution(haplotypes: dict[str, str], species_map: dict[str, str]) -> float:
    """Return the fraction of species whose haplotypes are not shared across species."""

    species_to_haplotypes: dict[str, set[str]] = defaultdict(set)
    haplotype_to_species: dict[str, set[str]] = defaultdict(set)
    for sample_name, haplotype_id in haplotypes.items():
        species = species_map.get(sample_name)
        if not species:
            continue
        species_to_haplotypes[species].add(haplotype_id)
        haplotype_to_species[haplotype_id].add(species)

    if not species_to_haplotypes:
        return 0.0
    if len(species_to_haplotypes) == 1:
        return 1.0

    resolved_species = 0
    for species, species_haplotypes in species_to_haplotypes.items():
        if all(len(haplotype_to_species[haplotype_id]) == 1 for haplotype_id in species_haplotypes):
            resolved_species += 1
    return resolved_species / len(species_to_haplotypes)


def species_specific_haplotype_ratio(haplotypes: dict[str, str], species_map: dict[str, str]) -> float:
    """Return the fraction of haplotypes observed in only one species."""

    haplotype_to_species: dict[str, set[str]] = defaultdict(set)
    for sample_name, haplotype_id in haplotypes.items():
        species = species_map.get(sample_name)
        if species:
            haplotype_to_species[haplotype_id].add(species)
    if not haplotype_to_species:
        return 0.0
    species_specific = sum(1 for species_set in haplotype_to_species.values() if len(species_set) == 1)
    return species_specific / len(haplotype_to_species)


def _p_distance(left: str, right: str) -> tuple[float | None, int]:
    valid_sites = 0
    differing_sites = 0
    for left_base, right_base in zip(left.upper(), right.upper(), strict=False):
        if left_base not in CANONICAL_BASES or right_base not in CANONICAL_BASES:
            continue
        valid_sites += 1
        if left_base != right_base:
            differing_sites += 1
    if valid_sites < 20:
        return None, valid_sites
    return differing_sites / valid_sites, valid_sites


def divergence_stats(aligned_block: dict[str, str], species_map: dict[str, str]) -> dict[str, float | None]:
    """Compute pairwise p-distance diagnostics among and within species."""

    _sequence_length(aligned_block)
    samples = list(aligned_block)
    species_counts = Counter(species_map.get(sample_name) for sample_name in samples if species_map.get(sample_name))
    has_within_species_pairs = any(count >= 2 for count in species_counts.values())
    pair_distances: dict[tuple[str, str], float] = {}
    interspecific: list[float] = []
    intraspecific: list[float] = []
    for left, right in combinations(samples, 2):
        distance, _valid_sites = _p_distance(aligned_block[left], aligned_block[right])
        if distance is None:
            continue
        pair_distances[(left, right)] = distance
        pair_distances[(right, left)] = distance
        if species_map.get(left) == species_map.get(right):
            intraspecific.append(distance)
        else:
            interspecific.append(distance)

    interspecific_mean = float(np.mean(interspecific)) if len(interspecific) >= 2 else None
    intraspecific_mean = float(np.mean(intraspecific)) if len(intraspecific) >= 2 else None
    if len(interspecific) >= 2 and len(intraspecific) >= 2:
        barcoding_gap = min(interspecific) - max(intraspecific)
    else:
        barcoding_gap = None

    nearest_neighbor_discrimination = (
        _nearest_neighbor_discrimination(samples, pair_distances, species_map)
        if has_within_species_pairs
        else None
    )
    misclassification_risk = (
        None if nearest_neighbor_discrimination is None else 1.0 - nearest_neighbor_discrimination
    )
    return {
        "interspecific": interspecific_mean,
        "intraspecific": intraspecific_mean,
        "nearest_neighbor_discrimination": nearest_neighbor_discrimination,
        "barcoding_gap": barcoding_gap,
        "misclassification_risk": misclassification_risk,
    }


def _nearest_neighbor_discrimination(
    samples: list[str],
    pair_distances: dict[tuple[str, str], float],
    species_map: dict[str, str],
) -> float | None:
    if len(samples) < 2:
        return None
    correct = 0
    for sample_name in samples:
        neighbours = [
            (distance, other_sample)
            for (left, other_sample), distance in pair_distances.items()
            if left == sample_name
        ]
        if not neighbours:
            continue
        nearest_distance = min(distance for distance, _sample in neighbours)
        nearest_samples = [other_sample for distance, other_sample in neighbours if distance == nearest_distance]
        if any(species_map.get(sample_name) == species_map.get(other_sample) for other_sample in nearest_samples):
            correct += 1
    return correct / len(samples)


def alignment_reliability_block(aligned_block: dict[str, str]) -> float:
    """Estimate the fraction of columns with low gap, ambiguity, and entropy."""

    length = _sequence_length(aligned_block)
    if length == 0:
        return 0.0
    sequences = [sequence.upper() for sequence in aligned_block.values()]
    reliable_columns: list[bool] = []
    max_entropy = math.log2(4) * 0.95
    sample_count = len(sequences)
    for column in zip(*(sequence for sequence in sequences), strict=False):
        counts = Counter(column)
        gap_ratio = counts["-"] / sample_count
        canonical_counts = np.array([counts[base] for base in "ACGT"], dtype=float)
        canonical_total = float(canonical_counts.sum())
        ambig_count = sample_count - counts["-"] - int(canonical_total)
        ambig_ratio = ambig_count / sample_count
        if canonical_total == 0:
            entropy = 0.0
        else:
            probabilities = canonical_counts[canonical_counts > 0] / canonical_total
            entropy = float(-(probabilities * np.log2(probabilities)).sum())
        reliable_columns.append(gap_ratio < 0.5 and ambig_ratio < 0.5 and entropy < max_entropy)
    return float(np.mean(reliable_columns))


def flanking_conservation(
    full_aligned: dict[str, str],
    region_start: int,
    region_end: int,
    max_scan: int = 300,
) -> tuple[int, int]:
    """Measure conserved contiguous flanks immediately left and right of a region."""

    genome_length = _sequence_length(full_aligned)
    if genome_length == 0 or max_scan <= 0:
        return 0, 0
    return flanking_conservation_from_mask(
        conserved_column_mask(full_aligned),
        region_start=region_start,
        region_end=region_end,
        genome_length=genome_length,
        max_scan=max_scan,
    )


def flanking_conservation_from_mask(
    conserved: list[bool],
    region_start: int,
    region_end: int,
    genome_length: int,
    max_scan: int = 300,
) -> tuple[int, int]:
    """Measure conserved flanks from a precomputed conserved-column mask."""

    if genome_length == 0 or max_scan <= 0:
        return 0, 0
    scan_limit = min(max_scan, genome_length, len(conserved))
    start = region_start % genome_length
    end = region_end % genome_length
    left_bp = 0
    for offset in range(1, scan_limit + 1):
        position = (start - offset) % genome_length
        if not conserved[position]:
            break
        left_bp += 1

    right_bp = 0
    for offset in range(scan_limit):
        position = (end + offset) % genome_length
        if not conserved[position]:
            break
        right_bp += 1
    return left_bp, right_bp


def conserved_column_mask(full_aligned: dict[str, str]) -> list[bool]:
    """Return a per-column mask where the dominant canonical base reaches 90%."""

    _sequence_length(full_aligned)
    sequences = [sequence.upper() for sequence in full_aligned.values()]
    conserved: list[bool] = []
    for column in zip(*(sequence for sequence in sequences), strict=False):
        canonical_counts = Counter(base for base in column if base in CANONICAL_BASES)
        valid_bases = sum(canonical_counts.values())
        if valid_bases == 0:
            conserved.append(False)
            continue
        conserved.append(max(canonical_counts.values()) / valid_bases >= 0.9)
    return conserved
