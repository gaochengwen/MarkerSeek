"""Primer design and in silico PCR helpers."""

from __future__ import annotations

from dataclasses import dataclass

CANONICAL_BASES = frozenset("ACGT")
DNA_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


@dataclass(frozen=True)
class FlankCandidate:
    start: int
    end: int
    length: int
    identity: float


@dataclass
class PrimerPair:
    fwd_seq: str
    rev_seq: str
    fwd_gc: float
    rev_gc: float
    fwd_tm: float
    rev_tm: float
    fwd_self_any_th: float
    rev_self_any_th: float
    primer3_penalty: float
    amplicon_min_len: int = 0
    amplicon_max_len: int = 0
    amplicon_mean_len: float = 0.0
    cross_species_success_rate: float = 0.0
    amplicon_variable_sites: int = 0
    amplicon_indel_sites: int = 0
    primer_score: float = 0.0
    rank: int = 0
    target_start: int = 0
    target_end: int = 0

    @property
    def fwd_len(self) -> int:
        return len(self.fwd_seq)

    @property
    def rev_len(self) -> int:
        return len(self.rev_seq)


@dataclass(frozen=True)
class AmpliconHit:
    sample_name: str
    amplicon_seq: str
    fwd_pos: int
    rev_pos: int
    length: int
    fwd_mismatches: int
    rev_mismatches: int
    multiple_hits: bool


def find_conserved_flanks(
    aligned_block: dict[str, str],
    conserved_mask: list[bool],
    target_start: int,
    target_end: int,
    *,
    max_scan: int = 300,
    min_len: int = 18,
    max_len: int = 24,
    identity_threshold: float = 0.95,
) -> tuple[FlankCandidate | None, FlankCandidate | None]:
    """Find conserved flanks immediately outside a target interval."""

    genome_length = _aligned_length(aligned_block)
    if (
        genome_length == 0
        or len(conserved_mask) < genome_length
        or max_scan <= 0
        or min_len <= 0
        or max_len < min_len
        or target_end <= target_start
    ):
        return None, None

    target_start = max(0, min(target_start, genome_length))
    target_end = max(0, min(target_end, genome_length))
    left = _scan_left_flank(
        aligned_block,
        conserved_mask,
        target_start,
        max_scan=max_scan,
        min_len=min_len,
        max_len=max_len,
        identity_threshold=identity_threshold,
    )
    right = _scan_right_flank(
        aligned_block,
        conserved_mask,
        target_end,
        genome_length,
        max_scan=max_scan,
        min_len=min_len,
        max_len=max_len,
        identity_threshold=identity_threshold,
    )
    return left, right


def design_primers_for_region(
    reference_seq: str,
    target_start: int,
    target_end: int,
    left_flank: FlankCandidate,
    right_flank: FlankCandidate,
    *,
    max_pairs: int = 5,
    primer_settings: dict | None = None,
) -> list[PrimerPair]:
    """Design primer pairs around a target region with primer3-py."""

    if target_end <= target_start or left_flank.end > target_start or right_flank.start < target_end:
        return []
    template_start = max(0, left_flank.start)
    template_end = min(len(reference_seq), right_flank.end)
    if template_end <= template_start:
        return []

    seq_args = {
        "SEQUENCE_ID": f"target_{target_start}_{target_end}",
        "SEQUENCE_TEMPLATE": reference_seq[template_start:template_end].upper(),
        "SEQUENCE_INCLUDED_REGION": [0, template_end - template_start],
        "SEQUENCE_TARGET": [target_start - template_start, target_end - target_start],
    }
    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MIN_TM": 55.0,
        "PRIMER_MAX_TM": 65.0,
        "PRIMER_OPT_GC_PERCENT": 50.0,
        "PRIMER_MIN_GC": 35.0,
        "PRIMER_MAX_GC": 65.0,
        "PRIMER_PRODUCT_SIZE_RANGE": [[80, 3000]],
        "PRIMER_NUM_RETURN": max_pairs,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
    }
    if primer_settings:
        global_args.update(primer_settings)

    primer3_bindings = _primer3_bindings()
    design_func = getattr(primer3_bindings, "design_primers", None) or getattr(
        primer3_bindings, "designPrimers"
    )
    raw_result = design_func(seq_args, global_args)
    returned = int(raw_result.get("PRIMER_PAIR_NUM_RETURNED", 0) or 0)

    pairs: list[PrimerPair] = []
    for index in range(returned):
        left_seq = str(raw_result.get(f"PRIMER_LEFT_{index}_SEQUENCE", "")).upper()
        right_seq = str(raw_result.get(f"PRIMER_RIGHT_{index}_SEQUENCE", "")).upper()
        if not left_seq or not right_seq:
            continue
        pairs.append(
            PrimerPair(
                fwd_seq=left_seq,
                rev_seq=right_seq,
                fwd_gc=round(float(raw_result.get(f"PRIMER_LEFT_{index}_GC_PERCENT", _gc_percent(left_seq))), 1),
                rev_gc=round(float(raw_result.get(f"PRIMER_RIGHT_{index}_GC_PERCENT", _gc_percent(right_seq))), 1),
                fwd_tm=round(float(raw_result.get(f"PRIMER_LEFT_{index}_TM", 0.0) or 0.0), 1),
                rev_tm=round(float(raw_result.get(f"PRIMER_RIGHT_{index}_TM", 0.0) or 0.0), 1),
                fwd_self_any_th=round(_primer3_metric(raw_result, f"PRIMER_LEFT_{index}_SELF_ANY"), 1),
                rev_self_any_th=round(_primer3_metric(raw_result, f"PRIMER_RIGHT_{index}_SELF_ANY"), 1),
                primer3_penalty=float(raw_result.get(f"PRIMER_PAIR_{index}_PENALTY", 0.0) or 0.0),
                target_start=target_start,
                target_end=target_end,
            )
        )
    return pairs


def fuzzy_search(
    primer_seq: str,
    target_seq: str,
    *,
    max_mismatch: int = 1,
    anchor_3prime: int = 5,
) -> list[tuple[int, int]]:
    """Return forward-strand matches allowing body mismatches outside the 3' anchor."""

    return _fuzzy_search(primer_seq, target_seq, max_mismatch=max_mismatch, anchor_3prime=anchor_3prime)


def in_silico_pcr(
    forward: str,
    reverse: str,
    sample_seqs: dict[str, str],
    *,
    min_amplicon: int = 80,
    max_amplicon: int = 3000,
    max_mismatch: int = 1,
    anchor_3prime: int = 5,
) -> dict[str, AmpliconHit | None]:
    """Search each sample for a single valid plus-strand amplicon."""

    results: dict[str, AmpliconHit | None] = {}
    reverse_plus = reverse_complement(reverse)
    reverse_len = len(reverse_plus)

    for sample_name, raw_sequence in sample_seqs.items():
        sequence = raw_sequence.upper().replace("-", "").replace(".", "")
        fwd_hits = fuzzy_search(
            forward,
            sequence,
            max_mismatch=max_mismatch,
            anchor_3prime=anchor_3prime,
        )
        rev_hits = _fuzzy_search(
            reverse_plus,
            sequence,
            max_mismatch=max_mismatch,
            anchor_3prime=anchor_3prime,
            anchor_side="left",
        )

        valid_pairs: list[tuple[int, int, int, int, int]] = []
        for fwd_pos, fwd_mismatches in fwd_hits:
            for rev_pos, rev_mismatches in rev_hits:
                if rev_pos <= fwd_pos:
                    continue
                amplicon_end = rev_pos + reverse_len
                length = amplicon_end - fwd_pos
                if min_amplicon <= length <= max_amplicon:
                    valid_pairs.append((length, fwd_pos, rev_pos, fwd_mismatches, rev_mismatches))

        if not valid_pairs:
            results[sample_name] = None
            continue

        valid_pairs.sort(key=lambda item: (item[0], item[1], item[2]))
        length, fwd_pos, rev_pos, fwd_mismatches, rev_mismatches = valid_pairs[0]
        results[sample_name] = AmpliconHit(
            sample_name=sample_name,
            amplicon_seq=sequence[fwd_pos : rev_pos + reverse_len],
            fwd_pos=fwd_pos,
            rev_pos=rev_pos,
            length=length,
            fwd_mismatches=fwd_mismatches,
            rev_mismatches=rev_mismatches,
            multiple_hits=len(valid_pairs) > 1,
        )
    return results


def score_primer_pair(
    pair: PrimerPair,
    amplicons: dict[str, AmpliconHit | None],
    total_samples: int,
) -> float:
    """Score a primer pair on primer3 quality, amplification rate, and information."""

    successful = [hit for hit in amplicons.values() if hit is not None and not hit.multiple_hits]
    cross_species = len(successful) / total_samples if total_samples else 0.0
    pair.cross_species_success_rate = cross_species
    if cross_species == 0:
        pair.primer_score = 0.0
        return 0.0

    penalty_term = 1 - _clip(pair.primer3_penalty / 5.0, 0, 1)
    var_density = pair.amplicon_variable_sites / max(pair.amplicon_mean_len, 1)
    indel_density = pair.amplicon_indel_sites / max(pair.amplicon_mean_len, 1)
    info_term = (0.4 * _clip(var_density / 0.10, 0, 1)) + (
        0.2 * _clip(indel_density / 0.05, 0, 1)
    )
    primer_score = 100 * ((0.4 * penalty_term) + (0.4 * cross_species) + (0.2 * info_term))
    pair.primer_score = round(_clip(primer_score, 0, 100), 1)
    return pair.primer_score


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(DNA_COMPLEMENT)[::-1]


def _primer3_bindings():
    from primer3 import bindings

    return bindings


def _primer3_metric(raw_result: dict, base_key: str) -> float:
    return float(
        raw_result.get(f"{base_key}_TH", raw_result.get(base_key, 0.0)) or 0.0
    )


def _aligned_length(aligned_block: dict[str, str]) -> int:
    lengths = {len(sequence) for sequence in aligned_block.values()}
    if not lengths:
        return 0
    if len(lengths) != 1:
        raise ValueError("Aligned sequences must have equal length.")
    return lengths.pop()


def _scan_left_flank(
    aligned_block: dict[str, str],
    conserved_mask: list[bool],
    target_start: int,
    *,
    max_scan: int,
    min_len: int,
    max_len: int,
    identity_threshold: float,
) -> FlankCandidate | None:
    scan_start = target_start - 1
    scan_stop = max(-1, target_start - max_scan - 1)
    first_conserved = None
    for position in range(scan_start, scan_stop, -1):
        if conserved_mask[position]:
            first_conserved = position
            break
    if first_conserved is None:
        return None

    start = first_conserved
    while start > 0 and conserved_mask[start - 1] and first_conserved - start + 1 < max_len:
        start -= 1
    end = first_conserved + 1
    return _make_flank(aligned_block, start, end, min_len=min_len, identity_threshold=identity_threshold)


def _scan_right_flank(
    aligned_block: dict[str, str],
    conserved_mask: list[bool],
    target_end: int,
    genome_length: int,
    *,
    max_scan: int,
    min_len: int,
    max_len: int,
    identity_threshold: float,
) -> FlankCandidate | None:
    scan_stop = min(genome_length, target_end + max_scan)
    first_conserved = None
    for position in range(target_end, scan_stop):
        if conserved_mask[position]:
            first_conserved = position
            break
    if first_conserved is None:
        return None

    end = first_conserved + 1
    while end < genome_length and conserved_mask[end] and end - first_conserved < max_len:
        end += 1
    return _make_flank(
        aligned_block,
        first_conserved,
        end,
        min_len=min_len,
        identity_threshold=identity_threshold,
    )


def _make_flank(
    aligned_block: dict[str, str],
    start: int,
    end: int,
    *,
    min_len: int,
    identity_threshold: float,
) -> FlankCandidate | None:
    length = end - start
    if length < min_len:
        return None
    identity = _mean_column_identity(aligned_block, start, end)
    if identity < identity_threshold:
        return None
    return FlankCandidate(start=start, end=end, length=length, identity=identity)


def _mean_column_identity(aligned_block: dict[str, str], start: int, end: int) -> float:
    identities: list[float] = []
    sequences = [sequence.upper() for sequence in aligned_block.values()]
    for position in range(start, end):
        counts: dict[str, int] = {}
        for sequence in sequences:
            base = sequence[position]
            if base in CANONICAL_BASES:
                counts[base] = counts.get(base, 0) + 1
        valid = sum(counts.values())
        if valid == 0:
            identities.append(0.0)
        else:
            identities.append(max(counts.values()) / valid)
    return sum(identities) / len(identities) if identities else 0.0


def _fuzzy_search(
    primer_seq: str,
    target_seq: str,
    *,
    max_mismatch: int,
    anchor_3prime: int,
    anchor_side: str = "right",
) -> list[tuple[int, int]]:
    query = primer_seq.upper().encode("ascii")
    target = target_seq.upper().encode("ascii")
    query_len = len(query)
    target_len = len(target)
    if query_len == 0 or target_len < query_len:
        return []

    anchor_len = min(max(anchor_3prime, 0), query_len)
    if anchor_side == "left":
        anchor_indices = range(0, anchor_len)
        body_indices = range(anchor_len, query_len)
    else:
        anchor_indices = range(query_len - anchor_len, query_len)
        body_indices = range(0, query_len - anchor_len)

    hits: list[tuple[int, int]] = []
    for start in range(target_len - query_len + 1):
        if not _anchor_matches(query, target, start, anchor_indices):
            continue
        mismatches = _body_mismatches(query, target, start, body_indices, max_mismatch)
        if mismatches <= max_mismatch:
            hits.append((start, mismatches))
    return hits


def _anchor_matches(query: bytes, target: bytes, start: int, indices: range) -> bool:
    for index in indices:
        target_base = target[start + index]
        if target_base != ord("N") and target_base != query[index]:
            return False
    return True


def _body_mismatches(
    query: bytes,
    target: bytes,
    start: int,
    indices: range,
    max_mismatch: int,
) -> int:
    mismatches = 0
    for index in indices:
        target_base = target[start + index]
        if target_base != ord("N") and target_base != query[index]:
            mismatches += 1
            if mismatches > max_mismatch:
                break
    return mismatches


def _gc_percent(sequence: str) -> float:
    if not sequence:
        return 0.0
    uppercase = sequence.upper()
    return 100.0 * sum(1 for base in uppercase if base in {"G", "C"}) / len(uppercase)


def _clip(value: float, lower: float, upper: float) -> float:
    return max(lower, min(value, upper))
