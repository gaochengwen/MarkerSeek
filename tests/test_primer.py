from __future__ import annotations

from markerseek.models import PrimerResult
from markerseek.primer import AmpliconHit, PrimerPair, fuzzy_search, in_silico_pcr, reverse_complement, score_primer_pair


def _mutate(sequence: str, position: int) -> str:
    replacement = "A" if sequence[position] != "A" else "C"
    return f"{sequence[:position]}{replacement}{sequence[position + 1:]}"


def _pair() -> PrimerPair:
    return PrimerPair(
        fwd_seq="ACGTTGCACTGATCGTAC",
        rev_seq="TTAACCGGATGCCGTAAC",
        fwd_gc=50.0,
        rev_gc=50.0,
        fwd_tm=60.0,
        rev_tm=60.0,
        fwd_self_any_th=0.0,
        rev_self_any_th=0.0,
        primer3_penalty=1.0,
        amplicon_mean_len=100.0,
        amplicon_variable_sites=10,
        amplicon_indel_sites=0,
    )


def _hit(sample_name: str) -> AmpliconHit:
    return AmpliconHit(
        sample_name=sample_name,
        amplicon_seq="A" * 100,
        fwd_pos=0,
        rev_pos=82,
        length=100,
        fwd_mismatches=0,
        rev_mismatches=0,
        multiple_hits=False,
    )


def test_fuzzy_search_anchors_three_prime_end() -> None:
    primer = "ACGTTGCACTGATCGTAC"
    template = f"AAAAAA{primer}TTTTTT"

    anchor_mutated = template
    for offset in range(13, 18):
        anchor_mutated = _mutate(anchor_mutated, 6 + offset)
    body_one_mutation = _mutate(template, 6 + 3)
    body_two_mutations = _mutate(_mutate(template, 6 + 3), 6 + 8)

    assert fuzzy_search(primer, anchor_mutated, max_mismatch=1, anchor_3prime=5) == []
    assert fuzzy_search(primer, body_one_mutation, max_mismatch=1, anchor_3prime=5) == [(6, 1)]
    assert fuzzy_search(primer, body_two_mutations, max_mismatch=1, anchor_3prime=5) == []


def test_in_silico_pcr_returns_amplicon_with_correct_length_distribution() -> None:
    forward = "ACGTTGCACTGATCGTAC"
    reverse = "TTAACCGGATGCCGTAAC"
    reverse_site = reverse_complement(reverse)
    sample_seqs = {
        "s1": f"GGGGG{forward}{'A' * 70}{reverse_site}CCCCC",
        "s2": f"GGGGG{forward}{'C' * 70}{reverse_site}CCCCC",
        "s3": f"GGGGG{forward}{'G' * 70}{reverse_site}CCCCC",
        "s4": f"GGGGG{forward}{'T' * 120}{reverse_site}CCCCC",
    }

    hits = in_silico_pcr(forward, reverse, sample_seqs, min_amplicon=80, max_amplicon=300)
    lengths = [hit.length for hit in hits.values() if hit is not None]

    assert len(lengths) == 4
    assert max(lengths) - min(lengths) == 50


def test_in_silico_pcr_no_match_returns_none() -> None:
    hits = in_silico_pcr(
        "ACGTTGCACTGATCGTAC",
        "TTAACCGGATGCCGTAAC",
        {"missing": "G" * 200},
        min_amplicon=80,
        max_amplicon=300,
    )

    assert hits["missing"] is None


def test_score_pair_zero_when_cross_species_zero() -> None:
    pair = _pair()

    assert score_primer_pair(pair, {}, total_samples=4) == 0.0


def test_score_pair_rewards_high_cross_species() -> None:
    low_pair = _pair()
    high_pair = _pair()
    low_amplicons = {f"s{index}": (_hit(f"s{index}") if index < 3 else None) for index in range(10)}
    high_amplicons = {f"s{index}": _hit(f"s{index}") for index in range(10)}

    low_score = score_primer_pair(low_pair, low_amplicons, total_samples=10)
    high_score = score_primer_pair(high_pair, high_amplicons, total_samples=10)

    assert high_score > low_score


def test_primer_id_format() -> None:
    result = PrimerResult(
        label_name="trnH-GUG-psbA",
        rank=2,
        fwd_seq="ACGTACGTACGTACGTAC",
        rev_seq="TGCATGCATGCATGCATG",
        fwd_len=18,
        rev_len=18,
        fwd_gc=50.0,
        rev_gc=50.0,
        fwd_tm=58.0,
        rev_tm=58.0,
        fwd_self_any_th=0.0,
        rev_self_any_th=0.0,
        primer3_penalty=0.1,
        target_start=1,
        target_end=200,
        amplicon_min_len=100,
        amplicon_max_len=120,
        amplicon_mean_len=110.0,
        cross_species_success_rate=1.0,
        amplicon_variable_sites=5,
        amplicon_indel_sites=0,
        primer_score=90.0,
    )

    assert result.primer_id == "trnH-GUG-psbA_p2"
