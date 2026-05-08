from __future__ import annotations

from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from markerseek import analysis
from markerseek.models import AnalysisResult, FeatureResult, PrimerResult
from markerseek.primer import (
    AmpliconHit,
    FlankCandidate,
    PrimerPair,
    fuzzy_search,
    in_silico_pcr,
    reverse_complement,
    score_primer_pair,
)


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


def test_primer_targets_keep_full_feature_first_for_short_hotspot() -> None:
    feature = FeatureResult(
        feature_id="IGS_001",
        feature_type="igs",
        parent_gene="rpl32|trnL-UAG",
        label_name="rpl32-trnL-UAG",
        start=1001,
        end=1650,
        strand=0,
        length_bp=650,
        region="SSC",
        pi=0.05,
        valid_sites=650,
        spans_origin=False,
        is_hotspot=True,
    )

    targets = analysis._primer_target_candidates(
        feature,
        [0.02] * 3000,
        max_amplicon=3000,
        primer_min_len=18,
    )

    assert targets[0].full_feature is True
    assert (targets[0].start, targets[0].end) == (1000, 1650)


def test_primer_pairs_fall_back_to_broad_template_without_conserved_flanks(monkeypatch) -> None:
    feature = FeatureResult(
        feature_id="IGS_001",
        feature_type="igs",
        parent_gene="rpl32|trnL-UAG",
        label_name="rpl32-trnL-UAG",
        start=1001,
        end=1650,
        strand=0,
        length_bp=650,
        region="SSC",
        pi=0.05,
        valid_sites=650,
        spans_origin=False,
        is_hotspot=True,
    )
    target = analysis.PrimerTargetCandidate(1000, 1650, 0.05, full_feature=True)
    template_calls: list[tuple[int, int, int, int]] = []

    monkeypatch.setattr(analysis, "find_conserved_flanks", lambda *args, **kwargs: (None, None))

    def fake_template_design(
        reference_seq,
        template_start,
        template_end,
        target_start,
        target_end,
        *,
        max_pairs=5,
        primer_settings=None,
    ):
        template_calls.append((template_start, template_end, target_start, target_end))
        return [
            PrimerPair(
                fwd_seq="ACGTTGCACTGATCGTAC",
                rev_seq="TTAACCGGATGCCGTAAC",
                fwd_gc=50.0,
                rev_gc=50.0,
                fwd_tm=60.0,
                rev_tm=60.0,
                fwd_self_any_th=0.0,
                rev_self_any_th=0.0,
                primer3_penalty=1.0,
                target_start=target_start,
                target_end=target_end,
            )
        ]

    monkeypatch.setattr(analysis, "design_primers_for_template", fake_template_design)

    pairs = analysis._primer_pairs_for_target(
        "A" * 3000,
        {"ref": "A" * 3000, "sample": "A" * 3000},
        [False] * 3000,
        feature,
        target,
        primer_min_len=18,
        flank_max_len=120,
        genome_length=3000,
        max_amplicon=3000,
        max_pairs=5,
        strict_primer_settings=analysis._screening_primer_settings({"PRIMER_NUM_RETURN": 5}, 5, relaxed=False),
        relaxed_primer_settings=analysis._screening_primer_settings({"PRIMER_NUM_RETURN": 5}, 5, relaxed=True),
    )

    assert pairs
    assert template_calls
    assert template_calls[0][0] < target.start
    assert template_calls[0][1] > target.end
    assert (pairs[0].target_start, pairs[0].target_end) == (target.start, target.end)


def test_primer_design_uses_high_pi_subregion_for_long_hotspot(monkeypatch) -> None:
    feature = FeatureResult(
        feature_id="ycf1",
        feature_type="gene",
        parent_gene="ycf1",
        label_name="ycf1",
        start=101,
        end=5100,
        strand=1,
        length_bp=5000,
        region="SSC",
        pi=0.01,
        valid_sites=5000,
        spans_origin=False,
        is_hotspot=True,
        hotspot_rank=1,
    )
    position_pi = [0.001] * 6000
    for position in range(2400, 2600):
        position_pi[position] = 0.2
    result = AnalysisResult(
        reference_name="ref",
        genome_length=6000,
        sample_count=2,
        regions=[],
        position_pi=position_pi,
        windows=[],
        features=[feature],
        aligned_sequences={"ref": "A" * 6000, "sample": "A" * 6000},
        sample_order=["ref", "sample"],
    )
    raw_records = [
        analysis.GenomeInput(Path("ref.gb"), "ref", SeqRecord(Seq("A" * 6000), id="ref")),
        analysis.GenomeInput(Path("sample.gb"), "sample", SeqRecord(Seq("A" * 6000), id="sample")),
    ]
    observed_targets: list[tuple[int, int]] = []

    def fake_find_flanks(aligned_sequences, conserved_mask, target_start, target_end, **kwargs):
        observed_targets.append((target_start, target_end))
        return (
            FlankCandidate(target_start - 24, target_start, 24, 1.0),
            FlankCandidate(target_end, target_end + 24, 24, 1.0),
        )

    def fake_design_primers(
        reference_seq,
        target_start,
        target_end,
        left_flank,
        right_flank,
        *,
        max_pairs=5,
        primer_settings=None,
    ):
        return [
            PrimerPair(
                fwd_seq="ACGTTGCACTGATCGTAC",
                rev_seq="TTAACCGGATGCCGTAAC",
                fwd_gc=50.0,
                rev_gc=50.0,
                fwd_tm=60.0,
                rev_tm=60.0,
                fwd_self_any_th=0.0,
                rev_self_any_th=0.0,
                primer3_penalty=1.0,
                target_start=target_start,
                target_end=target_end,
            )
        ]

    def fake_pcr(forward, reverse, sample_seqs, **kwargs):
        return {
            sample_name: AmpliconHit(sample_name, "A" * 300, 0, 282, 300, 0, 0, False)
            for sample_name in sample_seqs
        }

    monkeypatch.setattr(analysis, "find_conserved_flanks", fake_find_flanks)
    monkeypatch.setattr(analysis, "design_primers_for_region", fake_design_primers)
    monkeypatch.setattr(analysis, "in_silico_pcr", fake_pcr)

    analysis.design_primers_for_hotspots(
        result,
        raw_records,
        {"PRIMER_MIN_SIZE": 18, "PRIMER_NUM_RETURN": 1},
        {"min_amplicon": 80, "max_amplicon": 500},
        conserved_mask=[True] * 6000,
    )

    assert result.primers
    primer = result.primers[0]
    assert primer.target_end - primer.target_start + 1 <= 464
    assert primer.target_start <= 2500 <= primer.target_end
    assert all(target_end - target_start <= 464 for target_start, target_end in observed_targets)
