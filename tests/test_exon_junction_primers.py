from __future__ import annotations

import pytest
from Bio.Seq import Seq

import bioinformatics
from bioinformatics import (
    PrimerDesignParams,
    generate_exon_junction_primer_pairs,
    generate_primer_pairs,
)
from models import ExonInterval


REFERENCE = (
    "CACCGTGGAGCGACTACTTATCCGAGCCCT"
    "ATGGCACAGCAGCACTGTTGCTGTAGGCA"
    "ACGCTGGGGCATGAACACTACTCATCTTTAT"
)
EXONS = [
    ExonInterval(1, 30, "1"),
    ExonInterval(31, 60, "2"),
    ExonInterval(61, 90, "3"),
]


def _params(**overrides: object) -> PrimerDesignParams:
    values: dict[str, object] = {
        "identity_threshold": 1.0,
        "coverage_threshold": 1.0,
        "min_primer_len": 11,
        "max_primer_len": 11,
        "min_gc": 0.0,
        "max_gc": 100.0,
        "min_tm": 0.0,
        "max_tm": 100.0,
        "target_tm": 60.0,
        "min_amplicon": 22,
        "max_amplicon": 75,
        "target_amplicon": 45,
        "top_pairs": 5000,
    }
    values.update(overrides)
    return PrimerDesignParams(**values)


def _alignment(reference: str = REFERENCE) -> str:
    return f">ref transcript\n{reference}\n>other\n{reference}\n"


def _unique_candidates(pairs):
    candidates = {}
    for pair in pairs:
        for candidate in (pair.forward, pair.reverse):
            key = (
                candidate.orientation,
                candidate.start,
                candidate.end,
                candidate.sequence,
            )
            candidates.setdefault(key, candidate)
    return list(candidates.values())


def _matches(candidate, junction_position: int):
    return [
        match
        for match in candidate.junctions
        if match.junction_position == junction_position
    ]


def test_junction_primers_apply_5_and_3_prime_matches_by_orientation(monkeypatch):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)

    pairs = generate_exon_junction_primer_pairs(
        _alignment(),
        _params(),
        reference_id="ref",
        exons=EXONS,
        min_5_prime_match=7,
        min_3_prime_match=4,
    )

    assert pairs
    assert all(pair.reference_accession == "ref" for pair in pairs)
    assert all(
        pair.forward.spans_exon_junction or pair.reverse.spans_exon_junction
        for pair in pairs
    )

    candidates = _unique_candidates(pairs)
    assert all(
        candidate.spans_exon_junction == bool(candidate.junctions)
        for candidate in candidates
    )
    spanning = [candidate for candidate in candidates if candidate.spans_exon_junction]
    assert {candidate.orientation for candidate in spanning} == {"F", "R"}

    expected_exons = {30: ("1", "2"), 60: ("2", "3")}
    for candidate in spanning:
        reference_window = REFERENCE[candidate.start - 1 : candidate.end]
        expected_sequence = (
            reference_window
            if candidate.orientation == "F"
            else str(Seq(reference_window).reverse_complement())
        )
        assert candidate.sequence == expected_sequence

        for match in candidate.junctions:
            junction = match.junction_position
            assert candidate.start <= junction < candidate.end
            assert match.reference_accession == "ref"
            assert (
                match.left_exon_number,
                match.right_exon_number,
            ) == expected_exons[junction]

            left_bases = junction - candidate.start + 1
            right_bases = candidate.end - junction
            if candidate.orientation == "F":
                assert match.primer_5_prime_bases == left_bases
                assert match.primer_3_prime_bases == right_bases
            else:
                assert match.primer_5_prime_bases == right_bases
                assert match.primer_3_prime_bases == left_bases
            assert match.primer_5_prime_bases >= 7
            assert match.primer_3_prime_bases >= 4

    # With an 11-base primer and required matches of 7 + 4, there is exactly
    # one valid placement in each orientation at a given junction.
    forward_at_30 = {
        (candidate.start, candidate.end)
        for candidate in spanning
        if candidate.orientation == "F" and _matches(candidate, 30)
    }
    reverse_at_30 = {
        (candidate.start, candidate.end)
        for candidate in spanning
        if candidate.orientation == "R" and _matches(candidate, 30)
    }
    assert forward_at_30 == {(24, 34)}
    assert reverse_at_30 == {(27, 37)}


def test_primer_shorter_than_required_matches_cannot_span_a_junction(monkeypatch):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)

    pairs = generate_exon_junction_primer_pairs(
        _alignment(),
        _params(min_primer_len=10, max_primer_len=10),
        reference_id="ref",
        exons=EXONS,
        min_5_prime_match=7,
        min_3_prime_match=4,
    )

    assert pairs == []


def test_junction_filter_is_applied_before_top_pairs_limit(monkeypatch):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)

    pairs = generate_exon_junction_primer_pairs(
        _alignment(),
        _params(top_pairs=1),
        reference_id="ref",
        exons=EXONS,
    )

    assert len(pairs) == 1
    assert pairs[0].rank == 1
    assert (
        pairs[0].forward.spans_exon_junction
        or pairs[0].reverse.spans_exon_junction
    )


def test_generator_uses_ungapped_reference_coordinates_and_sequence(monkeypatch):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)
    complement = REFERENCE.translate(str.maketrans("ACGT", "TGCA"))
    aligned_reference = f"{REFERENCE[:10]}--{REFERENCE[10:]}"
    aligned_other = f"{complement[:10]}GG{complement[10:]}"
    alignment = (
        f">ref transcript\n{aligned_reference}\n"
        f">other-1\n{aligned_other}\n"
        f">other-2\n{aligned_other}\n"
    )

    pairs = generate_exon_junction_primer_pairs(
        alignment,
        _params(identity_threshold=0.0, coverage_threshold=0.0),
        reference_id="ref",
        exons=EXONS,
    )

    assert pairs
    candidates = _unique_candidates(pairs)
    for candidate in candidates:
        reference_window = REFERENCE[candidate.start - 1 : candidate.end]
        expected_sequence = (
            reference_window
            if candidate.orientation == "F"
            else str(Seq(reference_window).reverse_complement())
        )
        assert candidate.sequence == expected_sequence
        assert 1 <= candidate.start <= candidate.end <= len(REFERENCE)

    observed_junctions = {
        match.junction_position
        for candidate in candidates
        for match in candidate.junctions
    }
    assert observed_junctions == {30, 60}


def test_monoexonic_reference_returns_no_junction_pairs(monkeypatch):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)

    assert generate_exon_junction_primer_pairs(
        _alignment(),
        _params(),
        reference_id="ref",
        exons=[ExonInterval(1, len(REFERENCE), "1")],
    ) == []


def test_missing_reference_id_is_rejected(monkeypatch):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)

    with pytest.raises(ValueError):
        generate_exon_junction_primer_pairs(
            _alignment(),
            _params(),
            reference_id="not-in-alignment",
            exons=EXONS,
        )


@pytest.mark.parametrize(
    "raw_exons",
    [
        [(0, 30, "1"), (31, 90, "2")],
        [(1, 30, "1"), (31, 20, "2")],
        [(1, 30, "1"), (32, 90, "2")],
        [(1, 30, "1"), (30, 90, "2")],
        [(31, 60, "2"), (1, 30, "1"), (61, 90, "3")],
        [(1, 30, "1"), (31, 91, "2")],
    ],
    ids=[
        "zero-based-start",
        "reversed-interval",
        "non-contiguous-gap",
        "overlap",
        "unsorted",
        "outside-reference",
    ],
)
def test_invalid_or_non_contiguous_exons_are_rejected(monkeypatch, raw_exons):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)

    with pytest.raises(ValueError):
        exons = [ExonInterval(start, end, number) for start, end, number in raw_exons]
        generate_exon_junction_primer_pairs(
            _alignment(),
            _params(),
            reference_id="ref",
            exons=exons,
        )


def test_conventional_primer_generation_remains_backward_compatible(monkeypatch):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)

    pairs = generate_primer_pairs(_alignment(), _params(top_pairs=10))

    assert pairs
    assert [pair.rank for pair in pairs] == list(range(1, len(pairs) + 1))
    assert all(pair.forward.orientation == "F" for pair in pairs)
    assert all(pair.reverse.orientation == "R" for pair in pairs)
    assert all(
        not candidate.spans_exon_junction and candidate.junctions == []
        for candidate in _unique_candidates(pairs)
    )


def test_invalid_top_pairs_is_rejected():
    with pytest.raises(ValueError, match="maior que zero"):
        generate_primer_pairs(_alignment(), _params(top_pairs=0))
