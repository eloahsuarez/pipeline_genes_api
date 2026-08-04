from io import StringIO

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

import bioinformatics
from bioinformatics import (
    PrimerDesignParams,
    find_conserved_regions,
    generate_primer_pairs,
    generate_primer_pairs_for_sequence,
    pairwise_align_fasta,
    parse_alignment_fasta,
)


SEQUENCE = (
    "ACGTCGATGCTAGCTACGATCGTACGATGCTAGCTACGAT"
    "GTCAGTCCGATGACCTAGGCTAACGTTCGACATGCTAGCA"
)


def _permissive_params(**overrides: object) -> PrimerDesignParams:
    values: dict[str, object] = {
        "min_primer_len": 8,
        "max_primer_len": 8,
        "min_gc": 0,
        "max_gc": 100,
        "min_tm": 0,
        "max_tm": 100,
        "target_tm": 60,
        "min_amplicon": 20,
        "max_amplicon": 45,
        "target_amplicon": 30,
        "top_pairs": 300,
    }
    values.update(overrides)
    return PrimerDesignParams(**values)


def test_conserved_region_exact():
    alignment = ">a\nAACCGGTTAA\n>b\nAACCGGTTAA\n>c\nAACCGGTTAA\n"
    regions = find_conserved_regions(alignment, 1.0, 1.0, 4)
    assert len(regions) == 1
    assert regions[0].sequence == "AACCGGTTAA"
    assert regions[0].length == 10


def test_gap_breaks_exact_region():
    alignment = ">a\nAACCGGTTAA\n>b\nAACC-GTTAA\n>c\nAACCGGTTAA\n"
    regions = find_conserved_regions(alignment, 1.0, 1.0, 4)
    assert [r.sequence for r in regions] == ["AACC", "GTTAA"]


def test_primer_generation_returns_pairs():
    seq = "ACGTCGATGCTAGCTACGATCGTACGATGCTAGCTACGATCGTACGATGCTAGCTACGATCGTACGATGCTAGCTACGATCGTACGATGCTAGCTACGATCGT"
    alignment = f">a\n{seq}\n>b\n{seq}\n>c\n{seq}\n"
    params = PrimerDesignParams(
        min_primer_len=18, max_primer_len=22, min_gc=25, max_gc=75,
        min_tm=40, max_tm=80, min_amplicon=60, max_amplicon=110,
        target_amplicon=80, top_pairs=5,
    )
    pairs = generate_primer_pairs(alignment, params)
    assert pairs
    assert pairs[0].forward.orientation == "F"
    assert pairs[0].reverse.orientation == "R"


def test_pairwise_alignment_preserves_headers_and_ungapped_sequences():
    fasta = (
        ">first alpha description\nACGTTGCA\n"
        ">second beta description\nACTTGCA\n"
    )

    result = pairwise_align_fasta(fasta)
    records = list(SeqIO.parse(StringIO(result), "fasta"))

    assert [(record.id, record.description) for record in records] == [
        ("first", "first alpha description"),
        ("second", "second beta description"),
    ]
    assert len(records[0].seq) == len(records[1].seq)
    assert str(records[0].seq).replace("-", "") == "ACGTTGCA"
    assert str(records[1].seq).replace("-", "") == "ACTTGCA"


@pytest.mark.parametrize(
    "fasta",
    [
        ">one\nACGT\n",
        ">one\nACGT\n>two\nACGT\n>three\nACGT\n",
    ],
)
def test_pairwise_alignment_requires_exactly_two_records(fasta):
    with pytest.raises(ValueError, match="exatamente duas"):
        pairwise_align_fasta(fasta)


def test_pairwise_alignment_enforces_cell_limit(monkeypatch):
    monkeypatch.setattr(bioinformatics, "MAX_PAIRWISE_ALIGNMENT_CELLS", 20)

    with pytest.raises(ValueError, match="limite local"):
        pairwise_align_fasta(">one\nACGT\n>two\nACGT\n")


def test_alignment_parser_accepts_one_record_and_rejects_empty_input():
    assert parse_alignment_fasta(">only description\nAcgT\n") == [("only", "ACGT")]

    with pytest.raises(ValueError, match="não contém sequências"):
        parse_alignment_fasta("")
    with pytest.raises(ValueError, match="sequência vazia"):
        parse_alignment_fasta(">only\n")


def test_alignment_parser_still_rejects_different_lengths():
    with pytest.raises(ValueError, match="mesmo comprimento"):
        parse_alignment_fasta(">one\nACGT\n>two\nACG\n")


def test_direct_primer_design_normalizes_sequence_and_sets_reference(monkeypatch):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)
    raw_sequence = (
        f"  {SEQUENCE[:40].lower().replace('t', 'u')}\n"
        f"{SEQUENCE[40:].lower().replace('t', 'u')}  "
    )

    pairs = generate_primer_pairs_for_sequence(
        raw_sequence,
        _permissive_params(top_pairs=10),
        reference_accession=" NM_TEST ",
    )

    assert pairs
    assert all(pair.reference_accession == "NM_TEST" for pair in pairs)
    for pair in pairs:
        for candidate in (pair.forward, pair.reverse):
            window = SEQUENCE[candidate.start - 1 : candidate.end]
            expected = (
                window
                if candidate.orientation == "F"
                else str(Seq(window).reverse_complement())
            )
            assert candidate.sequence == expected


def test_direct_primer_design_preserves_n_coordinate_and_never_crosses_it(monkeypatch):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)
    n_position = 41
    sequence = f"{SEQUENCE[: n_position - 1]}N{SEQUENCE[n_position - 1 :]}"

    pairs = generate_primer_pairs_for_sequence(sequence, _permissive_params(), "NM_WITH_N")
    candidates = [candidate for pair in pairs for candidate in (pair.forward, pair.reverse)]

    assert pairs
    assert any(candidate.start > n_position for candidate in candidates)
    assert all(not (candidate.start <= n_position <= candidate.end) for candidate in candidates)
    for candidate in candidates:
        window = sequence[candidate.start - 1 : candidate.end]
        expected = (
            window
            if candidate.orientation == "F"
            else str(Seq(window).reverse_complement())
        )
        assert candidate.sequence == expected
        assert "N" not in candidate.sequence


def test_two_sequence_alignment_feeds_conservation_and_primer_pipeline(monkeypatch):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)
    sequence_with_deletion = f"{SEQUENCE[:35]}{SEQUENCE[36:]}"
    alignment = pairwise_align_fasta(
        f">reference full sequence\n{SEQUENCE}\n"
        f">variant one-base deletion\n{sequence_with_deletion}\n"
    )

    parsed = parse_alignment_fasta(alignment)
    pairs = generate_primer_pairs(alignment, _permissive_params(top_pairs=10))

    assert len(parsed) == 2
    assert len(parsed[0][1]) == len(parsed[1][1])
    assert pairs


def test_aligned_ambiguous_column_is_not_compressed_into_a_primer(monkeypatch):
    monkeypatch.setattr(bioinformatics, "calculate_tm", lambda sequence, params: 60.0)
    n_position = 41
    sequence = f"{SEQUENCE[: n_position - 1]}N{SEQUENCE[n_position - 1 :]}"
    alignment = f">one\n{sequence}\n>two\n{sequence}\n"

    pairs = generate_primer_pairs(alignment, _permissive_params())
    candidates = [candidate for pair in pairs for candidate in (pair.forward, pair.reverse)]

    assert pairs
    assert any(candidate.start > n_position for candidate in candidates)
    assert all(not (candidate.start <= n_position <= candidate.end) for candidate in candidates)
