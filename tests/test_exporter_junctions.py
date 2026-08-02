from __future__ import annotations

import csv

from openpyxl import load_workbook

from exporter import export_primers_xlsx, export_sequences_csv
from models import (
    ExonInterval,
    ExonJunctionMatch,
    PrimerCandidate,
    PrimerPair,
    SequenceRecord,
)


def _candidate(orientation: str, sequence: str) -> PrimerCandidate:
    return PrimerCandidate(
        orientation=orientation,
        start=10,
        end=29,
        sequence=sequence,
        length=20,
        gc_percent=50.0,
        tm_c=60.0,
        score=0.0,
    )


def test_xlsx_export_includes_exon_junction_metadata(tmp_path):
    forward = _candidate("F", "ACGTACGTACGTACGTACGT")
    forward.spans_exon_junction = True
    forward.junctions = [
        ExonJunctionMatch(
            reference_accession="NM_TEST.1",
            junction_position=20,
            left_exon_number="1",
            right_exon_number="2",
            primer_5_prime_bases=11,
            primer_3_prime_bases=9,
        )
    ]
    reverse = _candidate("R", "TGCATGCATGCATGCATGCA")
    pair = PrimerPair(
        rank=1,
        forward=forward,
        reverse=reverse,
        amplicon_start=10,
        amplicon_end=80,
        amplicon_length=71,
        score=1.0,
        reference_accession="NM_TEST.1",
    )

    path = tmp_path / "primers.xlsx"
    export_primers_xlsx(path, [pair])

    worksheet = load_workbook(path, read_only=True)["Pares de primers"]
    headers = [cell.value for cell in worksheet[1]]
    values = [cell.value for cell in worksheet[2]]
    row = dict(zip(headers, values))
    assert row["Referência da junção"] == "NM_TEST.1"
    assert row["F cruza junção"] is True
    assert "20 (E1-E2; 5'=11; 3'=9)" in row["F junções"]
    assert row["R cruza junção"] is False


def test_sequence_csv_includes_exon_annotations(tmp_path):
    record = SequenceRecord(
        uid="NM_TEST.1",
        accession="NM_TEST.1",
        definition="transcript",
        organism="Homo sapiens",
        sequence="ACGT",
        length=4,
        molecule_type="mRNA",
        exons=[ExonInterval(1, 2, "1"), ExonInterval(3, 4, "2")],
    )
    path = tmp_path / "sequences.csv"
    export_sequences_csv(path, [record])

    with path.open(encoding="utf-8-sig", newline="") as handle:
        row = next(csv.DictReader(handle))
    assert row["molecule_type"] == "mRNA"
    assert row["exons"] == "1:1-2; 2:3-4"
