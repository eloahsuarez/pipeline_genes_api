from __future__ import annotations

from models import ExonInterval, SequenceRecord
from ncbi_client import NcbiClient, NcbiSearchParams


SEQUENCE = "ACGT" * 15


def _params() -> NcbiSearchParams:
    return NcbiSearchParams(
        gene="HBB",
        organism="Homo sapiens",
        sequence_type="mRNA",
        refseq_only=False,
        exclude_predicted=False,
        exclude_partial=False,
        require_gene_feature=True,
        min_length=1,
        max_length=1000,
    )


def _exon_feature(start: int, end: int, number: str | None) -> str:
    number_qualifier = ""
    if number is not None:
        number_qualifier = f"""
          <GBQualifier>
            <GBQualifier_name>number</GBQualifier_name>
            <GBQualifier_value>{number}</GBQualifier_value>
          </GBQualifier>
        """
    return f"""
      <GBFeature>
        <GBFeature_key>exon</GBFeature_key>
        <GBFeature_location>{start}..{end}</GBFeature_location>
        <GBFeature_intervals>
          <GBInterval>
            <GBInterval_from>{start}</GBInterval_from>
            <GBInterval_to>{end}</GBInterval_to>
          </GBInterval>
        </GBFeature_intervals>
        <GBFeature_quals>
          {number_qualifier}
        </GBFeature_quals>
      </GBFeature>
    """


def _xml(exon_features: str = "") -> str:
    return f"""
    <GBSet>
      <GBSeq>
        <GBSeq_locus>TEST_TRANSCRIPT</GBSeq_locus>
        <GBSeq_length>{len(SEQUENCE)}</GBSeq_length>
        <GBSeq_moltype>mRNA</GBSeq_moltype>
        <GBSeq_definition>Homo sapiens HBB transcript</GBSeq_definition>
        <GBSeq_primary-accession>NM_TEST</GBSeq_primary-accession>
        <GBSeq_accession-version>NM_TEST.1</GBSeq_accession-version>
        <GBSeq_organism>Homo sapiens</GBSeq_organism>
        <GBSeq_feature-table>
          <GBFeature>
            <GBFeature_key>gene</GBFeature_key>
            <GBFeature_location>1..{len(SEQUENCE)}</GBFeature_location>
            <GBFeature_quals>
              <GBQualifier>
                <GBQualifier_name>gene</GBQualifier_name>
                <GBQualifier_value>HBB</GBQualifier_value>
              </GBQualifier>
            </GBFeature_quals>
          </GBFeature>
          {exon_features}
        </GBSeq_feature-table>
        <GBSeq_sequence>{SEQUENCE.lower()}</GBSeq_sequence>
      </GBSeq>
    </GBSet>
    """


def test_genbank_xml_parser_preserves_exons_numbers_and_molecule_type():
    xml = _xml(
        _exon_feature(1, 20, "1")
        + _exon_feature(21, 45, "2")
        + _exon_feature(46, 60, None)
    )

    records = NcbiClient._parse_genbank_xml(xml, _params())

    assert len(records) == 1
    record = records[0]
    assert record.accession == "NM_TEST.1"
    assert record.molecule_type == "mRNA"
    assert record.exons == [
        ExonInterval(1, 20, "1"),
        ExonInterval(21, 45, "2"),
        ExonInterval(46, 60, ""),
    ]
    assert "exon" in record.feature_keys

    serialized = record.to_dict()
    assert serialized["molecule_type"] == "mRNA"
    assert serialized["exons"] == [
        {"start": 1, "end": 20, "number": "1"},
        {"start": 21, "end": 45, "number": "2"},
        {"start": 46, "end": 60, "number": ""},
    ]


def test_genbank_xml_parser_does_not_infer_exons_when_features_are_absent():
    records = NcbiClient._parse_genbank_xml(_xml(), _params())

    assert len(records) == 1
    assert records[0].molecule_type == "mRNA"
    assert records[0].exons == []


def test_genbank_xml_parser_preserves_coordinates_for_iupac_bases():
    ambiguous = "ACGTRYSWKMBDHVN" * 4
    xml = _xml().replace(SEQUENCE.lower(), ambiguous.lower())

    records = NcbiClient._parse_genbank_xml(xml, _params())

    assert len(records) == 1
    assert len(records[0].sequence) == len(ambiguous)
    assert records[0].sequence == "".join(
        base if base in "ACGT" else "N" for base in ambiguous
    )


def test_sequence_record_new_fields_have_independent_backward_compatible_defaults():
    first = SequenceRecord(
        uid="one",
        accession="one",
        definition="first",
        organism="Homo sapiens",
        sequence="ACGT",
        length=4,
    )
    second = SequenceRecord(
        uid="two",
        accession="two",
        definition="second",
        organism="Homo sapiens",
        sequence="TGCA",
        length=4,
    )

    assert first.molecule_type == ""
    assert first.exons == []
    first.exons.append(ExonInterval(1, 4, "1"))
    assert second.exons == []
