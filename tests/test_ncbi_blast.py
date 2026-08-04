from __future__ import annotations

from types import SimpleNamespace

import pytest
import requests

from ncbi_blast import (
    BLAST_URL,
    ENTREZ_QUERY,
    REFSEQ_MRNA_DATABASE,
    REFSEQ_SELECT_DATABASE,
    BlastSubmission,
    NcbiBlastClient,
    NcbiBlastError,
    NcbiBlastParams,
    NcbiPrimerHit,
    PrimerPairQuery,
    correlate_primer_hits,
    parse_blast_xml2,
    primer_pairs_to_fasta,
)


def _pair(
    rank: int,
    forward: str = "ACGTACGTACGTACGTACGT",
    reverse: str = "TGCATGCATGCATGCATGCA",
):
    return SimpleNamespace(
        rank=rank,
        forward=SimpleNamespace(sequence=forward),
        reverse=SimpleNamespace(sequence=reverse),
    )


class FakeResponse:
    def __init__(self, text: str, status_code: int = 200) -> None:
        self.text = text
        self.status_code = status_code

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")


class FakeSession:
    def __init__(self, *, post=(), get=()) -> None:
        self.post_responses = list(post)
        self.get_responses = list(get)
        self.post_calls: list[tuple[str, dict]] = []
        self.get_calls: list[tuple[str, dict]] = []

    def post(self, url, **kwargs):
        self.post_calls.append((url, kwargs))
        if not self.post_responses:
            raise AssertionError("POST inesperado")
        response = self.post_responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response

    def get(self, url, **kwargs):
        self.get_calls.append((url, kwargs))
        if not self.get_responses:
            raise AssertionError("GET inesperado")
        response = self.get_responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response


def _hsp(
    *,
    query_start: int = 1,
    query_end: int = 20,
    subject_start: int,
    subject_end: int,
    identities: int = 20,
    alignment_length: int = 20,
    bit_score: float = 40.0,
    qseq: str = "",
    hseq: str = "",
) -> str:
    hit_frame = 1 if subject_start <= subject_end else -1
    hit_strand = "Plus" if hit_frame == 1 else "Minus"
    aligned_sequences = ""
    if qseq or hseq:
        aligned_sequences = f"<qseq>{qseq}</qseq><hseq>{hseq}</hseq>"
    return f"""
      <Hsp>
        <bit-score>{bit_score}</bit-score><evalue>1e-8</evalue>
        <query-from>{query_start}</query-from><query-to>{query_end}</query-to>
        <hit-from>{subject_start}</hit-from><hit-to>{subject_end}</hit-to>
        <query-frame>1</query-frame><query-strand>Plus</query-strand>
        <hit-frame>{hit_frame}</hit-frame><hit-strand>{hit_strand}</hit-strand>
        <identity>{identities}</identity><align-len>{alignment_length}</align-len>
        {aligned_sequences}
      </Hsp>
    """


def _hit(
    accession: str,
    title: str,
    hsp: str,
    *,
    taxid: int | None = 9606,
    scientific_name: str = "Homo sapiens",
) -> str:
    taxid_xml = f"<taxid>{taxid}</taxid>" if taxid is not None else ""
    return f"""
      <Hit>
        <description><HitDescr><id>ref|{accession}|</id>
          <accession>{accession}</accession><title>{title}</title>
          {taxid_xml}<sciname>{scientific_name}</sciname>
        </HitDescr></description>
        <hsps>{hsp}</hsps>
      </Hit>
    """


def _search(
    query_id: str,
    hits: str = "",
    query_length: int = 20,
    *,
    include_statistics: bool = True,
    db_num: int | None = 10,
    db_len: int | None = 1000,
) -> str:
    statistics = ""
    if include_statistics:
        db_num_xml = f"<db-num>{db_num}</db-num>" if db_num is not None else ""
        db_len_xml = f"<db-len>{db_len}</db-len>" if db_len is not None else ""
        statistics = (
            "<stat><Statistics>"
            f"{db_num_xml}{db_len_xml}"
            "</Statistics></stat>"
        )
    return f"""
      <Search>
        <query-id>Query_1</query-id><query-title>{query_id}</query-title>
        <query-len>{query_length}</query-len><hits>{hits}</hits>
        {statistics}
      </Search>
    """


def _xml(*searches: str) -> str:
    return (
        '<?xml version="1.0"?>'
        '<BlastXML2 xmlns="http://www.ncbi.nlm.nih.gov">'
        + "".join(searches)
        + "</BlastXML2>"
    )


def _xml_with_diagnostics(
    *searches: str,
    messages: tuple[str, ...] = (),
    error: str = "",
    db_num: int | None = None,
    db_len: int | None = None,
) -> str:
    diagnostics = "".join(f"<message>{message}</message>" for message in messages)
    if error:
        diagnostics += f"<error>{error}</error>"
    statistics = ""
    if db_num is not None or db_len is not None:
        db_num_xml = f"<db-num>{db_num}</db-num>" if db_num is not None else ""
        db_len_xml = f"<db-len>{db_len}</db-len>" if db_len is not None else ""
        statistics = f"<stat><Statistics>{db_num_xml}{db_len_xml}</Statistics></stat>"
    return (
        '<?xml version="1.0"?>'
        '<BlastXML2 xmlns="http://www.ncbi.nlm.nih.gov">'
        + diagnostics
        + statistics
        + "".join(searches)
        + "</BlastXML2>"
    )


def _client(session: FakeSession, sleeps: list[float] | None = None) -> NcbiBlastClient:
    sleep_calls = sleeps if sleeps is not None else []
    return NcbiBlastClient(
        "pesquisadora@example.org",
        tool="GeneConservadoGUI",
        session=session,
        sleeper=sleep_calls.append,
    )


def test_primer_pairs_to_fasta_uses_first_n_and_normalizes_sequences():
    fasta, queries = primer_pairs_to_fasta(
        [
            _pair(1, "acgt acgt", "tgca\ntgca"),
            _pair(2, "AAAA", "CCCC"),
        ],
        1,
    )

    assert fasta == ">pair_0001_F\nACGTACGT\n>pair_0001_R\nTGCATGCA\n"
    assert queries == [
        PrimerPairQuery(
            pair_id="pair_0001",
            pair_rank=1,
            forward_sequence="ACGTACGT",
            reverse_sequence="TGCATGCA",
        )
    ]


@pytest.mark.parametrize(
    "pairs, message",
    [
        ([], "Nenhum par"),
        ([_pair(1, forward="ACGTX")], "bases inválidas"),
        ([_pair(1), SimpleNamespace(rank=1, forward=None, reverse=None)], "ranks duplicados"),
        ([SimpleNamespace(rank=1, forward=None, reverse=None)], "Forward.*ausente"),
        ([SimpleNamespace(rank="x", forward=None, reverse=None)], "rank inteiro"),
    ],
)
def test_primer_pairs_to_fasta_rejects_invalid_pairs(pairs, message):
    with pytest.raises(NcbiBlastError, match=message):
        primer_pairs_to_fasta(pairs, 10)


def test_analyze_submits_exact_refseq_select_payload_and_gets_xml2_s():
    session = FakeSession(
        post=[FakeResponse("QBlastInfoBegin\nRID = RID123\nRTOE = 0\nQBlastInfoEnd")],
        get=[FakeResponse(_xml(_search("pair_0001_F"), _search("pair_0001_R")))],
    )
    params = NcbiBlastParams(
        top_pairs=1,
        hitlist_size=77,
        target_gene="BRCA1",
        database=REFSEQ_SELECT_DATABASE,
        megablast=True,
        short_query_adjust=True,
    )

    report = _client(session).analyze([_pair(1)], params)

    assert session.post_calls == [
        (
            BLAST_URL,
            {
                "data": {
                    "CMD": "Put",
                    "PROGRAM": "blastn",
                    "DATABASE": REFSEQ_SELECT_DATABASE,
                    "MEGABLAST": "on",
                    "ENTREZ_QUERY": ENTREZ_QUERY,
                    "SHORT_QUERY_ADJUST": "true",
                    "HITLIST_SIZE": "77",
                    "email": "pesquisadora@example.org",
                    "tool": "GeneConservadoGUI",
                    "QUERY": (
                        ">pair_0001_F\nACGTACGTACGTACGTACGT\n"
                        ">pair_0001_R\nTGCATGCATGCATGCATGCA\n"
                    ),
                },
                "timeout": 60.0,
            },
        )
    ]
    assert session.get_calls == [
        (
            BLAST_URL,
            {
                "params": {
                    "CMD": "Get",
                    "RID": "RID123",
                    "FORMAT_TYPE": "XML2_S",
                    "HITLIST_SIZE": "77",
                    "email": "pesquisadora@example.org",
                    "tool": "GeneConservadoGUI",
                },
                "timeout": 60.0,
            },
        )
    ]
    assert report.rid == "RID123"
    assert report.query_fasta.startswith(">pair_0001_F")
    assert (report.database, report.organism, report.taxid) == (
        REFSEQ_SELECT_DATABASE,
        "Homo sapiens",
        9606,
    )
    assert report.results[0].verdict == "Nenhum produto conjunto"
    assert "não cobre todas as isoformas" in report.warnings[0]


def test_submit_uses_sensitive_short_oligo_payload_by_default():
    session = FakeSession(
        post=[FakeResponse("RID = RID-SENSITIVE\nRTOE = 1")],
    )
    fasta = ">pair_0001_F\nACGT\n"

    submission = _client(session).submit(fasta, NcbiBlastParams())

    assert submission.rid == "RID-SENSITIVE"
    payload = session.post_calls[0][1]["data"]
    assert "MEGABLAST" not in payload
    assert payload == {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": REFSEQ_MRNA_DATABASE,
        "ENTREZ_QUERY": ENTREZ_QUERY,
        "HITLIST_SIZE": "50000",
        "email": "pesquisadora@example.org",
        "tool": "GeneConservadoGUI",
        "QUERY": fasta,
        "SHORT_QUERY_ADJUST": "false",
        "WORD_SIZE": "7",
        "EXPECT": "30000",
        "NUCL_REWARD": "1",
        "NUCL_PENALTY": "-3",
        "GAPCOSTS": "5 2",
        "FILTER": "F",
    }


def test_remote_database_constants_use_ncbi_aliases():
    assert REFSEQ_SELECT_DATABASE == "refseq_select_rna"
    assert REFSEQ_MRNA_DATABASE == "refseq_mrna"


@pytest.mark.parametrize(
    "database",
    [REFSEQ_SELECT_DATABASE, REFSEQ_MRNA_DATABASE],
)
def test_submit_sends_each_supported_database_exactly(database):
    session = FakeSession(
        post=[FakeResponse("RID = RID-DATABASE\nRTOE = 1")],
    )

    submission = _client(session).submit(
        ">pair_0001_F\nACGT\n",
        NcbiBlastParams(database=database),
    )

    assert submission.rid == "RID-DATABASE"
    assert session.post_calls[0][1]["data"]["DATABASE"] == database


def test_submit_rejects_an_unsupported_database_before_http_request():
    session = FakeSession()

    with pytest.raises(NcbiBlastError):
        _client(session).submit(
            ">pair_0001_F\nACGT\n",
            NcbiBlastParams(database="refseq_select"),
        )

    assert session.post_calls == []


def test_parser_rejects_ncbi_error_reported_in_message_element():
    xml_text = _xml_with_diagnostics(
        _search("pair_0001_F"),
        messages=(
            "Error: No alias or index file found for nucleotide database "
            "[refseq_select]",
        ),
        db_num=10,
        db_len=1000,
    )

    with pytest.raises(NcbiBlastError, match="No alias or index file"):
        parse_blast_xml2(xml_text)


@pytest.mark.parametrize(
    ("db_num", "db_len"),
    [(0, 1000), (10, 0), (0, 0)],
)
def test_parser_rejects_explicitly_empty_database_statistics(db_num, db_len):
    xml_text = _xml_with_diagnostics(
        _search("pair_0001_F"),
        db_num=db_num,
        db_len=db_len,
    )

    with pytest.raises(NcbiBlastError):
        parse_blast_xml2(xml_text)


def test_parser_keeps_benign_no_hits_message_with_positive_database_statistics():
    xml_text = _xml_with_diagnostics(
        _search("pair_0001_F"),
        messages=("No hits found",),
        db_num=10,
        db_len=1000,
    )

    parsed = parse_blast_xml2(xml_text)

    assert parsed.hits == []
    assert parsed.raw_hit_counts == {"pair_0001_F": 0}


def test_parser_rejects_result_xml_without_search_statistics():
    xml_text = _xml(
        _search("pair_0001_F", include_statistics=False),
    )

    with pytest.raises(NcbiBlastError, match=r"faltam estatísticas.*pair_0001_F"):
        parse_blast_xml2(xml_text)


def test_parser_requires_positive_statistics_for_each_search():
    xml_text = _xml(
        _search("pair_0001_F"),
        _search("pair_0001_R", db_len=None),
    )

    with pytest.raises(
        NcbiBlastError,
        match=r"db-len inválida.*pair_0001_R",
    ):
        parse_blast_xml2(xml_text)


def test_parser_detects_error_element_after_a_benign_message():
    xml_text = _xml_with_diagnostics(
        _search("pair_0001_F"),
        messages=("No hits found",),
        error="Database unavailable",
        db_num=10,
        db_len=1000,
    )

    with pytest.raises(NcbiBlastError, match="Database unavailable"):
        parse_blast_xml2(xml_text)


@pytest.mark.parametrize("returned_orientation", ["F", "R"])
def test_analyze_rejects_nonempty_xml_missing_a_submitted_query(returned_orientation):
    session = FakeSession(
        post=[FakeResponse("RID = RID-PARTIAL\nRTOE = 0")],
        get=[FakeResponse(_xml(_search(f"pair_0001_{returned_orientation}")))],
    )
    missing_orientation = "R" if returned_orientation == "F" else "F"

    with pytest.raises(
        NcbiBlastError,
        match=rf"incompleta.*pair_0001_{missing_orientation}",
    ):
        _client(session).analyze([_pair(1)], NcbiBlastParams(top_pairs=1))


def test_analyze_rejects_unknown_query_even_when_it_has_no_hits():
    xml_text = _xml(
        _search("pair_0001_F"),
        _search("pair_0001_R"),
        _search("pair_9999_F"),
    )
    session = FakeSession(
        post=[FakeResponse("RID = RID-UNKNOWN-EMPTY\nRTOE = 0")],
        get=[FakeResponse(xml_text)],
    )

    with pytest.raises(
        NcbiBlastError,
        match=r"não pertencem.*pair_9999_F",
    ):
        _client(session).analyze([_pair(1)], NcbiBlastParams(top_pairs=1))


def test_analyze_accepts_explicit_there_are_no_hits_response():
    session = FakeSession(
        post=[FakeResponse("RID = RID-NO-HITS\nRTOE = 0")],
        get=[FakeResponse("Status=READY\nThereAreHits=no")],
    )

    report = _client(session).analyze([_pair(1)], NcbiBlastParams(top_pairs=1))

    assert report.rid == "RID-NO-HITS"
    assert len(report.results) == 1
    assert report.results[0].products == []
    assert report.results[0].verdict == "Nenhum produto conjunto"


def test_analyze_does_not_hide_error_xml_behind_no_hits_header():
    invalid_xml = _xml_with_diagnostics(
        _search("pair_0001_F"),
        _search("pair_0001_R"),
        messages=(
            "Error: No alias or index file found for nucleotide database [invalid]",
        ),
        db_num=0,
        db_len=0,
    )
    session = FakeSession(
        post=[FakeResponse("RID = RID-INVALID-DATABASE\nRTOE = 0")],
        get=[FakeResponse("Status=READY\nThereAreHits=no\n" + invalid_xml)],
    )

    with pytest.raises(NcbiBlastError, match="No alias or index file"):
        _client(session).analyze([_pair(1)], NcbiBlastParams(top_pairs=1))


def test_xml2_namespace_combines_target_and_off_target_and_keeps_metrics():
    target_title = "Homo sapiens alvo (BRCA1), transcript variant 1, mRNA"
    other_title = "Homo sapiens outro gene (OTHER1), mRNA"
    xml_text = _xml(
        _search(
            "pair_0001_F",
            _hit(
                "NM_TARGET.7",
                target_title,
                _hsp(
                    query_start=3,
                    query_end=20,
                    subject_start=102,
                    subject_end=119,
                    identities=17,
                    alignment_length=18,
                    bit_score=42,
                ),
            )
            + _hit("NM_OTHER.1", other_title, _hsp(subject_start=500, subject_end=519)),
        ),
        _search(
            "pair_0001_R",
            _hit("NM_TARGET.7", target_title, _hsp(subject_start=250, subject_end=231))
            + _hit("NM_OTHER.1", other_title, _hsp(subject_start=680, subject_end=661)),
        ),
    )
    session = FakeSession(
        post=[FakeResponse("RID = RID-MULTI\nRTOE = 0")],
        get=[FakeResponse(xml_text)],
    )
    params = NcbiBlastParams(
        top_pairs=1,
        target_gene="BRCA1",
        target_accessions=("NM_TARGET",),
        max_amplicon=200,
        min_identity_pct=80,
        min_query_coverage_pct=80,
    )

    report = _client(session).analyze([_pair(1)], params)

    result = report.results[0]
    assert result.verdict == "Pode amplificar outro gene"
    assert [product.classification for product in result.products] == ["target", "off_target"]
    target, off_target = result.products
    assert (target.accession, target.gene, target.start, target.end, target.length) == (
        "NM_TARGET.7",
        "BRCA1",
        100,
        250,
        151,
    )
    assert target.title == target_title
    assert target.organism == "Homo sapiens"
    assert target.forward_hit.query_coverage_pct == pytest.approx(90.0)
    assert target.forward_hit.identity_pct == pytest.approx(100 * 17 / 18)
    assert target.reverse_hit.query_coverage_pct == 100.0
    assert off_target.gene == "OTHER1"
    assert report.to_dict()["results"][0]["products"][0]["forward_hit"][
        "identity_pct"
    ] == pytest.approx(100 * 17 / 18)


def _manual_hit(
    orientation: str,
    accession: str,
    start: int,
    end: int,
    *,
    title: str = "Homo sapiens gene (GENE1), mRNA",
    score: float = 40.0,
) -> NcbiPrimerHit:
    return NcbiPrimerHit(
        pair_id="pair_0001",
        pair_rank=1,
        primer_orientation=orientation,
        query_id=f"pair_0001_{orientation}",
        query_length=20,
        accession=accession,
        title=title,
        taxid=9606,
        scientific_name="Homo sapiens",
        identity_pct=100.0,
        query_coverage_pct=100.0,
        identities=20,
        alignment_length=20,
        query_start=1,
        query_end=20,
        subject_start=start,
        subject_end=end,
        query_strand="plus",
        subject_strand="plus" if start <= end else "minus",
        evalue=1e-8,
        bit_score=score,
    )


def test_pair_correlation_rejects_separate_accessions_same_strand_outward_and_distance():
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]
    hits = [
        # Acesso diferente.
        _manual_hit("F", "NM_A.1", 100, 119),
        _manual_hit("R", "NM_B.1", 250, 231),
        # Mesmo strand.
        _manual_hit("F", "NM_SAME.1", 100, 119),
        _manual_hit("R", "NM_SAME.1", 180, 199),
        # Oligos 3′ voltados para fora.
        _manual_hit("F", "NM_OUT.1", 119, 100),
        _manual_hit("R", "NM_OUT.1", 300, 319),
        # Distância acima do limite.
        _manual_hit("F", "NM_FAR.1", 1, 20),
        _manual_hit("R", "NM_FAR.1", 500, 481),
    ]

    result = correlate_primer_hits(
        hits,
        queries,
        NcbiBlastParams(target_gene="GENE1", max_amplicon=250),
    )[0]

    assert result.products == []
    assert result.verdict == "Nenhum produto conjunto"


def test_pair_correlation_accepts_reverse_plus_and_forward_minus_and_deduplicates():
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]
    hits = [
        _manual_hit("R", "NM_REV.1", 100, 119, score=20),
        _manual_hit("F", "NM_REV.1", 260, 241, score=20),
        # HSP duplicado com pontuação maior deve substituir o anterior.
        _manual_hit("R", "NM_REV.1", 100, 119, score=50),
    ]

    result = correlate_primer_hits(
        hits,
        queries,
        NcbiBlastParams(target_accessions=("NM_REV",), max_amplicon=200),
    )[0]

    assert result.verdict == "Específico no banco"
    assert len(result.products) == 1
    assert (result.products[0].start, result.products[0].end) == (100, 260)
    assert result.products[0].reverse_hit.bit_score == 50


@pytest.mark.parametrize("orientation", ["F", "R"])
def test_pair_correlation_accepts_inward_same_primer_products(orientation):
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]
    hits = [
        _manual_hit(orientation, "NM_HOMO.1", 50, 69),
        _manual_hit(orientation, "NM_HOMO.1", 180, 161),
    ]

    result = correlate_primer_hits(
        hits,
        queries,
        NcbiBlastParams(target_gene="CCR5", max_amplicon=1000),
    )[0]

    assert len(result.products) == 1
    product = result.products[0]
    assert product.primer_combination == f"{orientation}-{orientation}"
    assert (product.start, product.end, product.length) == (50, 180, 131)
    assert product.left_hit.subject_strand == "plus"
    assert product.right_hit.subject_strand == "minus"
    assert {
        id(product.left_hit),
        id(product.right_hit),
    } == {
        id(product.forward_hit),
        id(product.reverse_hit),
    }


def test_pair_correlation_finds_clk3_reverse_reverse_product_of_192_bp():
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]
    title = "Homo sapiens CDC like kinase 3 (CLK3), mRNA"
    hits = [
        _manual_hit("R", "NM_CLK3.1", 100, 119, title=title),
        _manual_hit("R", "NM_CLK3.1", 291, 272, title=title),
    ]

    result = correlate_primer_hits(
        hits,
        queries,
        NcbiBlastParams(target_gene="CCR5", max_amplicon=1000),
    )[0]

    assert result.verdict == "Apenas fora do alvo"
    assert len(result.products) == 1
    product = result.products[0]
    assert product.primer_combination == "R-R"
    assert product.gene == "CLK3"
    assert (product.start, product.end, product.length) == (100, 291, 192)
    assert product.left_hit is hits[0]
    assert product.right_hit is hits[1]


@pytest.mark.parametrize(
    ("max_amplicon", "expected_products"),
    [(2968, 0), (2969, 1), (4000, 1)],
)
def test_cntnap3b_reverse_reverse_product_obeys_inclusive_max_amplicon(
    max_amplicon,
    expected_products,
):
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]
    title = "Homo sapiens contactin associated protein family member 3B (CNTNAP3B), mRNA"
    hits = [
        _manual_hit("R", "NM_CNTNAP3B.1", 100, 119, title=title),
        _manual_hit("R", "NM_CNTNAP3B.1", 3068, 3049, title=title),
    ]

    result = correlate_primer_hits(
        hits,
        queries,
        NcbiBlastParams(target_gene="CCR5", max_amplicon=max_amplicon),
    )[0]

    assert len(result.products) == expected_products
    if expected_products:
        product = result.products[0]
        assert product.primer_combination == "R-R"
        assert product.gene == "CNTNAP3B"
        assert (product.start, product.end, product.length) == (100, 3068, 2969)


def test_pair_correlation_keeps_swapped_nkd1_forward_reverse_product_of_347_bp():
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]
    title = "Homo sapiens NKD inhibitor of WNT signaling pathway 1 (NKD1), mRNA"
    reverse_plus = _manual_hit("R", "NM_NKD1.1", 100, 119, title=title)
    forward_minus = _manual_hit("F", "NM_NKD1.1", 446, 427, title=title)

    result = correlate_primer_hits(
        [reverse_plus, forward_minus],
        queries,
        NcbiBlastParams(target_gene="CCR5", max_amplicon=1000),
    )[0]

    assert len(result.products) == 1
    product = result.products[0]
    assert product.primer_combination == "F-R"
    assert product.gene == "NKD1"
    assert (product.start, product.end, product.length) == (100, 446, 347)
    assert product.left_hit is reverse_plus
    assert product.right_hit is forward_minus
    assert product.forward_hit is forward_minus
    assert product.reverse_hit is reverse_plus


def test_pair_correlation_does_not_pair_the_same_hsp_with_itself():
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]
    hit = _manual_hit("R", "NM_SINGLE.1", 100, 119)

    result = correlate_primer_hits(
        [hit, hit],
        queries,
        NcbiBlastParams(target_gene="CCR5", max_amplicon=1000),
    )[0]

    assert result.products == []


def test_pair_correlation_preserves_distinct_combinations_at_same_coordinates():
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]
    hits = [
        _manual_hit("F", "NM_SHARED.1", 100, 119),
        _manual_hit("F", "NM_SHARED.1", 200, 181),
        _manual_hit("R", "NM_SHARED.1", 100, 119),
        _manual_hit("R", "NM_SHARED.1", 200, 181),
    ]

    result = correlate_primer_hits(
        hits,
        queries,
        NcbiBlastParams(target_gene="CCR5", max_amplicon=1000),
    )[0]

    products_by_combination = {
        product.primer_combination: product for product in result.products
    }
    assert set(products_by_combination) == {"F-F", "F-R", "R-R"}
    assert len(result.products) == 3
    assert {
        (product.start, product.end, product.length)
        for product in result.products
    } == {(100, 200, 101)}


def test_partial_hit_can_form_potential_product_in_sensitive_mode():
    xml_text = _xml(
        _search(
            "pair_0001_F",
            _hit(
                "NM_PARTIAL.1",
                "Homo sapiens gene (GENE1), mRNA",
                _hsp(
                    query_start=1,
                    query_end=18,
                    subject_start=100,
                    subject_end=117,
                    identities=18,
                    alignment_length=18,
                ),
            ),
        ),
        _search(
            "pair_0001_R",
            _hit(
                "NM_PARTIAL.1",
                "Homo sapiens gene (GENE1), mRNA",
                _hsp(subject_start=220, subject_end=201),
            ),
        ),
    )
    parsed = parse_blast_xml2(
        xml_text,
        min_identity_pct=80,
        min_query_coverage_pct=80,
    )
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]

    result = correlate_primer_hits(
        parsed.hits,
        queries,
        NcbiBlastParams(target_gene="GENE1", max_amplicon=250),
    )[0]

    assert len(parsed.hits) == 2
    assert parsed.hits[0].query_coverage_pct == 90.0
    assert parsed.hits[0].covers_query_three_prime is False
    assert parsed.hits[0].estimated_mismatches == 2
    assert parsed.hits[0].estimated_three_prime_mismatches == 2
    assert parsed.hits[0].subject_three_prime == 119
    assert len(result.products) == 1
    assert (
        result.products[0].start,
        result.products[0].end,
        result.products[0].length,
    ) == (100, 220, 121)


@pytest.mark.parametrize(
    ("qseq", "hseq", "identities", "alignment_length", "expected_three_prime"),
    [
        ("AAA-" + "A" * 17, "A" * 21, 20, 21, 0),
        ("A" * 17 + "-AAA", "A" * 21, 20, 21, 1),
        ("A" * 20, "AAA-" + "A" * 16, 19, 20, 0),
        ("A" * 20, "A" * 17 + "-AA", 19, 20, 1),
    ],
    ids=(
        "subject-insertion-outside-3prime",
        "subject-insertion-inside-3prime",
        "subject-deletion-outside-3prime",
        "subject-deletion-inside-3prime",
    ),
)
def test_internal_indels_are_positioned_in_three_prime_mismatch_window(
    qseq,
    hseq,
    identities,
    alignment_length,
    expected_three_prime,
):
    subject_bases = sum(base != "-" for base in hseq)
    xml_text = _xml(
        _search(
            "pair_0001_F",
            _hit(
                "NM_INDEL.1",
                "Homo sapiens indel test (INDEL), mRNA",
                _hsp(
                    subject_start=100,
                    subject_end=99 + subject_bases,
                    identities=identities,
                    alignment_length=alignment_length,
                    qseq=qseq,
                    hseq=hseq,
                ),
            ),
        )
    )

    parsed = parse_blast_xml2(xml_text, three_prime_window=5)

    assert len(parsed.hits) == 1
    assert parsed.hits[0].estimated_mismatches == 1
    assert (
        parsed.hits[0].estimated_three_prime_mismatches
        == expected_three_prime
    )


@pytest.fixture
def cntnap3b_partial_reverse_reverse_xml():
    title = (
        "Homo sapiens contactin associated protein family member 3B "
        "(CNTNAP3B), mRNA"
    )
    return _xml(
        _search(
            "pair_0001_R",
            _hit(
                "NM_CNTNAP3B.1",
                title,
                _hsp(
                    query_start=5,
                    query_end=16,
                    subject_start=638,
                    subject_end=649,
                    identities=12,
                    alignment_length=12,
                )
                + _hsp(
                    query_start=5,
                    query_end=17,
                    subject_start=3598,
                    subject_end=3586,
                    identities=12,
                    alignment_length=13,
                ),
            ),
            query_length=18,
        )
    )


def test_cntnap3b_real_partial_hits_form_2969_bp_reverse_reverse_product(
    cntnap3b_partial_reverse_reverse_xml,
):
    parsed = parse_blast_xml2(cntnap3b_partial_reverse_reverse_xml)
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 18, "T" * 18)]

    result = correlate_primer_hits(
        parsed.hits,
        queries,
        NcbiBlastParams(
            target_gene="CCR5",
            max_amplicon=4000,
            max_estimated_mismatches=6,
        ),
    )[0]

    assert len(parsed.hits) == 2
    assert [hit.estimated_mismatches for hit in parsed.hits] == [6, 6]
    assert len(result.products) == 1
    product = result.products[0]
    assert product.primer_combination == "R-R"
    assert product.gene == "CNTNAP3B"
    assert (product.start, product.end, product.length) == (634, 3602, 2969)


def test_cntnap3b_real_partial_hits_respect_estimated_mismatch_limit(
    cntnap3b_partial_reverse_reverse_xml,
):
    parsed = parse_blast_xml2(cntnap3b_partial_reverse_reverse_xml)
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 18, "T" * 18)]

    result = correlate_primer_hits(
        parsed.hits,
        queries,
        NcbiBlastParams(
            target_gene="CCR5",
            max_amplicon=4000,
            max_estimated_mismatches=5,
        ),
    )[0]

    assert [hit.estimated_mismatches for hit in parsed.hits] == [6, 6]
    assert result.products == []


def test_target_gene_is_a_token_not_a_substring():
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]
    hits = [
        _manual_hit(
            "F", "NM_BIND.1", 100, 119,
            title="Homo sapiens TP53BP1 (TP53BP1), mRNA",
        ),
        _manual_hit(
            "R", "NM_BIND.1", 200, 181,
            title="Homo sapiens TP53BP1 (TP53BP1), mRNA",
        ),
    ]

    result = correlate_primer_hits(
        hits, queries, NcbiBlastParams(target_gene="TP53")
    )[0]

    assert result.verdict == "Apenas fora do alvo"
    assert result.products[0].gene == "TP53BP1"


def test_inferred_off_target_gene_overrides_target_token_elsewhere_in_title():
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]
    title = (
        "Homo sapiens CCR5-associated transcript (OTHER1), "
        "transcript variant 1, mRNA"
    )
    hits = [
        _manual_hit("F", "NM_OTHER.1", 100, 119, title=title),
        _manual_hit("R", "NM_OTHER.1", 200, 181, title=title),
    ]

    result = correlate_primer_hits(
        hits,
        queries,
        NcbiBlastParams(target_gene="CCR5"),
    )[0]

    assert result.verdict == "Apenas fora do alvo"
    assert result.products[0].classification == "off_target"
    assert result.products[0].gene == "OTHER1"


def test_target_title_token_remains_fallback_when_gene_cannot_be_inferred():
    queries = [PrimerPairQuery("pair_0001", 1, "A" * 20, "T" * 20)]
    title = "Homo sapiens curated CCR5 transcript without gene annotation"
    hits = [
        _manual_hit("F", "NM_UNKNOWN.1", 100, 119, title=title),
        _manual_hit("R", "NM_UNKNOWN.1", 200, 181, title=title),
    ]

    result = correlate_primer_hits(
        hits,
        queries,
        NcbiBlastParams(target_gene="CCR5"),
    )[0]

    assert result.verdict == "Específico no banco"
    assert result.products[0].classification == "target"
    assert result.products[0].gene == "CCR5"


def test_parser_filters_thresholds_and_nonhuman_taxid_but_keeps_missing_taxid():
    xml_text = _xml(
        _search(
            "pair_0001_F",
            _hit(
                "NM_LOW_ID.1",
                "low",
                _hsp(subject_start=1, subject_end=20, identities=15),
            )
            + _hit(
                "NM_LOW_COV.1",
                "low",
                _hsp(
                    query_start=1,
                    query_end=10,
                    subject_start=30,
                    subject_end=39,
                    identities=10,
                    alignment_length=10,
                ),
            )
            + _hit(
                "NM_MOUSE.1",
                "mouse",
                _hsp(subject_start=50, subject_end=69),
                taxid=10090,
                scientific_name="Mus musculus",
            )
            + _hit(
                "NM_NO_TAXID.1",
                "unknown taxid",
                _hsp(subject_start=80, subject_end=99),
                taxid=None,
                scientific_name="",
            ),
        )
    )

    parsed = parse_blast_xml2(
        xml_text,
        min_identity_pct=80,
        min_query_coverage_pct=80,
    )

    assert [hit.accession for hit in parsed.hits] == ["NM_NO_TAXID.1"]
    assert any("Taxid ausente" in warning for warning in parsed.warnings)
    assert any("fora de Homo sapiens" in warning for warning in parsed.warnings)


def test_report_warns_when_hitlist_may_be_truncated():
    xml_text = _xml(
        _search(
            "pair_0001_F",
            _hit("NM_ONE.1", "gene", _hsp(subject_start=1, subject_end=20)),
        ),
        _search("pair_0001_R"),
    )
    session = FakeSession(
        post=[FakeResponse("RID = RID-LIMIT\nRTOE = 0")],
        get=[FakeResponse(xml_text)],
    )

    report = _client(session).analyze(
        [_pair(1)], NcbiBlastParams(hitlist_size=1, target_gene="GENE")
    )

    assert any("pode ter sido truncada" in warning for warning in report.warnings)
    assert report.results[0].verdict == "Inconclusivo: limite de hits atingido"


def test_poll_waits_for_rtoe_and_waiting_status():
    sleeps: list[float] = []
    session = FakeSession(
        get=[FakeResponse("Status=WAITING"), FakeResponse(_xml())]
    )
    params = NcbiBlastParams(
        poll_interval_seconds=2,
        poll_timeout_seconds=100,
    )

    result = _client(session, sleeps).poll(BlastSubmission("RID-WAIT", 3), params)

    assert result.startswith("<?xml")
    assert sleeps == [20.0, 2]


@pytest.mark.parametrize("status, message", [("FAILED", "falhou"), ("UNKNOWN", "RID")])
def test_poll_rejects_failed_and_unknown_status(status, message):
    session = FakeSession(get=[FakeResponse(f"Status={status}")])

    with pytest.raises(NcbiBlastError, match=message):
        _client(session).poll(BlastSubmission("RID-BAD", 0), NcbiBlastParams())


def test_poll_times_out_with_injected_sleeper():
    sleeps: list[float] = []
    session = FakeSession(
        get=[FakeResponse("Status=WAITING"), FakeResponse("Status=WAITING")]
    )
    params = NcbiBlastParams(poll_interval_seconds=1, poll_timeout_seconds=22)

    with pytest.raises(NcbiBlastError, match="tempo máximo"):
        _client(session, sleeps).poll(BlastSubmission("RID-SLOW", 0), params)

    assert sleeps == [20.0, 1, 1]


@pytest.mark.parametrize(
    "response, message",
    [
        (FakeResponse("RTOE = 1"), "RID válido"),
        (FakeResponse("RID = RID-NO-RTOE"), "RTOE válido"),
        (FakeResponse("erro", status_code=503), "erro HTTP"),
    ],
)
def test_submit_rejects_missing_rid_rtoe_and_http_errors(response, message):
    session = FakeSession(post=[response])

    with pytest.raises(NcbiBlastError, match=message):
        _client(session).submit(">pair_0001_F\nACGT\n", NcbiBlastParams())


def test_network_exception_is_wrapped():
    session = FakeSession(post=[requests.ConnectionError("offline")])

    with pytest.raises(NcbiBlastError, match="Falha HTTP"):
        _client(session).submit(">pair_0001_F\nACGT\n", NcbiBlastParams())
