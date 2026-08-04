from __future__ import annotations

import re
import time
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from itertools import combinations
from typing import Any, Callable, Iterable, Mapping, Sequence

import requests


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
REFSEQ_MRNA_DATABASE = "refseq_mrna"
REFSEQ_SELECT_DATABASE = "refseq_select_rna"
SUPPORTED_DATABASES = (REFSEQ_MRNA_DATABASE, REFSEQ_SELECT_DATABASE)
ORGANISM_NAME = "Homo sapiens"
ORGANISM_TAXID = 9606
ENTREZ_QUERY = "txid9606[ORGN]"

PAIR_QUERY_RE = re.compile(r"(?i)^(pair[_-]?(\d+))[_-]([FR])$")
IUPAC_DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.IGNORECASE)
ACCESSION_VERSION_RE = re.compile(r"\.\d+$")
STATUS_RE = re.compile(r"(?im)^\s*Status\s*=\s*([A-Z]+)\s*$")
RID_RE = re.compile(r"(?im)^\s*RID\s*=\s*(\S+)\s*$")
RTOE_RE = re.compile(r"(?im)^\s*RTOE\s*=\s*(\d+)\s*$")
BLAST_ERROR_MESSAGE_RE = re.compile(
    r"(?:^|\s)error\s*:|no alias or index file|"
    r"(?:blast\s+)?database.{0,120}(?:not found|does not exist|unavailable)|"
    r"\bfailed\b|\bfailure\b",
    re.IGNORECASE,
)


class NcbiBlastError(RuntimeError):
    """Erro esperado durante a análise remota de especificidade no NCBI."""


@dataclass(slots=True)
class NcbiBlastParams:
    top_pairs: int = 1
    hitlist_size: int = 50_000
    min_identity_pct: float = 0.0
    min_query_coverage_pct: float = 0.0
    max_amplicon: int = 4_000
    max_estimated_mismatches: int = 6
    three_prime_window: int = 5
    word_size: int = 7
    expect: float = 30_000.0
    target_gene: str = ""
    target_accessions: Sequence[str] = field(default_factory=tuple)
    database: str = REFSEQ_MRNA_DATABASE
    organism: str = ORGANISM_NAME
    taxid: int = ORGANISM_TAXID
    program: str = "blastn"
    megablast: bool = False
    short_query_adjust: bool = False
    poll_interval_seconds: float = 60.0
    poll_timeout_seconds: float = 900.0
    http_timeout_seconds: float = 60.0

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["target_accessions"] = [str(value) for value in self.target_accessions]
        return data


@dataclass(slots=True)
class PrimerPairQuery:
    pair_id: str
    pair_rank: int
    forward_sequence: str
    reverse_sequence: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class BlastSubmission:
    rid: str
    rtoe_seconds: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class NcbiPrimerHit:
    pair_id: str
    pair_rank: int
    primer_orientation: str
    query_id: str
    query_length: int
    accession: str
    title: str
    taxid: int | None
    scientific_name: str
    identity_pct: float
    query_coverage_pct: float
    identities: int
    alignment_length: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    query_strand: str
    subject_strand: str
    evalue: float
    bit_score: float
    aligned_query_sequence: str = ""
    aligned_subject_sequence: str = ""
    estimated_mismatches: int = 0
    estimated_three_prime_mismatches: int = 0

    @property
    def subject_three_prime(self) -> int:
        """Coordenada do 3′ do oligo na orientação indicada pelo HSP."""

        direction = 1 if self.subject_strand == "plus" else -1
        return self.subject_end + direction * (self.query_length - self.query_end)

    @property
    def subject_five_prime(self) -> int:
        direction = 1 if self.subject_strand == "plus" else -1
        return self.subject_start - direction * (self.query_start - 1)

    @property
    def covers_query_three_prime(self) -> bool:
        # O XML2 de blastn numera a query submetida no sentido 5′→3′. Para
        # inferir amplificação, a extremidade alinhada precisa alcançar a base
        # 3′ real do oligo; hits parciais continuam no relatório de hits, mas
        # não são promovidos a produto conjunto.
        return self.query_end == self.query_length

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class PairedAmplificationProduct:
    pair_id: str
    pair_rank: int
    accession: str
    gene: str
    description: str
    taxid: int | None
    scientific_name: str
    start: int
    end: int
    length: int
    classification: str
    primer_combination: str
    forward_hit: NcbiPrimerHit
    reverse_hit: NcbiPrimerHit
    left_hit: NcbiPrimerHit
    right_hit: NcbiPrimerHit

    @property
    def is_target(self) -> bool:
        return self.classification == "target"

    @property
    def title(self) -> str:
        return self.description

    @property
    def organism(self) -> str:
        return self.scientific_name

    def to_dict(self) -> dict[str, Any]:
        return {
            "pair_id": self.pair_id,
            "pair_rank": self.pair_rank,
            "accession": self.accession,
            "gene": self.gene,
            "title": self.title,
            "description": self.description,
            "taxid": self.taxid,
            "organism": self.organism,
            "scientific_name": self.scientific_name,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "classification": self.classification,
            "is_target": self.is_target,
            "primer_combination": self.primer_combination,
            "forward_hit": self.forward_hit.to_dict(),
            "reverse_hit": self.reverse_hit.to_dict(),
            "left_hit": self.left_hit.to_dict(),
            "right_hit": self.right_hit.to_dict(),
        }


@dataclass(slots=True)
class NcbiPairSpecificityResult:
    pair_id: str
    pair_rank: int
    forward_sequence: str
    reverse_sequence: str
    products: list[PairedAmplificationProduct] = field(default_factory=list)
    verdict: str = "Nenhum produto conjunto"

    @property
    def target_products(self) -> list[PairedAmplificationProduct]:
        return [product for product in self.products if product.is_target]

    @property
    def off_target_products(self) -> list[PairedAmplificationProduct]:
        return [product for product in self.products if not product.is_target]

    def to_dict(self) -> dict[str, Any]:
        return {
            "pair_id": self.pair_id,
            "pair_rank": self.pair_rank,
            "forward_sequence": self.forward_sequence,
            "reverse_sequence": self.reverse_sequence,
            "verdict": self.verdict,
            "product_count": len(self.products),
            "target_product_count": len(self.target_products),
            "off_target_product_count": len(self.off_target_products),
            "products": [product.to_dict() for product in self.products],
        }


@dataclass(slots=True)
class NcbiSpecificityReport:
    params: NcbiBlastParams
    submission: BlastSubmission
    submitted_fasta: str
    results: list[NcbiPairSpecificityResult]
    warnings: list[str] = field(default_factory=list)
    created_at_utc: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds")
    )

    @property
    def database(self) -> str:
        return self.params.database

    @property
    def organism(self) -> str:
        return self.params.organism

    @property
    def taxid(self) -> int:
        return self.params.taxid

    @property
    def query_fasta(self) -> str:
        return self.submitted_fasta

    @property
    def rid(self) -> str:
        return self.submission.rid

    def to_dict(self) -> dict[str, Any]:
        return {
            "created_at_utc": self.created_at_utc,
            "database": self.database,
            "organism": self.organism,
            "taxid": self.taxid,
            "rid": self.rid,
            "params": self.params.to_dict(),
            "submission": self.submission.to_dict(),
            "query_fasta": self.query_fasta,
            "submitted_fasta": self.submitted_fasta,
            "warnings": list(self.warnings),
            "results": [result.to_dict() for result in self.results],
        }


@dataclass(slots=True)
class ParsedBlastXml:
    hits: list[NcbiPrimerHit]
    raw_hit_counts: dict[str, int]
    warnings: list[str] = field(default_factory=list)


def _validate_params(params: NcbiBlastParams) -> None:
    if params.top_pairs < 1:
        raise NcbiBlastError("A quantidade de pares deve ser maior que zero.")
    if params.hitlist_size < 1:
        raise NcbiBlastError("A quantidade máxima de hits deve ser maior que zero.")
    if not 0.0 <= params.min_identity_pct <= 100.0:
        raise NcbiBlastError("A identidade BLAST deve estar entre 0% e 100%.")
    if not 0.0 <= params.min_query_coverage_pct <= 100.0:
        raise NcbiBlastError("A cobertura BLAST deve estar entre 0% e 100%.")
    if params.max_amplicon < 1:
        raise NcbiBlastError("O tamanho máximo do amplicon deve ser maior que zero.")
    if params.max_estimated_mismatches < 0:
        raise NcbiBlastError("O máximo de diferenças estimadas não pode ser negativo.")
    if params.three_prime_window < 1:
        raise NcbiBlastError("A janela 3′ deve ter pelo menos uma base.")
    if params.word_size not in {7, 11, 15}:
        raise NcbiBlastError("O word size remoto deve ser 7, 11 ou 15.")
    if params.expect <= 0:
        raise NcbiBlastError("O E-value máximo do BLAST deve ser maior que zero.")
    if params.poll_interval_seconds <= 0 or params.poll_timeout_seconds <= 0:
        raise NcbiBlastError("Os tempos de polling devem ser maiores que zero.")
    if params.http_timeout_seconds <= 0:
        raise NcbiBlastError("O tempo máximo HTTP deve ser maior que zero.")
    if params.database not in SUPPORTED_DATABASES:
        raise NcbiBlastError(
            "Banco BLAST não suportado. Use refseq_mrna ou refseq_select_rna."
        )
    if params.program.casefold() != "blastn":
        raise NcbiBlastError("Esta análise usa obrigatoriamente o programa blastn.")
    if params.organism.casefold() != ORGANISM_NAME.casefold() or params.taxid != ORGANISM_TAXID:
        raise NcbiBlastError("Esta análise está configurada para Homo sapiens (taxid 9606).")


def _primer_sequence(pair: object, side: str) -> str:
    candidate = getattr(pair, side, None)
    sequence = "".join(str(getattr(candidate, "sequence", "") or "").split()).upper()
    label = "Forward" if side == "forward" else "Reverse"
    if not sequence:
        raise NcbiBlastError(f"O primer {label} está ausente ou vazio.")
    if not IUPAC_DNA_RE.fullmatch(sequence):
        raise NcbiBlastError(f"O primer {label} contém bases inválidas.")
    return sequence


def primer_pairs_to_fasta(
    pairs: Iterable[object], top_pairs: int
) -> tuple[str, list[PrimerPairQuery]]:
    """Converte os primeiros pares ranqueados em consultas FASTA F/R."""

    if top_pairs < 1:
        raise NcbiBlastError("A quantidade de pares deve ser maior que zero.")
    selected = list(pairs)[:top_pairs]
    if not selected:
        raise NcbiBlastError("Nenhum par de primers foi informado para a análise.")

    queries: list[PrimerPairQuery] = []
    used_pair_ids: set[str] = set()
    chunks: list[str] = []
    for position, pair in enumerate(selected, start=1):
        try:
            rank = int(getattr(pair, "rank"))
        except (AttributeError, TypeError, ValueError) as exc:
            raise NcbiBlastError(
                f"O par na posição {position} não possui um rank inteiro válido."
            ) from exc
        if rank < 1:
            raise NcbiBlastError("O rank de cada par deve ser maior que zero.")
        pair_id = f"pair_{rank:04d}"
        if pair_id in used_pair_ids:
            raise NcbiBlastError("Os pares de primers possuem ranks duplicados.")
        used_pair_ids.add(pair_id)
        forward = _primer_sequence(pair, "forward")
        reverse = _primer_sequence(pair, "reverse")
        query = PrimerPairQuery(
            pair_id=pair_id,
            pair_rank=rank,
            forward_sequence=forward,
            reverse_sequence=reverse,
        )
        queries.append(query)
        chunks.extend(
            [
                f">{pair_id}_F\n{forward}\n",
                f">{pair_id}_R\n{reverse}\n",
            ]
        )
    return "".join(chunks), queries


# Nome mantido para consumidores da primeira versão isolada deste módulo.
build_primer_pairs_fasta = primer_pairs_to_fasta


def _local_name(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1].casefold()


def _direct_child(element: ET.Element, name: str) -> ET.Element | None:
    expected = name.casefold()
    return next((child for child in element if _local_name(child) == expected), None)


def _direct_text(element: ET.Element, name: str, default: str = "") -> str:
    child = _direct_child(element, name)
    if child is None or child.text is None:
        return default
    return child.text.strip()


def _first_descendant(element: ET.Element, name: str) -> ET.Element | None:
    expected = name.casefold()
    return next((child for child in element.iter() if _local_name(child) == expected), None)


def _descendant_text(element: ET.Element, name: str, default: str = "") -> str:
    child = _first_descendant(element, name)
    if child is None or child.text is None:
        return default
    return child.text.strip()


def _required_int(element: ET.Element, name: str, context: str) -> int:
    raw = _direct_text(element, name)
    try:
        return int(raw)
    except (TypeError, ValueError) as exc:
        raise NcbiBlastError(f"O XML2 contém {name} inválido em {context}.") from exc


def _required_float(element: ET.Element, name: str, context: str) -> float:
    raw = _direct_text(element, name)
    try:
        return float(raw)
    except (TypeError, ValueError) as exc:
        raise NcbiBlastError(f"O XML2 contém {name} inválido em {context}.") from exc


def _query_identity(search: ET.Element) -> tuple[str, int, str, int]:
    title = _direct_text(search, "query-title")
    query_id_field = _direct_text(search, "query-id")
    candidates = [title.split(maxsplit=1)[0] if title else "", query_id_field]
    match = None
    for candidate in candidates:
        if not candidate:
            continue
        match = PAIR_QUERY_RE.fullmatch(candidate)
        if match is not None:
            break
    if match is None:
        label = title or query_id_field or "sem identificador"
        raise NcbiBlastError(f"O XML2 retornou uma consulta desconhecida: {label}.")
    rank = int(match.group(2))
    pair_id = f"pair_{rank:04d}"
    orientation = match.group(3).upper()
    return f"{pair_id}_{orientation}", rank, orientation, _required_int(
        search, "query-len", label if (label := title or query_id_field) else "consulta"
    )


def _parse_taxid(value: str) -> int | None:
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _strand(explicit: str, frame: str, start: int, end: int) -> str:
    normalized = explicit.strip().casefold()
    if normalized in {"plus", "+", "+1", "1"}:
        return "plus"
    if normalized in {"minus", "-", "-1"}:
        return "minus"
    try:
        numeric_frame = int(frame)
    except (TypeError, ValueError):
        numeric_frame = 0
    if numeric_frame < 0:
        return "minus"
    if numeric_frame > 0:
        return "plus"
    return "plus" if start <= end else "minus"


def _hit_descriptors(hit_element: ET.Element) -> list[tuple[str, str, int | None, str]]:
    descriptor_elements = [
        element for element in hit_element.iter() if _local_name(element) == "hitdescr"
    ]
    if not descriptor_elements:
        descriptor_elements = [hit_element]
    descriptors: list[tuple[str, str, int | None, str]] = []
    for descriptor in descriptor_elements:
        accession = _descendant_text(descriptor, "accession")
        if not accession:
            accession = _descendant_text(descriptor, "id")
            if "|" in accession:
                tokens = [token for token in accession.split("|") if token]
                accession = tokens[-1] if tokens else accession
        title = _descendant_text(descriptor, "title")
        taxid = _parse_taxid(_descendant_text(descriptor, "taxid"))
        scientific_name = _descendant_text(descriptor, "sciname")
        if not accession:
            raise NcbiBlastError("O XML2 contém um hit sem accession.")
        descriptors.append((accession, title, taxid, scientific_name))
    return descriptors


def _element_text(element: ET.Element) -> str:
    return " ".join("".join(element.itertext()).split())


def _validate_blast_response(root: ET.Element) -> None:
    """Rejeita respostas que parecem XML válido, mas representam falha do BLAST."""

    for element in root.iter():
        local_name = _local_name(element)
        if local_name not in {"error", "message"}:
            continue
        message = _element_text(element)
        if local_name == "error" or (message and BLAST_ERROR_MESSAGE_RE.search(message)):
            details = message or "erro não detalhado"
            raise NcbiBlastError(f"O NCBI retornou erro no XML2: {details}")

    for field_name in ("db-num", "db-len"):
        for element in root.iter():
            if _local_name(element) != field_name:
                continue
            raw_value = _element_text(element)
            try:
                value = int(raw_value)
            except (TypeError, ValueError) as exc:
                raise NcbiBlastError(
                    f"O XML2 contém estatística de banco inválida ({field_name})."
                ) from exc
            if value <= 0:
                raise NcbiBlastError(
                    "O NCBI respondeu com um banco vazio ou indisponível "
                    f"({field_name}={value}); o resultado não é confiável."
                )


def _validate_search_statistics(search: ET.Element, query_id: str) -> None:
    """Exige evidência de que cada consulta pesquisou um banco disponível."""

    statistics = _first_descendant(search, "statistics")
    if statistics is None:
        raise NcbiBlastError(
            "A resposta XML2 do NCBI está incompleta: faltam estatísticas "
            f"do banco para a consulta {query_id}."
        )
    for field_name in ("db-num", "db-len"):
        raw_value = _descendant_text(statistics, field_name)
        try:
            value = int(raw_value)
        except (TypeError, ValueError) as exc:
            raise NcbiBlastError(
                "A resposta XML2 do NCBI está incompleta: estatística "
                f"{field_name} inválida para a consulta {query_id}."
            ) from exc
        if value <= 0:
            raise NcbiBlastError(
                "O NCBI respondeu com um banco vazio ou indisponível para a "
                f"consulta {query_id} ({field_name}={value}); o resultado não "
                "é confiável."
            )


def _estimated_mismatch_counts(
    *,
    query_length: int,
    query_start: int,
    query_end: int,
    aligned_query: str,
    aligned_subject: str,
    identities: int,
    alignment_length: int,
    three_prime_window: int,
) -> tuple[int, int]:
    """Estima diferenças no oligo inteiro a partir do HSP local do BLAST.

    Bases da query fora do HSP são tratadas conservadoramente como diferenças.
    Quando qseq/hseq estão presentes, substituições e gaps dentro do HSP são
    posicionados também na janela 3′. Isso não substitui o alinhamento global do
    Primer-BLAST, mas evita descartar sementes locais com pontas não alinhadas.
    """

    query_start = max(1, query_start)
    query_end = min(query_length, query_end)
    window_start = max(1, query_length - three_prime_window + 1)
    mismatch_positions = set(range(1, query_start))
    mismatch_positions.update(range(query_end + 1, query_length + 1))
    insertion_differences = 0
    insertion_three_prime_differences = 0

    if aligned_query and aligned_subject and len(aligned_query) == len(aligned_subject):
        query_position = query_start
        for query_base, subject_base in zip(aligned_query.upper(), aligned_subject.upper()):
            if query_base == "-":
                if subject_base != "-":
                    insertion_differences += 1
                    # Uma inserção no subject fica entre bases da query. Use a
                    # base adjacente seguinte e, quando o gap é terminal, a
                    # anterior; assim ela participa da janela 3′ sem inventar
                    # uma coordenada fora do oligo.
                    adjacent_query_position = min(
                        query_length,
                        max(1, query_position),
                    )
                    if adjacent_query_position >= window_start:
                        insertion_three_prime_differences += 1
                continue
            if subject_base == "-" or query_base != subject_base:
                mismatch_positions.add(query_position)
            query_position += 1
        internal_differences = sum(
            1
            for query_base, subject_base in zip(
                aligned_query.upper(), aligned_subject.upper()
            )
            if query_base != subject_base
        )
    else:
        internal_differences = max(0, alignment_length - identities)

    outside_hsp = (query_start - 1) + (query_length - query_end)
    total = max(
        len(mismatch_positions) + insertion_differences,
        outside_hsp + internal_differences,
    )
    three_prime = (
        sum(position >= window_start for position in mismatch_positions)
        + insertion_three_prime_differences
    )
    return total, three_prime


def parse_blast_xml2(
    xml_text: str,
    *,
    min_identity_pct: float = 0.0,
    min_query_coverage_pct: float = 0.0,
    expected_taxid: int = ORGANISM_TAXID,
    three_prime_window: int = 5,
) -> ParsedBlastXml:
    """Lê XML2/XML2_S do BLAST, inclusive quando há namespace padrão."""

    if not xml_text.strip():
        return ParsedBlastXml(hits=[], raw_hit_counts={})
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        raise NcbiBlastError("A resposta XML2 do NCBI é inválida.") from exc

    _validate_blast_response(root)

    hits: list[NcbiPrimerHit] = []
    raw_hit_counts: dict[str, int] = {}
    missing_taxid_accessions: set[str] = set()
    foreign_taxid_hits = 0
    searches = [
        element
        for element in root.iter()
        if _local_name(element) == "search"
        and _direct_child(element, "query-len") is not None
    ]
    if not searches:
        raise NcbiBlastError("A resposta XML2 do NCBI não contém consultas BLAST.")
    for search in searches:
        query_id, rank, orientation, query_length = _query_identity(search)
        if query_length < 1:
            raise NcbiBlastError(f"O XML2 retornou tamanho inválido para {query_id}.")
        _validate_search_statistics(search, query_id)
        if query_id in raw_hit_counts:
            raise NcbiBlastError(
                f"O XML2 retornou a consulta {query_id} mais de uma vez."
            )
        pair_id = f"pair_{rank:04d}"
        hit_elements = [
            element for element in search.iter() if _local_name(element) == "hit"
        ]
        raw_hit_counts[query_id] = len(hit_elements)
        for hit_element in hit_elements:
            hsp_elements = [
                element for element in hit_element.iter() if _local_name(element) == "hsp"
            ]
            for accession, title, taxid, scientific_name in _hit_descriptors(hit_element):
                if taxid is not None and taxid != expected_taxid:
                    foreign_taxid_hits += 1
                    continue
                if taxid is None:
                    missing_taxid_accessions.add(accession)
                for hsp_index, hsp in enumerate(hsp_elements, start=1):
                    context = f"{query_id}, hit {accession}, HSP {hsp_index}"
                    query_start = _required_int(hsp, "query-from", context)
                    query_end = _required_int(hsp, "query-to", context)
                    subject_start = _required_int(hsp, "hit-from", context)
                    subject_end = _required_int(hsp, "hit-to", context)
                    identities = _required_int(hsp, "identity", context)
                    alignment_length = _required_int(hsp, "align-len", context)
                    if alignment_length < 1:
                        raise NcbiBlastError(f"O XML2 contém alinhamento vazio em {context}.")
                    query_span = abs(query_end - query_start) + 1
                    query_coverage_pct = 100.0 * query_span / query_length
                    identity_pct = 100.0 * identities / alignment_length
                    aligned_query = _direct_text(hsp, "qseq")
                    aligned_subject = _direct_text(hsp, "hseq")
                    estimated_mismatches, estimated_three_prime_mismatches = (
                        _estimated_mismatch_counts(
                            query_length=query_length,
                            query_start=query_start,
                            query_end=query_end,
                            aligned_query=aligned_query,
                            aligned_subject=aligned_subject,
                            identities=identities,
                            alignment_length=alignment_length,
                            three_prime_window=three_prime_window,
                        )
                    )
                    hit = NcbiPrimerHit(
                        pair_id=pair_id,
                        pair_rank=rank,
                        primer_orientation=orientation,
                        query_id=query_id,
                        query_length=query_length,
                        accession=accession,
                        title=title,
                        taxid=taxid,
                        scientific_name=scientific_name,
                        identity_pct=identity_pct,
                        query_coverage_pct=query_coverage_pct,
                        identities=identities,
                        alignment_length=alignment_length,
                        query_start=query_start,
                        query_end=query_end,
                        subject_start=subject_start,
                        subject_end=subject_end,
                        query_strand=_strand(
                            _direct_text(hsp, "query-strand"),
                            _direct_text(hsp, "query-frame"),
                            query_start,
                            query_end,
                        ),
                        subject_strand=_strand(
                            _direct_text(hsp, "hit-strand"),
                            _direct_text(hsp, "hit-frame"),
                            subject_start,
                            subject_end,
                        ),
                        evalue=_required_float(hsp, "evalue", context),
                        bit_score=_required_float(hsp, "bit-score", context),
                        aligned_query_sequence=aligned_query,
                        aligned_subject_sequence=aligned_subject,
                        estimated_mismatches=estimated_mismatches,
                        estimated_three_prime_mismatches=(
                            estimated_three_prime_mismatches
                        ),
                    )
                    if (
                        identity_pct >= min_identity_pct
                        and query_coverage_pct >= min_query_coverage_pct
                    ):
                        hits.append(hit)

    warnings: list[str] = []
    if missing_taxid_accessions:
        accessions = ", ".join(sorted(missing_taxid_accessions))
        warnings.append(
            "Taxid ausente no XML2 para accession(s) mantidos na análise: " + accessions + "."
        )
    if foreign_taxid_hits:
        warnings.append(
            f"{foreign_taxid_hits} hit(s) fora de Homo sapiens foram descartados pelo taxid."
        )
    return ParsedBlastXml(hits=hits, raw_hit_counts=raw_hit_counts, warnings=warnings)


def _accession_without_version(accession: str) -> str:
    return ACCESSION_VERSION_RE.sub("", accession.strip()).upper()


def _gene_token_in_title(target_gene: str, title: str) -> bool:
    symbol = target_gene.strip()
    if not symbol or not title:
        return False
    pattern = re.compile(
        rf"(?<![A-Za-z0-9]){re.escape(symbol)}(?![A-Za-z0-9])", re.IGNORECASE
    )
    return pattern.search(title) is not None


def _infer_gene_symbol(title: str) -> str:
    patterns = (
        r"(?i)\[gene=([A-Za-z][A-Za-z0-9_.-]*)\]",
        r"(?i)\bgene\s*[:=]\s*([A-Za-z][A-Za-z0-9_.-]*)",
        r"\(([A-Za-z][A-Za-z0-9_.-]*)\)(?=,\s*(?:transcript|mRNA|RNA|ncRNA|rRNA))",
    )
    for pattern in patterns:
        match = re.search(pattern, title)
        if match:
            return match.group(1)
    return ""


def _product_is_target(
    accession: str,
    title: str,
    *,
    target_accessions: set[str],
    target_gene: str,
    inferred_gene: str = "",
) -> bool:
    if _accession_without_version(accession) in target_accessions:
        return True
    normalized_target = target_gene.strip().casefold()
    normalized_inferred = inferred_gene.strip().casefold()
    if normalized_inferred:
        return bool(normalized_target) and normalized_inferred == normalized_target
    return _gene_token_in_title(target_gene, title)


def _compatible_product_geometry(
    first_hit: NcbiPrimerHit,
    second_hit: NcbiPrimerHit,
    max_amplicon: int,
) -> tuple[int, int, int, NcbiPrimerHit, NcbiPrimerHit] | None:
    if first_hit is second_hit:
        return None
    if first_hit.accession.casefold() != second_hit.accession.casefold():
        return None
    if first_hit.subject_strand == second_hit.subject_strand:
        return None

    plus_hit = first_hit if first_hit.subject_strand == "plus" else second_hit
    minus_hit = first_hit if first_hit.subject_strand == "minus" else second_hit
    if plus_hit.subject_strand != "plus" or minus_hit.subject_strand != "minus":
        return None
    # O oligo no strand plus estende para coordenadas crescentes; no minus,
    # para coordenadas decrescentes. Seus 3′ precisam apontar para o interior.
    if plus_hit.subject_three_prime >= minus_hit.subject_three_prime:
        return None

    start = plus_hit.subject_five_prime
    end = minus_hit.subject_five_prime
    if start < 1 or end < start:
        return None
    length = end - start + 1
    if length > max_amplicon:
        return None
    return start, end, length, plus_hit, minus_hit


def _candidate_groups_by_accession(
    hits: Iterable[NcbiPrimerHit],
    pair_id: str,
) -> Iterable[
    tuple[str, Iterable[tuple[NcbiPrimerHit, NcbiPrimerHit]]]
]:
    """Gera combinações somente dentro do mesmo acesso.

    Além de ser biologicamente necessário, o agrupamento evita um produto
    cartesiano quadrático entre dezenas de milhares de hits de acessos distintos.
    """

    grouped: dict[str, dict[str, list[NcbiPrimerHit]]] = {}
    for hit in hits:
        if hit.pair_id != pair_id or hit.primer_orientation not in {"F", "R"}:
            continue
        orientations = grouped.setdefault(
            hit.accession.casefold(),
            {"F": [], "R": []},
        )
        orientations[hit.primer_orientation].append(hit)

    for orientations in grouped.values():
        forward_hits = orientations["F"]
        reverse_hits = orientations["R"]
        yield (
            "F-R",
            (
                (forward_hit, reverse_hit)
                for forward_hit in forward_hits
                for reverse_hit in reverse_hits
            ),
        )
        yield "F-F", combinations(forward_hits, 2)
        yield "R-R", combinations(reverse_hits, 2)


def correlate_primer_hits(
    hits: Iterable[NcbiPrimerHit],
    queries: Sequence[PrimerPairQuery],
    params: NcbiBlastParams,
) -> list[NcbiPairSpecificityResult]:
    """Cruza F-R, F-F e R-R no mesmo acesso em geometria amplificável."""

    hits_list = list(hits)
    target_accessions = {
        _accession_without_version(str(accession))
        for accession in params.target_accessions
        if str(accession).strip()
    }
    results: list[NcbiPairSpecificityResult] = []
    for query in queries:
        products_by_key: dict[
            tuple[str, int, int, str], PairedAmplificationProduct
        ] = {}
        for primer_combination, candidate_pairs in _candidate_groups_by_accession(
            hits_list,
            query.pair_id,
        ):
            for first_hit, second_hit in candidate_pairs:
                if first_hit is second_hit:
                    continue
                if (
                    first_hit.estimated_mismatches
                    > params.max_estimated_mismatches
                    or second_hit.estimated_mismatches
                    > params.max_estimated_mismatches
                ):
                    continue
                geometry = _compatible_product_geometry(
                    first_hit, second_hit, params.max_amplicon
                )
                if geometry is None:
                    continue
                start, end, length, left_hit, right_hit = geometry
                description = first_hit.title or second_hit.title
                gene = _infer_gene_symbol(description)
                is_target = _product_is_target(
                    first_hit.accession,
                    description,
                    target_accessions=target_accessions,
                    target_gene=params.target_gene,
                    inferred_gene=gene,
                )
                if is_target and not gene:
                    gene = params.target_gene.strip()
                if primer_combination == "F-R":
                    forward_hit = (
                        first_hit
                        if first_hit.primer_orientation == "F"
                        else second_hit
                    )
                    reverse_hit = (
                        first_hit
                        if first_hit.primer_orientation == "R"
                        else second_hit
                    )
                else:
                    # Compatibilidade com consumidores antigos: em F-F/R-R estes
                    # campos representam os dois sítios, enquanto left/right
                    # registram inequivocamente a geometria física do amplicon.
                    forward_hit, reverse_hit = first_hit, second_hit
                product = PairedAmplificationProduct(
                    pair_id=query.pair_id,
                    pair_rank=query.pair_rank,
                    accession=first_hit.accession,
                    gene=gene,
                    description=description,
                    taxid=(
                        first_hit.taxid
                        if first_hit.taxid is not None
                        else second_hit.taxid
                    ),
                    scientific_name=(
                        first_hit.scientific_name or second_hit.scientific_name
                    ),
                    start=start,
                    end=end,
                    length=length,
                    classification="target" if is_target else "off_target",
                    primer_combination=primer_combination,
                    forward_hit=forward_hit,
                    reverse_hit=reverse_hit,
                    left_hit=left_hit,
                    right_hit=right_hit,
                )
                key = (
                    first_hit.accession.casefold(),
                    start,
                    end,
                    primer_combination,
                )
                existing = products_by_key.get(key)
                if existing is None:
                    products_by_key[key] = product
                    continue
                existing_score = (
                    existing.left_hit.bit_score + existing.right_hit.bit_score
                )
                new_score = left_hit.bit_score + right_hit.bit_score
                existing_quality = (
                    min(
                        existing.left_hit.query_coverage_pct,
                        existing.right_hit.query_coverage_pct,
                    ),
                    min(existing.left_hit.identity_pct, existing.right_hit.identity_pct),
                    existing_score,
                    -max(existing.left_hit.evalue, existing.right_hit.evalue),
                )
                new_quality = (
                    min(left_hit.query_coverage_pct, right_hit.query_coverage_pct),
                    min(left_hit.identity_pct, right_hit.identity_pct),
                    new_score,
                    -max(left_hit.evalue, right_hit.evalue),
                )
                if new_quality > existing_quality:
                    products_by_key[key] = product

        products = sorted(
            products_by_key.values(),
            key=lambda product: (
                not product.is_target,
                product.accession.casefold(),
                product.start,
                product.end,
                product.primer_combination,
            ),
        )
        has_target = any(product.is_target for product in products)
        has_off_target = any(not product.is_target for product in products)
        if has_target and has_off_target:
            verdict = "Pode amplificar outro gene"
        elif has_target:
            verdict = "Específico no banco"
        elif has_off_target:
            verdict = "Apenas fora do alvo"
        else:
            verdict = "Nenhum produto conjunto"
        results.append(
            NcbiPairSpecificityResult(
                pair_id=query.pair_id,
                pair_rank=query.pair_rank,
                forward_sequence=query.forward_sequence,
                reverse_sequence=query.reverse_sequence,
                products=products,
                verdict=verdict,
            )
        )
    return results


Sleeper = Callable[[float], None]


class NcbiBlastClient:
    """Cliente do Common URL API para BLAST de oligos em RefSeq RNA humano."""

    def __init__(
        self,
        email: str,
        *,
        tool: str = "GeneConservadoGUI",
        session: requests.Session | None = None,
        sleeper: Sleeper = time.sleep,
        log: Callable[[str], None] | None = None,
    ) -> None:
        if not email.strip():
            raise ValueError("O NCBI exige um e-mail válido nos parâmetros da ferramenta.")
        if not tool.strip():
            raise ValueError("Informe o nome da ferramenta para o NCBI.")
        self.email = email.strip()
        self.tool = tool.strip()
        self.session = session or requests.Session()
        self.sleeper = sleeper
        self.log = log or (lambda _message: None)

    @staticmethod
    def _response_text(response: Any, operation: str) -> str:
        try:
            response.raise_for_status()
        except requests.RequestException as exc:
            raise NcbiBlastError(f"O NCBI retornou erro HTTP durante {operation}: {exc}") from exc
        except Exception as exc:
            raise NcbiBlastError(f"O NCBI retornou erro HTTP durante {operation}: {exc}") from exc
        return str(getattr(response, "text", "") or "")

    def submit(self, fasta: str, params: NcbiBlastParams) -> BlastSubmission:
        _validate_params(params)
        payload = {
            "CMD": "Put",
            "PROGRAM": "blastn",
            "DATABASE": params.database,
            "ENTREZ_QUERY": ENTREZ_QUERY,
            "HITLIST_SIZE": str(params.hitlist_size),
            "email": self.email,
            "tool": self.tool,
            "QUERY": fasta,
        }
        if params.megablast:
            payload.update(
                {
                    "MEGABLAST": "on",
                    "SHORT_QUERY_ADJUST": (
                        "true" if params.short_query_adjust else "false"
                    ),
                }
            )
            profile = "megablast rápido"
        else:
            payload.update(
                {
                    "SHORT_QUERY_ADJUST": "false",
                    "WORD_SIZE": str(params.word_size),
                    "EXPECT": f"{params.expect:g}",
                    "NUCL_REWARD": "1",
                    "NUCL_PENALTY": "-3",
                    "GAPCOSTS": "5 2",
                    "FILTER": "F",
                }
            )
            profile = "blastn sensível para oligos curtos"
        self.log(
            f"NCBI BLAST: submetendo os primers ao banco {params.database} humano "
            f"com {profile}."
        )
        try:
            response = self.session.post(
                BLAST_URL,
                data=payload,
                timeout=params.http_timeout_seconds,
            )
        except requests.RequestException as exc:
            raise NcbiBlastError(f"Falha HTTP ao submeter a busca ao NCBI: {exc}") from exc
        text = self._response_text(response, "a submissão")
        rid_match = RID_RE.search(text)
        rtoe_match = RTOE_RE.search(text)
        if rid_match is None:
            raise NcbiBlastError("O NCBI não retornou um RID válido para a busca.")
        if rtoe_match is None:
            raise NcbiBlastError("O NCBI não retornou um RTOE válido para a busca.")
        return BlastSubmission(rid=rid_match.group(1), rtoe_seconds=int(rtoe_match.group(1)))

    def poll(self, submission: BlastSubmission, params: NcbiBlastParams) -> str:
        _validate_params(params)
        elapsed = 0.0
        initial_delay = max(float(submission.rtoe_seconds), 20.0)
        if initial_delay >= params.poll_timeout_seconds:
            raise NcbiBlastError("A busca BLAST excedeu o tempo máximo de polling.")
        self.sleeper(initial_delay)
        elapsed += initial_delay

        payload = {
            "CMD": "Get",
            "RID": submission.rid,
            "FORMAT_TYPE": "XML2_S",
            "HITLIST_SIZE": str(params.hitlist_size),
            "email": self.email,
            "tool": self.tool,
        }
        while True:
            try:
                response = self.session.get(
                    BLAST_URL,
                    params=payload,
                    timeout=params.http_timeout_seconds,
                )
            except requests.RequestException as exc:
                raise NcbiBlastError(f"Falha HTTP ao consultar a busca no NCBI: {exc}") from exc
            text = self._response_text(response, "o polling")
            status_match = STATUS_RE.search(text)
            status = status_match.group(1).upper() if status_match else ""
            if status == "WAITING":
                remaining = params.poll_timeout_seconds - elapsed
                if remaining <= 0:
                    raise NcbiBlastError("A busca BLAST excedeu o tempo máximo de polling.")
                delay = min(params.poll_interval_seconds, remaining)
                self.sleeper(delay)
                elapsed += delay
                if elapsed >= params.poll_timeout_seconds:
                    raise NcbiBlastError("A busca BLAST excedeu o tempo máximo de polling.")
                continue
            if status == "FAILED":
                raise NcbiBlastError("O NCBI informou que a busca BLAST falhou.")
            if status == "UNKNOWN":
                raise NcbiBlastError("O NCBI informou um RID desconhecido ou expirado.")
            if status and status != "READY":
                raise NcbiBlastError(f"O NCBI retornou um status BLAST inesperado: {status}.")
            first_xml = text.find("<")
            if first_xml >= 0:
                return text[first_xml:]
            if "ThereAreHits=no" in text:
                return ""
            if status == "READY":
                # Algumas respostas informam READY antes de devolver o XML final.
                remaining = params.poll_timeout_seconds - elapsed
                if remaining <= 0:
                    raise NcbiBlastError("A busca BLAST excedeu o tempo máximo de polling.")
                delay = min(params.poll_interval_seconds, remaining)
                self.sleeper(delay)
                elapsed += delay
                continue
            raise NcbiBlastError("O NCBI retornou uma resposta de polling não reconhecida.")

    def analyze(
        self,
        pairs: Iterable[object],
        params: NcbiBlastParams,
    ) -> NcbiSpecificityReport:
        _validate_params(params)
        fasta, queries = primer_pairs_to_fasta(pairs, params.top_pairs)
        submission = self.submit(fasta, params)
        xml_text = self.poll(submission, params)
        parsed = parse_blast_xml2(
            xml_text,
            min_identity_pct=(params.min_identity_pct if params.megablast else 0.0),
            min_query_coverage_pct=(
                params.min_query_coverage_pct if params.megablast else 0.0
            ),
            expected_taxid=params.taxid,
            three_prime_window=params.three_prime_window,
        )
        known_queries = {
            f"{query.pair_id}_{orientation}"
            for query in queries
            for orientation in ("F", "R")
        }
        returned_queries = set(parsed.raw_hit_counts)
        unknown_queries = sorted(returned_queries - known_queries)
        if unknown_queries:
            raise NcbiBlastError(
                "O NCBI retornou consultas que não pertencem aos pares submetidos: "
                + ", ".join(unknown_queries)
                + "."
            )

        # `poll()` retorna texto vazio apenas quando o NCBI declara
        # `ThereAreHits=no`, que é uma conclusão válida. Quando há XML, porém,
        # todas as consultas submetidas devem estar presentes: uma resposta
        # parcial não pode ser interpretada silenciosamente como ausência de hit.
        if xml_text.strip():
            missing_queries = sorted(known_queries - returned_queries)
            if missing_queries:
                raise NcbiBlastError(
                    "A resposta do NCBI está incompleta; faltam as consultas: "
                    + ", ".join(missing_queries)
                    + "."
                )

        if params.database == REFSEQ_MRNA_DATABASE:
            database_warning = (
                "RefSeq mRNA inclui múltiplos transcritos e isoformas anotadas, mas "
                "não cobre todas as variantes nem possíveis produtos no genoma."
            )
        else:
            database_warning = (
                "RefSeq Select contém transcritos representativos e não cobre todas "
                "as isoformas nem possíveis produtos no genoma."
            )
        if params.megablast:
            profile_warning = (
                "O perfil megablast prioriza correspondências quase exatas e pode "
                "omitir off-targets com vários mismatches."
            )
        else:
            profile_warning = (
                "O perfil sensível usa word size 7 e extrapola as pontas ausentes "
                "do HSP como diferenças estimadas; ele verifica F-R, F-F e R-R, "
                "mas não substitui o alinhamento global do Primer-BLAST."
            )
        warnings = [database_warning, profile_warning, *parsed.warnings]
        truncated = sorted(
            query_id
            for query_id, count in parsed.raw_hit_counts.items()
            if count >= params.hitlist_size
        )
        if truncated:
            warnings.append(
                "A lista de hits pode ter sido truncada no limite configurado para: "
                + ", ".join(truncated)
                + "."
            )
        results = correlate_primer_hits(parsed.hits, queries, params)
        truncated_pairs = {
            query_id.rsplit("_", 1)[0]
            for query_id in truncated
        }
        for result in results:
            if result.pair_id in truncated_pairs and not result.off_target_products:
                result.verdict = "Inconclusivo: limite de hits atingido"
        self.log(
            f"NCBI BLAST: {len(results)} par(es) analisado(s) em {params.database}."
        )
        return NcbiSpecificityReport(
            params=params,
            submission=submission,
            submitted_fasta=fasta,
            results=results,
            warnings=warnings,
        )


__all__ = [
    "BLAST_URL",
    "ENTREZ_QUERY",
    "NcbiBlastClient",
    "NcbiBlastError",
    "NcbiBlastParams",
    "NcbiSpecificityReport",
    "NcbiPairSpecificityResult",
    "NcbiPrimerHit",
    "ORGANISM_NAME",
    "ORGANISM_TAXID",
    "PairedAmplificationProduct",
    "ParsedBlastXml",
    "PrimerPairQuery",
    "REFSEQ_MRNA_DATABASE",
    "REFSEQ_SELECT_DATABASE",
    "SUPPORTED_DATABASES",
    "BlastSubmission",
    "build_primer_pairs_fasta",
    "correlate_primer_hits",
    "parse_blast_xml2",
    "primer_pairs_to_fasta",
]


# Alias transitório para código que tenha importado o nome inicial.
NcbiBlastReport = NcbiSpecificityReport
