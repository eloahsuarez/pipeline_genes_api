from __future__ import annotations

import re
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Callable, Iterable

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from models import ExonInterval, SequenceRecord


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass(slots=True)
class NcbiSearchParams:
    gene: str
    organism: str
    database: str = "nuccore"
    max_records: int = 100
    sequence_type: str = "mRNA"
    refseq_only: bool = True
    exclude_predicted: bool = True
    exclude_partial: bool = True
    require_gene_feature: bool = True
    min_length: int = 100
    max_length: int = 100000
    extra_query: str = ""


class NcbiClient:
    def __init__(
        self,
        email: str,
        api_key: str = "",
        tool: str = "GeneConservadoGUI",
        timeout: int = 60,
        log: Callable[[str], None] | None = None,
    ) -> None:
        if not email.strip():
            raise ValueError("O NCBI exige um e-mail válido nos parâmetros da ferramenta.")
        self.email = email.strip()
        self.api_key = api_key.strip()
        self.tool = tool
        self.timeout = timeout
        self.log = log or (lambda _: None)
        self.session = requests.Session()
        retry = Retry(
            total=4,
            connect=4,
            read=4,
            backoff_factor=0.8,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset({"GET", "POST"}),
        )
        self.session.mount("https://", HTTPAdapter(max_retries=retry))

    def _common(self) -> dict[str, str]:
        values = {"email": self.email, "tool": self.tool}
        if self.api_key:
            values["api_key"] = self.api_key
        return values

    @staticmethod
    def build_query(params: NcbiSearchParams) -> str:
        gene = params.gene.strip()
        organism = params.organism.strip()
        if not gene or not organism:
            raise ValueError("Informe o gene e o organismo.")

        parts = [f'"{gene}"[Gene Name]', f'"{organism}"[Organism]']
        sequence_type = params.sequence_type.lower()
        if sequence_type == "mrna":
            parts.append("biomol_mrna[PROP]")
        elif sequence_type == "genômico" or sequence_type == "genomico":
            parts.append("biomol_genomic[PROP]")
        elif sequence_type == "rna":
            parts.append("biomol_rna[PROP]")

        if params.refseq_only:
            parts.append("refseq[filter]")
        if params.exclude_predicted:
            parts.append("NOT predicted[Title]")
        if params.exclude_partial:
            parts.append("NOT partial[Title]")
        if params.extra_query.strip():
            parts.append(f"({params.extra_query.strip()})")
        return " AND ".join(parts)

    def search_ids(self, params: NcbiSearchParams) -> tuple[list[str], str, int]:
        query = self.build_query(params)
        payload = {
            **self._common(),
            "db": params.database,
            "term": query,
            "retmode": "json",
            "retmax": str(params.max_records),
            "idtype": "acc",
        }
        self.log(f"NCBI: pesquisando {params.database} com a consulta: {query}")
        response = self.session.get(
            f"{EUTILS_BASE}/esearch.fcgi", params=payload, timeout=self.timeout
        )
        response.raise_for_status()
        result = response.json().get("esearchresult", {})
        ids = result.get("idlist", [])
        count = int(result.get("count", 0))
        self.log(f"NCBI: {count} registro(s) encontrado(s); baixando até {len(ids)}.")
        return ids, query, count

    def fetch_records(
        self,
        ids: Iterable[str],
        params: NcbiSearchParams,
        batch_size: int = 100,
    ) -> list[SequenceRecord]:
        ids_list = list(ids)
        output: list[SequenceRecord] = []
        for offset in range(0, len(ids_list), batch_size):
            batch = ids_list[offset : offset + batch_size]
            payload = {
                **self._common(),
                "db": params.database,
                "id": ",".join(batch),
                "rettype": "gb",
                "retmode": "xml",
            }
            response = self.session.post(
                f"{EUTILS_BASE}/efetch.fcgi", data=payload, timeout=self.timeout
            )
            response.raise_for_status()
            output.extend(self._parse_genbank_xml(response.text, params))
            if offset + batch_size < len(ids_list):
                time.sleep(0.12 if self.api_key else 0.36)
        self.log(f"NCBI: {len(output)} sequência(s) permaneceram após os filtros locais.")
        return output

    def search_and_fetch(self, params: NcbiSearchParams) -> tuple[list[SequenceRecord], str, int]:
        ids, query, total = self.search_ids(params)
        if not ids:
            return [], query, total
        return self.fetch_records(ids, params), query, total

    @staticmethod
    def _feature_intervals(feature: ET.Element) -> list[tuple[int, int]]:
        """Extrai intervalos 1-based inclusivos de uma feature do GBSeq XML."""

        intervals: list[tuple[int, int]] = []
        for interval in feature.findall("./GBFeature_intervals/GBInterval"):
            start_text = interval.findtext("GBInterval_from", default="").strip()
            end_text = interval.findtext("GBInterval_to", default="").strip()
            point_text = interval.findtext("GBInterval_point", default="").strip()
            try:
                if start_text and end_text:
                    first, last = int(start_text), int(end_text)
                elif point_text:
                    first = last = int(point_text)
                else:
                    continue
            except ValueError:
                continue
            intervals.append((min(first, last), max(first, last)))

        if intervals:
            return intervals

        # Compatibilidade com respostas que tragam apenas GBFeature_location.
        location = feature.findtext("GBFeature_location", default="")
        for first, last in re.findall(r"<?(\d+)\.\.>?(\d+)", location):
            start, end = int(first), int(last)
            intervals.append((min(start, end), max(start, end)))
        if not intervals:
            match = re.fullmatch(r"[<>]?(\d+)", location.strip())
            if match:
                point = int(match.group(1))
                intervals.append((point, point))
        return intervals

    @staticmethod
    def _parse_genbank_xml(xml_text: str, params: NcbiSearchParams) -> list[SequenceRecord]:
        root = ET.fromstring(xml_text)
        records: list[SequenceRecord] = []
        requested_gene = params.gene.strip().casefold()

        for gbseq in root.findall(".//GBSeq"):
            def text(name: str, default: str = "") -> str:
                node = gbseq.find(name)
                return (node.text or default).strip() if node is not None else default

            accession = text("GBSeq_accession-version") or text("GBSeq_primary-accession")
            definition = text("GBSeq_definition")
            organism = text("GBSeq_organism")
            molecule_type = text("GBSeq_moltype")
            raw_seq = text("GBSeq_sequence").upper()
            nucleotide_letters = re.sub(r"[^A-Z]", "", raw_seq).replace("U", "T")
            # Símbolos IUPAC ambíguos precisam ocupar uma posição; removê-los
            # deslocaria todas as coordenadas de features posteriores.
            sequence = "".join(
                base if base in {"A", "C", "G", "T"} else "N"
                for base in nucleotide_letters
            )
            if not sequence:
                continue
            length = len(sequence)
            declared_length = text("GBSeq_length")
            if declared_length.isdigit() and int(declared_length) != length:
                continue
            if length < params.min_length or length > params.max_length:
                continue

            lowered_definition = definition.casefold()
            if params.exclude_predicted and "predicted" in lowered_definition:
                continue
            if params.exclude_partial and "partial" in lowered_definition:
                continue

            feature_keys: list[str] = []
            genes: list[str] = []
            exons: list[ExonInterval] = []
            for feature in gbseq.findall(".//GBFeature"):
                key_node = feature.find("GBFeature_key")
                feature_key = key_node.text.strip() if key_node is not None and key_node.text else ""
                if feature_key:
                    feature_keys.append(feature_key)
                qualifiers: dict[str, list[str]] = {}
                for qualifier in feature.findall(".//GBQualifier"):
                    qname = qualifier.findtext("GBQualifier_name", default="").strip()
                    qvalue = qualifier.findtext("GBQualifier_value", default="").strip()
                    if qname and qvalue:
                        qualifiers.setdefault(qname, []).append(qvalue)
                    if qname in {"gene", "locus_tag"} and qvalue:
                        genes.append(qvalue)

                if feature_key.casefold() == "exon":
                    feature_genes = {
                        value.casefold()
                        for name in ("gene", "locus_tag")
                        for value in qualifiers.get(name, [])
                    }
                    if feature_genes and requested_gene not in feature_genes:
                        continue
                    exon_number = next(iter(qualifiers.get("number", [])), "")
                    for start, end in NcbiClient._feature_intervals(feature):
                        exons.append(ExonInterval(start=start, end=end, number=exon_number))

            if params.sequence_type.lower() == "cds" and "CDS" not in feature_keys:
                continue
            if params.require_gene_feature:
                normalized = {g.casefold() for g in genes}
                if requested_gene not in normalized:
                    continue

            records.append(
                SequenceRecord(
                    uid=accession,
                    accession=accession,
                    definition=definition,
                    organism=organism,
                    sequence=sequence,
                    length=length,
                    genes=sorted(set(genes)),
                    feature_keys=sorted(set(feature_keys)),
                    exons=sorted(
                        {
                            (exon.start, exon.end, exon.number): exon
                            for exon in exons
                        }.values(),
                        key=lambda exon: (exon.start, exon.end, exon.number),
                    ),
                    molecule_type=molecule_type,
                )
            )
        return records
