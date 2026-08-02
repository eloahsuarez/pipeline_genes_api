from __future__ import annotations

import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Callable

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"


@dataclass(slots=True)
class ClustalParams:
    email: str
    title: str = "Gene conserved regions"
    sequence_type: str = "dna"
    out_format: str = "fa"
    iterations: int = 0
    dealign: bool = False
    mbed: bool = True
    mbed_iteration: bool = True
    order: str = "aligned"
    poll_seconds: int = 3
    timeout_minutes: int = 20


class EbiClustalClient:
    def __init__(self, log: Callable[[str], None] | None = None, timeout: int = 90) -> None:
        self.log = log or (lambda _: None)
        self.timeout = timeout
        self.session = requests.Session()
        retry = Retry(
            total=4,
            connect=4,
            read=4,
            backoff_factor=1,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset({"GET", "POST"}),
        )
        self.session.mount("https://", HTTPAdapter(max_retries=retry))
        self.session.headers.update({"User-Agent": "GeneConservadoGUI/1.0"})

    def submit(self, fasta_text: str, params: ClustalParams) -> str:
        if not params.email.strip():
            raise ValueError("O EMBL-EBI exige um e-mail válido.")
        if fasta_text.count(">") < 3:
            raise ValueError("O Clustal Omega requer pelo menos três sequências para alinhamento múltiplo.")
        if len(fasta_text.encode("utf-8")) > 4_000_000:
            raise ValueError("O FASTA ultrapassa o limite público de 4 MB do Clustal Omega do EBI.")

        payload = {
            "email": params.email.strip(),
            "title": params.title,
            "stype": params.sequence_type,
            "sequence": fasta_text,
            "outfmt": params.out_format,
            "iterations": str(params.iterations),
            "dealign": str(params.dealign).lower(),
            "mbed": str(params.mbed).lower(),
            "mbediteration": str(params.mbed_iteration).lower(),
            "order": params.order,
            "addformats": "true",
        }
        self.log("EBI: enviando as sequências ao Clustal Omega.")
        response = self.session.post(f"{BASE_URL}/run/", data=payload, timeout=self.timeout)
        response.raise_for_status()
        job_id = response.text.strip()
        if not job_id or "ERROR" in job_id.upper():
            raise RuntimeError(f"O EBI não retornou um jobId válido: {job_id}")
        self.log(f"EBI: tarefa criada: {job_id}")
        return job_id

    def status(self, job_id: str) -> str:
        response = self.session.get(f"{BASE_URL}/status/{job_id}", timeout=self.timeout)
        response.raise_for_status()
        return response.text.strip()

    def wait(self, job_id: str, params: ClustalParams) -> str:
        deadline = time.monotonic() + params.timeout_minutes * 60
        last = ""
        while time.monotonic() < deadline:
            current = self.status(job_id)
            if current != last:
                self.log(f"EBI: estado da tarefa = {current}")
                last = current
            if current == "FINISHED":
                return current
            if current in {"ERROR", "FAILURE", "NOT_FOUND"}:
                raise RuntimeError(f"A tarefa do Clustal terminou com estado {current}.")
            time.sleep(max(1, params.poll_seconds))
        raise TimeoutError("O alinhamento excedeu o tempo máximo configurado.")

    def result_types(self, job_id: str) -> list[dict[str, str]]:
        response = self.session.get(f"{BASE_URL}/resulttypes/{job_id}", timeout=self.timeout)
        response.raise_for_status()
        root = ET.fromstring(response.text)
        result: list[dict[str, str]] = []
        for node in root.findall(".//type"):
            item: dict[str, str] = {}
            for child in list(node):
                item[child.tag] = (child.text or "").strip()
            result.append(item)
        return result

    def get_result(self, job_id: str, result_type: str) -> bytes:
        response = self.session.get(
            f"{BASE_URL}/result/{job_id}/{result_type}", timeout=self.timeout
        )
        response.raise_for_status()
        return response.content

    def run(self, fasta_text: str, params: ClustalParams) -> tuple[str, str, str]:
        job_id = self.submit(fasta_text, params)
        self.wait(job_id, params)
        types = self.result_types(job_id)
        identifiers = [item.get("identifier", "") for item in types]
        preferred = next(
            (identifier for identifier in identifiers if "aln" in identifier and "fasta" in identifier),
            None,
        )
        if preferred is None:
            preferred = next((identifier for identifier in identifiers if "aln" in identifier), None)
        if preferred is None and identifiers:
            preferred = identifiers[0]
        if not preferred:
            raise RuntimeError("O EBI não informou nenhum formato de resultado.")
        content = self.get_result(job_id, preferred)
        try:
            alignment = content.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise RuntimeError("O resultado escolhido pelo EBI não é um alinhamento textual.") from exc
        self.log(f"EBI: alinhamento baixado no formato {preferred}.")
        return job_id, preferred, alignment
