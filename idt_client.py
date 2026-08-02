from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


TOKEN_URL = "https://www.idtdna.com/Identityserver/connect/token"
OLIGO_BASE = "https://www.idtdna.com/Restapi/v1/OligoAnalyzer"


@dataclass(slots=True)
class IdtCredentials:
    client_id: str
    client_secret: str
    username: str
    password: str


@dataclass(slots=True)
class IdtConditions:
    na_mm: float = 50.0
    mg_mm: float = 3.0
    dntp_mm: float = 0.8
    oligo_um: float = 0.25
    folding_temp_c: float = 37.0
    nucleotide_type: str = "DNA"


class IdtClient:
    def __init__(
        self,
        credentials: IdtCredentials,
        log: Callable[[str], None] | None = None,
        timeout: int = 60,
    ) -> None:
        self.credentials = credentials
        self.log = log or (lambda _: None)
        self.timeout = timeout
        self.token = ""
        self.session = requests.Session()
        retry = Retry(
            total=3,
            connect=3,
            read=3,
            backoff_factor=0.8,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset({"POST"}),
        )
        self.session.mount("https://", HTTPAdapter(max_retries=retry))

    def authenticate(self) -> str:
        c = self.credentials
        if not all([c.client_id, c.client_secret, c.username, c.password]):
            raise ValueError("Preencha Client ID, Client Secret, usuário e senha da IDT.")
        payload = {
            "grant_type": "password",
            "scope": "test",
            "username": c.username,
            "password": c.password,
        }
        self.log("IDT: solicitando token OAuth.")
        response = self.session.post(
            TOKEN_URL,
            data=payload,
            auth=(c.client_id, c.client_secret),
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            timeout=self.timeout,
        )
        response.raise_for_status()
        body = response.json()
        self.token = body.get("access_token", "")
        if not self.token:
            raise RuntimeError(f"A IDT não retornou access_token: {body}")
        self.session.headers.update(
            {"Authorization": f"Bearer {self.token}", "Accept": "application/json"}
        )
        self.log("IDT: autenticação concluída.")
        return self.token

    def _ensure_token(self) -> None:
        if not self.token:
            self.authenticate()

    def _post_json(self, endpoint: str, payload: dict[str, Any]) -> Any:
        self._ensure_token()
        response = self.session.post(
            f"{OLIGO_BASE}/{endpoint}", json=payload, timeout=self.timeout
        )
        if response.status_code == 401:
            self.token = ""
            self._ensure_token()
            response = self.session.post(
                f"{OLIGO_BASE}/{endpoint}", json=payload, timeout=self.timeout
            )
        response.raise_for_status()
        return response.json()

    def _post_params(self, endpoint: str, params: dict[str, str]) -> Any:
        self._ensure_token()
        response = self.session.post(
            f"{OLIGO_BASE}/{endpoint}", params=params, timeout=self.timeout
        )
        if response.status_code == 401:
            self.token = ""
            self._ensure_token()
            response = self.session.post(
                f"{OLIGO_BASE}/{endpoint}", params=params, timeout=self.timeout
            )
        response.raise_for_status()
        return response.json()

    def analyze(self, sequence: str, conditions: IdtConditions) -> dict[str, Any]:
        payload = {
            "Sequence": sequence,
            "NaConc": conditions.na_mm,
            "MgConc": conditions.mg_mm,
            "dNTPsConc": conditions.dntp_mm,
            "OligoConc": conditions.oligo_um,
            "NucleotideType": conditions.nucleotide_type,
        }
        return self._post_json("Analyze", payload)

    def hairpin(self, sequence: str, conditions: IdtConditions) -> list[dict[str, Any]]:
        payload = {
            "Sequence": sequence,
            "NaConc": conditions.na_mm,
            "MgConc": conditions.mg_mm,
            "FoldingTemp": conditions.folding_temp_c,
            "NucleotideType": conditions.nucleotide_type,
        }
        result = self._post_json("Hairpin", payload)
        return result if isinstance(result, list) else [result]

    def self_dimer(self, sequence: str) -> list[dict[str, Any]]:
        result = self._post_params("SelfDimer", {"primary": sequence})
        return result if isinstance(result, list) else [result]

    def hetero_dimer(self, primary: str, secondary: str) -> list[dict[str, Any]]:
        result = self._post_params(
            "HeteroDimer", {"primary": primary, "secondary": secondary}
        )
        return result if isinstance(result, list) else [result]

    @staticmethod
    def strongest_structure(rows: list[dict[str, Any]], key: str = "DeltaG") -> dict[str, Any]:
        if not rows:
            return {}
        def extract(item: dict[str, Any]) -> float:
            for name in (key, key.lower(), "deltaG", "thermo"):
                value = item.get(name)
                if isinstance(value, (int, float)):
                    return float(value)
            return float("inf")
        return min(rows, key=extract)

    def analyze_pair(
        self, forward: str, reverse: str, conditions: IdtConditions
    ) -> dict[str, Any]:
        self.log(f"IDT: analisando par {forward[:8]}… / {reverse[:8]}…")
        f_analysis = self.analyze(forward, conditions)
        r_analysis = self.analyze(reverse, conditions)
        f_hairpins = self.hairpin(forward, conditions)
        r_hairpins = self.hairpin(reverse, conditions)
        f_self = self.self_dimer(forward)
        r_self = self.self_dimer(reverse)
        hetero = self.hetero_dimer(forward, reverse)
        return {
            "forward": {
                "analysis": f_analysis,
                "hairpins": f_hairpins,
                "strongest_hairpin": self.strongest_structure(f_hairpins, "deltaG"),
                "self_dimers": f_self,
                "strongest_self_dimer": self.strongest_structure(f_self),
            },
            "reverse": {
                "analysis": r_analysis,
                "hairpins": r_hairpins,
                "strongest_hairpin": self.strongest_structure(r_hairpins, "deltaG"),
                "self_dimers": r_self,
                "strongest_self_dimer": self.strongest_structure(r_self),
            },
            "hetero_dimers": hetero,
            "strongest_hetero_dimer": self.strongest_structure(hetero),
        }
