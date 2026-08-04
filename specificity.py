from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass, field, is_dataclass
from datetime import datetime, timezone
from itertools import chain
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence


BLAST_OUTFMT_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "sstrand",
)
BLAST_OUTFMT = "6 " + " ".join(BLAST_OUTFMT_FIELDS)
PAIR_ID_RE = re.compile(r"(?i)(?<![a-z0-9])pair[_-]?(\d+)(?!\d)")
QUERY_ID_RE = re.compile(r"(?i)^(pair[_-]?\d+)[_-]([FR])$")
IUPAC_DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.IGNORECASE)


class SpecificityError(RuntimeError):
    """Erro esperado durante a análise local de especificidade."""


def _json_value(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        return {key: _json_value(item) for key, item in asdict(value).items()}
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


@dataclass(slots=True)
class SpecificityExecutables:
    blastn: str = ""
    makeblastdb: str = ""
    mfeprimer: str = ""

    def to_dict(self) -> dict[str, str]:
        return {
            "blastn": str(self.blastn),
            "makeblastdb": str(self.makeblastdb),
            "mfeprimer": str(self.mfeprimer),
        }


@dataclass(slots=True)
class SpecificityParams:
    top_pairs: int = 10
    min_identity_pct: float = 80.0
    min_query_coverage_pct: float = 80.0
    max_target_seqs: int = 100
    cpu: int = 4
    min_amplicon: int = 80
    max_amplicon: int = 250
    tm_cutoff_c: float = 30.0
    mono_mm: float = 50.0
    diva_mm: float = 3.0
    dntp_mm: float = 0.8
    oligo_nm: float = 50.0
    timeout_seconds: int = 1800

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class BlastHit:
    pair_id: str
    pair_rank: int
    primer_orientation: str
    query_id: str
    subject_id: str
    identity_pct: float
    query_coverage_pct: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bit_score: float
    query_length: int
    subject_strand: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class MfeAmplicon:
    pair_id: str
    pair_rank: int
    amplicon_id: str = ""
    sequence_id: str = ""
    start: int | None = None
    end: int | None = None
    length: int | None = None
    forward_tm_c: float | None = None
    reverse_tm_c: float | None = None
    source_format: str = "json"
    raw: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return _json_value(self)


@dataclass(slots=True)
class PairSpecificityResult:
    pair_id: str
    pair_rank: int
    forward_sequence: str
    reverse_sequence: str
    blast_forward_hits: list[BlastHit] = field(default_factory=list)
    blast_reverse_hits: list[BlastHit] = field(default_factory=list)
    mfe_amplicons: list[MfeAmplicon] = field(default_factory=list)
    verdict: str = "Sem produto previsto"

    @property
    def blast_forward_hit_count(self) -> int:
        return len(self.blast_forward_hits)

    @property
    def blast_reverse_hit_count(self) -> int:
        return len(self.blast_reverse_hits)

    @property
    def mfe_amplicon_count(self) -> int:
        return len(self.mfe_amplicons)

    def to_dict(self) -> dict[str, Any]:
        return {
            "pair_id": self.pair_id,
            "pair_rank": self.pair_rank,
            "forward_sequence": self.forward_sequence,
            "reverse_sequence": self.reverse_sequence,
            "blast_forward_hit_count": self.blast_forward_hit_count,
            "blast_reverse_hit_count": self.blast_reverse_hit_count,
            "mfe_amplicon_count": self.mfe_amplicon_count,
            "verdict": self.verdict,
            "blast_forward_hits": [hit.to_dict() for hit in self.blast_forward_hits],
            "blast_reverse_hits": [hit.to_dict() for hit in self.blast_reverse_hits],
            "mfe_amplicons": [amplicon.to_dict() for amplicon in self.mfe_amplicons],
        }


@dataclass(slots=True)
class PreparedReference:
    source_fasta: Path
    cached_fasta: Path
    blast_db_prefix: Path
    mfe_index_path: Path
    fingerprint: str
    manifest_path: Path

    def to_dict(self) -> dict[str, str]:
        return {
            "source_fasta": str(self.source_fasta),
            "cached_fasta": str(self.cached_fasta),
            "blast_db_prefix": str(self.blast_db_prefix),
            "mfe_index_path": str(self.mfe_index_path),
            "fingerprint": self.fingerprint,
            "manifest_path": str(self.manifest_path),
        }


@dataclass(slots=True)
class SpecificityReport:
    reference: PreparedReference
    params: SpecificityParams
    executables: SpecificityExecutables
    tool_versions: dict[str, str]
    results: list[PairSpecificityResult]
    created_at_utc: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds")
    )
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "created_at_utc": self.created_at_utc,
            "reference": self.reference.to_dict(),
            "params": self.params.to_dict(),
            "executables": self.executables.to_dict(),
            "tool_versions": dict(self.tool_versions),
            "warnings": list(self.warnings),
            "results": [result.to_dict() for result in self.results],
        }


Runner = Callable[..., subprocess.CompletedProcess[str]]


class SpecificityService:
    """Coordena BLAST+ e MFEprimer sem depender da interface gráfica."""

    CACHE_MARKER = ".gene_conservado_specificity_cache"
    REFERENCES_DIR = "specificity-references-v1"
    WORK_DIR = "specificity-work-v1"
    DEFAULT_TIMEOUT_SECONDS = 1800

    def __init__(
        self,
        executables: SpecificityExecutables,
        cache_dir: str | os.PathLike[str],
        log: Callable[[str], None] | None = None,
        runner: Runner | None = None,
    ) -> None:
        self.executables = executables
        self.cache_dir = Path(cache_dir).expanduser().resolve(strict=False)
        self.log = log or (lambda _message: None)
        self.runner: Runner = runner or subprocess.run
        self.project_dir = Path(__file__).resolve().parent
        self._resolved_executables: SpecificityExecutables | None = None
        self._validate_cache_root()

    @property
    def resolved_executables(self) -> SpecificityExecutables:
        if self._resolved_executables is None:
            self._resolved_executables = self._resolve_all_executables()
        return self._resolved_executables

    def _validate_cache_root(self) -> None:
        root = Path(self.cache_dir.anchor).resolve()
        forbidden = {root, Path.home().resolve(), self.project_dir.resolve()}
        if self.cache_dir in forbidden:
            raise SpecificityError(
                "Escolha um subdiretório dedicado para o cache de especificidade."
            )

    def _ensure_cache_root(self) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        if self.cache_dir.is_symlink():
            raise SpecificityError("O diretório de cache não pode ser um link simbólico.")
        if os.name != "nt":
            self.cache_dir.chmod(0o700)
        marker = self.cache_dir / self.CACHE_MARKER
        if not marker.exists():
            marker.write_text("Gene Conservado - cache de especificidade\n", encoding="utf-8")
            if os.name != "nt":
                marker.chmod(0o600)

    def _managed_dir(self, name: str, *, create: bool = True) -> Path:
        self._ensure_cache_root()
        path = self.cache_dir / name
        if path.is_symlink():
            raise SpecificityError(f"O caminho de cache gerenciado {path} é um link simbólico.")
        if create:
            path.mkdir(parents=True, exist_ok=True)
            if os.name != "nt":
                path.chmod(0o700)
        resolved = path.resolve(strict=False)
        if not resolved.is_relative_to(self.cache_dir):
            raise SpecificityError("O caminho de cache saiu do diretório permitido.")
        return path

    @staticmethod
    def _executable_names(name: str) -> tuple[str, ...]:
        if os.name == "nt":
            return (f"{name}.exe", name)
        return (name, f"{name}.exe")

    def _resolve_executable(self, name: str, explicit: str) -> str:
        if explicit and str(explicit).strip():
            configured = Path(str(explicit).strip()).expanduser()
            if configured.parent != Path(".") or configured.is_absolute():
                candidate = configured.resolve(strict=False)
                if candidate.is_file() and os.access(candidate, os.X_OK):
                    return str(candidate)
                raise SpecificityError(
                    f"O executável configurado para {name} não existe ou não pode ser executado: "
                    f"{candidate}"
                )
            located = shutil.which(str(explicit).strip())
            if located:
                return str(Path(located).resolve())
            raise SpecificityError(f"O executável configurado para {name} não foi encontrado.")

        local_bin = self.project_dir / ".tools" / "bin"
        for executable_name in self._executable_names(name):
            candidate = local_bin / executable_name
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate.resolve())

        for executable_name in self._executable_names(name):
            located = shutil.which(executable_name)
            if located:
                return str(Path(located).resolve())
        raise SpecificityError(
            f"{name} não foi encontrado. Configure o caminho do executável ou instale-o no PATH."
        )

    def _resolve_all_executables(self) -> SpecificityExecutables:
        return SpecificityExecutables(
            blastn=self._resolve_executable("blastn", self.executables.blastn),
            makeblastdb=self._resolve_executable("makeblastdb", self.executables.makeblastdb),
            mfeprimer=self._resolve_executable("mfeprimer", self.executables.mfeprimer),
        )

    def _run(
        self,
        command: Sequence[str | os.PathLike[str]],
        *,
        timeout: int,
        cwd: Path | None = None,
    ) -> subprocess.CompletedProcess[str]:
        argv = [str(argument) for argument in command]
        self.log(f"Especificidade: executando {Path(argv[0]).name}.")
        try:
            completed = self.runner(
                argv,
                cwd=str(cwd) if cwd is not None else None,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
                shell=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise SpecificityError(
                f"{Path(argv[0]).name} excedeu o tempo máximo de {timeout} segundos."
            ) from exc
        except FileNotFoundError as exc:
            raise SpecificityError(f"O executável {argv[0]} não foi encontrado.") from exc
        except OSError as exc:
            raise SpecificityError(f"Não foi possível executar {argv[0]}: {exc}") from exc

        return_code = int(getattr(completed, "returncode", 0))
        if return_code != 0:
            stderr = str(getattr(completed, "stderr", "") or "").strip()
            stdout = str(getattr(completed, "stdout", "") or "").strip()
            detail = stderr or stdout or f"código de saída {return_code}"
            if len(detail) > 2000:
                detail = detail[-2000:]
            raise SpecificityError(f"{Path(argv[0]).name} falhou: {detail}")
        return completed

    @staticmethod
    def _version_line(completed: subprocess.CompletedProcess[str]) -> str:
        output = "\n".join(
            part for part in (completed.stdout or "", completed.stderr or "") if part
        )
        return next((line.strip() for line in output.splitlines() if line.strip()), "detectado")

    def probe_tools(self) -> dict[str, str]:
        tools = self.resolved_executables
        commands = {
            "blastn": [tools.blastn, "-version"],
            "makeblastdb": [tools.makeblastdb, "-version"],
            "mfeprimer": [tools.mfeprimer, "version"],
        }
        versions: dict[str, str] = {}
        for name, command in commands.items():
            completed = self._run(command, timeout=30)
            versions[name] = self._version_line(completed)
        return versions

    @staticmethod
    def _validate_reference_fasta(path: Path) -> None:
        if not path.is_file():
            raise SpecificityError(f"O FASTA de referência não foi encontrado: {path}")
        if path.stat().st_size < 1:
            raise SpecificityError("O FASTA de referência está vazio.")
        found_header = False
        found_sequence = False
        try:
            with path.open(encoding="utf-8") as handle:
                for line in handle:
                    stripped = line.strip()
                    if not stripped:
                        continue
                    if stripped.startswith(">"):
                        found_header = True
                    elif found_header:
                        found_sequence = True
                        break
                    else:
                        break
        except UnicodeDecodeError as exc:
            raise SpecificityError("O banco deve ser um FASTA textual não compactado.") from exc
        if not (found_header and found_sequence):
            raise SpecificityError("O arquivo de referência não parece ser um FASTA válido.")

    @staticmethod
    def _sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    @staticmethod
    def _blast_index_files(prefix: Path) -> list[Path]:
        return [path for path in prefix.parent.glob(f"{prefix.name}.*") if path.is_file()]

    @staticmethod
    def _find_mfe_index(cached_fasta: Path) -> Path | None:
        candidates = [
            Path(f"{cached_fasta}.primerqc.bin"),
            cached_fasta.with_suffix(".primerqc.bin"),
        ]
        candidates.extend(cached_fasta.parent.glob(f"{cached_fasta.name}*.primerqc.bin"))
        for candidate in candidates:
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
        return None

    def _prepared_from_directory(
        self, source_fasta: Path, fingerprint: str, directory: Path
    ) -> PreparedReference | None:
        manifest = directory / "manifest.json"
        cached_fasta = directory / "reference.fasta"
        blast_prefix = directory / "blast_db"
        if not manifest.is_file() or not cached_fasta.is_file():
            return None
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            return None
        if not isinstance(data, dict) or data.get("fingerprint") != fingerprint:
            return None
        if not self._blast_index_files(blast_prefix):
            return None
        mfe_index = self._find_mfe_index(cached_fasta)
        if mfe_index is None:
            return None
        return PreparedReference(
            source_fasta=source_fasta,
            cached_fasta=cached_fasta,
            blast_db_prefix=blast_prefix,
            mfe_index_path=mfe_index,
            fingerprint=fingerprint,
            manifest_path=manifest,
        )

    def _safe_remove_managed(self, path: Path) -> None:
        if not path.exists():
            return
        if path.is_symlink():
            raise SpecificityError(f"Recusando remover link simbólico do cache: {path}")
        resolved = path.resolve()
        if not resolved.is_relative_to(self.cache_dir) or resolved == self.cache_dir:
            raise SpecificityError("Recusando remover um caminho fora do cache gerenciado.")
        shutil.rmtree(path)

    def prepare_reference(self, reference_fasta: str | os.PathLike[str]) -> PreparedReference:
        source = Path(reference_fasta).expanduser().resolve()
        self._validate_reference_fasta(source)
        fingerprint = self._sha256(source)
        references_dir = self._managed_dir(self.REFERENCES_DIR)
        target_dir = references_dir / fingerprint

        prepared = self._prepared_from_directory(source, fingerprint, target_dir)
        if prepared is not None:
            self.log("Especificidade: índices de referência encontrados no cache.")
            return prepared
        if target_dir.exists():
            self._safe_remove_managed(target_dir)

        tools = self.resolved_executables
        staging = Path(
            tempfile.mkdtemp(prefix=f".{fingerprint}.", dir=str(references_dir))
        )
        try:
            cached_fasta = staging / "reference.fasta"
            shutil.copyfile(source, cached_fasta)
            if self._sha256(cached_fasta) != fingerprint:
                raise SpecificityError("O FASTA de referência mudou durante a cópia.")

            blast_prefix = staging / "blast_db"
            self._run(
                [
                    tools.makeblastdb,
                    "-in",
                    cached_fasta,
                    "-dbtype",
                    "nucl",
                    "-out",
                    blast_prefix,
                ],
                timeout=self.DEFAULT_TIMEOUT_SECONDS,
                cwd=staging,
            )
            if not self._blast_index_files(blast_prefix):
                raise SpecificityError("makeblastdb terminou sem criar os arquivos do banco.")

            self._run(
                [tools.mfeprimer, "index", "-i", cached_fasta, "-k", "9", "-c", "4"],
                timeout=self.DEFAULT_TIMEOUT_SECONDS,
                cwd=staging,
            )
            mfe_index = self._find_mfe_index(cached_fasta)
            if mfe_index is None:
                raise SpecificityError("MFEprimer terminou sem criar o índice .primerqc.bin.")

            manifest = staging / "manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "fingerprint": fingerprint,
                        "source_fasta": str(source),
                        "source_size": source.stat().st_size,
                        "cached_fasta": cached_fasta.name,
                        "blast_db_prefix": blast_prefix.name,
                        "mfe_index": mfe_index.name,
                        "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    },
                    ensure_ascii=False,
                    indent=2,
                ),
                encoding="utf-8",
            )
            os.replace(staging, target_dir)
        except Exception:
            if staging.exists():
                self._safe_remove_managed(staging)
            raise

        prepared = self._prepared_from_directory(source, fingerprint, target_dir)
        if prepared is None:
            raise SpecificityError("O cache da referência ficou incompleto após a indexação.")
        self.log("Especificidade: banco BLAST+ e índice MFEprimer preparados.")
        return prepared

    @staticmethod
    def _validate_params(params: SpecificityParams) -> None:
        if params.top_pairs < 1:
            raise SpecificityError("A quantidade de pares deve ser maior que zero.")
        if not 0.0 <= params.min_identity_pct <= 100.0:
            raise SpecificityError("A identidade BLAST deve estar entre 0% e 100%.")
        if not 0.0 <= params.min_query_coverage_pct <= 100.0:
            raise SpecificityError("A cobertura BLAST deve estar entre 0% e 100%.")
        if params.max_target_seqs < 1 or params.cpu < 1:
            raise SpecificityError("Máximo de alvos e CPUs devem ser maiores que zero.")
        if params.min_amplicon < 1 or params.max_amplicon < params.min_amplicon:
            raise SpecificityError("Os limites de tamanho do amplicon são inválidos.")
        if params.timeout_seconds < 1:
            raise SpecificityError("O tempo máximo deve ser maior que zero.")
        for name, value in (
            ("Tm mínima", params.tm_cutoff_c),
            ("concentração monovalente", params.mono_mm),
            ("concentração divalente", params.diva_mm),
            ("concentração de dNTP", params.dntp_mm),
            ("concentração de oligo", params.oligo_nm),
        ):
            if value < 0:
                raise SpecificityError(f"A {name} não pode ser negativa.")

    @staticmethod
    def _primer_sequence(pair: object, side: str) -> str:
        candidate = getattr(pair, side, None)
        sequence = "".join(
            str(getattr(candidate, "sequence", "") or "").upper().split()
        )
        if not sequence or not IUPAC_DNA_RE.fullmatch(sequence):
            label = "Forward" if side == "forward" else "Reverse"
            raise SpecificityError(f"O primer {label} contém uma sequência inválida.")
        return sequence

    @staticmethod
    def _normalize_pair_id(value: str) -> str | None:
        match = PAIR_ID_RE.search(value)
        if not match:
            return None
        return f"pair_{int(match.group(1)):04d}"

    def _pair_inputs(self, pairs: Iterable[object], top_pairs: int) -> list[dict[str, Any]]:
        selected = list(pairs)[:top_pairs]
        if not selected:
            raise SpecificityError("Nenhum par de primers foi informado para a análise.")
        inputs: list[dict[str, Any]] = []
        used_ids: set[str] = set()
        for position, pair in enumerate(selected, start=1):
            try:
                rank = int(getattr(pair, "rank", position))
            except (TypeError, ValueError):
                rank = position
            if rank < 1:
                rank = position
            pair_id = f"pair_{rank:04d}"
            if pair_id in used_ids:
                raise SpecificityError("Os pares de primers possuem ranks duplicados.")
            used_ids.add(pair_id)
            inputs.append(
                {
                    "pair": pair,
                    "pair_id": pair_id,
                    "rank": rank,
                    "forward": self._primer_sequence(pair, "forward"),
                    "reverse": self._primer_sequence(pair, "reverse"),
                }
            )
        return inputs

    @staticmethod
    def _parse_blast_output(
        output: str,
        known_pairs: Mapping[str, int],
        params: SpecificityParams,
    ) -> list[BlastHit]:
        hits: list[BlastHit] = []
        for line_number, line in enumerate(output.splitlines(), start=1):
            if not line.strip():
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != len(BLAST_OUTFMT_FIELDS):
                raise SpecificityError(
                    f"Saída BLAST inválida na linha {line_number}: esperadas "
                    f"{len(BLAST_OUTFMT_FIELDS)} colunas."
                )
            (
                query_id,
                subject_id,
                identity,
                alignment_length,
                mismatches,
                gap_opens,
                query_start,
                query_end,
                subject_start,
                subject_end,
                evalue,
                bit_score,
                query_length,
                subject_strand,
            ) = columns
            match = QUERY_ID_RE.fullmatch(query_id)
            if match is None:
                raise SpecificityError(f"O BLAST retornou um identificador desconhecido: {query_id}")
            pair_id = SpecificityService._normalize_pair_id(match.group(1))
            if pair_id is None or pair_id not in known_pairs:
                raise SpecificityError(f"O BLAST retornou um par desconhecido: {query_id}")
            try:
                identity_pct = float(identity)
                aln_length = int(alignment_length)
                q_start = int(query_start)
                q_end = int(query_end)
                q_length = int(query_length)
                if q_length < 1:
                    raise ValueError
                coverage = 100.0 * (abs(q_end - q_start) + 1) / q_length
                hit = BlastHit(
                    pair_id=pair_id,
                    pair_rank=known_pairs[pair_id],
                    primer_orientation=match.group(2).upper(),
                    query_id=query_id,
                    subject_id=subject_id,
                    identity_pct=identity_pct,
                    query_coverage_pct=coverage,
                    alignment_length=aln_length,
                    mismatches=int(mismatches),
                    gap_opens=int(gap_opens),
                    query_start=q_start,
                    query_end=q_end,
                    subject_start=int(subject_start),
                    subject_end=int(subject_end),
                    evalue=float(evalue),
                    bit_score=float(bit_score),
                    query_length=q_length,
                    subject_strand=subject_strand,
                )
            except ValueError as exc:
                raise SpecificityError(
                    f"A saída numérica do BLAST é inválida na linha {line_number}."
                ) from exc
            if (
                hit.identity_pct >= params.min_identity_pct
                and hit.query_coverage_pct >= params.min_query_coverage_pct
            ):
                hits.append(hit)
        return hits

    @staticmethod
    def _normalize_key(value: str) -> str:
        return re.sub(r"[^a-z0-9]", "", value.casefold())

    @classmethod
    def _mapping_value(cls, data: Mapping[str, Any], names: set[str]) -> Any:
        normalized_names = {cls._normalize_key(name) for name in names}
        for key, value in data.items():
            if cls._normalize_key(str(key)) in normalized_names and value not in (None, ""):
                return value
        for value in data.values():
            if isinstance(value, Mapping):
                nested = cls._mapping_value(value, names)
                if nested not in (None, ""):
                    return nested
        return None

    @classmethod
    def _direct_mapping_value(cls, data: Mapping[str, Any], names: set[str]) -> Any:
        normalized_names = {cls._normalize_key(name) for name in names}
        for key, value in data.items():
            if cls._normalize_key(str(key)) in normalized_names and value not in (None, ""):
                return value
        return None

    @classmethod
    def _child_mapping(cls, data: Mapping[str, Any], name: str) -> Mapping[str, Any] | None:
        normalized_name = cls._normalize_key(name)
        for key, value in data.items():
            if cls._normalize_key(str(key)) == normalized_name and isinstance(value, Mapping):
                return value
        return None

    @classmethod
    def _pair_id_from_record(cls, record: Any) -> str | None:
        if isinstance(record, Mapping):
            direct = cls._mapping_value(
                record, {"pair_id", "pairid", "pair", "primer_pair", "primerpair", "name"}
            )
            if direct is not None:
                normalized = cls._normalize_pair_id(str(direct))
                if normalized is not None:
                    return normalized
            for key, value in record.items():
                normalized = cls._normalize_pair_id(str(key))
                if normalized is not None:
                    return normalized
                normalized = cls._pair_id_from_record(value)
                if normalized is not None:
                    return normalized
        elif isinstance(record, (list, tuple)):
            for value in record:
                normalized = cls._pair_id_from_record(value)
                if normalized is not None:
                    return normalized
        elif isinstance(record, str):
            return cls._normalize_pair_id(record)
        return None

    @staticmethod
    def _optional_int(value: Any) -> int | None:
        if value in (None, ""):
            return None
        try:
            return int(float(str(value)))
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _optional_float(value: Any) -> float | None:
        if value in (None, ""):
            return None
        try:
            return float(str(value))
        except (TypeError, ValueError):
            return None

    @classmethod
    def _amplicon_from_mapping(
        cls,
        record: Mapping[str, Any],
        known_pairs: Mapping[str, int],
        *,
        source_format: str,
    ) -> MfeAmplicon | None:
        pair_id = cls._pair_id_from_record(record)
        if pair_id is None and len(known_pairs) == 1:
            pair_id = next(iter(known_pairs))
        if pair_id is None or pair_id not in known_pairs:
            return None

        start = cls._optional_int(
            cls._direct_mapping_value(
                record, {"start", "amp_start", "ampstart", "product_start"}
            )
        )
        end = cls._optional_int(
            cls._direct_mapping_value(
                record, {"end", "amp_end", "ampend", "product_end"}
            )
        )
        length = cls._optional_int(
            cls._direct_mapping_value(
                record,
                {
                    "size",
                    "length",
                    "amp_size",
                    "ampsize",
                    "amplicon_length",
                    "ampliconlength",
                    "product_size",
                },
            )
        )
        forward_data = cls._child_mapping(record, "F")
        reverse_data = cls._child_mapping(record, "R")
        if start is None and forward_data is not None:
            start = cls._optional_int(
                cls._direct_mapping_value(forward_data, {"binding_start", "bindingstart"})
            )
            if start is None:
                forward_positions = [
                    position
                    for position in (
                        cls._optional_int(cls._direct_mapping_value(forward_data, {"start"})),
                        cls._optional_int(cls._direct_mapping_value(forward_data, {"end"})),
                    )
                    if position is not None
                ]
                if forward_positions:
                    # O JSON do MFEprimer 4.x usa posição inicial zero-based.
                    start = min(forward_positions) + 1
        if end is None and reverse_data is not None:
            end = cls._optional_int(
                cls._direct_mapping_value(reverse_data, {"binding_start", "bindingstart"})
            )
            if end is None:
                reverse_positions = [
                    position
                    for position in (
                        cls._optional_int(cls._direct_mapping_value(reverse_data, {"start"})),
                        cls._optional_int(cls._direct_mapping_value(reverse_data, {"end"})),
                    )
                    if position is not None
                ]
                if reverse_positions:
                    end = max(reverse_positions)
        if length is None and start is not None and end is not None:
            length = abs(end - start) + 1
        forward_tm = cls._optional_float(
            cls._direct_mapping_value(
                record, {"forward_tm", "forwardtm", "f_tm", "ftm", "fp_tm", "fptm"}
            )
        )
        reverse_tm = cls._optional_float(
            cls._direct_mapping_value(
                record, {"reverse_tm", "reversetm", "r_tm", "rtm", "rp_tm", "rptm"}
            )
        )
        if forward_tm is None and forward_data is not None:
            forward_tm = cls._optional_float(cls._direct_mapping_value(forward_data, {"tm"}))
        if reverse_tm is None and reverse_data is not None:
            reverse_tm = cls._optional_float(cls._direct_mapping_value(reverse_data, {"tm"}))
        raw = {str(key): _json_value(value) for key, value in record.items()}
        return MfeAmplicon(
            pair_id=pair_id,
            pair_rank=known_pairs[pair_id],
            amplicon_id=str(
                cls._direct_mapping_value(
                    record,
                    {"amplicon_id", "ampliconid", "amp_id", "ampid", "id", "name"},
                )
                or ""
            ),
            sequence_id=str(
                cls._direct_mapping_value(
                    record,
                    {
                        "sequence_id",
                        "sequenceid",
                        "seq_id",
                        "seqid",
                        "subject_id",
                        "subjectid",
                        "hit_id",
                        "hitid",
                        "hid",
                        "target",
                        "chrom",
                        "chromosome",
                        "template",
                    },
                )
                or ""
            ),
            start=start,
            end=end,
            length=length,
            forward_tm_c=forward_tm,
            reverse_tm_c=reverse_tm,
            source_format=source_format,
            raw=raw,
        )

    @classmethod
    def _amplicon_records_from_json(cls, data: Any) -> list[Mapping[str, Any]]:
        records: list[Mapping[str, Any]] = []

        def flatten_amp_list(value: Any) -> None:
            if isinstance(value, Mapping):
                if cls._pair_id_from_record(value) is not None or any(
                    cls._normalize_key(str(key))
                    in {"size", "length", "ampliconlength", "ampstart", "ampend"}
                    for key in value
                ):
                    records.append(value)
                else:
                    for nested in value.values():
                        flatten_amp_list(nested)
            elif isinstance(value, list):
                for nested in value:
                    flatten_amp_list(nested)

        if isinstance(data, list):
            flatten_amp_list(data)
            return records
        if not isinstance(data, Mapping):
            return records
        for key, value in data.items():
            if cls._normalize_key(str(key)) in {"amplist", "amplicons", "products"}:
                flatten_amp_list(value)
        return records

    @classmethod
    def _parse_mfe_json(
        cls, path: Path, known_pairs: Mapping[str, int]
    ) -> list[MfeAmplicon]:
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError) as exc:
            raise SpecificityError(f"O JSON do MFEprimer é inválido: {path.name}") from exc
        result: list[MfeAmplicon] = []
        for record in cls._amplicon_records_from_json(data):
            amplicon = cls._amplicon_from_mapping(record, known_pairs, source_format="json")
            if amplicon is not None:
                result.append(amplicon)
        return result

    @classmethod
    def _parse_mfe_tsv(
        cls, path: Path, known_pairs: Mapping[str, int]
    ) -> list[MfeAmplicon]:
        try:
            handle = path.open(encoding="utf-8", newline="")
        except (OSError, UnicodeDecodeError) as exc:
            raise SpecificityError(f"Não foi possível ler a saída MFEprimer: {path.name}") from exc
        try:
            with handle:
                header: str | None = None
                for line in handle:
                    if not line.strip():
                        continue
                    if line.lstrip().startswith("#"):
                        possible_header = line.lstrip()[1:].strip()
                        if (
                            "\t" in possible_header
                            and PAIR_ID_RE.search(possible_header) is None
                        ):
                            header = possible_header
                            break
                        continue
                    header = line.strip()
                    break
                if header is None:
                    return []
                reader = csv.DictReader(chain([header], handle), delimiter="\t")
                if not reader.fieldnames:
                    return []
                result: list[MfeAmplicon] = []
                for row in reader:
                    normalized_row = {
                        str(key): value for key, value in row.items() if key is not None
                    }
                    amplicon = cls._amplicon_from_mapping(
                        normalized_row, known_pairs, source_format="spec.tsv"
                    )
                    if amplicon is not None:
                        result.append(amplicon)
                return result
        except UnicodeDecodeError as exc:
            raise SpecificityError(
                f"Não foi possível ler a saída MFEprimer: {path.name}"
            ) from exc

    @staticmethod
    def _mfe_output_candidates(base: Path) -> tuple[list[Path], list[Path]]:
        json_candidates = [Path(f"{base}.json"), base.with_suffix(".json")]
        tsv_candidates = [
            Path(f"{base}.spec.tsv"),
            base.with_suffix(".spec.tsv"),
            Path(f"{base}.tsv"),
        ]
        return json_candidates, tsv_candidates

    @staticmethod
    def _deduplicate_amplicons(amplicons: Iterable[MfeAmplicon]) -> list[MfeAmplicon]:
        unique: list[MfeAmplicon] = []
        seen: set[tuple[Any, ...]] = set()
        for amplicon in amplicons:
            key = (
                amplicon.pair_id,
                amplicon.amplicon_id,
                amplicon.sequence_id,
                amplicon.start,
                amplicon.end,
                amplicon.length,
            )
            if key in seen:
                continue
            seen.add(key)
            unique.append(amplicon)
        return unique

    def _read_mfe_outputs(
        self, output_base: Path, known_pairs: Mapping[str, int]
    ) -> tuple[list[MfeAmplicon], list[str]]:
        json_candidates, tsv_candidates = self._mfe_output_candidates(output_base)
        json_path = next(
            (path for path in json_candidates if path.is_file() and path.stat().st_size > 0),
            None,
        )
        tsv_path = next(
            (path for path in tsv_candidates if path.is_file() and path.stat().st_size > 0),
            None,
        )
        warnings: list[str] = []
        if tsv_path is not None:
            # O MFEprimer 4.x pode limitar o AmpList do JSON a --max-amp-count,
            # mas mantém todos os produtos no .spec.tsv. O TSV é, portanto, a
            # fonte canônica sempre que estiver disponível.
            return self._deduplicate_amplicons(
                self._parse_mfe_tsv(tsv_path, known_pairs)
            ), warnings
        if json_path is not None:
            return self._deduplicate_amplicons(
                self._parse_mfe_json(json_path, known_pairs)
            ), warnings
        raise SpecificityError("MFEprimer terminou sem criar uma saída JSON ou .spec.tsv.")

    @staticmethod
    def _verdict(
        forward_hits: int,
        reverse_hits: int,
        amplicons: int,
        max_target_seqs: int,
    ) -> str:
        if amplicons == 0:
            return "Sem produto previsto"
        if amplicons > 1:
            return "Múltiplos produtos previstos"
        blast_limit_reached = (
            forward_hits >= max_target_seqs or reverse_hits >= max_target_seqs
        )
        if forward_hits <= 1 and reverse_hits <= 1 and not blast_limit_reached:
            return "Produto único no banco"
        return "Produto único; revisar sítios BLAST"

    def analyze(
        self,
        reference_fasta: str | os.PathLike[str],
        pairs: Iterable[object],
        params: SpecificityParams,
    ) -> SpecificityReport:
        self._validate_params(params)
        pair_inputs = self._pair_inputs(pairs, params.top_pairs)
        versions = self.probe_tools()
        prepared = self.prepare_reference(reference_fasta)
        tools = self.resolved_executables
        known_pairs = {item["pair_id"]: item["rank"] for item in pair_inputs}

        work_root = self._managed_dir(self.WORK_DIR)
        warnings: list[str] = []
        with tempfile.TemporaryDirectory(prefix="analysis-", dir=str(work_root)) as temporary:
            work_dir = Path(temporary)
            blast_query = work_dir / "primers.fasta"
            mfe_input = work_dir / "primer_pairs.tsv"
            output_base = work_dir / "mfeprimer_result"

            fasta_chunks: list[str] = []
            tsv_lines: list[str] = []
            for item in pair_inputs:
                fasta_chunks.extend(
                    [
                        f">{item['pair_id']}_F\n{item['forward']}\n",
                        f">{item['pair_id']}_R\n{item['reverse']}\n",
                    ]
                )
                tsv_lines.append(
                    f"{item['pair_id']}\t{item['forward']}\t{item['reverse']}\n"
                )
            blast_query.write_text("".join(fasta_chunks), encoding="utf-8")
            mfe_input.write_text("".join(tsv_lines), encoding="utf-8")

            blast_completed = self._run(
                [
                    tools.blastn,
                    "-task",
                    "blastn-short",
                    "-query",
                    blast_query,
                    "-db",
                    prepared.blast_db_prefix,
                    "-strand",
                    "both",
                    "-dust",
                    "no",
                    "-evalue",
                    "1000",
                    "-max_target_seqs",
                    str(params.max_target_seqs),
                    "-num_threads",
                    str(params.cpu),
                    "-outfmt",
                    BLAST_OUTFMT,
                ],
                timeout=params.timeout_seconds,
                cwd=work_dir,
            )
            blast_hits = self._parse_blast_output(
                blast_completed.stdout or "", known_pairs, params
            )

            self._run(
                [
                    tools.mfeprimer,
                    "spec",
                    "-i",
                    mfe_input,
                    "-d",
                    prepared.cached_fasta,
                    "-o",
                    output_base,
                    "-j",
                    "-s",
                    str(params.min_amplicon),
                    "-S",
                    str(params.max_amplicon),
                    "-t",
                    str(params.tm_cutoff_c),
                    "-c",
                    str(params.cpu),
                    "--mono",
                    str(params.mono_mm),
                    "--diva",
                    str(params.diva_mm),
                    "--dntp",
                    str(params.dntp_mm),
                    "--oligo",
                    str(params.oligo_nm),
                ],
                timeout=params.timeout_seconds,
                cwd=work_dir,
            )
            amplicons, mfe_warnings = self._read_mfe_outputs(output_base, known_pairs)
            warnings.extend(mfe_warnings)

        results: list[PairSpecificityResult] = []
        for item in pair_inputs:
            pair_id = item["pair_id"]
            forward_hits = [
                hit
                for hit in blast_hits
                if hit.pair_id == pair_id and hit.primer_orientation == "F"
            ]
            reverse_hits = [
                hit
                for hit in blast_hits
                if hit.pair_id == pair_id and hit.primer_orientation == "R"
            ]
            pair_amplicons = [amplicon for amplicon in amplicons if amplicon.pair_id == pair_id]
            results.append(
                PairSpecificityResult(
                    pair_id=pair_id,
                    pair_rank=item["rank"],
                    forward_sequence=item["forward"],
                    reverse_sequence=item["reverse"],
                    blast_forward_hits=forward_hits,
                    blast_reverse_hits=reverse_hits,
                    mfe_amplicons=pair_amplicons,
                    verdict=self._verdict(
                        len(forward_hits),
                        len(reverse_hits),
                        len(pair_amplicons),
                        params.max_target_seqs,
                    ),
                )
            )

        self.log(f"Especificidade: {len(results)} pares analisados no banco selecionado.")
        return SpecificityReport(
            reference=prepared,
            params=params,
            executables=tools,
            tool_versions=versions,
            results=results,
            warnings=warnings,
        )

    def clear_cache(self) -> None:
        if not self.cache_dir.exists():
            return
        marker = self.cache_dir / self.CACHE_MARKER
        if not marker.is_file():
            raise SpecificityError(
                "O diretório não possui a marca de cache do Gene Conservado; nada foi removido."
            )
        for name in (self.REFERENCES_DIR, self.WORK_DIR):
            path = self.cache_dir / name
            if path.exists():
                self._safe_remove_managed(path)
        self.log("Especificidade: índices e arquivos temporários locais removidos.")


__all__ = [
    "BlastHit",
    "MfeAmplicon",
    "PairSpecificityResult",
    "PreparedReference",
    "SpecificityError",
    "SpecificityExecutables",
    "SpecificityParams",
    "SpecificityReport",
    "SpecificityService",
]
