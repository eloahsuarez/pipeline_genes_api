from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any


@dataclass(slots=True)
class ExonInterval:
    """Intervalo de exon nas coordenadas 1-based inclusivas do transcrito."""

    start: int
    end: int
    number: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class SequenceRecord:
    uid: str
    accession: str
    definition: str
    organism: str
    sequence: str
    length: int
    genes: list[str] = field(default_factory=list)
    feature_keys: list[str] = field(default_factory=list)
    selected: bool = True
    exons: list[ExonInterval] = field(default_factory=list)
    molecule_type: str = ""

    def fasta(self) -> str:
        safe_definition = " ".join(self.definition.split())
        return f">{self.accession} {safe_definition}\n{self.sequence.upper()}\n"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class ConservedRegion:
    alignment_start: int
    alignment_end: int
    consensus_start: int
    consensus_end: int
    length: int
    sequence: str
    mean_identity: float
    mean_coverage: float
    reference_accession: str = ""
    reference_start: int = 0
    reference_end: int = 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class ExonJunctionMatch:
    """Metadados de uma junção atravessada por um primer."""

    reference_accession: str
    junction_position: int
    left_exon_number: str
    right_exon_number: str
    primer_5_prime_bases: int
    primer_3_prime_bases: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class PrimerCandidate:
    orientation: str
    start: int
    end: int
    sequence: str
    length: int
    gc_percent: float
    tm_c: float
    score: float
    spans_exon_junction: bool = False
    junctions: list[ExonJunctionMatch] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class PrimerPair:
    rank: int
    forward: PrimerCandidate
    reverse: PrimerCandidate
    amplicon_start: int
    amplicon_end: int
    amplicon_length: int
    score: float
    idt: dict[str, Any] = field(default_factory=dict)
    reference_accession: str = ""

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        return data
