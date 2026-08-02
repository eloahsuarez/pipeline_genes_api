from __future__ import annotations

import csv
import json
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Iterable

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

from models import ConservedRegion, PrimerPair, SequenceRecord


def _serializable(value: Any) -> Any:
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, dict):
        return {key: _serializable(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_serializable(item) for item in value]
    return value


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(_serializable(data), ensure_ascii=False, indent=2), encoding="utf-8")


def export_sequences_csv(path: Path, records: Iterable[SequenceRecord]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "selected", "accession", "organism", "length", "definition", "genes",
                "molecule_type", "exons", "sequence",
            ],
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "selected": record.selected,
                    "accession": record.accession,
                    "organism": record.organism,
                    "length": record.length,
                    "definition": record.definition,
                    "genes": "; ".join(record.genes),
                    "molecule_type": record.molecule_type,
                    "exons": "; ".join(
                        f"{exon.number or '?'}:{exon.start}-{exon.end}"
                        for exon in record.exons
                    ),
                    "sequence": record.sequence,
                }
            )


def export_regions_csv(path: Path, regions: Iterable[ConservedRegion]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        fieldnames = list(ConservedRegion.__dataclass_fields__.keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for region in regions:
            writer.writerow(asdict(region))


def _style_sheet(ws) -> None:
    header_fill = PatternFill("solid", fgColor="D9EAF7")
    for cell in ws[1]:
        cell.font = Font(bold=True)
        cell.fill = header_fill
        cell.alignment = Alignment(horizontal="center")
    ws.freeze_panes = "A2"
    for column in ws.columns:
        values = [str(cell.value or "") for cell in column[:50]]
        width = min(max(len(value) for value in values) + 2, 60)
        ws.column_dimensions[get_column_letter(column[0].column)].width = max(width, 10)


def export_primers_xlsx(path: Path, pairs: list[PrimerPair]) -> None:
    wb = Workbook()
    ws = wb.active
    ws.title = "Pares de primers"
    headers = [
        "Rank", "Score", "Amplicon início", "Amplicon fim", "Amplicon bp",
        "Forward", "F início", "F fim", "F tamanho", "F GC %", "F Tm °C",
        "Reverse", "R início", "R fim", "R tamanho", "R GC %", "R Tm °C",
        "IDT F Tm", "IDT R Tm", "Hairpin F ΔG", "Hairpin R ΔG",
        "Self-dimer F ΔG", "Self-dimer R ΔG", "Heterodimer ΔG",
        "Referência da junção", "F cruza junção", "F junções",
        "R cruza junção", "R junções",
    ]
    ws.append(headers)
    for pair in pairs:
        idt = pair.idt or {}
        f = idt.get("forward", {})
        r = idt.get("reverse", {})
        hetero = idt.get("strongest_hetero_dimer", {})

        def junction_summary(candidate) -> str:
            return "; ".join(
                (
                    f"{match.junction_position} "
                    f"(E{match.left_exon_number or '?'}-E{match.right_exon_number or '?'}; "
                    f"5'={match.primer_5_prime_bases}; 3'={match.primer_3_prime_bases})"
                )
                for match in candidate.junctions
            )

        ws.append(
            [
                pair.rank, round(pair.score, 3), pair.amplicon_start, pair.amplicon_end,
                pair.amplicon_length, pair.forward.sequence, pair.forward.start,
                pair.forward.end, pair.forward.length, round(pair.forward.gc_percent, 2),
                round(pair.forward.tm_c, 2), pair.reverse.sequence, pair.reverse.start,
                pair.reverse.end, pair.reverse.length, round(pair.reverse.gc_percent, 2),
                round(pair.reverse.tm_c, 2),
                f.get("analysis", {}).get("MeltTemp"),
                r.get("analysis", {}).get("MeltTemp"),
                f.get("strongest_hairpin", {}).get("deltaG"),
                r.get("strongest_hairpin", {}).get("deltaG"),
                f.get("strongest_self_dimer", {}).get("DeltaG"),
                r.get("strongest_self_dimer", {}).get("DeltaG"),
                hetero.get("DeltaG"),
                pair.reference_accession,
                pair.forward.spans_exon_junction,
                junction_summary(pair.forward),
                pair.reverse.spans_exon_junction,
                junction_summary(pair.reverse),
            ]
        )
    _style_sheet(ws)

    raw = wb.create_sheet("IDT bruto")
    raw.append(["Rank", "JSON"])
    for pair in pairs:
        raw.append([pair.rank, json.dumps(pair.idt, ensure_ascii=False)])
    _style_sheet(raw)
    wb.save(path)
