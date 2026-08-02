from __future__ import annotations

from collections import Counter
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from io import StringIO
from typing import Iterable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from models import (
    ConservedRegion,
    ExonInterval,
    ExonJunctionMatch,
    PrimerCandidate,
    PrimerPair,
)


DNA = set("ACGT")


@dataclass(slots=True)
class PrimerDesignParams:
    identity_threshold: float = 1.0
    coverage_threshold: float = 1.0
    min_primer_len: int = 18
    max_primer_len: int = 25
    min_gc: float = 35.0
    max_gc: float = 65.0
    min_tm: float = 55.0
    max_tm: float = 65.0
    target_tm: float = 60.0
    min_amplicon: int = 80
    max_amplicon: int = 250
    target_amplicon: int = 150
    top_pairs: int = 50
    na_mm: float = 50.0
    mg_mm: float = 3.0
    dntp_mm: float = 0.8
    primer_conc_nm: float = 250.0


def _validate_primer_design_params(params: PrimerDesignParams) -> None:
    if not (0.0 <= params.identity_threshold <= 1.0):
        raise ValueError("A identidade mínima deve estar entre 0% e 100%.")
    if not (0.0 <= params.coverage_threshold <= 1.0):
        raise ValueError("A cobertura mínima deve estar entre 0% e 100%.")
    if params.min_primer_len < 1 or params.max_primer_len < params.min_primer_len:
        raise ValueError("Os comprimentos mínimo e máximo do primer são inválidos.")
    if not (0.0 <= params.min_gc <= params.max_gc <= 100.0):
        raise ValueError("A faixa de GC deve estar entre 0% e 100%.")
    if params.max_tm < params.min_tm:
        raise ValueError("A Tm máxima não pode ser menor que a Tm mínima.")
    if params.min_amplicon < 1 or params.max_amplicon < params.min_amplicon:
        raise ValueError("Os comprimentos mínimo e máximo do amplicon são inválidos.")
    if params.top_pairs < 1:
        raise ValueError("A quantidade de pares deve ser maior que zero.")


def parse_alignment_fasta(text: str) -> list[tuple[str, str]]:
    records = [(record.id, str(record.seq).upper()) for record in SeqIO.parse(StringIO(text), "fasta")]
    if len(records) < 2:
        raise ValueError("O alinhamento FASTA contém menos de duas sequências.")
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise ValueError("As sequências do alinhamento não possuem o mesmo comprimento.")
    return records


def _column_statistics(sequences: list[str], index: int) -> tuple[str, float, float]:
    chars = [sequence[index] for sequence in sequences]
    canonical = [char for char in chars if char in DNA]
    coverage = len(canonical) / len(chars)
    if not canonical:
        return "N", 0.0, coverage
    counts = Counter(canonical)
    base, count = counts.most_common(1)[0]
    identity = count / len(canonical)
    return base, identity, coverage


def build_consensus_and_mask(
    alignment_text: str,
    identity_threshold: float,
    coverage_threshold: float,
) -> tuple[str, list[bool], list[float], list[float], list[int]]:
    records = parse_alignment_fasta(alignment_text)
    sequences = [sequence for _, sequence in records]
    aln_len = len(sequences[0])
    ungapped_consensus: list[str] = []
    mask: list[bool] = []
    identities: list[float] = []
    coverages: list[float] = []
    alignment_positions: list[int] = []

    for index in range(aln_len):
        base, identity, coverage = _column_statistics(sequences, index)
        if base not in DNA:
            continue
        ungapped_consensus.append(base)
        mask.append(identity >= identity_threshold and coverage >= coverage_threshold)
        identities.append(identity)
        coverages.append(coverage)
        alignment_positions.append(index + 1)
    return "".join(ungapped_consensus), mask, identities, coverages, alignment_positions


def find_conserved_regions(
    alignment_text: str,
    identity_threshold: float = 1.0,
    coverage_threshold: float = 1.0,
    min_length: int = 18,
) -> list[ConservedRegion]:
    consensus, mask, identities, coverages, aln_positions = build_consensus_and_mask(
        alignment_text, identity_threshold, coverage_threshold
    )
    regions: list[ConservedRegion] = []
    start: int | None = None

    def close(end_exclusive: int) -> None:
        nonlocal start
        if start is None:
            return
        length = end_exclusive - start
        if length >= min_length:
            regions.append(
                ConservedRegion(
                    alignment_start=aln_positions[start],
                    alignment_end=aln_positions[end_exclusive - 1],
                    consensus_start=start + 1,
                    consensus_end=end_exclusive,
                    length=length,
                    sequence=consensus[start:end_exclusive],
                    mean_identity=sum(identities[start:end_exclusive]) / length,
                    mean_coverage=sum(coverages[start:end_exclusive]) / length,
                )
            )
        start = None

    for index, is_conserved in enumerate(mask):
        if is_conserved and start is None:
            start = index
        elif not is_conserved and start is not None:
            close(index)
    close(len(mask))
    return regions


def gc_percent(sequence: str) -> float:
    if not sequence:
        return 0.0
    return 100.0 * sum(base in "GC" for base in sequence.upper()) / len(sequence)


def calculate_tm(sequence: str, params: PrimerDesignParams) -> float:
    try:
        return float(
            mt.Tm_NN(
                sequence,
                dnac1=params.primer_conc_nm,
                dnac2=params.primer_conc_nm,
                Na=params.na_mm,
                Mg=params.mg_mm,
                dNTPs=params.dntp_mm,
                saltcorr=7,
            )
        )
    except Exception:
        return float(mt.Tm_Wallace(sequence))


def _homopolymer_penalty(sequence: str) -> float:
    penalty = 0.0
    for base in DNA:
        if base * 5 in sequence:
            penalty += 10.0
        elif base * 4 in sequence:
            penalty += 3.0
    return penalty


def _three_prime_penalty(sequence: str) -> float:
    tail = sequence[-5:]
    gc_tail = sum(base in "GC" for base in tail)
    if gc_tail == 0:
        return 4.0
    if gc_tail >= 4:
        return 5.0
    return 0.0


def _candidate(
    sequence: str,
    orientation: str,
    start: int,
    end: int,
    params: PrimerDesignParams,
) -> PrimerCandidate | None:
    oligo = sequence if orientation == "F" else str(Seq(sequence).reverse_complement())
    gc = gc_percent(oligo)
    if not (params.min_gc <= gc <= params.max_gc):
        return None
    tm_value = calculate_tm(oligo, params)
    if not (params.min_tm <= tm_value <= params.max_tm):
        return None
    score = (
        abs(tm_value - params.target_tm) * 2.0
        + abs(gc - 50.0) * 0.12
        + _homopolymer_penalty(oligo)
        + _three_prime_penalty(oligo)
    )
    return PrimerCandidate(
        orientation=orientation,
        start=start,
        end=end,
        sequence=oligo,
        length=len(oligo),
        gc_percent=gc,
        tm_c=tm_value,
        score=score,
    )


def _pair_primer_candidates(
    forward: list[PrimerCandidate],
    reverse: list[PrimerCandidate],
    params: PrimerDesignParams,
    *,
    require_exon_junction: bool = False,
    reference_accession: str = "",
) -> list[PrimerPair]:
    reverse_by_start = sorted(reverse, key=lambda item: item.start)
    reverse_starts = [item.start for item in reverse_by_start]
    junction_reverse = [item for item in reverse_by_start if item.spans_exon_junction]
    junction_reverse_starts = [item.start for item in junction_reverse]
    pairs: list[PrimerPair] = []
    for fwd in forward:
        if require_exon_junction and not fwd.spans_exon_junction:
            eligible_reverse = junction_reverse
            eligible_starts = junction_reverse_starts
        else:
            eligible_reverse = reverse_by_start
            eligible_starts = reverse_starts
        min_reverse_start = fwd.start + params.min_amplicon - params.max_primer_len
        max_reverse_start = fwd.start + params.max_amplicon
        lo = bisect_left(eligible_starts, min_reverse_start)
        hi = bisect_right(eligible_starts, max_reverse_start)
        for rev in eligible_reverse[lo:hi]:
            if rev.start <= fwd.end:
                continue
            if require_exon_junction and not (
                fwd.spans_exon_junction or rev.spans_exon_junction
            ):
                continue
            amplicon_length = rev.end - fwd.start + 1
            if not (params.min_amplicon <= amplicon_length <= params.max_amplicon):
                continue
            tm_difference = abs(fwd.tm_c - rev.tm_c)
            pair_score = (
                fwd.score
                + rev.score
                + tm_difference * 3.0
                + abs(amplicon_length - params.target_amplicon) * 0.03
            )
            pairs.append(
                PrimerPair(
                    rank=0,
                    forward=fwd,
                    reverse=rev,
                    amplicon_start=fwd.start,
                    amplicon_end=rev.end,
                    amplicon_length=amplicon_length,
                    score=pair_score,
                    reference_accession=reference_accession,
                )
            )

    pairs.sort(key=lambda pair: pair.score)
    unique: list[PrimerPair] = []
    seen: set[tuple[object, ...]] = set()
    for pair in pairs:
        if require_exon_junction:
            # Em regiões repetidas, oligos iguais podem representar amplicons e
            # junções diferentes; preserve cada mapeamento no resultado.
            key = (
                pair.forward.sequence,
                pair.reverse.sequence,
                pair.forward.start,
                pair.forward.end,
                pair.reverse.start,
                pair.reverse.end,
            )
        else:
            key = (pair.forward.sequence, pair.reverse.sequence)
        if key in seen:
            continue
        seen.add(key)
        pair.rank = len(unique) + 1
        unique.append(pair)
        if len(unique) >= params.top_pairs:
            break
    return unique


def generate_primer_pairs(alignment_text: str, params: PrimerDesignParams) -> list[PrimerPair]:
    _validate_primer_design_params(params)
    consensus, mask, _, _, _ = build_consensus_and_mask(
        alignment_text, params.identity_threshold, params.coverage_threshold
    )
    forward: list[PrimerCandidate] = []
    reverse: list[PrimerCandidate] = []
    n = len(consensus)

    for start0 in range(n):
        for length in range(params.min_primer_len, params.max_primer_len + 1):
            end0 = start0 + length
            if end0 > n:
                break
            if not all(mask[start0:end0]):
                continue
            window = consensus[start0:end0]
            f = _candidate(window, "F", start0 + 1, end0, params)
            r = _candidate(window, "R", start0 + 1, end0, params)
            if f:
                forward.append(f)
            if r:
                reverse.append(r)

    return _pair_primer_candidates(forward, reverse, params)


def _find_reference_alignment(
    alignment_records: list[tuple[str, str]], reference_id: str
) -> str:
    exact = [sequence for record_id, sequence in alignment_records if record_id == reference_id]
    if len(exact) == 1:
        return exact[0]
    if len(exact) > 1:
        raise ValueError(
            f"O transcrito de referência {reference_id!r} aparece mais de uma vez no alinhamento."
        )

    # Alguns formatadores FASTA envolvem o acesso em identificadores com barras.
    aliases = [
        sequence
        for record_id, sequence in alignment_records
        if reference_id in record_id.split("|")
    ]
    if len(aliases) == 1:
        return aliases[0]
    if len(aliases) > 1:
        raise ValueError(
            f"O transcrito de referência {reference_id!r} é ambíguo no alinhamento."
        )
    raise ValueError(
        f"O transcrito de referência {reference_id!r} não foi encontrado no alinhamento."
    )


def _validated_transcript_exons(
    exons: Iterable[ExonInterval], transcript_length: int
) -> list[ExonInterval]:
    ordered = list(exons)
    if ordered != sorted(ordered, key=lambda exon: (exon.start, exon.end)):
        raise ValueError("Os éxons não estão na ordem 5′→3′ do transcrito.")
    for exon in ordered:
        if exon.start < 1 or exon.end < exon.start or exon.end > transcript_length:
            raise ValueError(
                "As coordenadas dos éxons são inválidas para o transcrito de referência."
            )
    for left, right in zip(ordered, ordered[1:]):
        if right.start != left.end + 1:
            raise ValueError(
                "Os éxons anotados não são contíguos nas coordenadas do transcrito. "
                "Use um registro de mRNA/RNA processado, não uma sequência genômica."
            )
    return ordered


def _junction_matches(
    *,
    start: int,
    end: int,
    orientation: str,
    reference_accession: str,
    exon_pairs: list[tuple[ExonInterval, ExonInterval]],
    min_5_prime_match: int,
    min_3_prime_match: int,
) -> list[ExonJunctionMatch]:
    matches: list[ExonJunctionMatch] = []
    for left_exon, right_exon in exon_pairs:
        junction_position = left_exon.end
        if not (start <= junction_position < end):
            continue
        left_bases = junction_position - start + 1
        right_bases = end - junction_position
        if orientation == "F":
            five_prime_bases, three_prime_bases = left_bases, right_bases
        else:
            five_prime_bases, three_prime_bases = right_bases, left_bases
        if (
            five_prime_bases < min_5_prime_match
            or three_prime_bases < min_3_prime_match
        ):
            continue
        matches.append(
            ExonJunctionMatch(
                reference_accession=reference_accession,
                junction_position=junction_position,
                left_exon_number=left_exon.number,
                right_exon_number=right_exon.number,
                primer_5_prime_bases=five_prime_bases,
                primer_3_prime_bases=three_prime_bases,
            )
        )
    return matches


def generate_exon_junction_primer_pairs(
    alignment_text: str,
    params: PrimerDesignParams,
    reference_id: str,
    exons: Iterable[ExonInterval],
    min_5_prime_match: int = 7,
    min_3_prime_match: int = 4,
) -> list[PrimerPair]:
    """Gera pares em que ao menos um primer cruza uma junção éxon–éxon.

    As coordenadas de ``exons`` são 1-based inclusivas sobre o transcrito
    processado escolhido como referência. Os candidatos continuam sujeitos aos
    mesmos filtros de conservação, GC, Tm e tamanho de amplicon do modo comum.
    """

    reference_id = reference_id.strip()
    _validate_primer_design_params(params)
    if not reference_id:
        raise ValueError("Informe o transcrito de referência para o desenho por junção.")
    if min_5_prime_match < 1 or min_3_prime_match < 1:
        raise ValueError("As ancoragens mínimas 5′ e 3′ devem ser maiores que zero.")

    alignment_records = parse_alignment_fasta(alignment_text)
    reference_alignment = _find_reference_alignment(alignment_records, reference_id)
    reference_columns = [
        column
        for column, base in enumerate(reference_alignment)
        if base not in {"-", "."}
    ]
    transcript_length = len(reference_columns)
    ordered_exons = _validated_transcript_exons(exons, transcript_length)
    if len(ordered_exons) < 2:
        return []

    exon_pairs = list(zip(ordered_exons, ordered_exons[1:]))
    junction_positions = [left_exon.end for left_exon, _ in exon_pairs]
    reference_sequence = "".join(reference_alignment[column] for column in reference_columns)
    aligned_sequences = [sequence for _, sequence in alignment_records]
    conserved: list[bool] = []
    for position, column in enumerate(reference_columns):
        reference_base = reference_sequence[position]
        canonical = [sequence[column] for sequence in aligned_sequences if sequence[column] in DNA]
        coverage = len(canonical) / len(aligned_sequences)
        identity = canonical.count(reference_base) / len(canonical) if canonical else 0.0
        conserved.append(
            reference_base in DNA
            and identity >= params.identity_threshold
            and coverage >= params.coverage_threshold
        )

    forward: list[PrimerCandidate] = []
    reverse: list[PrimerCandidate] = []
    relevant_starts: set[int] = set()
    for junction_position in junction_positions:
        first = max(0, junction_position - params.max_amplicon)
        last = min(transcript_length, junction_position + params.max_amplicon)
        relevant_starts.update(range(first, last))

    for start0 in sorted(relevant_starts):
        for length in range(params.min_primer_len, params.max_primer_len + 1):
            end0 = start0 + length
            if end0 > transcript_length:
                break
            if not all(conserved[start0:end0]):
                continue

            # Uma inserção alinhada dentro da janela quebraria a continuidade do
            # oligo nas outras sequências; não a comprima silenciosamente.
            first_column = reference_columns[start0]
            last_column = reference_columns[end0 - 1]
            if last_column - first_column + 1 != length:
                continue

            start, end = start0 + 1, end0
            window = reference_sequence[start0:end0]
            for orientation, output in (("F", forward), ("R", reverse)):
                candidate = _candidate(window, orientation, start, end, params)
                if candidate is None:
                    continue
                candidate.junctions = _junction_matches(
                    start=start,
                    end=end,
                    orientation=orientation,
                    reference_accession=reference_id,
                    exon_pairs=exon_pairs,
                    min_5_prime_match=min_5_prime_match,
                    min_3_prime_match=min_3_prime_match,
                )
                candidate.spans_exon_junction = bool(candidate.junctions)
                output.append(candidate)

    return _pair_primer_candidates(
        forward,
        reverse,
        params,
        require_exon_junction=True,
        reference_accession=reference_id,
    )


def records_to_fasta(records: Iterable[object]) -> str:
    chunks: list[str] = []
    for record in records:
        if getattr(record, "selected", True):
            chunks.append(record.fasta())
    return "".join(chunks)
