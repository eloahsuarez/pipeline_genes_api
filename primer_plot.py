from __future__ import annotations

import math
from dataclasses import dataclass
from numbers import Real
from typing import Iterable, Literal

from models import PrimerPair


Direction = Literal["left", "right"]
Orientation = Literal["F", "R"]


@dataclass(frozen=True, slots=True)
class PixelSpan:
    """Intervalo horizontal fechado, já convertido para coordenadas do canvas."""

    start_x: float
    end_x: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.start_x) or not math.isfinite(self.end_x):
            raise ValueError("As coordenadas em pixels devem ser finitas.")
        if self.end_x < self.start_x:
            raise ValueError("O fim do intervalo em pixels não pode preceder o início.")

    @property
    def left_x(self) -> float:
        return self.start_x

    @property
    def right_x(self) -> float:
        return self.end_x

    @property
    def width(self) -> float:
        return self.end_x - self.start_x


@dataclass(frozen=True, slots=True)
class PrimerArrowGeometry:
    """Intervalo de um primer e sua direção 5′→3′ no canvas."""

    orientation: Orientation
    span: PixelSpan

    def __post_init__(self) -> None:
        if self.orientation not in {"F", "R"}:
            raise ValueError("A orientação do primer deve ser 'F' ou 'R'.")

    @property
    def direction(self) -> Direction:
        return "right" if self.orientation == "F" else "left"

    @property
    def points_right(self) -> bool:
        return self.orientation == "F"

    @property
    def points_left(self) -> bool:
        return self.orientation == "R"

    @property
    def tail_x(self) -> float:
        return self.span.start_x if self.points_right else self.span.end_x

    @property
    def tip_x(self) -> float:
        return self.span.end_x if self.points_right else self.span.start_x

    @property
    def start_x(self) -> float:
        return self.span.start_x

    @property
    def end_x(self) -> float:
        return self.span.end_x


@dataclass(frozen=True, slots=True)
class PrimerPairGeometry:
    """Geometria horizontal e posição vertical de um par de primers."""

    pair: PrimerPair
    rank: int
    y: float
    forward: PrimerArrowGeometry
    reverse: PrimerArrowGeometry
    amplicon: PixelSpan


@dataclass(frozen=True, slots=True)
class PrimerMapLayout:
    """Resultado imutável e independente de Tkinter para renderizar o mapa."""

    pairs: tuple[PrimerPairGeometry, ...]
    height: float
    target_length: int
    plot_left: float
    plot_right: float
    first_row_y: float
    row_height: float
    bottom_padding: float


def _positive_integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} deve ser um número inteiro.")
    if value < 1:
        raise ValueError(f"{name} deve ser maior que zero.")
    return value


def _finite_number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} deve ser um número real.")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} deve ser finito.")
    return number


def _validated_scale(
    target_length: int, plot_left: float, plot_right: float
) -> tuple[int, float, float]:
    length = _positive_integer(target_length, "target_length")
    left = _finite_number(plot_left, "plot_left")
    right = _finite_number(plot_right, "plot_right")
    if right <= left:
        raise ValueError("plot_right deve ser maior que plot_left.")
    return length, left, right


def _validated_sequence_position(position: int, target_length: int, name: str) -> int:
    value = _positive_integer(position, name)
    if value > target_length:
        raise ValueError(
            f"{name}={value} excede o comprimento da sequência ({target_length})."
        )
    return value


def sequence_position_to_pixel(
    position: int,
    target_length: int,
    plot_left: float,
    plot_right: float,
) -> float:
    """Mapeia uma posição 1-based para o centro correspondente no eixo horizontal.

    As posições 1 e ``target_length`` coincidem com ``plot_left`` e
    ``plot_right``. Para uma sequência de uma base, a única posição ocupa o meio
    do intervalo disponível.
    """

    length, left, right = _validated_scale(target_length, plot_left, plot_right)
    value = _validated_sequence_position(position, length, "position")
    if length == 1:
        return (left + right) / 2.0
    return left + ((value - 1) / (length - 1)) * (right - left)


def sequence_span_to_pixels(
    start: int,
    end: int,
    target_length: int,
    plot_left: float,
    plot_right: float,
) -> PixelSpan:
    """Converte um intervalo 1-based inclusivo em fronteiras de bases no canvas.

    Diferentemente das posições pontuais da régua, um intervalo ocupa a largura
    integral de suas bases: a borda esquerda é ``(start - 1) / N`` e a direita
    é ``end / N``. Assim, até um alvo de uma única base ocupa todo o eixo.
    """

    length, left, right = _validated_scale(target_length, plot_left, plot_right)
    first = _validated_sequence_position(start, length, "start")
    last = _validated_sequence_position(end, length, "end")
    if last < first:
        raise ValueError("end não pode ser menor que start.")
    plot_width = right - left
    return PixelSpan(
        left + ((first - 1) / length) * plot_width,
        left + (last / length) * plot_width,
    )


def _nice_step_ceiling(raw_step: float) -> int:
    exponent = math.floor(math.log10(raw_step))
    magnitude = 10**exponent
    fraction = raw_step / magnitude
    for preferred in (1, 2, 5, 10):
        if fraction <= preferred:
            return max(1, int(preferred * magnitude))
    raise AssertionError("Não foi possível calcular o intervalo da régua.")


def _ticks_for_step(target_length: int, step: int) -> tuple[int, ...]:
    ticks = {1, target_length}
    ticks.update(range(step, target_length, step))
    return tuple(sorted(ticks))


def nice_tick_positions(target_length: int, max_ticks: int = 9) -> tuple[int, ...]:
    """Retorna marcas legíveis da régua, sempre incluindo os extremos.

    Os intervalos intermediários usam a progressão convencional 1, 2, 5 × 10ⁿ.
    O total retornado nunca excede ``max_ticks``.
    """

    length = _positive_integer(target_length, "target_length")
    maximum = _positive_integer(max_ticks, "max_ticks")
    if maximum < 2:
        raise ValueError("max_ticks deve ser pelo menos 2.")
    if length == 1:
        return (1,)
    if length <= maximum:
        return tuple(range(1, length + 1))

    raw_step = (length - 1) / (maximum - 1)
    step = _nice_step_ceiling(raw_step)
    ticks = _ticks_for_step(length, step)
    while len(ticks) > maximum:
        step = _nice_step_ceiling(step * 1.0000001)
        ticks = _ticks_for_step(length, step)
    return ticks


def _pair_rank(pair: PrimerPair) -> int:
    try:
        return _positive_integer(pair.rank, "pair.rank")
    except AttributeError as exc:
        raise ValueError("Todo par deve informar rank.") from exc


def _candidate_span(
    pair: PrimerPair,
    candidate_name: Literal["forward", "reverse"],
    expected_orientation: Orientation,
    target_length: int,
    plot_left: float,
    plot_right: float,
) -> tuple[int, int, PrimerArrowGeometry]:
    try:
        candidate = getattr(pair, candidate_name)
        orientation = candidate.orientation
        start = candidate.start
        end = candidate.end
    except AttributeError as exc:
        raise ValueError(f"O par não contém um primer {candidate_name} completo.") from exc

    if orientation != expected_orientation:
        raise ValueError(
            f"O primer {candidate_name} deve ter orientação {expected_orientation}."
        )
    span = sequence_span_to_pixels(start, end, target_length, plot_left, plot_right)
    return start, end, PrimerArrowGeometry(expected_orientation, span)


def _validate_pair_and_build_geometry(
    pair: PrimerPair,
    rank: int,
    y: float,
    target_length: int,
    plot_left: float,
    plot_right: float,
) -> PrimerPairGeometry:
    forward_start, forward_end, forward = _candidate_span(
        pair, "forward", "F", target_length, plot_left, plot_right
    )
    reverse_start, reverse_end, reverse = _candidate_span(
        pair, "reverse", "R", target_length, plot_left, plot_right
    )
    try:
        amplicon_start = pair.amplicon_start
        amplicon_end = pair.amplicon_end
        amplicon_length = pair.amplicon_length
    except AttributeError as exc:
        raise ValueError("O par não contém coordenadas completas do amplicon.") from exc

    amplicon = sequence_span_to_pixels(
        amplicon_start, amplicon_end, target_length, plot_left, plot_right
    )
    expected_length = amplicon_end - amplicon_start + 1
    actual_length = _positive_integer(amplicon_length, "pair.amplicon_length")
    if actual_length != expected_length:
        raise ValueError(
            "pair.amplicon_length não corresponde às coordenadas inclusivas do amplicon."
        )
    if not (
        amplicon_start <= forward_start <= forward_end <= amplicon_end
        and amplicon_start <= reverse_start <= reverse_end <= amplicon_end
    ):
        raise ValueError("O amplicon deve conter integralmente os dois primers.")

    return PrimerPairGeometry(
        pair=pair,
        rank=rank,
        y=y,
        forward=forward,
        reverse=reverse,
        amplicon=amplicon,
    )


def build_primer_map_layout(
    pairs: Iterable[PrimerPair],
    target_length: int,
    plot_left: float,
    plot_right: float,
    limit: int | None = None,
    first_row_y: float = 120.0,
    row_height: float = 54.0,
    bottom_padding: float = 24.0,
) -> PrimerMapLayout:
    """Cria a geometria dos melhores pares, ordenada por rank.

    ``first_row_y`` é o centro vertical da primeira linha. Cada linha seguinte
    avança exatamente ``row_height``. A altura inclui meia linha após o último
    centro e então ``bottom_padding``; sem pares, reserva apenas o cabeçalho e o
    preenchimento inferior.
    """

    length, left, right = _validated_scale(target_length, plot_left, plot_right)
    first_y = _finite_number(first_row_y, "first_row_y")
    row = _finite_number(row_height, "row_height")
    bottom = _finite_number(bottom_padding, "bottom_padding")
    if first_y < 0:
        raise ValueError("first_row_y não pode ser negativo.")
    if row <= 0:
        raise ValueError("row_height deve ser maior que zero.")
    if bottom < 0:
        raise ValueError("bottom_padding não pode ser negativo.")

    maximum: int | None = None
    if limit is not None:
        maximum = _positive_integer(limit, "limit")

    try:
        supplied_pairs = list(pairs)
    except TypeError as exc:
        raise TypeError("pairs deve ser iterável.") from exc

    ranked = [(_pair_rank(pair), pair) for pair in supplied_pairs]
    ranks = [rank for rank, _ in ranked]
    if len(set(ranks)) != len(ranks):
        raise ValueError("Os ranks dos pares devem ser únicos.")
    ranked.sort(key=lambda item: item[0])

    # Valide também pares que ficariam fora do limite. Assim, o layout nunca
    # mascara dados inconsistentes apenas porque uma linha deixou de ser exibida.
    validated: list[PrimerPairGeometry] = []
    for index, (rank, pair) in enumerate(ranked):
        validated.append(
            _validate_pair_and_build_geometry(
                pair,
                rank,
                first_y + index * row,
                length,
                left,
                right,
            )
        )

    if maximum is not None:
        validated = validated[:maximum]

    # Recalcule Y depois do corte para manter o contrato independente da entrada.
    geometries = tuple(
        PrimerPairGeometry(
            pair=geometry.pair,
            rank=geometry.rank,
            y=first_y + index * row,
            forward=geometry.forward,
            reverse=geometry.reverse,
            amplicon=geometry.amplicon,
        )
        for index, geometry in enumerate(validated)
    )
    if geometries:
        height = geometries[-1].y + row / 2.0 + bottom
    else:
        height = first_y + bottom

    return PrimerMapLayout(
        pairs=geometries,
        height=height,
        target_length=length,
        plot_left=left,
        plot_right=right,
        first_row_y=first_y,
        row_height=row,
        bottom_padding=bottom,
    )


__all__ = [
    "PixelSpan",
    "PrimerArrowGeometry",
    "PrimerMapLayout",
    "PrimerPairGeometry",
    "build_primer_map_layout",
    "nice_tick_positions",
    "sequence_position_to_pixel",
    "sequence_span_to_pixels",
]
