from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest

from models import PrimerCandidate, PrimerPair
from primer_plot import (
    PixelSpan,
    build_primer_map_layout,
    nice_tick_positions,
    sequence_position_to_pixel,
    sequence_span_to_pixels,
)


def _candidate(orientation: str, start: int, end: int) -> PrimerCandidate:
    return PrimerCandidate(
        orientation=orientation,
        start=start,
        end=end,
        sequence="A" * (end - start + 1),
        length=end - start + 1,
        gc_percent=50.0,
        tm_c=60.0,
        score=1.0,
    )


def _pair(
    rank: int,
    forward: tuple[int, int],
    reverse: tuple[int, int],
    *,
    amplicon: tuple[int, int] | None = None,
) -> PrimerPair:
    amplicon_start, amplicon_end = amplicon or (forward[0], reverse[1])
    return PrimerPair(
        rank=rank,
        forward=_candidate("F", *forward),
        reverse=_candidate("R", *reverse),
        amplicon_start=amplicon_start,
        amplicon_end=amplicon_end,
        amplicon_length=amplicon_end - amplicon_start + 1,
        score=float(rank),
    )


def test_sequence_positions_map_extremes_and_middle_to_pixels():
    assert sequence_position_to_pixel(1, 101, 10, 210) == pytest.approx(10)
    assert sequence_position_to_pixel(51, 101, 10, 210) == pytest.approx(110)
    assert sequence_position_to_pixel(101, 101, 10, 210) == pytest.approx(210)


def test_one_base_target_maps_its_only_position_to_the_plot_middle():
    assert sequence_position_to_pixel(1, 1, 20, 80) == pytest.approx(50)
    assert sequence_span_to_pixels(1, 1, 1, 20, 80) == PixelSpan(20, 80)


def test_inclusive_sequence_span_maps_outer_base_boundaries():
    span = sequence_span_to_pixels(21, 41, 101, 10, 210)

    assert span.start_x == pytest.approx(10 + 20 / 101 * 200)
    assert span.end_x == pytest.approx(10 + 41 / 101 * 200)
    assert span.left_x == span.start_x
    assert span.right_x == span.end_x
    assert span.width == pytest.approx(21 / 101 * 200)


@pytest.mark.parametrize("target_length", [0, -1, 1.5, True])
def test_mapping_rejects_invalid_target_lengths(target_length):
    with pytest.raises((TypeError, ValueError)):
        sequence_position_to_pixel(1, target_length, 0, 100)


@pytest.mark.parametrize("position", [0, -1, 11, 1.5, True])
def test_mapping_rejects_positions_outside_the_1_based_target(position):
    with pytest.raises((TypeError, ValueError)):
        sequence_position_to_pixel(position, 10, 0, 100)


@pytest.mark.parametrize(
    ("left", "right"),
    [(10, 10), (20, 10), (float("nan"), 10), (0, float("inf")), (False, 10)],
)
def test_mapping_rejects_invalid_pixel_ranges(left, right):
    with pytest.raises((TypeError, ValueError)):
        sequence_position_to_pixel(1, 10, left, right)


@pytest.mark.parametrize(
    ("start", "end"),
    [(0, 5), (1, 11), (8, 7), (1.0, 5), (1, True)],
)
def test_span_mapping_rejects_invalid_inclusive_coordinates(start, end):
    with pytest.raises((TypeError, ValueError)):
        sequence_span_to_pixels(start, end, 10, 0, 100)


def test_nice_ticks_include_extremes_and_use_readable_intervals():
    assert nice_tick_positions(3358, max_ticks=9) == (
        1,
        500,
        1000,
        1500,
        2000,
        2500,
        3000,
        3358,
    )


def test_nice_ticks_return_every_position_for_a_short_target():
    assert nice_tick_positions(5, max_ticks=9) == (1, 2, 3, 4, 5)
    assert nice_tick_positions(1, max_ticks=9) == (1,)


@pytest.mark.parametrize("length", [11, 21, 99, 100, 101, 10_001])
def test_nice_ticks_never_exceed_the_requested_count(length):
    ticks = nice_tick_positions(length, max_ticks=5)

    assert ticks[0] == 1
    assert ticks[-1] == length
    assert len(ticks) <= 5
    assert tuple(sorted(set(ticks))) == ticks


@pytest.mark.parametrize(
    ("target_length", "max_ticks"),
    [(0, 9), (100, 1), (100, 0), (100, 2.5), (100, True)],
)
def test_nice_ticks_validate_length_and_tick_limit(target_length, max_ticks):
    with pytest.raises((TypeError, ValueError)):
        nice_tick_positions(target_length, max_ticks)


def test_layout_orders_pairs_by_rank_applies_limit_and_uses_deterministic_rows():
    pairs = [
        _pair(3, (30, 39), (70, 79)),
        _pair(1, (10, 19), (50, 59)),
        _pair(2, (20, 29), (60, 69)),
    ]

    layout = build_primer_map_layout(
        pairs,
        target_length=100,
        plot_left=0,
        plot_right=990,
        limit=2,
        first_row_y=120,
        row_height=54,
        bottom_padding=24,
    )

    assert [geometry.rank for geometry in layout.pairs] == [1, 2]
    assert [geometry.pair for geometry in layout.pairs] == [pairs[1], pairs[2]]
    assert [geometry.y for geometry in layout.pairs] == [120, 174]
    assert layout.height == pytest.approx(225)
    assert layout.target_length == 100
    assert layout.plot_left == 0
    assert layout.plot_right == 990


def test_layout_builds_forward_right_and_reverse_left_arrow_geometry():
    pair = _pair(1, (10, 29), (60, 80))

    geometry = build_primer_map_layout([pair], 100, 0, 990).pairs[0]

    assert geometry.forward.direction == "right"
    assert geometry.forward.points_right is True
    assert geometry.forward.points_left is False
    assert geometry.forward.tail_x == pytest.approx(89.1)
    assert geometry.forward.tip_x == pytest.approx(287.1)
    assert geometry.forward.start_x == pytest.approx(89.1)
    assert geometry.forward.end_x == pytest.approx(287.1)

    assert geometry.reverse.direction == "left"
    assert geometry.reverse.points_left is True
    assert geometry.reverse.points_right is False
    assert geometry.reverse.tail_x == pytest.approx(792)
    assert geometry.reverse.tip_x == pytest.approx(584.1)
    assert geometry.reverse.start_x == pytest.approx(584.1)
    assert geometry.reverse.end_x == pytest.approx(792)

    assert geometry.amplicon.start_x == pytest.approx(89.1)
    assert geometry.amplicon.end_x == pytest.approx(792)


def test_empty_layout_has_no_rows_and_reserves_header_plus_bottom_padding():
    layout = build_primer_map_layout(
        [], 50, 15, 500, first_row_y=100, row_height=40, bottom_padding=12
    )

    assert layout.pairs == ()
    assert layout.height == pytest.approx(112)


def test_layout_dataclasses_are_frozen_and_rows_are_a_tuple():
    layout = build_primer_map_layout([_pair(1, (1, 10), (20, 30))], 30, 0, 290)

    assert isinstance(layout.pairs, tuple)
    with pytest.raises(FrozenInstanceError):
        layout.height = 999
    with pytest.raises(FrozenInstanceError):
        layout.pairs[0].forward.span.start_x = 999


@pytest.mark.parametrize(
    ("name", "value"),
    [
        ("limit", 0),
        ("limit", True),
        ("first_row_y", -1),
        ("first_row_y", float("nan")),
        ("row_height", 0),
        ("row_height", -1),
        ("bottom_padding", -1),
        ("bottom_padding", float("inf")),
    ],
)
def test_layout_validates_limit_and_vertical_dimensions(name, value):
    kwargs = {name: value}
    with pytest.raises((TypeError, ValueError)):
        build_primer_map_layout([_pair(1, (1, 10), (20, 30))], 30, 0, 290, **kwargs)


def test_layout_rejects_duplicate_ranks():
    with pytest.raises(ValueError, match="ranks.*únicos"):
        build_primer_map_layout(
            [_pair(1, (1, 10), (20, 30)), _pair(1, (2, 11), (21, 31))],
            40,
            0,
            390,
        )


def test_layout_validates_pairs_even_when_they_fall_outside_the_limit():
    invalid_hidden_pair = _pair(2, (40, 49), (80, 110))

    with pytest.raises(ValueError, match="excede o comprimento"):
        build_primer_map_layout(
            [_pair(1, (1, 10), (20, 30)), invalid_hidden_pair],
            100,
            0,
            990,
            limit=1,
        )


def test_layout_rejects_wrong_primer_orientation():
    pair = _pair(1, (1, 10), (20, 30))
    pair.forward.orientation = "R"

    with pytest.raises(ValueError, match="forward.*orientação F"):
        build_primer_map_layout([pair], 30, 0, 290)


def test_layout_rejects_amplicon_length_inconsistent_with_inclusive_coordinates():
    pair = _pair(1, (1, 10), (20, 30))
    pair.amplicon_length = 29

    with pytest.raises(ValueError, match="amplicon_length"):
        build_primer_map_layout([pair], 30, 0, 290)


def test_layout_rejects_amplicon_that_does_not_contain_both_primers():
    pair = _pair(1, (1, 10), (20, 30), amplicon=(2, 30))

    with pytest.raises(ValueError, match="conter integralmente"):
        build_primer_map_layout([pair], 30, 0, 290)
