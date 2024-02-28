from typing import Tuple

from .helper_functions import IR


def ir_has_valid_gap_size(ir: IR):
    left_strand_end_idx: int = ir[0][1]
    right_strand_start_idx: int = ir[1][0]

    return right_strand_start_idx - left_strand_end_idx - 1 >= 3


def ir_pair_co_located(ir_a: IR, ir_b: IR) -> bool:
    ir_a_left_strand, ir_a_right_strand = ir_a[0], ir_a[1]

    # Accumulate paired base indices of ir_a
    paired_base_idxs_a = [
        idx for idx in range(ir_a_left_strand[0], ir_a_left_strand[1] + 1)
    ] + [idx for idx in range(ir_a_right_strand[0], ir_a_right_strand[1] + 1)]

    # Accumulate paired base indices of ir_b
    ir_b_left_strand, ir_b_right_strand = ir_b[0], ir_b[1]
    paired_base_idxs_b = [
        idx for idx in range(ir_b_left_strand[0], ir_b_left_strand[1] + 1)
    ] + [idx for idx in range(ir_b_right_strand[0], ir_b_right_strand[1] + 1)]

    # Return true if the intersection of the two set of accumulated indices is > 0
    return any([ir_b_bases in paired_base_idxs_a for ir_b_bases in paired_base_idxs_b])


def ir_pair_not_nested(ir_a: IR, ir_b: IR) -> bool:
    # ir_a comes entirely before ir_b
    ir_a_strictly_before_ir_b = ir_a[1][1] < ir_b[0][0]

    # ir_b comes entirely before ir_a
    ir_b_strictly_before_ir_a = ir_b[1][1] < ir_a[0][0]

    return ir_a_strictly_before_ir_b or ir_b_strictly_before_ir_a


def ir_pair_wholly_nested(ir_a: IR, ir_b: IR) -> bool:
    ir_a_left_strand_start_idx: int = ir_a[0][0]
    ir_a_left_strand_end_idx: int = ir_a[0][1]
    ir_a_right_strand_start_idx: int = ir_a[1][0]
    ir_a_right_strand_end_idx: int = ir_a[1][1]

    ir_b_left_strand_start_idx: int = ir_b[0][0]
    ir_b_left_strand_end_idx: int = ir_b[0][1]
    ir_b_right_strand_start_idx: int = ir_b[1][0]
    ir_b_right_strand_end_idx: int = ir_b[1][1]

    ir_a_in_ir_b: bool = (
        ir_b_left_strand_end_idx < ir_a_left_strand_start_idx
        and ir_b_right_strand_start_idx > ir_a_right_strand_end_idx
    )
    ir_b_in_ir_a: bool = (
        ir_a_left_strand_end_idx < ir_b_left_strand_start_idx
        and ir_a_right_strand_start_idx > ir_b_right_strand_end_idx
    )

    return ir_b_in_ir_a or ir_a_in_ir_b


def ir_pair_partially_nested(ir_a: IR, ir_b: IR) -> bool:
    """Returns true if ir_a is partially nested with ir_b."""
    ir_a_gap_first_idx: int = ir_a[0][1] + 1
    ir_a_gap_last_idx: int = ir_a[1][0] - 1
    ir_a_gap: Tuple[int, int] = (ir_a_gap_first_idx, ir_a_gap_last_idx)

    ir_b_gap_first_idx: int = ir_b[0][1] + 1
    ir_b_gap_last_idx: int = ir_b[1][0] - 1
    ir_b_gap: Tuple[int, int] = (ir_b_gap_first_idx, ir_b_gap_last_idx)

    def strand_in_gap(strand: Tuple[int, int], gap: Tuple[int, int]) -> bool:
        gap_start: int = gap[0]
        gap_stop: int = gap[1]

        if strand[0] >= gap_start and strand[1] <= gap_stop:
            return True
        return False

    # Check if one of ir_a's strands is contained in ir_b's gap and the other is not
    ir_a_left_strand: Tuple[int, int] = ir_a[0]
    ir_a_right_strand: Tuple[int, int] = ir_a[1]
    if (
        (
            strand_in_gap(ir_a_left_strand, ir_b_gap)
            and not strand_in_gap(ir_a_right_strand, ir_b_gap)
        )
    ) or (
        (
            strand_in_gap(ir_a_right_strand, ir_b_gap)
            and not strand_in_gap(ir_a_left_strand, ir_b_gap)
        )
    ):
        return True

    # Check if one of ir_b's strands is contained in ir_a's gap and the other is not
    ir_b_left_strand: Tuple[int, int] = ir_b[0]
    ir_b_right_strand: Tuple[int, int] = ir_b[1]
    if (
        (
            strand_in_gap(ir_b_left_strand, ir_a_gap)
            and not strand_in_gap(ir_b_right_strand, ir_a_gap)
        )
    ) or (
        (
            strand_in_gap(ir_b_right_strand, ir_a_gap)
            and not strand_in_gap(ir_b_left_strand, ir_a_gap)
        )
    ):
        return True

    return False


def ir_pair_invalid_relative_pos(ir_a: IR, ir_b: IR) -> bool:
    if ir_pair_co_located(ir_a, ir_b) or ir_pair_partially_nested(ir_a, ir_b):
        return True
