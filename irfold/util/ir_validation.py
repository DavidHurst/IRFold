from typing import Tuple

from .helper_functions import IR


def ir_pair_match_same_bases(ir_a: IR, ir_b: IR) -> bool:
    # Check if IRs match the same bases
    ir_a_left_strand, ir_a_right_strand = ir_a[0], ir_a[1]
    paired_base_idxs_a = [
        idx for idx in range(ir_a_left_strand[0], ir_a_left_strand[1] + 1)
    ] + [idx for idx in range(ir_a_right_strand[0], ir_a_right_strand[1] + 1)]

    ir_b_left_strand, ir_b_right_strand = ir_b[0], ir_b[1]
    paired_base_idxs_b = [
        idx for idx in range(ir_b_left_strand[0], ir_b_left_strand[1] + 1)
    ] + [idx for idx in range(ir_b_right_strand[0], ir_b_right_strand[1] + 1)]

    return any([ir_b_bases in paired_base_idxs_a for ir_b_bases in paired_base_idxs_b])


def ir_pair_disjoint(ir_a: IR, ir_b: IR) -> bool:
    # ir_a comes entirely before ir_b
    ir_a_strictly_before_ir_b = ir_a[1][1] < ir_b[0][0]

    # ir_b comes entirely before ir_a
    ir_b_strictly_before_ir_a = ir_b[1][1] < ir_a[0][0]

    return ir_a_strictly_before_ir_b or ir_b_strictly_before_ir_a


def ir_pair_wholly_nested(outer_ir: IR, nested_ir: IR) -> bool:
    # ToDo: Make this return true for nesting of A in B or B in A
    return outer_ir[0][0] < nested_ir[0][0] and nested_ir[1][1] < outer_ir[1][0]


def ir_pair_forms_valid_loop(ir_a: IR, ir_b: IR) -> bool:
    if (
        ir_pair_wholly_nested(ir_a, ir_b)
        or ir_pair_wholly_nested(ir_b, ir_a)
        or ir_pair_disjoint(ir_a, ir_b)
    ):
        return True

    ir_a_left_strand: Tuple[int, int] = ir_a[0]
    ir_b_left_strand: Tuple[int, int] = ir_b[0]

    latest_left_string_base_idx: int = (
        ir_a_left_strand[1]
        if ir_a_left_strand[1] > ir_b_left_strand[1]
        else ir_b_left_strand[1]
    )
    # ToDo: Introduce variables here for clarity
    earliest_right_string_base_idx: int = (
        ir_a[1][0] if ir_a[1][0] < ir_b[1][0] else ir_b[1][0]
    )

    # Invalid loop
    bases_inbetween_latest_left_and_earliest_right_bases: int = (
        earliest_right_string_base_idx - latest_left_string_base_idx - 1
    )
    if bases_inbetween_latest_left_and_earliest_right_bases < 3:
        return False

    return True

def ir_has_valid_gap_size(ir):
    left_strand_end_idx: int = ir[0][1]
    right_strand_start_idx: int = ir[1][0]

    return right_strand_start_idx - left_strand_end_idx - 1 >= 3
