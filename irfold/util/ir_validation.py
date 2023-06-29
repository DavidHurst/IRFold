import itertools
from typing import Tuple, List

from .helper_functions import IR


def ir_has_valid_gap_size(ir):
    left_strand_end_idx: int = ir[0][1]
    right_strand_start_idx: int = ir[1][0]

    return right_strand_start_idx - left_strand_end_idx - 1 >= 3


def ir_pair_match_same_bases(ir_a: IR, ir_b: IR) -> bool:
    ir_a_left_strand, ir_a_right_strand = ir_a[0], ir_a[1]
    paired_base_idxs_a = [
        idx for idx in range(ir_a_left_strand[0], ir_a_left_strand[1] + 1)
    ] + [idx for idx in range(ir_a_right_strand[0], ir_a_right_strand[1] + 1)]

    ir_b_left_strand, ir_b_right_strand = ir_b[0], ir_b[1]
    paired_base_idxs_b = [
        idx for idx in range(ir_b_left_strand[0], ir_b_left_strand[1] + 1)
    ] + [idx for idx in range(ir_b_right_strand[0], ir_b_right_strand[1] + 1)]

    return any([ir_b_bases in paired_base_idxs_a for ir_b_bases in paired_base_idxs_b])


def irs_match_same_bases(ir_list: List[IR]) -> bool:
    """Returns true if any IRs in the given list match the same bases"""
    return any(
        [
            ir_pair_match_same_bases(ir_a, ir_b)
            for ir_a, ir_b in list(itertools.combinations(ir_list, 2))
        ]
    )


# =================== Disjoint checks ===================


def ir_pair_disjoint(ir_a: IR, ir_b: IR) -> bool:
    # ir_a comes entirely before ir_b
    ir_a_strictly_before_ir_b = ir_a[1][1] < ir_b[0][0]

    # ir_b comes entirely before ir_a
    ir_b_strictly_before_ir_a = ir_b[1][1] < ir_a[0][0]

    return ir_a_strictly_before_ir_b or ir_b_strictly_before_ir_a


def ir_a_precedes_ir_b(ir_a: IR, ir_b: IR) -> bool:
    return ir_a[1][1] < ir_b[0][0]  # ir_a comes entirely before ir_b


def irs_disjoint(ir_list: List[IR]) -> bool:
    # Sort IRs by left strand start
    sorted_irs: List[IR] = sorted(ir_list, key=lambda ir: ir[0][0])

    # Check that first IR starts and ends before second IR which starts and ends before third IR...
    # ToDo: List comprehension-ify this
    for i in range(len(ir_list) - 1):
        preceding_ir: IR = sorted_irs[i]
        succeeding_ir: IR = sorted_irs[i + 1]

        if not ir_a_precedes_ir_b(preceding_ir, succeeding_ir):
            return False

    return True


# =================== Wholly nested checks ===================


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


def ir_a_subsumes_ir_b(ir_a: IR, ir_b: IR) -> bool:
    ir_a_left_strand_end_idx: int = ir_a[0][1]
    ir_a_right_strand_start_idx: int = ir_a[1][0]

    ir_b_left_strand_start_idx: int = ir_b[0][0]
    ir_b_right_strand_end_idx: int = ir_b[1][1]

    return (
        ir_a_left_strand_end_idx < ir_b_left_strand_start_idx
        and ir_a_right_strand_start_idx > ir_b_right_strand_end_idx
    )


def irs_wholly_nested(ir_list: List[IR]) -> bool:
    """Checks for one IR at each nesting level only i.e. first IR subsumes second which subsumes third IR etc."""
    # Sort IRs by left strand start
    sorted_irs: List[IR] = sorted(ir_list, key=lambda ir: ir[0][0])

    # Check that first IR subsumes second IR which subsumes third IR which subsumes fourth IR...
    # ToDo: List comprehension-ify this
    for i in range(len(ir_list) - 1):
        preceding_ir: IR = sorted_irs[i]
        succeeding_ir: IR = sorted_irs[i + 1]

        if not ir_a_subsumes_ir_b(preceding_ir, succeeding_ir):
            return False

    return True


# =================== Bases in intersection(s) checks ===================


def ir_pair_intersection_has_valid_base_count(ir_a: IR, ir_b: IR) -> bool:
    ir_a_left_strand: Tuple[int, int] = ir_a[0]
    ir_b_left_strand: Tuple[int, int] = ir_b[0]

    latest_left_string_base_idx: int = (
        ir_a_left_strand[1]
        if ir_a_left_strand[1] > ir_b_left_strand[1]
        else ir_b_left_strand[1]
    )

    earliest_right_string_base_idx: int = (
        ir_a[1][0] if ir_a[1][0] < ir_b[1][0] else ir_b[1][0]
    )

    num_bases_inbetween_latest_left_and_earliest_right_bases: int = (
        earliest_right_string_base_idx - latest_left_string_base_idx - 1
    )

    return num_bases_inbetween_latest_left_and_earliest_right_bases >= 3


# =================== Valid loop checks ===================


def ir_pair_forms_valid_loop(ir_a: IR, ir_b: IR) -> bool:
    return (
        ir_pair_wholly_nested(ir_a, ir_b)
        or ir_pair_disjoint(ir_a, ir_b)
        or ir_pair_intersection_has_valid_base_count(ir_a, ir_b)
    )


def irs_form_valid_loop(ir_list: List[IR]) -> bool:
    return (
        irs_disjoint(ir_list)
        or irs_wholly_nested(ir_list)
        or all(
            [
                ir_pair_intersection_has_valid_base_count(ir_a, ir_b)
                for ir_a, ir_b in list(itertools.combinations(ir_list, 2))
            ]
        )
    )
