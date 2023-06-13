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
