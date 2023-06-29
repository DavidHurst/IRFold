import sys
import itertools


from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))


from irfold.util import (
    irs_to_dot_bracket,
    calc_free_energy,
    ir_pair_match_same_bases,
    ir_pair_forms_valid_loop,
    ir_has_valid_gap_size,
    irs_match_same_bases,
)


def label_irs_in_db_repr(db_repr, irs):
    """Assigns the index value of IRs in the list provided to each base that
    IRs pair respectively"""
    ids = ["A", "B", "C", "D", "E", "F", "G", "H"]
    db_repr_id_assigned = list(db_repr)

    for i, ir in enumerate(irs):
        n_base_pairs: int = ir[0][1] - ir[0][0] + 1

        lhs, rhs = ir[0], ir[1]

        db_repr_id_assigned[lhs[0] : lhs[1] + 1] = [ids[i] for _ in range(n_base_pairs)]
        db_repr_id_assigned[rhs[0] : rhs[1] + 1] = [ids[i] for _ in range(n_base_pairs)]

    return "".join(db_repr_id_assigned)


def ir_pair_free_energy_calculation_variations(ir_a, ir_b, seq_len, seq, out_dir):
    # FE(IR_a) + FE(IR_b)
    ir_a_db_repr = irs_to_dot_bracket([ir_a], seq_len)
    ir_b_db_repr = irs_to_dot_bracket([ir_b], seq_len)

    ir_a_fe = round(calc_free_energy(ir_a_db_repr, seq, out_dir), 4)
    ir_b_fe = round(calc_free_energy(ir_b_db_repr, seq, out_dir), 4)
    fe_sum = round(ir_a_fe + ir_b_fe, 4)

    # FE(IR_a U IR_b)
    ir_a_b_db_repr = irs_to_dot_bracket([ir_a, ir_b], seq_len)
    fe_union = round(calc_free_energy(ir_a_b_db_repr, seq, out_dir), 4)

    return ir_a_fe, ir_b_fe, fe_union, fe_sum


def get_all_ir_pairs_not_matching_same_bases_valid_gap_sz(all_irs):
    valid_irs = [ir for ir in all_irs if ir_has_valid_gap_size(ir)]

    unique_idx_pairs = list(
        itertools.combinations([i for i in range(len(valid_irs))], 2)
    )
    unique_ir_pairs = [(valid_irs[i], valid_irs[j]) for i, j in unique_idx_pairs]

    return [
        ir_pair
        for ir_pair in unique_ir_pairs
        if not ir_pair_match_same_bases(ir_pair[0], ir_pair[1])
    ]


def eval_ir_pair_structure_and_mfe(pair_idx, ir_a, ir_b, seq_len, seq, out_dir, print_out):
    (
        ir_a_fe,
        ir_b_fe,
        pairs_fe_union,
        pairs_fe_sum,
    ) = ir_pair_free_energy_calculation_variations(ir_a, ir_b, seq_len, seq, out_dir)

    additivity_assumption_held = pairs_fe_sum == pairs_fe_union
    valid_loop_formed = ir_pair_forms_valid_loop(ir_a, ir_b)

    db_repr = irs_to_dot_bracket([ir_a, ir_b], seq_len)
    db_repr_id_assigned = label_irs_in_db_repr(db_repr, [ir_a, ir_b])

    if print_out:
        ir_a_fmt = f"[{str(ir_a[0][0]).zfill(2)}->{str(ir_a[0][1]).zfill(2)} {str(ir_a[1][0]).zfill(2)}->{str(ir_a[1][1]).zfill(2)}]"
        ir_b_fmt = f"[{str(ir_b[0][0]).zfill(2)}->{str(ir_b[0][1]).zfill(2)} {str(ir_b[1][0]).zfill(2)}->{str(ir_b[1][1]).zfill(2)}]"

        print_space = " " * 4
        print(f"Pair #{pair_idx}: A = {ir_a_fmt}, B = {ir_b_fmt}")
        print(print_space, f"A DB    : {irs_to_dot_bracket([ir_a], seq_len)}")
        print(print_space, f"B DB    : {irs_to_dot_bracket([ir_b], seq_len)}")
        print(print_space, f"A U B DB: {db_repr}")
        print(print_space, f"A U B DB: {db_repr_id_assigned}")
        print()
        print(print_space, f"FE(A)        : {ir_a_fe:.4f}")
        print(print_space, f"FE(B)        : {ir_b_fe:.4f}")
        print(print_space, f"FE(A) + FE(B): {pairs_fe_sum:.4f}")
        print(print_space, f"FE(A U B)    : {pairs_fe_union:.4f}")
        print()
        print(print_space, f"Valid loop      : {valid_loop_formed}")
        print(print_space, f"Assumption holds: {additivity_assumption_held}")

    return (
        additivity_assumption_held,
        valid_loop_formed,
        ir_a_fe,
        ir_b_fe,
        pairs_fe_sum,
        pairs_fe_union,
    )


def ir_triplet_free_energy_calculation_variations(
    ir_a, ir_b, ir_c, seq_len, seq, out_dir
):
    # FE(IR_a) + FE(IR_b) + FE(IR_c)
    ir_a_db_repr = irs_to_dot_bracket([ir_a], seq_len)
    ir_b_db_repr = irs_to_dot_bracket([ir_b], seq_len)
    ir_c_db_repr = irs_to_dot_bracket([ir_c], seq_len)

    ir_a_fe = round(calc_free_energy(ir_a_db_repr, seq, out_dir), 4)
    ir_b_fe = round(calc_free_energy(ir_b_db_repr, seq, out_dir), 4)
    ir_c_fe = round(calc_free_energy(ir_c_db_repr, seq, out_dir), 4)
    fe_sum = round(ir_a_fe + ir_b_fe + ir_c_fe, 4)

    # FE(IR_a U IR_b U IR_c)
    ir_a_b_c_db_repr = irs_to_dot_bracket([ir_a, ir_b, ir_c], seq_len)
    fe_union = round(calc_free_energy(ir_a_b_c_db_repr, seq, out_dir), 4)

    return ir_a_fe, ir_b_fe, ir_c_fe, fe_union, fe_sum


def get_all_ir_triplets_not_matching_same_bases_valid_gap_sz(all_irs):
    valid_irs = [ir for ir in all_irs if ir_has_valid_gap_size(ir)]

    unique_idx_triplets = list(
        itertools.combinations([i for i in range(len(valid_irs))], 3)
    )
    unique_ir_triplets = [
        [valid_irs[i], valid_irs[j], valid_irs[k]] for i, j, k in unique_idx_triplets
    ]

    return [
        ir_triplet
        for ir_triplet in unique_ir_triplets
        if not irs_match_same_bases(ir_triplet)
    ]


def eval_ir_triplet_structure_and_mfe(
    pair_idx, ir_a, ir_b, ir_c, seq_len, seq, out_dir, print_out
):
    (
        ir_a_fe,
        ir_b_fe,
        ir_c_fe,
        triplet_fe_union,
        triplet_fe_sum,
    ) = ir_triplet_free_energy_calculation_variations(
        ir_a, ir_b, ir_c, seq_len, seq, out_dir
    )

    additivity_assumption_held = triplet_fe_sum == triplet_fe_union

    db_repr = irs_to_dot_bracket([ir_a, ir_b, ir_c], seq_len)
    db_repr_ir_labelled = label_irs_in_db_repr(db_repr, [ir_a, ir_b, ir_c])

    num_invalid_loop_forming_pairs = len(
        [
            True
            for ir_i, ir_j in itertools.combinations([ir_a, ir_b, ir_c], 2)
            if not ir_pair_forms_valid_loop(ir_i, ir_j)
        ]
    )

    if print_out:
        ir_a_fmt = (
            f"[{str(ir_a[0][0]).zfill(2)}->{str(ir_a[0][1]).zfill(2)} "
            f"{str(ir_a[1][0]).zfill(2)}->{str(ir_a[1][1]).zfill(2)}]"
        )
        ir_b_fmt = (
            f"[{str(ir_b[0][0]).zfill(2)}->{str(ir_b[0][1]).zfill(2)} "
            f"{str(ir_b[1][0]).zfill(2)}->{str(ir_b[1][1]).zfill(2)}]"
        )
        ir_c_fmt = (
            f"[{str(ir_c[0][0]).zfill(2)}->{str(ir_c[0][1]).zfill(2)} "
            f"{str(ir_c[1][0]).zfill(2)}->{str(ir_c[1][1]).zfill(2)}]"
        )

        print_space = " " * 4
        print(f"Triplet #{pair_idx}: A = {ir_a_fmt}, B = {ir_b_fmt}, C = {ir_c_fmt}")
        print(print_space, f"Dot bracket reprs:")
        print(
            print_space,
            print_space,
            f"A        : {irs_to_dot_bracket([ir_a], seq_len)}",
        )
        print(
            print_space,
            print_space,
            f"B        : {irs_to_dot_bracket([ir_b], seq_len)}",
        )
        print(
            print_space,
            print_space,
            f"C        : {irs_to_dot_bracket([ir_c], seq_len)}",
        )
        print(print_space, print_space, f"A U B U C: {db_repr}")
        print(print_space, print_space, f"A U B U C: {db_repr_ir_labelled}")
        print()
        print(print_space, f"FE(A)        : {ir_a_fe:.4f}")
        print(print_space, f"FE(B)        : {ir_b_fe:.4f}")
        print(print_space, f"FE(C)        : {ir_c_fe:.4f}")
        print(print_space, f"FE(A) + FE(B) + FE(C): {triplet_fe_sum:.4f}")
        print(print_space, f"FE(A U B U C)        : {triplet_fe_union:.4f}")
        print()
        print(print_space, f"Assumption holds: {additivity_assumption_held}")
        print(print_space, f"Num. invalid IR pairs: {num_invalid_loop_forming_pairs}")

    return (
        additivity_assumption_held,
        num_invalid_loop_forming_pairs,
        ir_a_fe,
        ir_b_fe,
        ir_c_fe,
        triplet_fe_sum,
        triplet_fe_union,
    )


def ir_quadruplet_free_energy_calculation_variations(
    ir_a, ir_b, ir_c, ir_d, sequence_len, seq, out_dir
):
    # FE(IR_a) + FE(IR_b) + FE(IR_c) + FE(IR_d)
    ir_a_db_repr = irs_to_dot_bracket([ir_a], sequence_len)
    ir_b_db_repr = irs_to_dot_bracket([ir_b], sequence_len)
    ir_c_db_repr = irs_to_dot_bracket([ir_c], sequence_len)
    ir_d_db_repr = irs_to_dot_bracket([ir_d], sequence_len)

    ir_a_fe = round(calc_free_energy(ir_a_db_repr, seq, out_dir), 4)
    ir_b_fe = round(calc_free_energy(ir_b_db_repr, seq, out_dir), 4)
    ir_c_fe = round(calc_free_energy(ir_c_db_repr, seq, out_dir), 4)
    ir_d_fe = round(calc_free_energy(ir_d_db_repr, seq, out_dir), 4)
    fe_sum = round(ir_a_fe + ir_b_fe + ir_c_fe, 4)

    # FE(IR_a U IR_b U IR_c U IR_d)
    ir_a_b_c_d_db_repr = irs_to_dot_bracket([ir_a, ir_b, ir_c, ir_d], sequence_len)
    fe_union = round(calc_free_energy(ir_a_b_c_d_db_repr, seq, out_dir), 4)

    return ir_a_fe, ir_b_fe, ir_c_fe, ir_d_fe, fe_union, fe_sum


def get_all_ir_quadruplets_not_matching_same_bases_valid_gap_sz(all_irs):
    valid_irs = [ir for ir in all_irs if ir_has_valid_gap_size(ir)]

    unique_idx_quadruplets = list(
        itertools.combinations([i for i in range(len(valid_irs))], 4)
    )
    unique_ir_quadruplets = [
        [valid_irs[i], valid_irs[j], valid_irs[k], valid_irs[l]]
        for i, j, k, l in unique_idx_quadruplets
    ]

    return [
        ir_quadruplet
        for ir_quadruplet in unique_ir_quadruplets
        if not irs_match_same_bases(ir_quadruplet)
    ]


def eval_ir_quadruplet_structure_and_mfe(
    pair_idx, ir_a, ir_b, ir_c, ir_d, seq_len, seq, out_dir, print_out
):
    (
        ir_a_fe,
        ir_b_fe,
        ir_c_fe,
        ir_d_fe,
        quadruplet_fe_union,
        quadruplet_fe_sum,
    ) = ir_quadruplet_free_energy_calculation_variations(
        ir_a, ir_b, ir_c, ir_d, seq_len, seq, out_dir
    )

    additivity_assumption_held = quadruplet_fe_sum == quadruplet_fe_union

    db_repr = irs_to_dot_bracket([ir_a, ir_b, ir_c, ir_d], seq_len)
    db_repr_ir_labelled = label_irs_in_db_repr(db_repr, [ir_a, ir_b, ir_c, ir_d])

    num_invalid_loop_forming_pairs = len(
        [
            True
            for ir_i, ir_j in itertools.combinations([ir_a, ir_b, ir_c, ir_d], 2)
            if not ir_pair_forms_valid_loop(ir_i, ir_j)
        ]
    )

    if print_out:
        ir_a_fmt = (
            f"[{str(ir_a[0][0]).zfill(2)}->{str(ir_a[0][1]).zfill(2)} "
            f"{str(ir_a[1][0]).zfill(2)}->{str(ir_a[1][1]).zfill(2)}]"
        )
        ir_b_fmt = (
            f"[{str(ir_b[0][0]).zfill(2)}->{str(ir_b[0][1]).zfill(2)} "
            f"{str(ir_b[1][0]).zfill(2)}->{str(ir_b[1][1]).zfill(2)}]"
        )
        ir_c_fmt = (
            f"[{str(ir_c[0][0]).zfill(2)}->{str(ir_c[0][1]).zfill(2)} "
            f"{str(ir_c[1][0]).zfill(2)}->{str(ir_c[1][1]).zfill(2)}]"
        )
        ir_d_fmt = (
            f"[{str(ir_d[0][0]).zfill(2)}->{str(ir_d[0][1]).zfill(2)} "
            f"{str(ir_d[1][0]).zfill(2)}->{str(ir_d[1][1]).zfill(2)}]"
        )

        print_space = " " * 4
        print(
            f"Quadruplet #{pair_idx}: A = {ir_a_fmt}, B = {ir_b_fmt}, C = {ir_c_fmt}, D = {ir_d_fmt}"
        )
        print(print_space, f"Dot bracket reprs:")
        print(
            print_space,
            print_space,
            f"A            : {irs_to_dot_bracket([ir_a], seq_len)}",
        )
        print(
            print_space,
            print_space,
            f"B            : {irs_to_dot_bracket([ir_b], seq_len)}",
        )
        print(
            print_space,
            print_space,
            f"C            : {irs_to_dot_bracket([ir_c], seq_len)}",
        )
        print(
            print_space,
            print_space,
            f"D            : {irs_to_dot_bracket([ir_d], seq_len)}",
        )
        print(print_space, print_space, f"A U B U C U D: {db_repr}")
        print(print_space, print_space, f"A U B U C U D: {db_repr_ir_labelled}")
        print()
        print(print_space, f"FE(A)        : {ir_a_fe:.4f}")
        print(print_space, f"FE(B)        : {ir_b_fe:.4f}")
        print(print_space, f"FE(C)        : {ir_c_fe:.4f}")
        print(print_space, f"FE(D)        : {ir_d_fe:.4f}")
        print(print_space, f"FE(A) + FE(B) + FE(C) + FE(D): {quadruplet_fe_sum:.4f}")
        print(print_space, f"FE(A U B U C U D)            : {quadruplet_fe_union:.4f}")
        print()
        print(print_space, f"Assumption holds: {additivity_assumption_held}")
        print(print_space, f"Num. invalid IR pairs: {num_invalid_loop_forming_pairs}")

    return (
        additivity_assumption_held,
        num_invalid_loop_forming_pairs,
        ir_a_fe,
        ir_b_fe,
        ir_c_fe,
        ir_d_fe,
        quadruplet_fe_sum,
        quadruplet_fe_union,
    )
