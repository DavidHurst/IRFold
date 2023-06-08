import functools
import random
import sys
import itertools
import warnings

warnings.filterwarnings("error")

from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from src import IRFold


def irs_nested(outer_ir, nested_ir):
    return outer_ir[0][0] < nested_ir[0][0] and nested_ir[1][1] < outer_ir[1][0]


def irs_fe_union_fe_sum(ir_a, ir_b):
    # FE(IR_a) + FE(IR_b)
    ir_a_db_repr = IRFold.irs_to_dot_bracket([ir_a], seq_len)
    ir_b_db_repr = IRFold.irs_to_dot_bracket([ir_b], seq_len)

    ir_a_fe = round(IRFold.calc_free_energy(ir_a_db_repr, seq, data_dir), 4)
    ir_b_fe = round(IRFold.calc_free_energy(ir_b_db_repr, seq, data_dir), 4)
    fe_sum = round(ir_a_fe + ir_b_fe, 4)

    # FE(IR_a U IR_b)
    ir_a_b_db_repr = IRFold.irs_to_dot_bracket([ir_a, ir_b], seq_len)
    fe_union = round(IRFold.calc_free_energy(ir_a_b_db_repr, seq, data_dir), 4)

    return (ir_a_fe, ir_b_fe), fe_union, fe_sum


def id_base_pairs(db_repr, irs):
    """Assigns the index value of IRs in the list provided to each base that
    IRs pair respectively"""
    ids = ["A", "B", "C", "D", "E"]
    db_repr_id_assigned = list(db_repr)

    for i, ir in enumerate(irs):
        n_base_pairs: int = ir[0][1] - ir[0][0] + 1

        lhs, rhs = ir[0], ir[1]

        # IUPACpal returns base pairings using 1-based indexing
        db_repr_id_assigned[lhs[0] - 1 : lhs[1]] = [ids[i] for _ in range(n_base_pairs)]
        db_repr_id_assigned[rhs[0] - 1 : rhs[1]] = [ids[i] for _ in range(n_base_pairs)]

    return "".join(db_repr_id_assigned)


if __name__ == "__main__":
    verbose = True
    random.seed(4)

    data_dir = "."
    seq_len = 35
    seq = "".join(random.choice("ACGU") for _ in range(seq_len))

    print(f"Seq. length: {seq_len}")
    print(f"Seq.       : {seq}")

    find_irs_kwargs = {
        "sequence": seq,
        "min_len": 2,
        "max_len": seq_len,
        "max_gap": seq_len - 1,
        "mismatches": 0,
        "out_dir": data_dir,
    }
    found_irs = IRFold.find_irs(**find_irs_kwargs)
    n_irs = len(found_irs)

    # Find all compatible IR pairs
    print("Compatible IR Pairs".center(50, "="))
    # ToDo: Also evaluate triples, quadruples etc.
    unique_idx_pairs = list(itertools.combinations([i for i in range(n_irs)], 2))
    unique_ir_pairs = [(found_irs[i], found_irs[j]) for i, j in unique_idx_pairs]
    compatible_ir_pairs = [
        ir_pair
        for ir_pair in unique_ir_pairs
        if not IRFold.irs_incompatible(ir_pair[0], ir_pair[1])
    ]

    print(f"IUPACpal Parameters:")
    for kw, kwarg in find_irs_kwargs.items():
        if kw.startswith("m"):
            print(f"  - {kw} = {kwarg}")
    print(f"Num. IRs found         : {n_irs}")
    print(f"Num. unique pairs      : {len(unique_ir_pairs)}")
    print(f"Num. compatible pairs  : {len(compatible_ir_pairs)}")
    print(f"Num. incompatible pairs: {len(unique_ir_pairs) - len(compatible_ir_pairs)}")

    # Classify compatible IR pairs as nested or not nested
    print("Nested Compatible IR Pairs".center(50, "="))
    nested_ir_pairs, non_nested_ir_pairs = [], []
    for p in compatible_ir_pairs:
        ir_i, ir_j = p[0], p[1]

        if irs_nested(ir_i, ir_j) or irs_nested(ir_j, ir_i):
            nested_ir_pairs.append(p)
        else:
            non_nested_ir_pairs.append(p)

    n_nested_pairs = len(nested_ir_pairs)
    n_non_nested_pairs = len(non_nested_ir_pairs)
    print(f"Num. nested pairs    : {n_nested_pairs}")
    print(f"Num. non-nested pairs: {n_non_nested_pairs}")

    print("Nested Compatible IR Pair Free Energies".center(50, "="))
    assumption_holds_count_nested_compatible = []
    for p in nested_ir_pairs:

        ir_i, ir_j = p[0], p[1]

        (ir_i_fe, ir_j_fe), fe_union, fe_sum = irs_fe_union_fe_sum(ir_i, ir_j)

        assumption_held = fe_sum == fe_union
        assumption_holds_count_nested_compatible.append(assumption_held)

        if verbose:
            ir_i_fmt = f"[{str(ir_i[0][0]).zfill(2)}->{str(ir_i[0][1]).zfill(2)} {str(ir_i[1][0]).zfill(2)}->{str(ir_i[1][1]).zfill(2)}]"
            ir_j_fmt = f"[{str(ir_j[0][0]).zfill(2)}->{str(ir_j[0][1]).zfill(2)} {str(ir_j[1][0]).zfill(2)}->{str(ir_j[1][1]).zfill(2)}]"

            db_repr = IRFold.irs_to_dot_bracket([ir_i, ir_j], seq_len)
            db_repr_id_assigned = id_base_pairs(db_repr, [ir_i, ir_j])

            print(f"Pair:")
            print(f"  A = {ir_i_fmt}")
            print(f"  B = {ir_j_fmt}")
            if irs_nested(ir_i, ir_j):
                print(f"  B is nested in A:")
            else:
                print(f"  A is nested in B:")
            print(f"    A U B = {db_repr}")
            print(f"    A U B = {db_repr_id_assigned}")
            print(f"  FE(A) = {ir_i_fe:.4f}")
            if ir_i_fe > 9000:
                print(f"    [i] A would never be selected")
            print(f"  FE(B) = {ir_j_fe:.4f}")
            if ir_j_fe > 9000:
                print(f"    [i] B would never be selected")
            print(f"  FE(A) + FE(B): {fe_sum:.4f}")
            print(f"  FE(A U B)    : {fe_union:.4f}")
            print(f"  Assumption holds: {assumption_held}")
            # print(f"  FE(A) + FE(B) == FE(A U B): {assumption_held}")

    print("Non-Nested Compatible IR Pair Free Energies".center(50, "="))
    assumption_holds_count_non_nested_compatible = []
    for p in non_nested_ir_pairs:
        ir_i, ir_j = p[0], p[1]

        (ir_i_fe, ir_j_fe), fe_union, fe_sum = irs_fe_union_fe_sum(ir_i, ir_j)

        assumption_held = fe_sum == fe_union
        assumption_holds_count_non_nested_compatible.append(assumption_held)

        if verbose:
            ir_i_fmt = f"[{str(ir_i[0][0]).zfill(2)}->{str(ir_i[0][1]).zfill(2)} {str(ir_i[1][0]).zfill(2)}->{str(ir_i[1][1]).zfill(2)}]"
            ir_j_fmt = f"[{str(ir_j[0][0]).zfill(2)}->{str(ir_j[0][1]).zfill(2)} {str(ir_j[1][0]).zfill(2)}->{str(ir_j[1][1]).zfill(2)}]"

            db_repr = IRFold.irs_to_dot_bracket([ir_i, ir_j], seq_len)
            db_repr_id_assigned = id_base_pairs(db_repr, [ir_i, ir_j])

            print(f"Pair:")
            print(f"  A = {ir_i_fmt}")
            print(f"  B = {ir_j_fmt}")
            print(f"  A U B = {db_repr}")
            print(f"  A U B = {db_repr_id_assigned}")
            print(f"  FE(A) = {ir_i_fe:.4f}")
            if ir_i_fe > 9000:
                print(f"    [i] A would never be selected")
            print(f"  FE(B) = {ir_j_fe:.4f}")
            if ir_j_fe > 9000:
                print(f"    [i] B would never be selected")
            print(f"  FE(A) + FE(B): {fe_sum:.4f}")
            print(f"  FE(A U B)    : {fe_union:.4f}")
            print(f"  Assumption holds: {assumption_held}")
            # print(f"  FE(A) + FE(B) == FE(A U B): {assumption_held}")

    print("Assumption Holding Summary Stats.".center(50, "="))
    print("Nested:")
    print(f"  Num. pairs: {n_nested_pairs}")
    times_assumption_held_nested = len(
        [res for res in assumption_holds_count_nested_compatible if res == True]
    )
    print(f"  Assumption Held: {times_assumption_held_nested}/{n_nested_pairs} times")

    print("Non-Nested:")
    print(f"  Num. pairs: {n_non_nested_pairs}")
    times_assumption_held_not_nested = len(
        [res for res in assumption_holds_count_non_nested_compatible if res == True]
    )
    print(
        f"  Assumption Held: {times_assumption_held_not_nested}/{n_non_nested_pairs} times"
    )

    print("Overall:")
    print(f"  Num. pairs: {len(compatible_ir_pairs)}")
    total_times_assumption_held = (
        times_assumption_held_nested + times_assumption_held_not_nested
    )
    print(
        f"  Assumption Held: {total_times_assumption_held}/{len(compatible_ir_pairs)} times"
    )
