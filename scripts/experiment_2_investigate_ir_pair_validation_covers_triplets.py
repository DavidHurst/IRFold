import random
import sys
import itertools
import pandas as pd
from math import comb

from pathlib import Path


sys.path.append(str(Path(__file__).resolve().parents[1]))

from irfold import IRFoldBase
from irfold.util import (
    irs_to_dot_bracket,
    calc_free_energy,
    ir_pair_match_same_bases,
    ir_pair_forms_valid_loop,
    ir_pair_wholly_nested,
    ir_pair_disjoint,
    ir_has_valid_gap_size,
    irs_match_same_bases,
    irs_form_valid_loop,
)

DATA_DIR = str(Path(__file__).parent.parent / "data")


def ir_triplet_free_energy_calculation_variations(
    ir_a, ir_b, ir_c, sequence_len, out_dir
):
    # FE(IR_a) + FE(IR_b) + FE(IR_c)
    ir_a_db_repr = irs_to_dot_bracket([ir_a], sequence_len)
    ir_b_db_repr = irs_to_dot_bracket([ir_b], sequence_len)
    ir_c_db_repr = irs_to_dot_bracket([ir_c], sequence_len)

    ir_a_fe = round(calc_free_energy(ir_a_db_repr, seq, out_dir), 4)
    ir_b_fe = round(calc_free_energy(ir_b_db_repr, seq, out_dir), 4)
    ir_c_fe = round(calc_free_energy(ir_c_db_repr, seq, out_dir), 4)
    fe_sum = round(ir_a_fe + ir_b_fe + ir_c_fe, 4)

    # FE(IR_a U IR_b U IR_c)
    ir_a_b_c_db_repr = irs_to_dot_bracket([ir_a, ir_b, ir_c], sequence_len)
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


def eval_ir_triplet_structure_and_mfe(pair_idx, ir_a, ir_b, ir_c, seq_len, print_out):
    (
        ir_a_fe,
        ir_b_fe,
        ir_c_fe,
        triplet_fe_union,
        triplet_fe_sum,
    ) = ir_triplet_free_energy_calculation_variations(
        ir_a, ir_b, ir_c, seq_len, DATA_DIR
    )

    additivity_assumption_held = triplet_fe_sum == triplet_fe_union

    db_repr = irs_to_dot_bracket([ir_a, ir_b, ir_c], seq_len)
    db_repr_ir_labelled = label_irs_in_db_repr(db_repr, [ir_a, ir_b, ir_c])

    contains_invalid_loop_pair = any(
        [
            not ir_pair_forms_valid_loop(ir_i, ir_j)
            for ir_i, ir_j in itertools.combinations([ir_a, ir_b, ir_c], 2)
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
        print(
            print_space, f"Triplet contains invalid pair: {contains_invalid_loop_pair}"
        )

    return (
        additivity_assumption_held,
        contains_invalid_loop_pair,
        ir_a_fe,
        ir_b_fe,
        ir_c_fe,
        triplet_fe_sum,
        triplet_fe_union,
    )


if __name__ == "__main__":
    """Evaluate the IR pairs found in a single sequence"""

    random.seed(223332)
    experiment_results = {
        "seq": [],
        "seq_len": [],
        "ir_triplet_index": [],
        "ir_a": [],
        "ir_b": [],
        "ir_c": [],
        "ir_pairs": [],
        "ir_a_fe": [],
        "ir_b_fe": [],
        "ir_c_fe": [],
        "ir_triplet_fe_sum": [],
        "ir_triplet_fe_union": [],
        "fe_additivity_assumption_held": [],
        "ir_triplet_contains_invalid_loop_forming_pair": [],
    }

    sequence_len = 30
    seq = "UGAUGACAAAUGCUUAACCCAAGCACGGCA"  # "".join(random.choice("ACGU") for _ in range(sequence_len))

    print(f"Seq. length: {sequence_len}")
    print(f"Seq.       : {seq}")

    find_irs_kwargs = {
        "sequence": seq,
        "min_len": 2,
        "max_len": sequence_len,
        "max_gap": sequence_len - 1,
        "mismatches": 0,
        "out_dir": DATA_DIR,
    }
    found_irs = IRFoldBase.find_irs(**find_irs_kwargs)
    n_irs = len(found_irs)

    # Find all compatible IR triplets
    print("Compatible IR Triplets".center(60, "="))

    ir_triplets_not_matching_same_bases = (
        get_all_ir_triplets_not_matching_same_bases_valid_gap_sz(found_irs)
    )

    for i, (ir_i, ir_j, ir_k) in enumerate(ir_triplets_not_matching_same_bases):
        (
            assumption_held,
            contains_invalid_loop_pair,
            ir_a_fe,
            ir_b_fe,
            ir_c_fe,
            fe_sum,
            fe_union,
        ) = eval_ir_triplet_structure_and_mfe(
            i, ir_i, ir_j, ir_k, sequence_len, print_out=True
        )
        experiment_results["seq"].append(seq)
        experiment_results["seq_len"].append(sequence_len)
        experiment_results["ir_triplet_index"].append(i)
        experiment_results["ir_a"].append(ir_i)
        experiment_results["ir_b"].append(ir_j)
        experiment_results["ir_c"].append(ir_j)
        experiment_results["ir_pairs"].append(
            list(itertools.combinations([ir_i, ir_j, ir_k], 2))
        )
        experiment_results["ir_a_fe"].append(ir_a_fe)
        experiment_results["ir_b_fe"].append(ir_b_fe)
        experiment_results["ir_c_fe"].append(ir_c_fe)
        experiment_results["ir_triplet_fe_sum"].append(fe_sum)
        experiment_results["ir_triplet_fe_union"].append(fe_union)
        experiment_results["fe_additivity_assumption_held"].append(assumption_held)
        experiment_results["ir_triplet_contains_invalid_loop_forming_pair"].append(
            contains_invalid_loop_pair
        )

    print("=" * 60)
    print(f"IUPACpal Parameters:")
    for kw, kwarg in find_irs_kwargs.items():
        if kw.startswith("m"):
            print(f"  - {kw} = {kwarg}")
    print(f"Num. IRs found                       : {n_irs}")
    print(f"Num. unique IR triplets              : {comb(n_irs, 3)}")
    print(
        f"Num. triplets not matching same bases: {len(ir_triplets_not_matching_same_bases)}"
    )
    print(
        f'Num triplets containing invalid loop forming pairs: {len([val for val in experiment_results["ir_triplet_contains_invalid_loop_forming_pair"] if val is True])}'
    )

    # Save results
    df = pd.DataFrame(experiment_results)
    df.to_csv(
        f"{DATA_DIR}/experiment_2_results_ir_pair_validation_validating_triplets.csv"
    )
