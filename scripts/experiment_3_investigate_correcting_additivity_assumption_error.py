import itertools
import random
import sys
from pathlib import Path

import pandas as pd

from helper_functions import (
    get_all_ir_pairs_not_matching_same_bases_valid_gap_sz,
    ir_pair_free_energy_calculation_variations,
    get_all_ir_triplets_not_matching_same_bases_valid_gap_sz,
    eval_ir_triplet_structure_and_mfe,
)

sys.path.append(str(Path(__file__).resolve().parents[1]))

from irfold import IRFoldBase


DATA_DIR = str(Path(__file__).parent.parent / "data")


if __name__ == "__main__":
    # Get all valid gap IR pairs which don't match same bases
    # random.seed(223332)

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

    print("Compatible IR Pairs FEs and Correction Variables".center(60, "="))

    ir_pairs_not_matching_same_bases = (
        get_all_ir_pairs_not_matching_same_bases_valid_gap_sz(found_irs)
    )

    # Calculate variable for each IR pair which corrects additivity assumption for that pair
    ir_pair_correction_variables = {}
    for i, (ir_a, ir_b) in enumerate(ir_pairs_not_matching_same_bases):
        (
            ir_a_fe,
            ir_b_fe,
            pairs_fe_union,
            pairs_fe_sum,
        ) = ir_pair_free_energy_calculation_variations(
            ir_a, ir_b, sequence_len, seq, DATA_DIR
        )

        free_energy_difference = pairs_fe_sum - pairs_fe_union

        # Calculate required correction variable's value
        correction_var = 0
        if abs(free_energy_difference) > 0:
            if pairs_fe_sum > pairs_fe_union:
                correction_var = -abs(free_energy_difference)
            else:
                correction_var = abs(free_energy_difference)
        else:
            correction_var = 0
        ir_pair_correction_variables[str(ir_a) + str(ir_b)] = correction_var

        ir_a_fmt = f"[{str(ir_a[0][0]).zfill(2)}->{str(ir_a[0][1]).zfill(2)} {str(ir_a[1][0]).zfill(2)}->{str(ir_a[1][1]).zfill(2)}]"
        ir_b_fmt = f"[{str(ir_b[0][0]).zfill(2)}->{str(ir_b[0][1]).zfill(2)} {str(ir_b[1][0]).zfill(2)}->{str(ir_b[1][1]).zfill(2)}]"

        # print_space = " " * 4
        # print(f"Pair #{i}: A = {ir_a_fmt}, B = {ir_b_fmt}")
        # print(print_space, f"FE(A)        : {ir_a_fe:.4f}")
        # print(print_space, f"FE(B)        : {ir_b_fe:.4f}")
        # print(print_space, f"FE(A) + FE(B): {pairs_fe_sum:.4f}")
        # print(print_space, f"FE(A U B)    : {pairs_fe_union:.4f}")
        # print(print_space, f"Correction   : {correction_var:.4f}")

    print("IR Triplets Corrected FEs".center(60, "="))

    experiment_results = {
        "seq": [],
        "seq_len": [],
        "ir_triplet_index": [],
        "triplet_additive_fe": [],
        "triplet_union_fe": [],
        "corrected_additive_fe_first_ir_pair": [],
        "corrected_additive_fe_second_ir_pair": [],
        "corrected_additive_fe_third_ir_pair": [],
        "corrected_additive_fe_all_pairs": [],
    }

    # Get all valid gap IR triplets which don't match same bases
    ir_triplets_not_matching_same_bases = (
        get_all_ir_triplets_not_matching_same_bases_valid_gap_sz(found_irs)
    )

    # Calculate and compare free energies of IR triplets with IR pair correction variables active and inactive
    for i, (ir_i, ir_j, ir_k) in enumerate(ir_triplets_not_matching_same_bases):
        (_, _, _, _, _, fe_sum, fe_union,) = eval_ir_triplet_structure_and_mfe(
            i, ir_i, ir_j, ir_k, sequence_len, seq, DATA_DIR, print_out=False
        )

        ir_a_fmt = (
            f"[{str(ir_i[0][0]).zfill(2)}->{str(ir_i[0][1]).zfill(2)} "
            f"{str(ir_i[1][0]).zfill(2)}->{str(ir_i[1][1]).zfill(2)}]"
        )
        ir_b_fmt = (
            f"[{str(ir_j[0][0]).zfill(2)}->{str(ir_j[0][1]).zfill(2)} "
            f"{str(ir_j[1][0]).zfill(2)}->{str(ir_j[1][1]).zfill(2)}]"
        )
        ir_c_fmt = (
            f"[{str(ir_k[0][0]).zfill(2)}->{str(ir_k[0][1]).zfill(2)} "
            f"{str(ir_k[1][0]).zfill(2)}->{str(ir_k[1][1]).zfill(2)}]"
        )

        correction_variable_values = []
        for j, ir_pair in enumerate(itertools.combinations([ir_i, ir_j, ir_k], 2)):
            correction_variable_values.append(
                ir_pair_correction_variables[str(ir_pair[0]) + str(ir_pair[1])]
            )

        print_space = " " * 4
        print(f"Triplet #{i}: A = {ir_a_fmt}, B = {ir_b_fmt}, C = {ir_c_fmt}")
        print(print_space, f"Uncorrected FE          : {fe_sum:.4f}")
        print(print_space, f"Expected FE             : {fe_union:.4f}")
        for j, val in enumerate(correction_variable_values):
            print(print_space, f'Corrected w/ - Pair {j}   : {fe_sum + val:.4f}')
        print(print_space, f'Corrected w/ - All pairs: {fe_sum + sum(correction_variable_values):.4f}')

        experiment_results["seq"].append(seq)
        experiment_results["seq_len"].append(sequence_len)
        experiment_results["ir_triplet_index"].append(i)
        experiment_results["triplet_additive_fe"].append(fe_sum)
        experiment_results["triplet_union_fe"].append(fe_union)
        experiment_results["corrected_additive_fe_first_ir_pair"].append(fe_sum + correction_variable_values[0])
        experiment_results["corrected_additive_fe_second_ir_pair"].append(fe_sum + correction_variable_values[1])
        experiment_results["corrected_additive_fe_third_ir_pair"].append(fe_sum + correction_variable_values[2])
        experiment_results["corrected_additive_fe_all_pairs"].append(fe_sum + sum(correction_variable_values))

    df = pd.DataFrame(experiment_results)
    df.to_csv(f"{DATA_DIR}/experiment_3/experiment_3_results.csv")