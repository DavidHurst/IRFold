import random
import sys
import itertools
import pandas as pd
from math import comb

from pathlib import Path

from helper_functions import (
    eval_ir_quadruplet_structure_and_mfe,
    get_all_ir_quadruplets_not_matching_same_bases_valid_gap_sz,
)

sys.path.append(str(Path(__file__).resolve().parents[1]))

from irfold import IRFoldBase


DATA_DIR = (Path(__file__).parent.parent / "data").resolve()
EXPERIMENT_2_DATA_DIR = (DATA_DIR / "experiment_2").resolve()


if __name__ == "__main__":
    """Evaluate the IR pairs found in a single sequence"""

    random.seed(223332)
    experiment_results = {
        "seq": [],
        "seq_len": [],
        "ir_quadruplet_index": [],
        "ir_a": [],
        "ir_b": [],
        "ir_c": [],
        "ir_d": [],
        "ir_pairs": [],
        "ir_a_fe": [],
        "ir_b_fe": [],
        "ir_c_fe": [],
        "ir_d_fe": [],
        "ir_quadruplet_fe_sum": [],
        "ir_quadruplet_fe_union": [],
        "fe_additivity_assumption_held": [],
        "num_invalid_loop_forming_pairs_in_quadruplet": [],
    }

    sequence_len = 30
    seq = "UGAUGACAAAUGCUUAACCCAAGCACGGCA"  # "".join(random.choice("ACGU") for _ in range(sequence_len))

    print(f"Seq. length: {sequence_len}")
    print(f"Seq.       : {seq}")

    found_irs = IRFoldBase.find_irs(seq, out_dir=str(EXPERIMENT_2_DATA_DIR))
    n_irs = len(found_irs)

    # Find all compatible IR quadruplets
    print("Compatible IR Quadruplets".center(60, "="))

    ir_quadruplets_not_matching_same_bases = (
        get_all_ir_quadruplets_not_matching_same_bases_valid_gap_sz(found_irs)
    )

    for i, (ir_i, ir_j, ir_k, ir_l) in enumerate(
        ir_quadruplets_not_matching_same_bases
    ):
        (
            assumption_held,
            num_invalid_loop_forming_pairs,
            ir_a_fe,
            ir_b_fe,
            ir_c_fe,
            ir_d_fe,
            fe_sum,
            fe_union,
        ) = eval_ir_quadruplet_structure_and_mfe(
            i, ir_i, ir_j, ir_k, ir_l, sequence_len, seq, DATA_DIR, print_out=True
        )
        experiment_results["seq"].append(seq)
        experiment_results["seq_len"].append(sequence_len)
        experiment_results["ir_quadruplet_index"].append(i)
        experiment_results["ir_a"].append(ir_i)
        experiment_results["ir_b"].append(ir_j)
        experiment_results["ir_c"].append(ir_k)
        experiment_results["ir_d"].append(ir_l)
        experiment_results["ir_pairs"].append(
            list(itertools.combinations([ir_i, ir_j, ir_k, ir_l], 2))
        )
        experiment_results["ir_a_fe"].append(ir_a_fe)
        experiment_results["ir_b_fe"].append(ir_b_fe)
        experiment_results["ir_c_fe"].append(ir_c_fe)
        experiment_results["ir_d_fe"].append(ir_d_fe)
        experiment_results["ir_quadruplet_fe_sum"].append(fe_sum)
        experiment_results["ir_quadruplet_fe_union"].append(fe_union)
        experiment_results["fe_additivity_assumption_held"].append(assumption_held)
        experiment_results["num_invalid_loop_forming_pairs_in_quadruplet"].append(
            num_invalid_loop_forming_pairs
        )

    print("Summary Stats.".center(60, "="))
    print(f"Num. IRs found                                       : {n_irs}")
    print(f"Num. unique IR quadruplets                           : {comb(n_irs, 4)}")
    print(
        f"Num. quadruplets not matching same bases             : {len(ir_quadruplets_not_matching_same_bases)}"
    )
    print(
        f'Num quadruplets containing invalid loop forming pairs: {len([val for val in experiment_results["num_invalid_loop_forming_pairs_in_quadruplet"] if val > 0])}'
    )

    # Save results
    df = pd.DataFrame(experiment_results)
    df.to_csv(f"{EXPERIMENT_2_DATA_DIR}/results_quadruplets.csv")
