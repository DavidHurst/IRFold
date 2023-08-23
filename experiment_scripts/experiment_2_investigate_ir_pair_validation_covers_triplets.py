import sys
import itertools
import pandas as pd
from math import comb

from pathlib import Path

from helper_functions import (
    get_all_ir_triplets_not_matching_same_bases_valid_gap_sz,
    eval_ir_triplet_structure_and_mfe,
)

sys.path.append(str(Path(__file__).resolve().parents[1]))

from irfold import IRFoldBase

DATA_DIR = (Path(__file__).parent.parent / "data").resolve()
EXPERIMENT_2_DATA_DIR = (DATA_DIR / "experiment_2").resolve()


if __name__ == "__main__":
    """Evaluate the IR pairs found in a single sequence"""

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
        "num_invalid_loop_forming_pairs_in_triplet": [],
    }

    sequence_len = 30
    seq = "UGAUGACAAAUGCUUAACCCAAGCACGGCA"

    print(f"Seq. length: {sequence_len}")
    print(f"Seq.       : {seq}")

    found_irs = IRFoldBase._find_irs(seq, out_dir=str(EXPERIMENT_2_DATA_DIR))
    n_irs = len(found_irs)

    # Find all compatible IR triplets
    print("Compatible IR Triplets".center(60, "="))

    ir_triplets_not_matching_same_bases = (
        get_all_ir_triplets_not_matching_same_bases_valid_gap_sz(found_irs)
    )

    for i, (ir_i, ir_j, ir_k) in enumerate(ir_triplets_not_matching_same_bases):
        (
            assumption_held,
            num_invalid_loop_forming_pairs,
            ir_a_fe,
            ir_b_fe,
            ir_c_fe,
            fe_sum,
            fe_union,
        ) = eval_ir_triplet_structure_and_mfe(
            i, ir_i, ir_j, ir_k, sequence_len, seq, DATA_DIR, print_out=True
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
        experiment_results["num_invalid_loop_forming_pairs_in_triplet"].append(
            num_invalid_loop_forming_pairs
        )

    print("=" * 60)
    print(f"Num. IRs found                       : {n_irs}")
    print(f"Num. unique IR triplets              : {comb(n_irs, 3)}")
    print(
        f"Num. triplets not matching same bases: {len(ir_triplets_not_matching_same_bases)}"
    )
    print(
        f'Num triplets containing invalid loop forming pairs: {len([val for val in experiment_results["num_invalid_loop_forming_pairs_in_triplet"] if val > 0])}'
    )

    # Save results
    df = pd.DataFrame(experiment_results)
    df.to_csv(f"{EXPERIMENT_2_DATA_DIR}/results_triplets.csv")
