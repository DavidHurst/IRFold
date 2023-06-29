import random
import sys
import pandas as pd
from math import comb

from pathlib import Path

from helper_functions import (
    get_all_ir_pairs_not_matching_same_bases_valid_gap_sz,
    eval_ir_pair_structure_and_mfe,
)

sys.path.append(str(Path(__file__).resolve().parents[1]))

from irfold import IRFoldBase
from irfold.util import (
    ir_pair_wholly_nested,
    ir_pair_disjoint,
    ir_has_valid_gap_size,
)

DATA_DIR = str(Path(__file__).parent.parent / "data")


if __name__ == "__main__":
    """Evaluate the IR pairs found in a single sequence"""

    random.seed(223332)
    experiment_results = {
        "seq": [],
        "seq_len": [],
        "ir_pair_index": [],
        "ir_a": [],
        "ir_b": [],
        "ir_a_fe": [],
        "ir_b_fe": [],
        "ir_pair_fe_sum": [],
        "ir_pair_fe_union": [],
        "fe_additivity_assumption_held": [],
        "ir_pair_forms_valid_loop": [],
        "ir_pair_wholly_nested": [],
        "ir_pair_partially_nested": [],
        "ir_pair_disjoint": [],
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

    # Find all compatible IR pairs
    print("Compatible IR Pairs".center(60, "="))

    ir_pairs_not_matching_same_bases = (
        get_all_ir_pairs_not_matching_same_bases_valid_gap_sz(found_irs)
    )
    wholly_nested_ir_pairs = []
    disjoint_ir_pairs = []
    partially_nested_ir_pairs = []

    # Evaluate IR pairs and group pairs into three groups: wholly nested, partially nested or disjoint
    for i, (ir_i, ir_j) in enumerate(ir_pairs_not_matching_same_bases):
        (
            assumption_held,
            valid_loop_formed,
            ir_a_fe,
            ir_b_fe,
            fe_sum,
            fe_union,
        ) = eval_ir_pair_structure_and_mfe(
            i, ir_i, ir_j, sequence_len, seq, DATA_DIR, print_out=True
        )
        experiment_results["seq"].append(seq)
        experiment_results["seq_len"].append(sequence_len)
        experiment_results["ir_pair_index"].append(i)
        experiment_results["ir_a"].append(ir_i)
        experiment_results["ir_b"].append(ir_j)
        experiment_results["ir_a_fe"].append(ir_a_fe)
        experiment_results["ir_b_fe"].append(ir_b_fe)
        experiment_results["ir_pair_fe_sum"].append(fe_sum)
        experiment_results["ir_pair_fe_union"].append(fe_union)
        experiment_results["fe_additivity_assumption_held"].append(assumption_held)
        experiment_results["ir_pair_forms_valid_loop"].append(valid_loop_formed)

        pair_type_indicator_binary_repr = [
            0,
            0,
            0,
        ]  # Wholly nested, partially nested, disjoint
        if ir_pair_wholly_nested(ir_i, ir_j):
            pair_type_indicator_binary_repr[0] = 1
            wholly_nested_ir_pairs.append((ir_i, ir_j))
        if not ir_pair_wholly_nested(ir_i, ir_j) and not ir_pair_disjoint(ir_i, ir_j):
            pair_type_indicator_binary_repr[1] = 1
            partially_nested_ir_pairs.append((ir_i, ir_j))
        if ir_pair_disjoint(ir_i, ir_j):
            pair_type_indicator_binary_repr[2] = 1
            disjoint_ir_pairs.append((ir_i, ir_j))

        experiment_results["ir_pair_wholly_nested"].append(
            pair_type_indicator_binary_repr[0]
        )
        experiment_results["ir_pair_partially_nested"].append(
            pair_type_indicator_binary_repr[1]
        )
        experiment_results["ir_pair_disjoint"].append(
            pair_type_indicator_binary_repr[2]
        )

    print("=" * 60)
    print(f"IUPACpal Parameters:")
    for kw, kwarg in find_irs_kwargs.items():
        if kw.startswith("m"):
            print(f"  - {kw} = {kwarg}")
    print(f"Num. IRs found                    : {n_irs}")
    print(f"Num. unique IR pairs              : {comb(n_irs, 2)}")
    print(
        f"Num. pairs not matching same bases: {len(ir_pairs_not_matching_same_bases)}"
    )
    print(f"Num. wholly nested pairs          : {len(wholly_nested_ir_pairs)}")
    print(f"Num. partially nested pairs       : {len(partially_nested_ir_pairs)}")
    print(f"Num. disjoint pairs               : {len(disjoint_ir_pairs)}")

    print("All IRs".center(60, "="))
    for i, ir in enumerate(found_irs):
        print(
            f'IR#{str(i).zfill(2)}: {str(ir).ljust(20, " ")}, Valid Gap Sz.: {ir_has_valid_gap_size(ir)}'
        )

    # Save results
    df = pd.DataFrame(experiment_results)
    df.to_csv(f"{DATA_DIR}/experiment_1/experiment_1_results.csv")
