import random
import sys
import itertools
import pandas as pd
import matplotlib.pyplot as plt
from math import comb

from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from ir_fold import IRFold0

DATA_DIR = str(Path(__file__).parent.parent / "data")


def irs_wholly_nested(outer_ir, nested_ir):
    # ToDo: Make this return true for nesting of A in B or B in A
    return outer_ir[0][0] < nested_ir[0][0] and nested_ir[1][1] < outer_ir[1][0]


def irs_disjoint(ir_a, ir_b):
    # ir_a comes entirely before ir_b
    ir_a_strictly_before_ir_b = ir_a[1][1] < ir_b[0][0]

    # ir_b comes entirely before ir_a
    ir_b_strictly_before_ir_a = ir_b[1][1] < ir_a[0][0]

    return ir_a_strictly_before_ir_b or ir_b_strictly_before_ir_a


def irs_fe_union_fe_sum(ir_a, ir_b, sequence_len, out_dir):
    # FE(IR_a) + FE(IR_b)
    ir_a_db_repr = IRFold0.irs_to_dot_bracket([ir_a], sequence_len)
    ir_b_db_repr = IRFold0.irs_to_dot_bracket([ir_b], sequence_len)

    ir_a_fe = round(IRFold0.calc_free_energy(ir_a_db_repr, seq, out_dir), 4)
    ir_b_fe = round(IRFold0.calc_free_energy(ir_b_db_repr, seq, out_dir), 4)
    fe_sum = round(ir_a_fe + ir_b_fe, 4)

    # FE(IR_a U IR_b)
    ir_a_b_db_repr = IRFold0.irs_to_dot_bracket([ir_a, ir_b], sequence_len)
    fe_union = round(IRFold0.calc_free_energy(ir_a_b_db_repr, seq, out_dir), 4)

    return (ir_a_fe, ir_b_fe), fe_union, fe_sum


def irs_form_valid_loop(ir_a, ir_b):
    """This won't scale up"""

    latest_left_string_base_idx = ir_a[0][1] if ir_a[0][1] > ir_b[0][1] else ir_b[0][1]
    earliest_right_string_base_idx = (
        ir_a[1][0] if ir_a[1][0] < ir_b[1][0] else ir_b[1][0]
    )

    # Case 1: Disjoint
    if irs_disjoint(ir_a, ir_b):
        return True

    # Case 2: Invalid hairpin
    # ToDo: Test for if the bases between being paired or not makes a difference
    bases_inbetween_parens = (
        earliest_right_string_base_idx - latest_left_string_base_idx - 1
    )
    if bases_inbetween_parens < 3:
        return False

    return True


def get_all_ir_pairs_not_matching_same_bases_valid_gap_sz(all_irs):
    valid_irs = [ir for ir in all_irs if ir[1][0] - ir[0][1] - 1 >= 3]

    unique_idx_pairs = list(
        itertools.combinations([i for i in range(len(valid_irs))], 2)
    )
    unique_ir_pairs = [(valid_irs[i], valid_irs[j]) for i, j in unique_idx_pairs]

    return [
        ir_pair
        for ir_pair in unique_ir_pairs
        if not IRFold0.ir_pair_match_same_bases(ir_pair[0], ir_pair[1])
    ]


def id_base_pairs(db_repr, irs):
    """Assigns the index value of IRs in the list provided to each base that
    IRs pair respectively"""
    ids = ["A", "B", "C", "D", "E"]
    db_repr_id_assigned = list(db_repr)

    for i, ir in enumerate(irs):
        n_base_pairs: int = ir[0][1] - ir[0][0] + 1

        lhs, rhs = ir[0], ir[1]

        db_repr_id_assigned[lhs[0] : lhs[1] + 1] = [ids[i] for _ in range(n_base_pairs)]
        db_repr_id_assigned[rhs[0] : rhs[1] + 1] = [ids[i] for _ in range(n_base_pairs)]

    return "".join(db_repr_id_assigned)


def eval_ir_pair_structure_and_mfe(ir_a, ir_b, seq_len, print_out=True):
    (ir_a_fe, ir_b_fe), fe_union, fe_sum = irs_fe_union_fe_sum(
        ir_a, ir_b, seq_len, DATA_DIR
    )

    additivity_assumption_held = fe_sum == fe_union
    valid_loop_formed = irs_form_valid_loop(ir_a, ir_b)

    db_repr = IRFold0.irs_to_dot_bracket([ir_a, ir_b], seq_len)
    db_repr_id_assigned = id_base_pairs(db_repr, [ir_a, ir_b])

    if print_out:
        ir_a_fmt = f"[{str(ir_a[0][0]).zfill(2)}->{str(ir_a[0][1]).zfill(2)} {str(ir_a[1][0]).zfill(2)}->{str(ir_a[1][1]).zfill(2)}]"
        ir_b_fmt = f"[{str(ir_b[0][0]).zfill(2)}->{str(ir_b[0][1]).zfill(2)} {str(ir_b[1][0]).zfill(2)}->{str(ir_b[1][1]).zfill(2)}]"

        print_space = " " * 4
        print(f"A = {ir_a_fmt}, B = {ir_b_fmt}")
        print(print_space, f"A    : {IRFold0.irs_to_dot_bracket([ir_a], seq_len)}")
        print(print_space, f"B    : {IRFold0.irs_to_dot_bracket([ir_b], seq_len)}")
        print(print_space, f"A U B: {db_repr}")
        print(print_space, f"A U B: {db_repr_id_assigned}")
        print()
        print(print_space, f"FE(A)        : {ir_a_fe:.4f}")
        print(print_space, f"FE(B)        : {ir_b_fe:.4f}")
        print(print_space, f"FE(A) + FE(B): {fe_sum:.4f}")
        print(print_space, f"FE(A U B)    : {fe_union:.4f}")
        print()
        print(print_space, f"Valid loop      : {valid_loop_formed}")
        print(print_space, f"Assumption holds: {additivity_assumption_held}")

    return (
        additivity_assumption_held,
        valid_loop_formed,
        ir_a,
        ir_b,
        ir_a_fe,
        ir_b_fe,
        fe_sum,
        fe_union,
    )


if __name__ == "__main__":
    # random.seed(1232)
    experiment_results = {
        "seq": [],
        "seq_len": [],
        "num_valid_ir_gap_valid_base_pairing_pairs": [],
        "ir_pair_index": [],
        "ir_a": [],
        "ir_b": [],
        "ir_a_fe": [],
        "ir_b_fe": [],
        "ir_pair_fe_sum": [],
        "ir_pair_fe_union": [],
        "irs_wholly_nested": [],
        "irs_partially_nested": [],
        "irs_disjoint": [],
        "additivity_assumption_held": [],
        "irs_form_valid_loop": [],
    }

    for _ in range(2):
        sequence_len = 35
        seq = "".join(random.choice("ACGU") for _ in range(sequence_len))

        print(f"Seq. length: {sequence_len}")
        print(f"Seq.       : {seq}")

        summary_results = {"Wholly Nested": [], "Partially Nested": [], "Disjoint": []}

        find_irs_kwargs = {
            "sequence": seq,
            "min_len": 2,
            "max_len": sequence_len,
            "max_gap": sequence_len - 1,
            "mismatches": 0,
            "out_dir": DATA_DIR,
        }
        found_irs = IRFold0.find_irs(**find_irs_kwargs)
        n_irs = len(found_irs)

        # Find all compatible IR pairs
        print("Compatible IR Pairs".center(60, "="))
        # ToDo: Also evaluate triples, quadruples etc.

        separate_base_matching_ir_pairs = (
            get_all_ir_pairs_not_matching_same_bases_valid_gap_sz(found_irs)
        )

        # Classify IR pairs into three groups: wholly nested, partially nested or disjoint
        wholly_nested_ir_pairs = [
            pair
            for pair in separate_base_matching_ir_pairs
            if irs_wholly_nested(pair[0], pair[1]) or irs_wholly_nested(pair[0], pair[1])
        ]
        disjoint_ir_pairs = [
            pair
            for pair in separate_base_matching_ir_pairs
            if irs_disjoint(pair[0], pair[1])
        ]
        partially_nested_ir_pairs = [
            pair
            for pair in separate_base_matching_ir_pairs
            if pair not in wholly_nested_ir_pairs and pair not in disjoint_ir_pairs
        ]

        print(f"IUPACpal Parameters:")
        for kw, kwarg in find_irs_kwargs.items():
            if kw.startswith("m"):
                print(f"  - {kw} = {kwarg}")
        print(f"IRs found             : {n_irs}")
        print(f"Unique pairs          : {comb(n_irs, 2)}")
        print(f"Not over-pairing bases: {len(separate_base_matching_ir_pairs)}")
        print(f"Wholly nested         : {len(wholly_nested_ir_pairs)}")
        print(f"Partially nested      : {len(partially_nested_ir_pairs)}")
        print(f"Disjoint              : {len(disjoint_ir_pairs)}")

        pair_index = 0
        for pair_type_bin_str, pair_type, pair_list in zip(
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            ["Wholly Nested", "Partially Nested", "Disjoint"],
            [wholly_nested_ir_pairs, partially_nested_ir_pairs, disjoint_ir_pairs],
        ):
            print(f"{pair_type} Pairs".center(60, "="))

            for (ir_i, ir_j) in pair_list:
                (
                    assumption_held,
                    valid_loop_formed,
                    ir_a,
                    ir_b,
                    ir_a_fe,
                    ir_b_fe,
                    fe_sum,
                    fe_union,
                ) = eval_ir_pair_structure_and_mfe(ir_i, ir_j, sequence_len)
                summary_results[pair_type].append((assumption_held, valid_loop_formed))
                experiment_results["seq"].append(seq)
                experiment_results["seq_len"].append(sequence_len)
                experiment_results["num_valid_ir_gap_valid_base_pairing_pairs"].append(len(separate_base_matching_ir_pairs))
                experiment_results['ir_pair_index'].append(pair_index)
                experiment_results["ir_a"].append(ir_i)
                experiment_results["ir_b"].append(ir_j)
                experiment_results["ir_a_fe"].append(ir_a_fe)
                experiment_results["ir_b_fe"].append(ir_b_fe)
                experiment_results["ir_pair_fe_sum"].append(fe_sum)
                experiment_results["ir_pair_fe_union"].append(fe_union)
                experiment_results["irs_wholly_nested"].append(pair_type_bin_str[0])
                experiment_results["irs_partially_nested"].append(pair_type_bin_str[1])
                experiment_results["irs_disjoint"].append(pair_type_bin_str[2])
                experiment_results["additivity_assumption_held"].append(assumption_held)
                experiment_results["irs_form_valid_loop"].append(valid_loop_formed)

                pair_index += 1

        print("Assumption Holding Summative Stats.".center(60, "="))
        for pair_type, r in summary_results.items():
            num_assumption_held_valid_loop = len([p for p in r if p[0] and p[1]])
            num_assumption_held_invalid_loop = len([p for p in r if p[0] and not p[1]])
            num_pairs = len(r)
            print(f"{pair_type}:")
            print(f"   Valid Loop  : {num_assumption_held_valid_loop}/{num_pairs} times")
            print(f"   Invalid Loop: {num_assumption_held_invalid_loop}/{num_pairs} times")
            print(
                f"   Overall Acc.: {num_assumption_held_valid_loop + num_assumption_held_invalid_loop}/{num_pairs} times"
            )

    # Plot results
    df = pd.DataFrame(experiment_results)
    df.to_csv(f'{DATA_DIR}/additivity_experiment_results.csv')

    # Plot disjoint sum fe vs union fe
    # disjoint_pairs_df = df[df.irs_disjoint == 1]
    # plt.scatter(disjoint_pairs_df.ir_pair_index, disjoint_pairs_df.ir_pair_fe_sum, label='$\sum FE(IR_i) + FE(IR_j)$', alpha=0.5)
    # plt.scatter(disjoint_pairs_df.ir_pair_index, disjoint_pairs_df.ir_pair_fe_union, label='$FE(IR_i \cup IR_j)$', alpha=0.5, marker='x')
    # plt.title('MFE Disjoint IR Pairs')
    # plt.xlabel('IR Pair')
    # plt.ylabel('MFE')
    # plt.legend()
    # plt.show()

    # Plot wholly nested sum fe vs union fe
    nested_df = df[df.irs_wholly_nested == 1]
    nested_df.fe_union_zeroed = nested_df.ir_pair_fe_union - nested_df.ir_pair_fe_union
    nested_df.fe_sum_corrected = nested_df.ir_pair_fe_sum - nested_df.ir_pair_fe_union
    plt.scatter(nested_df.ir_pair_index, nested_df.fe_union_zeroed, label='MFE GT Zeroed',
                alpha=0.5, marker='^')
    plt.scatter(nested_df.ir_pair_index, nested_df.fe_sum_corrected, label='Our MFE Corrected',
                alpha=0.5, marker='+')

    # Draw arrow showing direction of incorrectness from assumption value
    for arrow_base_x, arrow_base_y, arrow_tip_x, arrow_tip_y in zip(nested_df.ir_pair_index, nested_df.fe_union_zeroed, nested_df.ir_pair_index, nested_df.fe_sum_corrected):
        arrow_dx = arrow_base_x - arrow_tip_x
        arrow_dy = arrow_base_y - arrow_tip_y

        plt.arrow(arrow_base_x, arrow_base_y, -arrow_dx, -arrow_dy, head_width=0.5, head_length=0.5, length_includes_head=True)


    plt.title('MFE Wholly Nested IR Pairs')
    plt.xlabel('IR Pair')
    plt.ylabel('MFE')
    plt.legend()
    plt.show()


