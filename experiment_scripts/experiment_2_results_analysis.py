import pandas as pd

from pathlib import Path

DATA_DIR = (Path(__file__).parent.parent / "data").resolve()
EXPERIMENT_2_DATA_DIR = (DATA_DIR / "experiment_2").resolve()


if __name__ == "__main__":
    trips_res_df = pd.read_csv(f"{EXPERIMENT_2_DATA_DIR}/results_triplets.csv")
    quads_res_df = pd.read_csv(f"{EXPERIMENT_2_DATA_DIR}/results_quadruplets.csv")

    # Show that only IR triplets containing invalid loop forming IR pairs are themselves invalid
    print("IR Triplets".center(70, "="))

    valid_triplets = trips_res_df[trips_res_df.ir_triplet_fe_union < 90_000]
    valid_ir_pair_count_invalid_triplets = (
        valid_triplets.num_invalid_loop_forming_pairs_in_triplet.sum()
    )
    print(
        f"Type".center(20, " ")
        + "|"
        + f"Count".center(20, " ")
        + "|"
        + f"Invalid IR Pair Count".center(30, " ")
    )
    print("-" * 70)
    print(
        f"Valid Triplets".ljust(20, " ")
        + "|"
        + f"{len(valid_triplets)}".center(20, " ")
        + "|"
        + f"{valid_ir_pair_count_invalid_triplets}".center(30, " ")
    )

    invalid_triplets = trips_res_df[trips_res_df.ir_triplet_fe_union >= 90_000]
    invalid_ir_pair_count_invalid_triplets = (
        invalid_triplets.num_invalid_loop_forming_pairs_in_triplet.sum()
    )
    print(
        f"Invalid Triplets".ljust(20, " ")
        + "|"
        + f"{len(invalid_triplets)}".center(20, " ")
        + "|"
        + f"{invalid_ir_pair_count_invalid_triplets}".center(30, " ")
    )

    print()

    # Show that only IR quadruplets containing invalid loop forming IR pairs are themselves invalid
    print("IR Quadruplets".center(70, "="))

    valid_quadruplets = quads_res_df[quads_res_df.ir_quadruplet_fe_union < 90_000]
    valid_ir_pair_count_valid_quadruplets = (
        valid_quadruplets.num_invalid_loop_forming_pairs_in_quadruplet.sum()
    )
    print(
        f"Type".center(20, " ")
        + "|"
        + f"Count".center(20, " ")
        + "|"
        + f"Invalid IR Pair Count".center(30, " ")
    )
    print("-" * 70)
    print(
        f"Valid Quadruplets".ljust(20, " ")
        + "|"
        + f"{len(valid_quadruplets)}".center(20, " ")
        + "|"
        + f"{valid_ir_pair_count_valid_quadruplets}".center(30, " ")
    )

    invalid_quadruplets = quads_res_df[quads_res_df.ir_quadruplet_fe_union >= 90_000]
    invalid_ir_pair_count_invalid_quadruplets = (
        invalid_quadruplets.num_invalid_loop_forming_pairs_in_quadruplet.sum()
    )
    print(
        f"Invalid Quadruplets".ljust(20, " ")
        + "|"
        + f"{len(invalid_quadruplets)}".center(20, " ")
        + "|"
        + f"{invalid_ir_pair_count_invalid_quadruplets}".center(30, " ")
    )
