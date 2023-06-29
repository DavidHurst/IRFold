import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

DATA_DIR = str(Path(__file__).parent.parent / "data")

plt.rcParams["figure.figsize"] = (6, 4)

if __name__ == "__main__":
    results_df = pd.read_csv(f"{DATA_DIR}/experiment_3/results.csv")

    valid_loop_triplets_df = results_df[results_df["triplet_union_fe"] < 90_000]
    valid_loop_triplets_df = valid_loop_triplets_df.reindex(
        [i for i in range(len(valid_loop_triplets_df))]
    )

    # Standardise values so that 0 is the true free energy (triplet_union_fe)
    for col in [
        "corrected_additive_fe_first_ir_pair",
        "corrected_additive_fe_second_ir_pair",
        "corrected_additive_fe_third_ir_pair",
        "corrected_additive_fe_all_pairs",
    ]:
        valid_loop_triplets_df[col] = (
            valid_loop_triplets_df[col] - valid_loop_triplets_df.triplet_union_fe
        )
    valid_loop_triplets_df.triplet_union_fe = (
        valid_loop_triplets_df.triplet_union_fe
        - valid_loop_triplets_df.triplet_union_fe
    )

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    # Plot triplets additive free energy with first IR pair corrected
    # ax1.scatter(valid_loop_triplets_df.index, valid_loop_triplets_df.triplet_union_fe)
    ax1.bar(valid_loop_triplets_df.index,
        valid_loop_triplets_df.index,
        valid_loop_triplets_df.corrected_additive_fe_first_ir_pair,
    )
    ax1.set_title("Triplets' Free Energy w/ First IR Pair's Free Energy Corrected")
    ax1.set_xlabel("IR Triplet")
    ax1.set_ylabel("Diff. from Triplet\'s True Free Energy")

    # Plot triplets additive free energy with second IR pair corrected
    ax2.bar(
        valid_loop_triplets_df.index,
        valid_loop_triplets_df.corrected_additive_fe_second_ir_pair,
    )
    ax2.set_title("Triplets' Free Energy w/ Second IR Pair's Free Energy Corrected")
    ax2.set_xlabel("IR Triplet")
    ax2.set_ylabel("Diff. from Triplet\'s True Free Energy")

    # Plot triplets additive free energy with third IR pair corrected
    ax3.bar(
        valid_loop_triplets_df.index,
        valid_loop_triplets_df.corrected_additive_fe_third_ir_pair,
    )
    ax3.set_title("Triplets' Free Energy w/ Third IR Pair's Free Energy Corrected")
    ax3.set_xlabel("IR Triplet")
    ax3.set_ylabel("Diff. from Triplet\'s True Free Energy")

    # Plot triplets additive free energy with all pairs corrected
    ax4.bar(
        valid_loop_triplets_df.index,
        valid_loop_triplets_df.corrected_additive_fe_all_pairs,
    )
    ax4.set_title("Triplets' Free Energy w/ All IR Pairs' Free Energies Corrected")
    ax4.set_xlabel("IR Triplet")
    ax4.set_ylabel("Diff. from Triplet\'s True Free Energy")

    plt.tight_layout()
    plt.legend()
    plt.show()

    invalid_loop_triplets_df = results_df[
        results_df["triplet_union_fe"] > 90_000
    ].reindex()
