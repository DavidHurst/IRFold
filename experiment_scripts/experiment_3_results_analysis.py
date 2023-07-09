import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

DATA_DIR = (Path(__file__).parent.parent / "data").resolve()
EXPERIMENT_3_DATA_DIR = (DATA_DIR / "experiment_3").resolve()

plt.rcParams["figure.figsize"] = (14, 12)


def plot_analysis(df, name):
    df = df.reindex([i for i in range(len(df))])

    # Standardise values so that 0 is the true free energy (triplet_union_fe)
    for col in [
        "corrected_additive_fe_first_ir_pair",
        "corrected_additive_fe_second_ir_pair",
        "corrected_additive_fe_third_ir_pair",
        "corrected_additive_fe_all_pairs",
    ]:
        df[col] = df[col] - df.triplet_union_fe
    df.triplet_union_fe = df.triplet_union_fe - df.triplet_union_fe

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    # Plot triplets additive free energy with first IR pair corrected
    ax1.bar(df.index, df.corrected_additive_fe_first_ir_pair)
    ax1.set_title("First IR Pair's Free Energy Corrected")
    ax1.set_xlabel("IR Triplet")
    ax1.set_ylabel("Free Energy - Difference from True")

    # Plot triplets additive free energy with second IR pair corrected
    ax2.bar(df.index, df.corrected_additive_fe_second_ir_pair)
    ax2.set_title("Second IR Pair's Free Energy Corrected")
    ax2.set_xlabel("IR Triplet")
    ax2.set_ylabel("Free Energy - Difference from True")

    # Plot triplets additive free energy with third IR pair corrected
    ax3.bar(df.index, df.corrected_additive_fe_third_ir_pair)
    ax3.set_title("Third IR Pair's Free Energy Corrected")
    ax3.set_xlabel("IR Triplet")
    ax3.set_ylabel("Free Energy - Difference from True")

    # Plot triplets additive free energy with all pairs corrected
    ax4.bar(df.index, df.corrected_additive_fe_all_pairs)
    ax4.set_title("All IR Pairs' Free Energies Corrected")
    ax4.set_xlabel("IR Triplet")
    ax4.set_ylabel("Free Energy - Difference from True")

    plt.tight_layout()
    plt.savefig(
        f"{EXPERIMENT_3_DATA_DIR}/ir_triplets_fe_corrections_{name.lower()}.png"
    )
    plt.show()


if __name__ == "__main__":
    results_df = pd.read_csv(f"{EXPERIMENT_3_DATA_DIR}/results.csv")

    valid_loop_triplets_df = results_df[results_df["triplet_union_fe"] < 90_000]
    invalid_loop_triplets_df = results_df[results_df["triplet_union_fe"] > 90_000]

    plot_analysis(valid_loop_triplets_df, "Valid")
    plot_analysis(invalid_loop_triplets_df, "Invalid")
