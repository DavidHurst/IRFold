import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

DATA_DIR = str(Path(__file__).parent.parent / "data")

plt.rcParams["figure.figsize"] = (6, 4)


def plot_pair_analysis(pairs_df, pair_type):
    pairs_df = pairs_df.assign(
        assumption_error=(pairs_df.ir_pair_fe_sum - pairs_df.ir_pair_fe_union).abs()
    )
    pairs_df = pairs_df.sort_values("assumption_error")

    bar_plot_ax = pairs_df.plot(
        x="ir_pair_index",
        y="assumption_error",
        kind="bar",
        title=f"MFE {pair_type} Pairs \n Difference From Added MFE to Pair's Union MFE",
        xticks=[],
        xlabel="IR Pair",
        ylabel="Sorted, Absolute FE Difference",
        legend=False,
    )
    for container in bar_plot_ax.containers:
        bar_plot_ax.bar_label(
            container,
            labels=[
                int(val) if val <= 0 else "" for val in pairs_df["assumption_error"]
            ],
            padding=8,
        )
    plt.savefig(f"{DATA_DIR}/experiment_1_{pair_type}.png")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    results_df = pd.read_csv(f"{DATA_DIR}/experiment_1_results.csv")

    # Show how the additivity assumption holds for pairs that form different structures e.g. wholly nested, disjoint...
    pairs_dfs = [
        results_df[results_df.ir_pair_disjoint == 1],
        results_df[results_df.ir_pair_wholly_nested == 1],
        results_df[results_df.ir_pair_partially_nested == 1],
    ]
    pair_type_names = ["Disjoint", "Wholly Nested", "Partially Nested"]
    for pairs_df, pair_name in zip(pairs_dfs, pair_type_names):
        plot_pair_analysis(pairs_df, pair_name)
