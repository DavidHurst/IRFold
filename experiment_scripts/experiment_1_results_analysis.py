import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

DATA_DIR = (Path(__file__).parent.parent / "data").resolve()
EXPERIMENT_1_DATA_DIR = (DATA_DIR / 'experiment_1').resolve()

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
        title=f"MFE {pair_type.replace('_', ' ')} Pairs \n Difference From Added MFE to Pair's Union MFE",
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
    plt.tight_layout()
    plt.savefig(
        f"{EXPERIMENT_1_DATA_DIR}/ir_pairs_fe_sorted_abs_diff_{pair_type.lower()}.png"
    )
    plt.show()


if __name__ == "__main__":
    results_df = pd.read_csv(f"{EXPERIMENT_1_DATA_DIR}/results.csv")

    # Show free energy value of pairs, highlight pairs that form valid loops and those that don't
    form_valid_loop_highlighted_df = results_df.assign(
        colour=[
            "green" if val else "red" for val in results_df.ir_pair_forms_valid_loop
        ]
    )
    form_valid_loop_highlighted_df.ir_pair_fe_union = (
        form_valid_loop_highlighted_df.ir_pair_fe_union.apply(
            lambda x: x if x < 50_000 else 100
        )
    )
    bar_plot_ax = form_valid_loop_highlighted_df.plot.bar(
        x="ir_pair_index",
        y="ir_pair_fe_union",
        color=form_valid_loop_highlighted_df.colour,
        title=f"Free Energy Values of IR Pairs",
        xticks=[],
        xlabel="IR Pair",
        ylabel="Pair Union Free Energy",
        legend=False,  # ToDo: Add legend
    )
    plt.tight_layout()
    plt.savefig(
        f"{EXPERIMENT_1_DATA_DIR}/ir_pairs_forming_valid_loops_vs_free_energy.png"
    )
    plt.show()

    # Show how the additivity assumption holds for pairs that form different structures e.g. wholly nested, disjoint...
    pairs_dfs = [
        results_df[results_df.ir_pair_disjoint == 1],
        results_df[results_df.ir_pair_wholly_nested == 1],
        results_df[results_df.ir_pair_partially_nested == 1],
    ]
    pair_type_names = ["Disjoint", "Wholly_Nested", "Partially_Nested"]
    for pairs_df, pair_name in zip(pairs_dfs, pair_type_names):
        plot_pair_analysis(pairs_df, pair_name)
