from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

DATA_DIR = Path(__file__).parent.parent / "data"
EXPERIMENT_5_DIR_PATH = DATA_DIR / "experiment_5"

if __name__ == "__main__":
    rnafold_res_df = pd.read_csv(
        EXPERIMENT_5_DIR_PATH / "RNAfold_benchmarking_results.csv"
    )
    rnastructure_res_df = pd.read_csv(
        EXPERIMENT_5_DIR_PATH / "RNAstructure_benchmarking_results.csv"
    )
    ipknot_res_df = pd.read_csv(
        EXPERIMENT_5_DIR_PATH / "IPknot_benchmarking_results.csv"
    )
    irfold_val2_res_df = pd.read_csv(
        EXPERIMENT_5_DIR_PATH / "IRFoldVal2_benchmarking_results.csv"
    )
    irfold_corx2_res_df = pd.read_csv(
        EXPERIMENT_5_DIR_PATH / "IRFoldCorX2_benchmarking_results.csv"
    )

    print(f"Num. Samples: {len(irfold_corx2_res_df)}")

    # Print average F1 scores and running times for each model
    f1_scores = []
    running_times = []
    for res_df, df_name in zip(
        [
            rnafold_res_df,
            rnastructure_res_df,
            ipknot_res_df,
            irfold_val2_res_df,
            irfold_corx2_res_df,
        ],
        ["RNAfold", "RNAstructure", "IPknot", "IRFoldVal2", "IRFoldCorX2"],
    ):

        f1_scores.append(res_df["f1"].mean())
        running_times.append(res_df["execution_wall_time_secs"].mean())

    print(
        f"Model".ljust(15) + "Mean F1".center(9) + "Mean Running Time (secs)".center(21)
    )
    print("-" * (15 + 9 + 21))
    for model, f1, running_time in zip(
        ["RNAfold", "RNAstructure", "IPknot", "IRFoldVal2", "IRFoldCorX2"],
        f1_scores,
        running_times,
    ):
        print(
            f"{model}".ljust(15)
            + f"{f1:.2f}".center(9)
            + f"{running_time:.4f}".center(21)
        )

    # Print and plot performance over sequence length intervals
    seq_interval_sz = 10
    avg_df = rnafold_res_df.groupby("sequence_length").mean(numeric_only=True)
    seq_intervals = [
        [start, start + seq_interval_sz - 1]
        for start in range(
            rnafold_res_df.sequence_length.unique().min(),
            rnafold_res_df.sequence_length.unique().max(),
            seq_interval_sz,
        )
    ]

    interval_counts = [
        len(
            ipknot_res_df[
                (ipknot_res_df.sequence_length >= start)
                & (ipknot_res_df.sequence_length <= end)
            ]
        )
        for start, end in seq_intervals
    ]
    normalised_interval_counts = [
        (cnt - min(interval_counts)) / (max(interval_counts) - min(interval_counts))
        for cnt in interval_counts
    ]

    seq_intervals[-1][-1] += 1
    col_width = 15
    n_intervals = len(seq_intervals)

    print(
        "".join(
            ["              ".ljust(col_width)]
            + [f"F1 Mean".center(col_width * n_intervals)]
        )
    )
    print("-" * (col_width + (col_width * n_intervals)))

    print(
        "".join(
            ["Model\Interval".ljust(col_width)]
            + [f"{i[0]}-{i[1]}".center(col_width) for i in seq_intervals]
        )
    )
    print("-" * (col_width * (n_intervals + 1)))

    model_seq_len_interval_f1_means = {
        "RNAfold": [],
        "RNAstructure": [],
        "IPknot": [],
        "IRFoldVal2": [],
        "IRFoldCorX2": [],
    }
    for res_df, df_name in zip(
        [
            rnafold_res_df,
            rnastructure_res_df,
            ipknot_res_df,
            irfold_val2_res_df,
            irfold_corx2_res_df,
        ],
        ["RNAfold", "RNAstructure", "IPknot", "IRFoldVal2", "IRFoldCorX2"],
    ):
        avg_df = res_df.groupby("sequence_length").mean(numeric_only=True)

        seq_interval_means = []
        for start, end in seq_intervals:
            f1_mean = avg_df[(avg_df.index >= start) & (avg_df.index <= end)].f1.mean()
            seq_interval_means.append(f1_mean)

        model_seq_len_interval_f1_means[df_name] = seq_interval_means

        print(
            "".join(
                [df_name.ljust(col_width)]
                + [f" {f1:.2f} ".center(col_width) for f1 in seq_interval_means]
            )
        )

    plt.rcParams["figure.figsize"] = (8, 6)

    for model_name, f1_inverval_means in model_seq_len_interval_f1_means.items():
        plt.plot(
            [i for i in range(len(f1_inverval_means))],
            f1_inverval_means,
            label=model_name,
        )

    plt.bar(
        [i for i in range(n_intervals)],
        normalised_interval_counts,
        label="Sample Dist.",
        alpha=0.3,
    )

    plt.xticks(
        [i for i in range(n_intervals)],
        [f"{i[0]}-{i[1]}".center(9) for i in seq_intervals],
    )
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean F1")
    plt.legend(loc="upper center")
    plt.tight_layout()
    plt.savefig(str(EXPERIMENT_5_DIR_PATH / "models_performance_over_seq_len.png"))
    plt.show()
