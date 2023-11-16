from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA_DIR = Path(__file__).parent.parent / "data"
EXPERIMENT_5_DIR_PATH = DATA_DIR / "experiment_5"
BP_RNA_1_M_ST_FILES_DIR = DATA_DIR / "bpRNA-1m" / "stFiles"

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

    plt.rcParams["figure.figsize"] = (3.5, 3)

    print(f"Num. Samples: {len(irfold_val2_res_df)}")

    # Print average F1 scores and running times for each model
    f1_scores = []
    running_times = []
    for res_df, df_name in zip(
        [
            rnafold_res_df,
            rnastructure_res_df,
            ipknot_res_df,
            irfold_val2_res_df,
        ],
        ["RNAfold", "RNAstructure", "IPknot", "IRFold"],
    ):
        f1_scores.append(res_df["f1"].mean())
        running_times.append(res_df["execution_wall_time_secs"].mean())

    print(
        f"Model".ljust(15) + "Mean F1".center(9) + "Mean Running Time (secs)".center(31)
    )
    print("-" * (15 + 9 + 31))
    for model, f1, running_time in zip(
        ["RNAfold", "RNAstructure", "IPknot", "IRFold"],
        f1_scores,
        running_times,
    ):
        print(
            f"{model}".ljust(15)
            + f"{f1:.2f}".center(9)
            + f"{running_time:.4f}".center(31)
        )

    # Get models' results DFs for pseudoknotted and non-pseudoknotted structures
    sample_db_names = []
    sample_db_nums = []
    samples_contain_pk = []
    samples_contain_bulge = []
    for st_file_path in BP_RNA_1_M_ST_FILES_DIR.glob("*.st"):
        sample_id = str(st_file_path).split("/")[-1].split(".")[0]
        sample_db_name = sample_id.split("_")[1]
        sample_db_num = sample_id.split("_")[-1]

        contains_pk = False
        contains_bulge = False

        with open(st_file_path, "r") as st_file:
            lines = st_file.readlines()
            for line in lines[7:]:
                if "PK" in line:
                    contains_pk = True
                if "B" in line:
                    contains_bulge = True

        sample_db_names.append(sample_db_name)
        sample_db_nums.append(int(sample_db_num))
        samples_contain_pk.append(contains_pk)
        samples_contain_bulge.append(contains_bulge)
    samples_pk_status_df = pd.DataFrame(
        {
            "database_name": sample_db_names,
            "sequence_database_number": pd.Series(sample_db_nums, dtype="int64"),
            "contains_pk": samples_contain_pk,
            "contains_bulge": samples_contain_bulge,
        }
    )

    # Inner join results DFs with pk_status=True DF to get results of pseudoknotted sequences
    pk_samples = samples_pk_status_df[samples_pk_status_df.contains_pk == True]
    rnafold_res_df_pk = pd.merge(
        rnafold_res_df, pk_samples, on=["database_name", "sequence_database_number"]
    )
    rnastructure_res_df_pk = pd.merge(
        rnastructure_res_df,
        pk_samples,
        on=["database_name", "sequence_database_number"],
    )
    ipknot_res_df_pk = pd.merge(
        ipknot_res_df, pk_samples, on=["database_name", "sequence_database_number"]
    )
    irfold_val2_res_df_pk = pd.merge(
        irfold_val2_res_df, pk_samples, on=["database_name", "sequence_database_number"]
    )

    # print(rnastructure_res_df_pk.sample(3))

    # Inner join results DFs with pk_status=False DF to get results of pseudoknotted sequences
    no_pk_samples = samples_pk_status_df[samples_pk_status_df.contains_pk == False]
    rnafold_res_df_no_pk = pd.merge(
        rnafold_res_df, no_pk_samples, on=["database_name", "sequence_database_number"]
    )
    rnastructure_res_df_no_pk = pd.merge(
        rnastructure_res_df,
        no_pk_samples,
        on=["database_name", "sequence_database_number"],
    )
    ipknot_res_df_no_pk = pd.merge(
        ipknot_res_df, no_pk_samples, on=["database_name", "sequence_database_number"]
    )
    irfold_val2_res_df_no_pk = pd.merge(
        irfold_val2_res_df,
        no_pk_samples,
        on=["database_name", "sequence_database_number"],
    )

    print(f"Num no PK seqs: {len(rnastructure_res_df_no_pk)}")
    print(f"Num PK seqs: {len(rnastructure_res_df_pk)}")

    # print(rnastructure_res_df_no_pk.sample(3))

    markers = {
        "RNAfold": "o",
        "RNAstructure": "x",
        "IPknot": "s",
        "IRFold": "^",
        "IRFold Solver": "^",
    }

    # Plot PPV against Sens for all seqs
    ppv_scores = []
    sens_scores = []
    for res_df, df_name in zip(
        [
            rnafold_res_df,
            rnastructure_res_df,
            ipknot_res_df,
            irfold_val2_res_df,
        ],
        ["RNAfold", "RNAstructure", "IPknot", "IRFold"],
    ):
        ppv_scores.append(res_df["sensitivity"].mean())
        sens_scores.append(res_df["ppv"].mean())

    for model_name, mean_ppv, mean_sens in zip(
        ["RNAfold", "RNAstructure", "IPknot", "IRFold"], ppv_scores, sens_scores
    ):
        plt.scatter(
            mean_ppv, mean_sens, label=model_name, marker=markers[model_name], s=50
        )

    plt.legend()
    plt.xlabel("Positive Predictive Value")
    plt.ylabel("Sensitivity")
    plt.ylim(0.7, 0.9)
    plt.xlim(0.4, 0.8)
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_5_DIR_PATH}/ppv_vs_sens_all.png")
    plt.show()

    # Plot PPV against Sens for non-pseudoknotted structures
    ppv_scores = []
    sens_scores = []
    for res_df, df_name in zip(
        [
            rnafold_res_df_no_pk,
            rnastructure_res_df_no_pk,
            ipknot_res_df_no_pk,
            irfold_val2_res_df_no_pk,
        ],
        ["RNAfold", "RNAstructure", "IPknot", "IRFold"],
    ):
        ppv_scores.append(res_df["sensitivity"].mean())
        sens_scores.append(res_df["ppv"].mean())

    for model_name, mean_ppv, mean_sens in zip(
        ["RNAfold", "RNAstructure", "IPknot", "IRFold"], ppv_scores, sens_scores
    ):
        plt.scatter(
            mean_ppv, mean_sens, label=model_name, marker=markers[model_name], s=50
        )

    plt.legend()
    plt.xlabel("Positive Predictive Value")
    plt.ylabel("Sensitivity")
    plt.ylim(0.7, 0.9)
    plt.xlim(0.4, 0.8)
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_5_DIR_PATH}/ppv_vs_sens_no_pk.png")
    plt.show()

    # Plot PPV against Sens for pseudoknotted structures
    ppv_scores = []
    sens_scores = []
    for res_df, df_name in zip(
        [
            rnafold_res_df_pk,
            rnastructure_res_df_pk,
            ipknot_res_df_pk,
            irfold_val2_res_df_pk,
        ],
        ["RNAfold", "RNAstructure", "IPknot", "IRFold"],
    ):
        ppv_scores.append(res_df["sensitivity"].mean())
        sens_scores.append(res_df["ppv"].mean())

    for model_name, mean_ppv, mean_sens in zip(
        ["RNAfold", "RNAstructure", "IPknot", "IRFold"], ppv_scores, sens_scores
    ):
        plt.scatter(
            mean_ppv, mean_sens, label=model_name, marker=markers[model_name], s=50
        )

    plt.legend(loc="lower right")
    plt.xlabel("Positive Predictive Value")
    plt.ylabel("Sensitivity")
    plt.ylim(0.7, 0.9)
    plt.xlim(0.4, 0.8)
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_5_DIR_PATH}/ppv_vs_sens_pk.png")
    plt.show()

    # Plot PPV vs sens for IRFold on samples with and without bulge loops and without pseudoknots
    no_bulge_samples = samples_pk_status_df[
        samples_pk_status_df.contains_bulge == False
    ]
    bulge_samples = samples_pk_status_df[samples_pk_status_df.contains_bulge == True]
    irfold_val2_res_df_no_bulge = pd.merge(
        irfold_val2_res_df,
        no_bulge_samples,
        on=["database_name", "sequence_database_number"],
    )
    irfold_val2_res_df_bulge = pd.merge(
        irfold_val2_res_df,
        bulge_samples,
        on=["database_name", "sequence_database_number"],
    )

    print(f"Num no bulge: {len(irfold_val2_res_df_bulge)}")
    print(f"Num bulge: {len(irfold_val2_res_df_no_bulge)}")

    bulge_sens = irfold_val2_res_df_bulge["sensitivity"].mean()
    bulge_ppv = irfold_val2_res_df_bulge["ppv"].mean()
    no_bulge_sens = irfold_val2_res_df_no_bulge["sensitivity"].mean()
    no_bulge_ppv = irfold_val2_res_df_no_bulge["ppv"].mean()

    print(f'No bulge mean F1: {irfold_val2_res_df_no_bulge["f1"].mean():.2f}')
    print(f'Bulge mean F1: {irfold_val2_res_df_bulge["f1"].mean():.2f}')

    plt.scatter(bulge_ppv, bulge_sens, label="Bulge Loops", marker="s", s=50)
    plt.scatter(no_bulge_ppv, no_bulge_sens, label="No Bulge Loops", marker="^", s=50)
    plt.legend()
    plt.xlabel("Positive Predictive Value")
    plt.ylabel("Sensitivity")
    # plt.ylim(0.7, 0.9)
    # plt.xlim(0.4, 0.8)
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_5_DIR_PATH}/ppv_vs_sens_bulge.png")
    plt.show()

    plt.rcParams["figure.figsize"] = (8, 5)

    # Plot mean F1 as sequence length increases for sequence length intervals
    min = 5  # rnafold_res_df.sequence_length.min()
    max = rnafold_res_df.sequence_length.max() + 5
    interval_sz = 10
    vals = np.arange(min, max + interval_sz, interval_sz)
    intervals = [[vals[i], vals[i + 1]] for i in range(len(vals) - 1)]
    n_intervals = len(intervals)
    for range in intervals[:-1]:
        range[1] = range[1] - 1

    model_seq_len_interval_f1_means = {
        "RNAfold": [],
        "RNAstructure": [],
        "IPknot": [],
        "IRFold": [],
    }
    model_seq_len_interval_seq_counts = {
        "RNAfold": [],
        "RNAstructure": [],
        "IPknot": [],
        "IRFold": [],
    }
    model_seq_len_interval_mean_running_times = {
        "RNAfold": [],
        "RNAstructure": [],
        "IPknot": [],
        "IRFold": [],
    }
    for res_df, df_name in zip(
        [
            rnafold_res_df,
            rnastructure_res_df,
            ipknot_res_df,
            irfold_val2_res_df,
        ],
        ["RNAfold", "RNAstructure", "IPknot", "IRFold"],
    ):
        avg_df = res_df.groupby("sequence_length").mean(numeric_only=True)

        seq_interval_f1_means = []
        seq_interval_seq_counts = []
        seq_interval_running_times = []

        for start, end in intervals:
            f1_mean = avg_df[(avg_df.index >= start) & (avg_df.index <= end)].f1.mean()
            seq_count = len(
                res_df[
                    (res_df.sequence_length >= start) & (res_df.sequence_length <= end)
                ]
            )
            running_time_mean = avg_df[
                (avg_df.index >= start) & (avg_df.index <= end)
            ].execution_wall_time_secs.mean()

            seq_interval_f1_means.append(f1_mean)
            seq_interval_seq_counts.append(seq_count)
            seq_interval_running_times.append(running_time_mean)

        model_seq_len_interval_f1_means[df_name] = seq_interval_f1_means
        model_seq_len_interval_seq_counts[df_name] = seq_interval_seq_counts
        model_seq_len_interval_mean_running_times[df_name] = seq_interval_running_times

    fig, (ax1, ax2) = plt.subplots(
        2,
        1,
        sharey="row",
        sharex="col",
        gridspec_kw={"height_ratios": [4, 1]},
    )
    fig.subplots_adjust(wspace=0.025, bottom=0.2)
    fig.text(0.5, 0.03, "Sequence Length (bases)", ha="center")

    for model_name, means in model_seq_len_interval_f1_means.items():
        ax1.plot(
            np.arange(0, len(intervals)),
            means,
            label=model_name,
            marker=markers[model_name],
        )

    for model_name, counts in model_seq_len_interval_seq_counts.items():
        ax2bar = ax2.bar(np.arange(0, len(intervals)), counts, color="grey")

        ax2.bar_label(
            ax2bar,
            labels=[int(val) if val <= 5 else "" for val in counts],
        )

    ticks = [f"{i[0]}-{i[1]}" for i in intervals]
    ax1.set_xticks(
        np.arange(0, n_intervals),
        ticks,
    )
    ax2.set_xticks(
        np.arange(0, n_intervals),
        ticks,
    )

    ax1.legend()
    ax1.tick_params(axis="x", rotation=90)
    ax1.set_ylabel("Mean F1")
    ax1.grid(axis="y")

    ax2.tick_params(axis="x", rotation=90)
    ax2.set_ylabel("# Sequences")

    plt.savefig(str(EXPERIMENT_5_DIR_PATH / "mean_f1_vs_seq_len.png"))
    plt.show()

    # Plot running time as sequence length increases
    plt.rcParams["figure.figsize"] = (8, 3)
    for model_name, means in model_seq_len_interval_mean_running_times.items():
        plt.plot(
            np.arange(0, len(intervals)),
            means,
            label=model_name,
            marker=markers[model_name],
        )

    ticks = [f"{i[0]}-{i[1]}" for i in intervals]
    plt.xticks(
        np.arange(0, n_intervals),
        ticks,
    )
    plt.legend()
    plt.tick_params(axis="x", rotation=90)
    plt.ylabel("Mean Running \nTime (secs)")
    plt.xlabel("Sequence Length (bases)")
    plt.grid(axis="y")
    plt.tight_layout()
    plt.savefig(str(EXPERIMENT_5_DIR_PATH / "mean_running_time_vs_seq_len.png"))
    plt.show()

    # Plot running time as sequence length increases with IRFold solver
    irfold_solver_perf_df = pd.read_csv(
        EXPERIMENT_5_DIR_PATH / "IRFoldVal2_solver_performance.csv"
    )
    avg_df = irfold_solver_perf_df.groupby("seq_len").mean(numeric_only=True)

    seq_interval_running_times = []

    for start, end in intervals:
        running_time_mean = avg_df[
            (avg_df.index >= start) & (avg_df.index <= end)
        ].solver_solve_time.mean()
        seq_interval_running_times.append(running_time_mean)

    model_seq_len_interval_mean_running_times[
        "IRFold Solver"
    ] = seq_interval_running_times
    del model_seq_len_interval_mean_running_times["IRFold"]

    for model_name, means in model_seq_len_interval_mean_running_times.items():
        plt.plot(
            np.arange(0, len(intervals)),
            means,
            label=model_name,
            marker=markers[model_name],
        )

    ticks = [f"{i[0]}-{i[1]}" for i in intervals]
    plt.xticks(
        np.arange(0, n_intervals),
        ticks,
    )
    plt.legend()
    plt.tick_params(axis="x", rotation=90)
    plt.ylabel("Mean Running \nTime (secs)")
    plt.xlabel("Sequence Length (bases)")
    plt.grid(axis="y")
    plt.tight_layout()
    plt.savefig(
        str(EXPERIMENT_5_DIR_PATH / "mean_running_time_vs_seq_len_irfold_solver.png")
    )
    plt.show()
