from pathlib import Path

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

    plt.rcParams["figure.figsize"] = (5, 3.5)

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
        ["RNAfold", "RNAstructure", "IPknot", "IRFoldVal2"],
    ):
        f1_scores.append(res_df["f1"].mean())
        running_times.append(res_df["execution_wall_time_secs"].mean())

    print(
        f"Model".ljust(15) + "Mean F1".center(9) + "Mean Running Time (secs)".center(31)
    )
    print("-" * (15 + 9 + 31))
    for model, f1, running_time in zip(
        ["RNAfold", "RNAstructure", "IPknot", "IRFoldVal2"],
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

    print(rnastructure_res_df_pk.sample(3))

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

    print(rnastructure_res_df_no_pk.sample(3))

    markers = {"RNAfold": "o", "RNAstructure": "+", "IPknot": "s", "IRFoldVal2": "^"}

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
        ["RNAfold", "RNAstructure", "IPknot", "IRFoldVal2"],
    ):
        ppv_scores.append(res_df["sensitivity"].mean())
        sens_scores.append(res_df["ppv"].mean())

    for model_name, mean_ppv, mean_sens in zip(
        ["RNAfold", "RNAstructure", "IPknot", "IRFoldVal2"], ppv_scores, sens_scores
    ):
        plt.scatter(
            mean_ppv, mean_sens, label=model_name, marker=markers[model_name], s=50
        )

    plt.legend()
    plt.xlabel("Positive Predictive Value")
    plt.ylabel("Sensitivity")
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
        ["RNAfold", "RNAstructure", "IPknot", "IRFoldVal2"],
    ):
        ppv_scores.append(res_df["sensitivity"].mean())
        sens_scores.append(res_df["ppv"].mean())

    for model_name, mean_ppv, mean_sens in zip(
        ["RNAfold", "RNAstructure", "IPknot", "IRFoldVal2"], ppv_scores, sens_scores
    ):
        plt.scatter(
            mean_ppv, mean_sens, label=model_name, marker=markers[model_name], s=50
        )

    plt.legend()
    plt.xlabel("Positive Predictive Value")
    plt.ylabel("Sensitivity")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_5_DIR_PATH}/ppv_vs_sens_pk.png")
    plt.show()

    # Plot PPV vs sens for IRFold on samples with and without bulge loops
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

    bulge_sens = irfold_val2_res_df_bulge["sensitivity"].mean()
    bulge_ppv = irfold_val2_res_df_bulge["ppv"].mean()
    no_bulge_sens = irfold_val2_res_df_no_bulge["sensitivity"].mean()
    no_bulge_ppv = irfold_val2_res_df_no_bulge["ppv"].mean()

    print(f'No bulge mean F1: {irfold_val2_res_df_no_bulge["f1"].mean():.4f}')
    print(f'Bulge mean F1: {irfold_val2_res_df_bulge["f1"].mean():.4f}')

    plt.scatter(bulge_ppv, bulge_sens, label="Bulge Loops", marker="s", s=50)
    plt.scatter(no_bulge_ppv, no_bulge_sens, label="No Bulge Loops", marker="^", s=50)
    plt.legend()
    plt.xlabel("Positive Predictive Value")
    plt.ylabel("Sensitivity")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_5_DIR_PATH}/ppv_vs_sens_bulge.png")
    plt.show()

    # Plot mean F1 as sequence length increases
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

    model_seq_len_interval_f1_means = {
        "RNAfold": [],
        "RNAstructure": [],
        "IPknot": [],
        "IRFoldVal2": [],
    }
    for res_df, df_name in zip(
        [
            rnafold_res_df,
            rnastructure_res_df,
            ipknot_res_df,
            irfold_val2_res_df,
        ],
        ["RNAfold", "RNAstructure", "IPknot", "IRFoldVal2"],
    ):
        avg_df = res_df.groupby("sequence_length").mean(numeric_only=True)

        seq_interval_means = []
        for start, end in seq_intervals:
            f1_mean = avg_df[(avg_df.index >= start) & (avg_df.index <= end)].f1.mean()
            seq_interval_means.append(f1_mean)

        model_seq_len_interval_f1_means[df_name] = seq_interval_means

        # print(
        #     "".join(
        #         [df_name.ljust(col_width)]
        #         + [f" {f1:.2f} ".center(col_width) for f1 in seq_interval_means]
        #     )
        # )

    for model_name, f1_inverval_means in model_seq_len_interval_f1_means.items():
        plt.plot(
            [i for i in range(len(f1_inverval_means))],
            f1_inverval_means,
            label=model_name,
        )

    # plt.bar(
    #     [i for i in range(n_intervals)],
    #     normalised_interval_counts,
    #     label="Sample Dist.",
    #     alpha=0.3,
    # )

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
