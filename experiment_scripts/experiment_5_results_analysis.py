from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).parent.parent / "data"
EXPERIMENT_5_DIR_PATH = DATA_DIR / "experiment_5"

if __name__ == "__main__":
    rnafold_res_df = pd.read_csv(
        EXPERIMENT_5_DIR_PATH / "RNAfold_performance_metrics.csv"
    )
    rnastructure_res_df = pd.read_csv(
        EXPERIMENT_5_DIR_PATH / "RNAstructure_performance_metrics.csv"
    )
    irfold_val2_res_df = pd.read_csv(
        EXPERIMENT_5_DIR_PATH / "IRFoldVal2_performance_metrics.csv"
    )
    irfold_corx2_res_df = pd.read_csv(
        EXPERIMENT_5_DIR_PATH / "IRFoldCorX2_performance_metrics.csv"
    )
    irfold_corx3_res_df = pd.read_csv(
        EXPERIMENT_5_DIR_PATH / "IRFoldCorX3_performance_metrics.csv"
    )

    for res_df, df_name in zip(
        [
            rnafold_res_df,
            rnastructure_res_df,
            irfold_val2_res_df,
            irfold_corx2_res_df,
            irfold_corx3_res_df,
        ],
        ["RNAfold", "RNAstructure", "IRFoldVal2", "IRFoldCorX2", "IRFoldCorX3"],
    ):
        print(f"Model: {df_name}")
        for col_name in [
            "sensitivity",
            "ppv",
            "f1",
            "execution_wall_time_secs",
        ]:
            print(f"   Mean {col_name}: {res_df[col_name].mean():.2f}")

    # Compare RNA SSP programs' performances at sequence length intervals
    #     rnalib_avg_df = rnalib_res_df.groupby("seq_len").mean(numeric_only=True)
    #
    #     seq_interval_sz = 10
    #     seq_intervals = [
    #         [start, start + seq_interval_sz - 1]
    #         for start in range(
    #             rnalib_res_df.seq_len.unique().min(),
    #             rnalib_res_df.seq_len.unique().max(),
    #             seq_interval_sz,
    #         )
    #     ]
    #     seq_intervals_means = {
    #         "IRFoldVal2": [],
    #         "IRFoldCor2": [],
    #         "IRFoldCor3": [],
    #         "RNALib": [],
    #     }
    #     for start, stop in seq_intervals:
    #         irfold_val2_interval_mean = irfold_val2_avg_df[
    #             (irfold_val2_avg_df.index >= start) & (irfold_val2_avg_df.index <= stop)
    #         ].dot_bracket_repr_mfe.mean()
    #         irfold_cor2_interval_mean = irfold_cor2_avg_df[
    #             (irfold_cor2_avg_df.index >= start) & (irfold_cor2_avg_df.index <= stop)
    #         ].dot_bracket_repr_mfe.mean()
    #         irfold_cor3_interval_mean = irfold_cor3_avg_df[
    #             (irfold_cor3_avg_df.index >= start) & (irfold_cor3_avg_df.index <= stop)
    #         ].dot_bracket_repr_mfe.mean()
    #         rnalib_interval_mean = rnalib_avg_df[
    #             (rnalib_avg_df.index >= start) & (rnalib_avg_df.index <= stop)
    #         ].solution_mfe.mean()
    #
    #         seq_intervals_means["IRFoldVal2"].append(irfold_val2_interval_mean)
    #         seq_intervals_means["IRFoldCor2"].append(irfold_cor2_interval_mean)
    #         seq_intervals_means["IRFoldCor3"].append(irfold_cor3_interval_mean)
    #         seq_intervals_means["RNALib"].append(rnalib_interval_mean)
    #
    #     # See below for horrific (but dynamic!) table printing to stdout lol
    #     means_table_print_width = (
    #         (15 * len(seq_intervals_means["RNALib"]))
    #         + 15
    #         + len(seq_intervals_means["RNALib"])
    #         + 2
    #     )
    #     print("-" * means_table_print_width)
    #     print(
    #         "|"
    #         + " " * 15
    #         + "|"
    #         + "Free Energy of Solution - Mean".center(
    #             15 * len(seq_intervals_means["RNALib"]) + 4
    #         )
    #         + "|"
    #     )
    #     print("-" * means_table_print_width)
    #     print(
    #         "|"
    #         + f"Seq. Length->".center(15)
    #         + "|"
    #         + "".join([str(interval).center(15) + "|" for interval in seq_intervals])
    #     )
    #     print("-" * means_table_print_width)
    #     for prog_name, means in seq_intervals_means.items():
    #         print(
    #             "|"
    #             + f"{prog_name}".ljust(15)
    #             + "|"
    #             + "".join(
    #                 [str(round(mean, 4)).ljust(7, "0").center(15) + "|" for mean in means]
    #             )
    #         )
    #     print("-" * means_table_print_width)
