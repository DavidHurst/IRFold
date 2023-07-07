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
