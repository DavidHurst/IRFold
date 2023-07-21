import sys
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

DATA_DIR = (Path(__file__).parent.parent / "data").resolve()
EXPERIMENT_4_DATA_DIR = (DATA_DIR / "experiment_4").resolve()

if __name__ == "__main__":
    show = True

    irfold_base_res_df = pd.read_csv(
        EXPERIMENT_4_DATA_DIR / "IRFoldBase_solver_performance.csv"
    )
    irfold_val1_res_df = pd.read_csv(
        EXPERIMENT_4_DATA_DIR / "IRFoldVal1_solver_performance.csv"
    )
    irfold_val2_res_df = pd.read_csv(
        EXPERIMENT_4_DATA_DIR / "IRFoldVal2_solver_performance.csv"
    )
    irfold_corx2_res_df = pd.read_csv(
        EXPERIMENT_4_DATA_DIR / "IRFoldCorX2_solver_performance.csv"
    )
    irfold_corx3_res_df = pd.read_csv(
        EXPERIMENT_4_DATA_DIR / "IRFoldCorX3_solver_performance.csv"
    )
    irfold_corx4_res_df = pd.read_csv(
        EXPERIMENT_4_DATA_DIR / "IRFoldCorX4_solver_performance.csv"
    )


    irfold_base_avg_df = irfold_base_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold_val1_avg_df = irfold_val1_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold_val2_avg_df = irfold_val2_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold_corx2_avg_df = irfold_corx2_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold_corx3_avg_df = irfold_corx3_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold_corx4_avg_df = irfold_corx4_res_df.groupby("seq_len").mean(numeric_only=True)

    res_dfs = [
        irfold_base_res_df,
        irfold_val1_res_df,
        irfold_val2_res_df,
        irfold_corx2_res_df,
        irfold_corx3_res_df,
        irfold_corx4_res_df,
    ]
    avg_res_dfs = [
        irfold_base_avg_df,
        irfold_val1_avg_df,
        irfold_val2_avg_df,
        irfold_corx2_avg_df,
        irfold_corx3_avg_df,
        irfold_corx4_avg_df,
    ]
    df_model_names = [
        "IRFoldBase",
        "IRFoldVal1",
        "IRFoldVal2",
        "IRFoldCorX2",
        "IRFoldCorX3",
        "IRFoldCorX4",
    ]

    plt.rcParams["figure.figsize"] = (6, 4)

    # Compare num solver variables
    for res_df, avg_df, df_model_name in zip(res_dfs, avg_res_dfs, df_model_names):
        plt.plot(
            res_df.seq_len.unique(),
            avg_df.solver_num_booleans,
            label=df_model_name,
        )

    plt.legend()
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("# Solver Variables")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/num_solver_variables_comparison.png")
    if show:
        plt.show()

    # Compare solver wall time
    for res_df, avg_df, df_model_name in zip(res_dfs, avg_res_dfs, df_model_names):
        plt.plot(
            res_df.seq_len.unique(),
            avg_df.solver_solve_time,
            label=df_model_name,
        )

    plt.legend()
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean Solver Wall Time \n (Milliseconds)")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/solver_wall_time_comparison.png")
    if show:
        plt.show()

    # Compare solver iterations
    for res_df, avg_df, df_model_name in zip(res_dfs, avg_res_dfs, df_model_names):
        plt.plot(
            res_df.seq_len.unique(),
            avg_df.solver_iterations,
            label=df_model_name,
        )

    plt.legend()
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean # Solver Iterations")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/solver_iterations_comparison.png")
    if show:
        plt.show()

    # Compare solver branches explored
    for res_df, avg_df, df_model_name in zip(res_dfs, avg_res_dfs, df_model_names):
        plt.plot(
            res_df.seq_len.unique(),
            avg_df.solver_num_branches_explored,
            label=df_model_name,
        )

    plt.legend()
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean # Solver Branches Explored")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/solver_branches_explored_comparison.png")
    if show:
        plt.show()

    # Compare IRFold versions' objective function error as sequence length increases

    # IRFoldBase
    plt.plot(
        irfold_base_res_df.seq_len.unique(),
        (
            irfold_base_avg_df.obj_fn_final_value
            - irfold_base_avg_df.dot_bracket_repr_mfe
        ).abs(),
    )
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean Abs. Objective Function Error (kcal/mol)")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/IRFoldBase_obj_fn_err.png")
    if show:
        plt.show()

    # IRFoldVal1
    plt.plot(
        irfold_val1_res_df.seq_len.unique(),
        (
            irfold_val1_avg_df.obj_fn_final_value
            - irfold_val1_avg_df.dot_bracket_repr_mfe
        ).abs(),
    )
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean Abs. Objective Function Error (kcal/mol)")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/IRFoldVal1_obj_fn_err.png")
    if show:
        plt.show()

    plt.plot(
        irfold_val2_res_df.seq_len.unique(),
        (
            irfold_val2_avg_df.obj_fn_final_value
            - irfold_val2_avg_df.dot_bracket_repr_mfe
        ).abs(),
    )
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean Abs. Objective Function Error (kcal/mol)")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/IRFoldVal2_obj_fn_err.png")
    if show:
        plt.show()

    # IRFoldCor models
    for res_df, avg_df, model_name in zip(
        [
            irfold_corx2_res_df,
            irfold_corx3_res_df,
            irfold_corx4_res_df,
        ],
        [
            irfold_val2_avg_df,
            irfold_corx2_avg_df,
            irfold_corx3_avg_df,
            irfold_corx4_avg_df,
        ],
        [
            "IRFoldCorX2",
            "IRFoldCorX3",
            "IRFoldCorX4",
        ],
    ):
        plt.plot(
            res_df.seq_len.unique(),
            (avg_df.obj_fn_final_value - avg_df.dot_bracket_repr_mfe).abs(),
            label=model_name,
        )
    plt.legend()
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean Abs. Objective Function Error (kcal/mol)")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/IRFoldCorX_obj_fn_err.png")
    if show:
        plt.show()
