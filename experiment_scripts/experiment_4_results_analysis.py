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

    line_markers = {
        "IRFoldVal2": "x",
        "IRFoldCor2": "o",
        "IRFoldCor3": "s",
        "IRFoldCor4": "^",
    }

    plt.rcParams["figure.figsize"] = (5, 3.5)

    plt.plot(
        irfold_base_res_df.seq_len.unique(),
        irfold_base_avg_df.n_irs_found,
    )
    plt.grid()
    plt.ylabel("Mean # IRs Found")
    plt.xlabel("Sequence Length")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/num_irs_found.png")
    if show:
        plt.show()

    # Compare IRFold versions' prediction's free energy values
    # IRFoldBase
    plt.plot(
        irfold_base_res_df.seq_len.unique(),
        irfold_base_avg_df.dot_bracket_repr_mfe,
    )
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean Free Energy (kcal/mol)")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/IRFoldBase_preds_fe.png")
    if show:
        plt.show()

    # IRFoldVal1
    plt.plot(
        irfold_val1_res_df.seq_len.unique(),
        irfold_val1_avg_df.dot_bracket_repr_mfe,
    )
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean Free Energy (kcal/mol)")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/IRFoldVal1_preds_fe.png")
    if show:
        plt.show()

    plt.plot(
        irfold_val2_res_df.seq_len.unique(),
        irfold_val2_avg_df.dot_bracket_repr_mfe,
    )
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean Free Energy (kcal/mol)")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/IRFoldVal2_preds_fe.png")
    if show:
        plt.show()

    # Compare IRFold versions' objective function error as sequence length increases

    # IRFoldVal2
    plt.plot(
        irfold_val2_res_df.seq_len.unique(),
        (
            irfold_val2_avg_df.obj_fn_final_value
            - irfold_val2_avg_df.dot_bracket_repr_mfe
        ).abs(),
        marker=line_markers["IRFoldVal2"],
    )
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean Objective Function \n Error (kcal/mol)")
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
            irfold_corx2_avg_df,
            irfold_corx3_avg_df,
            irfold_corx4_avg_df,
        ],
        [
            "IRFoldCor2",
            "IRFoldCor3",
            "IRFoldCor4",
        ],
    ):
        plt.plot(
            res_df.seq_len.unique(),
            (avg_df.obj_fn_final_value - avg_df.dot_bracket_repr_mfe).abs(),
            label=model_name,
            marker=line_markers[model_name],
        )
        plt.plot(
            irfold_val2_res_df.seq_len.unique(),
            (
                irfold_val2_avg_df.obj_fn_final_value
                - irfold_val2_avg_df.dot_bracket_repr_mfe
            ).abs(),
            label="IRFoldVal2",
            marker=line_markers["IRFoldVal2"],
        )
        ax = plt.gca()
        ax.set_ylim([0, 30])
        plt.grid()
        plt.legend()
        plt.xlabel("Sequence Length")
        plt.ylabel("Mean Objective Function \n Error (kcal/mol)")
        plt.tight_layout()
        plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/{model_name}_obj_fn_err.png")
        if show:
            plt.show()

    # Compare num solver variables

    res_dfs = [
        irfold_val2_res_df,
        irfold_corx2_res_df,
        irfold_corx3_res_df,
        irfold_corx4_res_df,
    ]
    avg_res_dfs = [
        irfold_val2_avg_df,
        irfold_corx2_avg_df,
        irfold_corx3_avg_df,
        irfold_corx4_avg_df,
    ]
    model_names = [
        "IRFoldVal2",
        "IRFoldCor2",
        "IRFoldCor3",
        "IRFoldCor4",
    ]

    for res_df, avg_df, model_name in zip(res_dfs, avg_res_dfs, model_names):
        plt.plot(
            res_df.seq_len.unique(),
            avg_df.solver_num_booleans,
            label=model_name,
            marker=line_markers[model_name],
        )
    plt.legend()
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean # \n Objective Function Variables")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/solver_variable_count_comparison.png")
    if show:
        plt.show()

    # Compare solver wall time
    for res_df, avg_df, model_name in zip(res_dfs, avg_res_dfs, model_names):
        plt.plot(
            res_df.seq_len.unique(),
            avg_df.solver_solve_time,
            label=model_name,
            marker=line_markers[model_name],
        )

    plt.legend()
    plt.grid()
    plt.xlabel("Sequence Length")
    plt.ylabel("Mean Solver Wall Time \n (Milliseconds)")
    plt.tight_layout()
    plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/solver_wall_time_comparison.png")
    if show:
        plt.show()

    # # Compare solver iterations
    # for res_df, avg_df, df_model_name in zip(res_dfs, avg_res_dfs, df_model_names):
    #     plt.plot(
    #         res_df.seq_len.unique(),
    #         avg_df.solver_iterations,
    #         label=df_model_name,
    #     )
    #
    # plt.legend()
    # plt.grid()
    # plt.xlabel("Sequence Length")
    # plt.ylabel("Mean # Solver Iterations")
    # plt.tight_layout()
    # plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/solver_iterations_comparison.png")
    # if show:
    #     plt.show()
    #
    # # Compare solver branches explored
    # for res_df, avg_df, df_model_name in zip(res_dfs, avg_res_dfs, df_model_names):
    #     plt.plot(
    #         res_df.seq_len.unique(),
    #         avg_df.solver_num_branches_explored,
    #         label=df_model_name,
    #     )
    #
    # plt.legend()
    # plt.grid()
    # plt.xlabel("Sequence Length")
    # plt.ylabel("Mean # Solver Branches Explored")
    # plt.tight_layout()
    # plt.savefig(f"{EXPERIMENT_4_DATA_DIR}/solver_branches_explored_comparison.png")
    # if show:
    #     plt.show()
