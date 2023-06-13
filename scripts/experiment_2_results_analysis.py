import sys
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

DATA_DIR = str(Path(__file__).parent.parent / "data")

if __name__ == "__main__":
    irfold0_res_df = pd.read_csv(Path(DATA_DIR).resolve() / "IRFold0_performance.csv")
    irfold1_res_df = pd.read_csv(Path(DATA_DIR).resolve() / "IRFold1_performance.csv")
    irfold2_res_df = pd.read_csv(Path(DATA_DIR).resolve() / "IRFold2_performance.csv")
    rnalib_res_df = pd.read_csv(Path(DATA_DIR).resolve() / "RNAlib_performance.csv")

    plt.rcParams["figure.figsize"] = (12, 6)

    irfold0_avg_df = irfold0_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold1_avg_df = irfold1_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold2_avg_df = irfold2_res_df.groupby("seq_len").mean(numeric_only=True)

    # Solver performance comparison
    num_bool_vars_ax = plt.subplot(221)
    num_bool_vars_ax.plot(
        irfold1_res_df.seq_len.unique(),
        irfold1_avg_df.solver_num_booleans,
        label="IRFold1",
    )
    num_bool_vars_ax.plot(
        irfold2_res_df.seq_len.unique(),
        irfold2_avg_df.solver_num_booleans,
        label="IRFold2",
    )
    num_bool_vars_ax.title.set_text("Num. Booleans Solver Contains")
    num_bool_vars_ax.set_ylabel("Mean Num. Booleans")
    num_bool_vars_ax.set_xlabel("Sequence Length")

    solve_time_ax = plt.subplot(222)
    solve_time_ax.plot(
        irfold1_res_df.seq_len.unique(),
        irfold1_avg_df.solver_solve_time,
        label="IRFold1",
    )
    solve_time_ax.plot(
        irfold2_res_df.seq_len.unique(),
        irfold2_avg_df.solver_solve_time,
        label="IRFold2",
    )
    solve_time_ax.title.set_text("Solve Time")
    solve_time_ax.set_ylabel("Mean Solve Time (Milli. Secs.)")
    solve_time_ax.set_xlabel("Sequence Length")

    num_iterations_ax = plt.subplot(223)
    num_iterations_ax.plot(
        irfold1_res_df.seq_len.unique(),
        irfold1_avg_df.solver_iterations,
        label="IRFold1",
    )
    num_iterations_ax.plot(
        irfold2_res_df.seq_len.unique(),
        irfold2_avg_df.solver_iterations,
        label="IRFold2",
    )
    num_iterations_ax.title.set_text("Num. Solver Iterations")
    num_iterations_ax.set_ylabel("Mean Iterations")
    num_iterations_ax.set_xlabel("Sequence Length")

    num_branches_expored_ax = plt.subplot(224)
    num_branches_expored_ax.plot(
        irfold1_res_df.seq_len.unique(),
        irfold1_avg_df.solver_num_branches_explored,
        label="IRFold1",
    )
    num_branches_expored_ax.plot(
        irfold2_res_df.seq_len.unique(),
        irfold2_avg_df.solver_num_branches_explored,
        label="IRFold2",
    )
    num_branches_expored_ax.title.set_text("Num. Branches Explored")
    num_branches_expored_ax.set_ylabel("Mean Num. Branches")
    num_branches_expored_ax.set_xlabel("Sequence Length")

    for ax in [
        num_bool_vars_ax,
        solve_time_ax,
        num_iterations_ax,
        num_branches_expored_ax,
    ]:
        ax.legend()
        ax.grid()

    plt.suptitle("Solver Performance Comparison")
    plt.tight_layout()
    plt.savefig(f"{DATA_DIR}/experiment_2_solver_performance_comparison.png")
    plt.show()

    plt.rcParams["figure.figsize"] = (12, 4)

    # Compare IRFold versions: True MFE of IRFold solutions vs. MFE IRFold returns
    # This is to show the effect of additivity assumption inconsistency
    irfold0_ax = plt.subplot(131)
    irfold0_ax.plot(
        irfold0_res_df.seq_len.unique(),
        irfold0_avg_df.obj_fn_final_value,
        label="Objective Function Value",
    )
    irfold0_ax.plot(
        irfold0_res_df.seq_len.unique(),
        irfold0_avg_df.dot_bracket_repr_mfe,
        label="Solution MFE",
    )
    irfold0_ax.title.set_text("IRFold0")

    irfold1_ax = plt.subplot(132)
    irfold1_ax.plot(
        irfold1_res_df.seq_len.unique(),
        irfold1_avg_df.obj_fn_final_value,
        label="Objective Function Value",
    )
    irfold1_ax.plot(
        irfold1_res_df.seq_len.unique(),
        irfold1_avg_df.dot_bracket_repr_mfe,
        label="Solution MFE",
    )
    irfold1_ax.title.set_text("IRFold1")

    irfold2_ax = plt.subplot(133, sharex=irfold1_ax)
    irfold2_ax.plot(
        irfold2_res_df.seq_len.unique(),
        irfold2_avg_df.obj_fn_final_value,
        label="Objective Function Value",
    )
    irfold2_ax.plot(
        irfold2_res_df.seq_len.unique(),
        irfold2_avg_df.dot_bracket_repr_mfe,
        label="Solution MFE",
    )
    irfold2_ax.title.set_text("IRFold2")

    for ax in [irfold1_ax, irfold2_ax]:
        ax.legend()
        ax.grid()
        ax.set_ylabel("Mean Value (kcal/mol)")
        ax.set_xlabel("Sequence Length")
    plt.suptitle("Solver Objective Function Final Value vs. Solution MFE")
    plt.tight_layout()
    plt.savefig(f"{DATA_DIR}/experiment_2_mfe_vs_objective_value_comparison.png")
    plt.show()

    plt.rcParams["figure.figsize"] = (6, 4)

    # MFEs of final solutions
    mfe_data = [
        irfold1_res_df.dot_bracket_repr_mfe,
        irfold2_res_df.dot_bracket_repr_mfe,
        rnalib_res_df.solution_mfe,
    ]
    fig, ax = plt.subplots()

    ax.boxplot(mfe_data, labels=["IRFold1", "IRFold2", "RNAlib"])

    ax.title.set_text("MFE of Solution")
    plt.ylabel("MFE (kcal/mol)")
    plt.suptitle("SSP Program Solution MFE Comparison")
    plt.tight_layout()
    plt.savefig(f"{DATA_DIR}/experiment_2_mfe_comparison.png")
    plt.show()

    # Compare IRFold versions' MFE's as sequence length increases
    plt.plot(
        irfold1_res_df.seq_len.unique(),
        irfold1_avg_df.dot_bracket_repr_mfe,
        label="IRFold1",
    )
    plt.plot(
        irfold2_res_df.seq_len.unique(),
        irfold2_avg_df.dot_bracket_repr_mfe,
        label="IRFold2",
    )

    plt.legend()
    plt.grid()
    plt.ylabel("Mean MFE of Solution")
    plt.xlabel("Sequence Length")
    plt.suptitle("IRFold Variant's MFEs over Time")
    plt.tight_layout()
    plt.savefig(f"{DATA_DIR}/experiment_2_mfe_over_time_comparison.png")
    plt.show()
