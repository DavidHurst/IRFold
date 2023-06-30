import sys
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

DATA_DIR = str(Path(__file__).parent.parent / "data")

if __name__ == "__main__":
    irfold_base_res_df = pd.read_csv(
        Path(DATA_DIR).resolve() / "IRFoldBase_performance.csv"
    )
    irfold_val1_res_df = pd.read_csv(
        Path(DATA_DIR).resolve() / "IRFoldVal1_performance.csv"
    )
    irfold_val2_res_df = pd.read_csv(
        Path(DATA_DIR).resolve() / "IRFoldVal2_performance.csv"
    )
    irfold_cor2_res_df = pd.read_csv(
        Path(DATA_DIR).resolve() / "IRFoldCor2_performance.csv"
    )
    irfold_cor3_res_df = pd.read_csv(
        Path(DATA_DIR).resolve() / "IRFoldCor3_performance.csv"
    )
    rnalib_res_df = pd.read_csv(Path(DATA_DIR).resolve() / "RNAlib_performance.csv")

    show = True

    plt.rcParams["figure.figsize"] = (12, 6)

    irfold_base_avg_df = irfold_base_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold_val1_avg_df = irfold_val1_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold_val2_avg_df = irfold_val2_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold_cor2_avg_df = irfold_cor2_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold_cor3_avg_df = irfold_cor3_res_df.groupby("seq_len").mean(numeric_only=True)

    # Solver performance comparison
    num_bool_vars_ax = plt.subplot(221)
    num_bool_vars_ax.plot(
        irfold_val1_res_df.seq_len.unique(),
        irfold_val1_avg_df.solver_num_booleans,
        label="IRFoldVal1",
    )
    num_bool_vars_ax.plot(
        irfold_val2_res_df.seq_len.unique(),
        irfold_val2_avg_df.solver_num_booleans,
        label="IRFoldVal2",
    )
    num_bool_vars_ax.plot(
        irfold_cor2_res_df.seq_len.unique(),
        irfold_cor2_avg_df.solver_num_booleans,
        label="IRFoldCor2",
    )
    num_bool_vars_ax.plot(
        irfold_cor3_res_df.seq_len.unique(),
        irfold_cor3_avg_df.solver_num_booleans,
        label="IRFoldCor3",
    )
    num_bool_vars_ax.title.set_text("Num. Booleans Solver Contains")
    num_bool_vars_ax.set_ylabel("Mean Num. Booleans")
    num_bool_vars_ax.set_xlabel("Sequence Length")

    solve_time_ax = plt.subplot(222)
    solve_time_ax.plot(
        irfold_val1_res_df.seq_len.unique(),
        irfold_val1_avg_df.solver_solve_time,
        label="IRFoldVal1",
    )
    solve_time_ax.plot(
        irfold_val2_res_df.seq_len.unique(),
        irfold_val2_avg_df.solver_solve_time,
        label="IRFoldVal2",
    )
    solve_time_ax.plot(
        irfold_cor2_res_df.seq_len.unique(),
        irfold_cor2_avg_df.solver_solve_time,
        label="IRFoldCor2",
    )
    solve_time_ax.plot(
        irfold_cor3_res_df.seq_len.unique(),
        irfold_cor3_avg_df.solver_solve_time,
        label="IRFoldCor3",
    )
    solve_time_ax.title.set_text("Solve Time")
    solve_time_ax.set_ylabel("Mean Solve Time (Milli. Secs.)")
    solve_time_ax.set_xlabel("Sequence Length")

    num_iterations_ax = plt.subplot(223)
    num_iterations_ax.plot(
        irfold_val1_res_df.seq_len.unique(),
        irfold_val1_avg_df.solver_iterations,
        label="IRFoldVal1",
    )
    num_iterations_ax.plot(
        irfold_val2_res_df.seq_len.unique(),
        irfold_val2_avg_df.solver_iterations,
        label="IRFoldVal2",
    )
    num_iterations_ax.plot(
        irfold_cor2_res_df.seq_len.unique(),
        irfold_cor2_avg_df.solver_iterations,
        label="IRFoldCor2",
    )
    num_iterations_ax.plot(
        irfold_cor3_res_df.seq_len.unique(),
        irfold_cor3_avg_df.solver_iterations,
        label="IRFoldCor3",
    )
    num_iterations_ax.title.set_text("Num. Solver Iterations")
    num_iterations_ax.set_ylabel("Mean Iterations")
    num_iterations_ax.set_xlabel("Sequence Length")

    num_branches_expored_ax = plt.subplot(224)
    num_branches_expored_ax.plot(
        irfold_val1_res_df.seq_len.unique(),
        irfold_val1_avg_df.solver_num_branches_explored,
        label="IRFoldVal1",
    )
    num_branches_expored_ax.plot(
        irfold_val2_res_df.seq_len.unique(),
        irfold_val2_avg_df.solver_num_branches_explored,
        label="IRFoldVal2",
    )
    num_branches_expored_ax.plot(
        irfold_cor2_res_df.seq_len.unique(),
        irfold_cor2_avg_df.solver_num_branches_explored,
        label="IRFoldCor2",
    )
    num_branches_expored_ax.plot(
        irfold_cor3_res_df.seq_len.unique(),
        irfold_cor3_avg_df.solver_num_branches_explored,
        label="IRFoldCor3",
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
    plt.savefig(f"{DATA_DIR}/experiment_4/solver_performance_comparison.png")
    if show:
        plt.show()

    plt.rcParams["figure.figsize"] = (16, 10)

    # Compare IRFold versions' objective function error as sequence length increases
    # ToDo: Set x axis to be true free energy and y to be absolute value of objective function
    #           to show deviation from true and compounding of error over time

    # IRFoldBase
    irfold_base_ax = plt.subplot(331)
    irfold_base_ax.plot(
        irfold_base_res_df.seq_len.unique(),
        (
            irfold_base_avg_df.obj_fn_final_value
            - irfold_base_avg_df.dot_bracket_repr_mfe
        ).abs(),
    )
    irfold_base_ax.title.set_text("IRFoldBase")

    # IRFoldVal1
    irfold_val1_ax = plt.subplot(334)
    irfold_val1_ax.plot(
        irfold_val1_res_df.seq_len.unique(),
        (
            irfold_val1_avg_df.obj_fn_final_value
            - irfold_val1_avg_df.dot_bracket_repr_mfe
        ).abs(),
    )
    irfold_val1_ax.title.set_text("IRFoldVal1")

    # IRFoldVal2
    irfold_val2_ax = plt.subplot(335, sharex=irfold_val1_ax)
    irfold_val2_ax.plot(
        irfold_val2_res_df.seq_len.unique(),
        (
            irfold_val2_avg_df.obj_fn_final_value
            - irfold_val2_avg_df.dot_bracket_repr_mfe
        ).abs(),
    )
    irfold_val2_ax.title.set_text("IRFoldVal2")

    # IRFoldCor2
    irfold_cor2_ax = plt.subplot(337)
    irfold_cor2_ax.plot(
        irfold_cor2_res_df.seq_len.unique(),
        (
            irfold_cor2_avg_df.obj_fn_final_value
            - irfold_cor2_avg_df.dot_bracket_repr_mfe
        ).abs(),
    )
    irfold_cor2_ax.title.set_text("IRFoldCor2")

    # IRFoldCor3
    irfold_cor3_ax = plt.subplot(338)
    irfold_cor3_ax.plot(
        irfold_cor3_res_df.seq_len.unique(),
        (
            irfold_cor3_avg_df.obj_fn_final_value
            - irfold_cor3_avg_df.dot_bracket_repr_mfe
        ).abs(),
    )
    irfold_cor3_ax.title.set_text("IRFoldCor3")

    for ax in [
        irfold_base_ax,
        irfold_val1_ax,
        irfold_val2_ax,
        irfold_cor2_ax,
        irfold_cor3_ax,
    ]:
        ax.grid()
        ax.set_ylabel("Absolute Objective Function Error (kcal/mol)")
        ax.set_xlabel("Sequence Length")
    plt.suptitle("Solver Objective Function Error")
    plt.tight_layout()
    plt.savefig(f"{DATA_DIR}/experiment_4/objective_function_error_comparison.png")
    if show:
        plt.show()

    plt.rcParams["figure.figsize"] = (6, 4)

    # MFEs of final solutions
    mfe_data = [
        irfold_val2_res_df.dot_bracket_repr_mfe,
        irfold_cor2_res_df.dot_bracket_repr_mfe,
        irfold_cor3_res_df.dot_bracket_repr_mfe,
        rnalib_res_df.solution_mfe,
    ]
    fig, ax = plt.subplots()

    ax.boxplot(mfe_data, labels=["IRFoldVal2", "IRFoldCor2", "IRFoldCor3", "RNAlib"])

    ax.title.set_text("MFE of Solution")
    plt.ylabel("MFE (kcal/mol)")
    plt.suptitle("SSP Program Solution MFE Comparison")
    plt.tight_layout()
    plt.savefig(f"{DATA_DIR}/experiment_4/ssp_programs_mfe_comparison.png")
    if show:
        plt.show()

    # Compare RNA SSP programs' performances at sequence length intervals
    rnalib_avg_df = rnalib_res_df.groupby("seq_len").mean(numeric_only=True)

    seq_interval_sz = 10
    seq_intervals = [
        [start, start + seq_interval_sz - 1]
        for start in range(
            rnalib_res_df.seq_len.unique().min(),
            rnalib_res_df.seq_len.unique().max(),
            seq_interval_sz,
        )
    ]
    seq_intervals_means = {
        "IRFoldVal2": [],
        "IRFoldCor2": [],
        "IRFoldCor3": [],
        "RNALib": [],
    }
    for start, stop in seq_intervals:
        irfold_val2_interval_mean = irfold_val2_avg_df[
            (irfold_val2_avg_df.index >= start) & (irfold_val2_avg_df.index <= stop)
        ].dot_bracket_repr_mfe.mean()
        irfold_cor2_interval_mean = irfold_cor2_avg_df[
            (irfold_cor2_avg_df.index >= start) & (irfold_cor2_avg_df.index <= stop)
        ].dot_bracket_repr_mfe.mean()
        irfold_cor3_interval_mean = irfold_cor3_avg_df[
            (irfold_cor3_avg_df.index >= start) & (irfold_cor3_avg_df.index <= stop)
        ].dot_bracket_repr_mfe.mean()
        rnalib_interval_mean = rnalib_avg_df[
            (rnalib_avg_df.index >= start) & (rnalib_avg_df.index <= stop)
        ].solution_mfe.mean()

        seq_intervals_means["IRFoldVal2"].append(irfold_val2_interval_mean)
        seq_intervals_means["IRFoldCor2"].append(irfold_cor2_interval_mean)
        seq_intervals_means["IRFoldCor3"].append(irfold_cor3_interval_mean)
        seq_intervals_means["RNALib"].append(rnalib_interval_mean)

    # See below for horrific (but dynamic!) table printing to stdout lol
    means_table_print_width = (
        (15 * len(seq_intervals_means["RNALib"]))
        + 15
        + len(seq_intervals_means["RNALib"])
        + 2
    )
    print("-" * means_table_print_width)
    print(
        "|"
        + " " * 15
        + "|"
        + "Free Energy of Solution - Mean".center(
            15 * len(seq_intervals_means["RNALib"]) + 4
        )
        + "|"
    )
    print("-" * means_table_print_width)
    print(
        "|"
        + f"Seq. Length->".center(15)
        + "|"
        + "".join([str(interval).center(15) + "|" for interval in seq_intervals])
    )
    print("-" * means_table_print_width)
    for prog_name, means in seq_intervals_means.items():
        print(
            "|"
            + f"{prog_name}".ljust(15)
            + "|"
            + "".join(
                [str(round(mean, 4)).ljust(7, "0").center(15) + "|" for mean in means]
            )
        )
    print("-" * means_table_print_width)
