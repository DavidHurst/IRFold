import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

DATA_DIR = str(Path(__file__).parent.parent / "data")

if __name__ == "__main__":
    irfold0_res_df = pd.read_csv(Path(DATA_DIR).resolve() / "IRFold0_performance.csv")
    irfold1_res_df = pd.read_csv(Path(DATA_DIR).resolve() / "IRFold1_performance.csv")
    rnalib_res_df = pd.read_csv(Path(DATA_DIR).resolve() / "RNAlib_performance.csv")

    plt.rcParams["figure.figsize"] = (12, 4)

    irfold0_avg_df = irfold0_res_df.groupby("seq_len").mean(numeric_only=True)
    irfold1_avg_df = irfold1_res_df.groupby("seq_len").mean(numeric_only=True)

    # Solver performance comparison
    plt.subplot(1, 2, 1)
    plt.plot(
        irfold0_res_df.seq_len.unique(),
        irfold0_avg_df.solver_num_variables,
        label="IRFold0",
    )
    plt.plot(
        irfold1_res_df.seq_len.unique(),
        irfold1_avg_df.solver_num_variables,
        label="IRFold1",
    )
    plt.legend()
    plt.grid()
    plt.ylabel("Mean Num. Solver Variables")
    plt.xlabel("Sequence Length")

    plt.subplot(1, 2, 2)
    plt.plot(
        irfold0_res_df.seq_len.unique(),
        irfold0_avg_df.solver_num_constraints,
        label="IRFold0",
    )
    plt.plot(
        irfold1_res_df.seq_len.unique(),
        irfold1_avg_df.solver_num_constraints,
        label="IRFold1",
    )

    plt.legend()
    plt.grid()
    plt.ylabel("Mean Num. Solver Constraints")
    plt.xlabel("Sequence Length")
    plt.tight_layout()
    plt.savefig(f"{DATA_DIR}/experiment_2_solver_performance_comparison.png")
    plt.show()

    plt.rcParams["figure.figsize"] = (6, 4)

    # MFEs of final solutions
    mfe_data = [
        irfold0_res_df.solution_mfe,
        irfold1_res_df.solution_mfe,
        rnalib_res_df.solution_mfe,
    ]
    fig, ax = plt.subplots()

    ax.boxplot(mfe_data, labels=["IRFold0", "IRFold1", "RNAlib"])

    ax.title.set_text("MFE of Final Solution")
    plt.ylabel("MFE (kcal/mol)")
    plt.tight_layout()
    plt.savefig(f"{DATA_DIR}/experiment_2_mfe_comparison.png")
    plt.show()

    # Compare IRFold versions' MFE's as sequence length increases
    plt.plot(
        irfold0_res_df.seq_len.unique(),
        irfold0_avg_df.solution_mfe,
        label="IRFold0",
    )
    plt.plot(
        irfold1_res_df.seq_len.unique(),
        irfold1_avg_df.solution_mfe,
        label="IRFold1",
    )

    plt.legend()
    plt.grid()
    plt.ylabel("Mean MFE of Final Solution")
    plt.xlabel("Sequence Length")
    plt.tight_layout()
    plt.savefig(f"{DATA_DIR}/experiment_2_mfe_over_time_comparison.png")
    plt.show()