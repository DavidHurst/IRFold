import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

DATA_DIR = str(Path(__file__).parent.parent / "data")

plt.rcParams["figure.figsize"] = (6, 4)

if __name__ == "__main__":
    triplets_results_df = pd.read_csv(
        f"{DATA_DIR}/experiment_2/experiment_2_results_ir_pair_validation_validating_quadruplets.csv"
    )
    quadruplets_results_df = pd.read_csv(
        f"{DATA_DIR}/experiment_2/experiment_2_results_ir_pair_validation_validating_quadruplets.csv"
    )

    # Show that only IR triplets and quadruplets containing invalid loop forming IR pairs are themselves invalid
