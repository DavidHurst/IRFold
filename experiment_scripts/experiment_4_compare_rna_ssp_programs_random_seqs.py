import csv
import random
import RNA

import sys

from pathlib import Path
from tqdm import tqdm

sys.path.append(str(Path(__file__).resolve().parents[1]))
DATA_DIR = str(Path(__file__).parent.parent / "data")

from irfold import IRFoldBase, IRFoldVal1, IRFoldCor2, IRFoldVal2, IRFoldCor3

if __name__ == "__main__":
    rnalib_performance_file_path = (Path(DATA_DIR) / "RNAlib_performance.csv").resolve()
    with open(str(rnalib_performance_file_path), "w") as perf_file:
        writer = csv.writer(perf_file)
        writer.writerow(["dot_bracket_repr", "solution_mfe", "seq_len"])
    # random.seed(182)

    n_runs_per_seq_length = 10
    for seq_len in tqdm(range(7, 50), desc=f"Running trials"):
        for _ in range(n_runs_per_seq_length):
            seq = "".join(random.choice("ACGU") for _ in range(seq_len))
            seq_name = "random_seq_for_ssp_program_comparison"

            fold_params = {
                "sequence": seq,
                "out_dir": DATA_DIR,
                "seq_name": seq_name,
                "save_performance": True,
            }

            irfold_base_secondary_structure, irfold_base_obj_fn_val = IRFoldBase.fold(
                **fold_params
            )

            irfold_val1_secondary_structure, irfold_val1_obj_fn_val = IRFoldVal1.fold(
                **fold_params
            )

            irfold_val2_secondary_structure, irfold_val2_obj_fn_val = IRFoldVal2.fold(
                **fold_params
            )

            irfold_cor2_secondary_structure, irfold_cor2_obj_fn_val = IRFoldCor2.fold(
                **fold_params
            )

            irfold_cor3_secondary_structure, irfold_cor3_obj_fn_val = IRFoldCor3.fold(
                **fold_params
            )

            # ToDo: Also compare with RNAstrucutre and IRFoldCorX variations

            rnalib_secondary_structure, rnalib_mfe = RNA.fold(seq, "")

            with open(str(rnalib_performance_file_path), "a") as perf_file:
                writer = csv.writer(perf_file)
                writer.writerow(
                    [
                        rnalib_secondary_structure,
                        rnalib_mfe,
                        seq_len,
                    ]
                )
