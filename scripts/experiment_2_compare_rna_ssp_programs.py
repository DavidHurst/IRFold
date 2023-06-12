import csv
import random
import RNA
import sys

from pathlib import Path
from tqdm import tqdm

sys.path.append(str(Path(__file__).resolve().parents[1]))
DATA_DIR = str(Path(__file__).parent.parent / "data")

from irfold import IRFold0, IRFold1, IRFold2

if __name__ == "__main__":
    rnalib_performance_file_path = (Path(DATA_DIR) / "RNAlib_performance.csv").resolve()
    with open(str(rnalib_performance_file_path), "w") as perf_file:
        writer = csv.writer(perf_file)
        writer.writerow(["dot_bracket_repr", "solution_mfe", "seq_len"])
    random.seed(182)

    n_runs_per_seq_length = 10
    for seq_len in tqdm(range(10, 100), desc=f"Running trials"):
        for _ in range(n_runs_per_seq_length):
            seq = "".join(random.choice("ACGU") for _ in range(seq_len))
            seq_name = "random_seq_for_ssp_program_comparison"

            # print(f"Seq. length: {seq_len}")
            # print(f"Seq.       : {seq}")

            fold_params = {
                "sequence": seq,
                "min_len": 2,
                "max_len": seq_len,
                "max_gap": seq_len - 1,
                "mismatches": 0,
                "out_dir": DATA_DIR,
                "seq_name": seq_name,
                "save_performance": True,
            }

            # irfold0_secondary_structure, irfold0_mfe = IRFold0.fold(**fold_params)

            irfold1_secondary_structure, irfold1_mfe = IRFold1.fold(**fold_params)

            # irfold2_secondary_structure, irfold2_mfe = IRFold2.fold(**fold_params)

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

            # print(f"IRFold0 Solution".center(50, "="))
            # print(f"Dot Bracket: {irfold0_secondary_structure}")
            # print(f"MFE        : {irfold0_mfe:.4f}")
            #
            # print(f"IRFold1 Solution".center(50, "="))
            # print(f"Dot Bracket: {irfold1_secondary_structure}")
            # print(f"MFE        : {irfold1_mfe:.4f}")
            #
            # print("RNAlib Solution".center(50, "="))
            # print(f"Dot Bracket: {rnalib_secondary_structure}")
            # print(f"MFE        : {rnalib_mfe:.4f}")
