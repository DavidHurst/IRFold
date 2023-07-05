import re
import sys
from pathlib import Path

import RNA

sys.path.append(str(Path(__file__).resolve().parents[1]))

from irfold import IRFoldVal2, IRFoldCorX

DATA_DIR = Path(__file__).parent.parent / "data"
BP_RNA_1_M_BP_SEQ_FILES_DIR = DATA_DIR / "bpRNA-1m" / "bpseqFiles"
BP_RNA_1_M_FASTA_FILES_DIR = DATA_DIR / "bpRNA-1m" / "fastaFiles"

UNPAIRED_FLAG = (
    -1
)  # BP Seq file format uses 1-based indexing for sequences, we use 0-based so -1 means unpaired here not 0 as in file


def evaluate_prediction(true_idx_pairs, true_unpaired_idxs, pred):
    true_positives = 0
    false_negatives = 0
    false_positives = 0

    for idx in true_idx_pairs:
        # True positive: pair is in the in prediction and in the true secondary structure
        if (pred[idx[0]] == "(" and pred[idx[1]] == ")") or (
            pred[idx[0]] == ")" and pred[idx[1]] == "("
        ):
            true_positives += 1
        else:
            # False negative pair is not in the prediction is in the true secondary structure
            false_negatives += 1

    # False positive: pair is in the prediction but is not in the true secondary structure
    for idx in true_unpaired_idxs:
        if pred[idx] == "(" or pred[idx] == ")":
            false_positives += 1

    return true_positives, false_positives, false_negatives


def calc_metrics(true_positives, false_positives, false_negatives):
    sensitivity = true_positives / (true_positives + false_negatives)
    ppv = true_positives / (true_positives + false_positives)
    f1 = (2 * sensitivity * ppv) / (sensitivity + ppv)

    return sensitivity, ppv, f1


if __name__ == "__main__":
    for fasta_file_path in BP_RNA_1_M_FASTA_FILES_DIR.glob("*.fasta"):
        if "CRW" in str(fasta_file_path) and "30278" in str(fasta_file_path):
            with open(fasta_file_path, "r") as fasta_file:
                seq = fasta_file.readlines()[1].strip()
                seq_len = len(seq)
                # Might be tildes in sequences

            seq_bp_seq_file = list(BP_RNA_1_M_BP_SEQ_FILES_DIR.glob("*CRW*30278*"))[0]

            with open(seq_bp_seq_file, "r") as bp_seq_file:
                pair_lines = [line.strip() for line in bp_seq_file.readlines()][2:]

            pairs_idxs = [
                (
                    int(re.findall(r"-?\d+\.?\d*", pair_line)[0]) - 1,
                    int(re.findall(r"-?\d+\.?\d*", pair_line)[-1]) - 1,
                )
                for pair_line in pair_lines
            ]
            pairs_idxs = [sorted(p) for p in pairs_idxs]
            sorted_pair_idxs = sorted(pairs_idxs, key=lambda tup: tup[0])

            pairs_idxs_no_dups = []
            unpaired_indices = []
            for true_pair in sorted_pair_idxs:
                if (
                    true_pair not in pairs_idxs_no_dups
                    and UNPAIRED_FLAG not in true_pair
                ):
                    pairs_idxs_no_dups.append(true_pair)
                if UNPAIRED_FLAG in true_pair:
                    unpaired_indices.append(int(true_pair[1]))

            true_num_pairs = len(pairs_idxs_no_dups)

            fold_params = {
                "sequence": seq,
                "min_len": 2,
                "max_len": seq_len,
                "max_gap": seq_len - 1,
                "mismatches": 0,
                "out_dir": DATA_DIR,
                "save_performance": False,
            }

            # Get predicted secondary structures from models
            irfold_val2_pred, _ = IRFoldVal2.fold(**fold_params)

            fold_params.update({"max_n_tuple_sz_to_correct": 2})
            irfold_cor2_pred, _ = IRFoldCorX.fold(**fold_params)

            fold_params.update({"max_n_tuple_sz_to_correct": 3})
            irfold_cor3_pred, _ = IRFoldCorX.fold(**fold_params)

            rnalib_pred, _ = RNA.fold(seq, "")

            print(f"Seq. len         : {seq_len}")
            print(f"True num. pairs  : {true_num_pairs}")

            for model_pred, model_name in zip(
                [irfold_val2_pred, irfold_cor2_pred, irfold_cor3_pred, rnalib_pred],
                ["IRFoldVal2", "IRFoldCor2", "IRFoldCor3", "RNAlib"],
            ):
                print("-" * 40)
                print(f"{model_name}:")
                print(f"  Pred     : {model_pred}")

                TPs, FPs, FNs = evaluate_prediction(
                    pairs_idxs_no_dups, unpaired_indices, model_pred
                )
                sens, ppv, f1 = calc_metrics(TPs, FPs, FNs)
                for metric, metric_name in zip(
                    [sens, ppv, f1], ["Sensitivity", "PPV", "F1"]
                ):
                    print(f"  {metric_name.ljust(15)}: {metric:.2f}")
