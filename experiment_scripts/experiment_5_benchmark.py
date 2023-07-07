import csv
import re
import sys
import subprocess
import time

from pathlib import Path


import RNA

sys.path.append(str(Path(__file__).resolve().parents[1]))  # Add root
RNA_STRUCTURE_DATA_TABLES_DIR_PATH = str(
    (Path(__file__).parent / "RNAstructure" / "data_tables").resolve()
)
RNA_STRUCTURE_FOLD_PATH = str(
    (Path(__file__).parent / "RNAstructure" / "Fold").resolve()
)
RNA_STRUCTURE_CT2DOT_PATH = str(
    (Path(__file__).parent / "RNAstructure" / "ct2dot").resolve()
)
TEMP_DOT_BRACKET_FILE_PATH = str(
    (Path(__file__).parent / "RNAstructure" / "temp_dot_bracket.txt").resolve()
)
TEMP_CT_FILE_PATH = str(
    (Path(__file__).parent / "RNAstructure" / "temp_ct.ct").resolve()
)


from irfold import IRFoldVal2, IRFoldCorX


DATA_DIR = Path(__file__).parent.parent / "data"
EXPERIMENT_5_DIR_PATH = DATA_DIR / "experiment_5"
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
        # ToDo: Add square brackets as TPs as RNAstructure uses them for pseudoknot identification
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
    if true_positives + false_negatives == 0:
        sensitivity = 0
    else:
        sensitivity = true_positives / (true_positives + false_negatives)

    if true_positives + false_positives == 0:
        ppv = 0
    else:
        ppv = true_positives / (true_positives + false_positives)

    if ppv == 0 and sensitivity == 0:
        f1 = 0
    else:
        f1 = (2 * sensitivity * ppv) / (sensitivity + ppv)

    return sensitivity, ppv, f1


def append_performance_to_file(
    pred,
    pred_wall_time,
    true_pair_idxs,
    true_unpaired_idxs,
    seq_db_name,
    seq_db_num,
    seq_info,
    seq_len,
    performance_file_path,
):
    TPs, FPs, FNs = evaluate_prediction(true_pair_idxs, true_unpaired_idxs, pred)
    sensitivity, ppv, f1 = calc_metrics(TPs, FPs, FNs)
    for metric, metric_name in zip(
        [sensitivity, ppv, f1], ["Sensitivity", "PPV", "F1"]
    ):
        print(f"  {metric_name.ljust(15)}: {metric:.2f}")

    column_entries = [
        seq_db_name,
        seq_db_num,
        seq_info,
        seq_len,
        sensitivity,
        ppv,
        f1,
        pred_wall_time,
    ]

    with open(str(performance_file_path), "a") as perf_file:
        writer = csv.writer(perf_file)
        writer.writerow(column_entries)


def run_rnastructure(seq_fasta_file_path):
    # Clear temp files
    with open(TEMP_DOT_BRACKET_FILE_PATH) as _:
        pass
    with open(TEMP_CT_FILE_PATH) as _:
        pass

    start_time = time.monotonic()
    subprocess.run([RNA_STRUCTURE_FOLD_PATH, seq_fasta_file_path, TEMP_CT_FILE_PATH])
    wall_time = time.monotonic() - start_time

    subprocess.run(
        [
            RNA_STRUCTURE_CT2DOT_PATH,
            TEMP_CT_FILE_PATH,
            str(1),
            TEMP_DOT_BRACKET_FILE_PATH,
        ]
    )

    with open(TEMP_DOT_BRACKET_FILE_PATH, "r") as db_file:
        lines = db_file.readlines()
        pred = lines[2]

    return pred, wall_time


def run_rnafold(seq):
    start_time = time.monotonic()
    pred, _ = RNA.fold(seq, "")
    wall_time = time.monotonic() - start_time

    return pred, wall_time


def run_irfold_val2(seq):
    start_time = time.monotonic()
    pred, _ = IRFoldVal2.fold(seq, out_dir=str(EXPERIMENT_5_DIR_PATH))
    wall_time = time.monotonic() - start_time

    return pred, wall_time


def run_irfold_corx2(seq):
    start_time = time.monotonic()
    pred, _ = IRFoldCorX.fold(
        seq, max_n_tuple_sz_to_correct=2, out_dir=str(EXPERIMENT_5_DIR_PATH)
    )
    wall_time = time.monotonic() - start_time

    return pred, wall_time


def run_irfold_corx3(seq):
    start_time = time.monotonic()
    pred, _ = IRFoldCorX.fold(
        seq, max_n_tuple_sz_to_correct=3, out_dir=str(EXPERIMENT_5_DIR_PATH)
    )
    wall_time = time.monotonic() - start_time

    return pred, wall_time


if __name__ == "__main__":
    results_file_column_names = [
        "database",
        "sequence_database_number",
        "sequence_info",
        "sequence_length",
        "sensitivity",
        "ppv",
        "f1",
        "execution_wall_time_secs",
    ]
    rnafold_performance_file_path = (
        EXPERIMENT_5_DIR_PATH / "RNAfold_performance_metrics.csv"
    ).resolve()
    rnastructure_performance_file_path = (
        EXPERIMENT_5_DIR_PATH / "RNAstructure_performance_metrics.csv"
    ).resolve()
    irfold_val2_performance_file_path = (
        EXPERIMENT_5_DIR_PATH / "IRFoldVal2_performance_metrics.csv"
    ).resolve()
    irfold_corx2_performance_file_path = (
        EXPERIMENT_5_DIR_PATH / "IRFoldCorX2_performance_metrics.csv"
    ).resolve()
    irfold_corx3_performance_file_path = (
        EXPERIMENT_5_DIR_PATH / "IRFoldCorX3_performance_metrics.csv"
    ).resolve()

    for file_path in [
        rnafold_performance_file_path,
        rnastructure_performance_file_path,
        irfold_val2_performance_file_path,
        irfold_corx2_performance_file_path,
        irfold_corx3_performance_file_path,
    ]:
        with open(str(file_path), "w") as perf_file:
            writer = csv.writer(perf_file)
            writer.writerow(results_file_column_names)

    for test_seq_idx, fasta_file_path in enumerate(
        BP_RNA_1_M_FASTA_FILES_DIR.glob("*.fasta")
    ):
        fasta_file_name = str(fasta_file_path).split("/")[-1]
        database_name = fasta_file_name.split("_")[1]
        sequence_number = fasta_file_name.split("_")[2].split(".")[0]

        print(
            f"Test Sequence #{test_seq_idx}: {database_name}-{sequence_number}".center(
                70, "="
            )
        )

        if database_name not in ["SPR", "SRP"]:
            print("Not SPR or SRP databases.")
            continue

        with open(fasta_file_path, "r") as fasta_file:
            lines = fasta_file.readlines()

            if len(lines) < 2:
                print("No sequence found in fasta.")
                continue

            seq_info = lines[0].strip().split("|")[0][1:]
            seq = lines[1].strip()
            seq_len = len(seq)

        if "." in seq or "_" in seq:
            print("Unsupported notation for RNAstructure")
            continue

        if seq_len > 80:
            print("Sequence too long")
            continue

        # Get corresponding bpSeq file
        bp_seq_file_search_res = list(
            BP_RNA_1_M_BP_SEQ_FILES_DIR.glob(
                f"bpRNA_{database_name}_{sequence_number}.bpseq"
            )
        )
        if len(bp_seq_file_search_res) == 0:
            print("Could not find corresponding bpSeq file.")
            continue

        seq_bp_seq_file_path = str(bp_seq_file_search_res[0])
        bp_seq_file_name = str(seq_bp_seq_file_path).split("/")[-1]

        with open(seq_bp_seq_file_path, "r") as bp_seq_file:
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
            if true_pair not in pairs_idxs_no_dups and UNPAIRED_FLAG not in true_pair:
                pairs_idxs_no_dups.append(true_pair)
            if UNPAIRED_FLAG in true_pair:
                unpaired_indices.append(int(true_pair[1]))

        max_index_in_bp_seq_file = max(
            [max(tup) for tup in pairs_idxs_no_dups] + unpaired_indices
        )

        if max_index_in_bp_seq_file != seq_len - 1:
            print("bpSeq file and fasta file do not align.")
            continue

        print(f"Fasta file   : {fasta_file_name}")
        print(f"BPSeq file   : {bp_seq_file_name}")
        print(f"DB name      : {database_name}")
        print(f"Seq's DB num.: {sequence_number}")
        print(f"Seq. info    : {seq_info}")
        print(f"Seq. len.    : {seq_len}")
        print(f"BPSeq max idx: {max_index_in_bp_seq_file}")

        # Get predicted secondary structures from models
        for run_fold_fn, perf_file in zip(
            [
                run_rnastructure,
                run_rnafold,
                run_irfold_val2,
                run_irfold_corx2,
                run_irfold_corx3,
            ],
            [
                rnastructure_performance_file_path,
                rnafold_performance_file_path,
                irfold_val2_performance_file_path,
                irfold_corx2_performance_file_path,
                irfold_corx3_performance_file_path,
            ],
        ):
            print("-" * 70)
            fn_input = seq
            if "rnastructure" in run_fold_fn.__name__:
                fn_input = fasta_file_path

            pred, pred_wall_time = run_fold_fn(fn_input)
            append_performance_to_file(
                pred,
                pred_wall_time,
                pairs_idxs_no_dups,
                unpaired_indices,
                database_name,
                sequence_number,
                seq_info,
                seq_len,
                perf_file,
            )
