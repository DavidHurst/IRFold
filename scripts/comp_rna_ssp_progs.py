import random
import RNA
import sys

from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
DATA_DIR = str(Path(__file__).parent.parent / 'data')

from ir_fold import IRFold, IRFoldA

if __name__ == "__main__":
    # random.seed(4)

    seq_len = 35
    seq = "".join(random.choice("ACGU") for _ in range(seq_len))

    print(f"Seq. length: {seq_len}")
    print(f"Seq.       : {seq}")

    irfold0_secondary_structure, irfold0_mfe = IRFold.fold(
        sequence=seq,
        min_len=2,
        max_len=seq_len,
        max_gap=seq_len - 1,
        mismatches=0,
        out_dir=DATA_DIR,
    )

    irfold1_secondary_structure, irfold1_mfe = IRFoldA.fold(
        sequence=seq,
        min_len=2,
        max_len=seq_len,
        max_gap=seq_len - 1,
        mismatches=0,
        out_dir=DATA_DIR,
    )

    print(f"{IRFold.__name__} Solution".center(50, "="))
    print(f"Dot Bracket: {irfold0_secondary_structure}")
    print(f"MFE        : {irfold0_mfe:.4f}")

    print(f"{IRFoldA.__name__} Solution".center(50, "="))
    print(f"Dot Bracket: {irfold1_secondary_structure}")
    print(f"MFE        : {irfold1_mfe:.4f}")

    rna_secondary_structure, rna_mfe = RNA.fold(seq, "")
    print("RNAlib Solution".center(50, "="))
    print(f"Dot Bracket: {rna_secondary_structure}")
    print(f"MFE        : {rna_mfe:.4f}")
