import random
import RNA
import sys

from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
DATA_DIR = str(Path(__file__).parent.parent / "data")

from irfold import IRFold0, IRFold1

if __name__ == "__main__":
    # random.seed(4)

    seq_len = 55
    seq = "".join(random.choice("ACGU") for _ in range(seq_len))
    seq_name = "random"

    print(f"Seq. length: {seq_len}")
    print(f"Seq.       : {seq}")

    irfold0_secondary_structure, irfold0_mfe = IRFold0.fold(
        sequence=seq,
        min_len=2,
        max_len=seq_len,
        max_gap=seq_len - 1,
        mismatches=0,
        out_dir=DATA_DIR,
        seq_name=seq_name,
    )

    irfold1_secondary_structure, irfold1_mfe = IRFold1.fold(
        sequence=seq,
        min_len=2,
        max_len=seq_len,
        max_gap=seq_len - 1,
        mismatches=0,
        out_dir=DATA_DIR,
        seq_name=seq_name,
    )

    rnalib_secondary_structure, rnalib_mfe = RNA.fold(seq, "")

    print(f"{IRFold0.__name__} Solution".center(50, "="))
    print(f"Dot Bracket: {irfold0_secondary_structure}")
    print(f"MFE        : {irfold0_mfe:.4f}")

    print(f"{IRFold1.__name__} Solution".center(50, "="))
    print(f"Dot Bracket: {irfold1_secondary_structure}")
    print(f"MFE        : {irfold1_mfe:.4f}")

    print("RNAlib Solution".center(50, "="))
    print(f"Dot Bracket: {rnalib_secondary_structure}")
    print(f"MFE        : {rnalib_mfe:.4f}")
