from irfold import IRFold
import random
import RNA

if __name__ == "__main__":
    # random.seed(4)

    data_dir = "../data"
    seq_len = 15
    seq = "".join(random.choice("ACGU") for _ in range(seq_len))

    print(f"Seq. length: {seq_len}")
    print(f"Seq.       : {seq}")

    # ir_fold = (data_dir)
    our_secondary_structure, our_mfe = IRFold.fold(
        sequence=seq,
        min_len=2,
        max_len=seq_len,
        max_gap=seq_len - 1,
        mismatches=0,
    )

    print("Our Solution".center(50, "="))
    print(f"Dot Bracket: {our_secondary_structure}")
    print(f"MFE        : {our_mfe:.4f}")

    rna_secondary_structure, rna_mfe = RNA.fold(seq, "")
    print("RNAlib Solution".center(50, "="))
    print(f"Dot Bracket: {rna_secondary_structure}")
    print(f"MFE        : {rna_mfe:.4f}")
