from irfold import IRFold
import random
import RNA

# random.seed(4)

data_dir = "./data"
seq_len = 35
seq = "".join(random.choice("ACGU") for _ in range(seq_len))
seq_name = "test_seq"
seq_energies_file = f"{seq_name}_energies.txt"

print(f"Seq.        : {seq}")
print(f"Seq. length : {seq_len}")

ir_fold = IRFold(data_dir)
secondary_structure, mfe = ir_fold.fold(
    sequence=seq,
    seq_name=seq_name,
    min_len=2,
    max_len=seq_len,
    max_gap=seq_len - 1,
    mismatches=0,
)

print('Our Solution'.center(50, '='))
print(f"Dot Bracket: {secondary_structure}")
print(f"MFE        : {mfe:.4f}")

out = RNA.fold(seq, '')
print("RNAlib Solution".center(50, "="))
print(f"Dot Bracket: {out[0]}")
print(f"MFE        : {out[1]:.4f}")
