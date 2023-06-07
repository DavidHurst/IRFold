from irfold import IRFold
import random
import RNA

if __name__ == "__main__":
    # random.seed(4)

    data_dir = "../data"
    seq_len = 25
    seq = "".join(random.choice("ACGU") for _ in range(seq_len))

    print(f"Seq. length: {seq_len}")
    print(f"Seq.       : {seq}")

    # Check if FE(IR_i) + FE(IR_j) == FE(IR_i + IR_j)
    find_irs_params = {
        "sequence": seq,
        "min_len": 2,
        "max_len": seq_len,
        "max_gap": seq_len - 1,
        "mismatches": 0,
        "out_dir": data_dir,
    }
    found_irs = IRFold.find_irs(**find_irs_params)

    compatible_irs = []
    for ir_1, ir_2 in zip(found_irs[1:], found_irs[:-1]):
        if not IRFold.irs_incompatible(ir_1, ir_2):
            compatible_irs.append(ir_1)
            compatible_irs.append(ir_2)

        if len(compatible_irs) >= 2:
            break

    ir_1, ir_2 = compatible_irs[0], compatible_irs[1]

    ir_1_db_repr = IRFold.irs_to_dot_bracket([ir_1], seq_len)
    ir_2_db_repr = IRFold.irs_to_dot_bracket([ir_2], seq_len)
    ir_1_2_db_repr = IRFold.irs_to_dot_bracket([ir_1, ir_2], seq_len)

    ir_1_fe = IRFold.calc_free_energy(ir_1_db_repr, seq, data_dir)
    ir_2_fe = IRFold.calc_free_energy(ir_2_db_repr, seq, data_dir)
    ir_1_2_fe = IRFold.calc_free_energy(ir_1_2_db_repr, seq, data_dir)

    print(f'FE(IR#1)           : {ir_1_fe:.3f}')
    print(f'FE(IR#2)           : {ir_2_fe:.3f}')
    print(f'FE(IR#1) + FE(IR#2): {ir_1_fe + ir_2_fe:.3f}')
    print(f'FE(IR#1 u IR#2)    : {ir_1_2_fe:.3f}')

    # ir_fold = (data_dir)
    # our_secondary_structure, our_mfe = IRFold.fold(
    #     sequence=seq,
    #     min_len=2,
    #     max_len=seq_len,
    #     max_gap=seq_len - 1,
    #     mismatches=0,
    #     out_dir="../data",
    # )
    #
    # print("Our Solution".center(50, "="))
    # print(f"Dot Bracket: {our_secondary_structure}")
    # print(f"MFE        : {our_mfe:.4f}")
    #
    # rna_secondary_structure, rna_mfe = RNA.fold(seq, "")
    # print("RNAlib Solution".center(50, "="))
    # print(f"Dot Bracket: {rna_secondary_structure}")
    # print(f"MFE        : {rna_mfe:.4f}")
