from irfold import IRFoldVal2


if __name__ == "__main__":
    seq_len = 50
    seq = 'GGGUGGCGCCAAUGUGAGACCCCGUACAGUCGCUUGAACAUGAGGUCGAA'

    folded_sequence, obj_fn_value = IRFoldVal2.fold(seq)

    print(f'Sequence: {seq}')
    print(f'Sequence length: {seq_len}')

    print(f'Folded sequence (Dot Bracket) : {folded_sequence}')
    print(f'Objective Function Final Value: {obj_fn_value}')