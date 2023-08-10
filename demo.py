import random
import sys

from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from irfold import (
    IRFoldVal2,
)

if __name__ == "__main__":
    seq_len = 50
    seq = "".join(random.choice("ACGU") for _ in range(seq_len))

    folded_sequence, obj_fn_value = IRFoldVal2.fold(seq)

    print(f'Sequence: {seq}')
    print(f'Sequence length: {seq_len}')

    print(f'Folded sequence (Dot Bracket) : {folded_sequence}')
    print(f'Objective Function Final Value: {obj_fn_value}')