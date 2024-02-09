from irfold import IRfold
from irfold.util.helper_functions import calc_free_energy

if __name__ == "__main__":

    seq = "UGAUGACUUAUGCUUAACCAAAGCACGGCA"

    folded_sequence, obj_fn_value = IRfold.fold(seq, show_prog=True)

    print(f"Sequence length               : {len(seq)}")
    print(f"Sequence                      : {seq}")
    print(f"Folded sequence (Dot bracket) : {folded_sequence}")
    print(f"Objective function final value: {obj_fn_value:.4f}")
    print(
        f"Pred. free energy (ViennaRNA) : {calc_free_energy(folded_sequence, seq, './'):.4f}"
    )
