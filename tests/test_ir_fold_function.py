import random

import pytest

from irfold import IRFold0


@pytest.mark.xfail
def test_ir_free_energy_additivity_irfold0(data_dir):
    assumption_held_count = 0
    seq_lengths = [i for i in range(10, 50)]
    n_seq_lengths = len(seq_lengths)
    n_runs_per_seq_length = 10
    for seq_len in seq_lengths:
        for _ in range(n_runs_per_seq_length):
            seq = "".join(random.choice("ACGU") for _ in range(seq_len))
            seq_name = "random_seq_for_ssp_program_comparison"
        seq_len = 60
        seq = "".join(random.choice("ACGU") for _ in range(seq_len))
        seq_name = "rna_seq_15_bases_3_irs"

        fold_params = {
            "sequence": seq,
            "min_len": 2,
            "max_len": seq_len,
            "max_gap": seq_len - 1,
            "mismatches": 0,
            "out_dir": data_dir,
            "seq_name": seq_name,
            "save_performance": True,
        }
        irfold0_secondary_structure, irfold0_mfe = IRFold0.fold(**fold_params)
        mfe_of_dot_bracket = IRFold0.calc_free_energy(
            irfold0_secondary_structure, seq, data_dir, seq_name
        )

        assumption_held_count = (
            assumption_held_count + 1
            if mfe_of_dot_bracket == irfold0_mfe
            else assumption_held_count
        )

    n_assumption_holding_evaluations = n_seq_lengths * n_runs_per_seq_length
    assert assumption_held_count == n_assumption_holding_evaluations


# Test that when a subchild calls its own function, that is called and not a parent's
