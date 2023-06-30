import pytest
from irfold import IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2, IRFoldCor3
from irfold.util import calc_free_energy


@pytest.mark.parametrize("irfold", [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2])
def test_objective_function_does_not_exactly_match_true_free_energy(
    irfold, find_irs_params, rna_seq_30_bases_19_irs, data_dir
):
    seq = rna_seq_30_bases_19_irs[0]
    find_irs_params["sequence"] = seq
    secondary_struct_pred, obj_fn_value = irfold.fold(**find_irs_params)

    true_free_energy = calc_free_energy(
        secondary_struct_pred, seq, data_dir, "test_irfold_obj_fn"
    )

    assert true_free_energy != obj_fn_value


def test_objective_function_more_correct_with_increasing_model_version(
    find_irs_params, rna_seq_30_bases_19_irs, data_dir
):
    seq = rna_seq_30_bases_19_irs[0]
    find_irs_params["sequence"] = seq

    model_outputs = [
        model.fold(**find_irs_params)
        for model in [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2, IRFoldCor3]
    ]

    true_free_energies = [
        calc_free_energy(struct_pred, seq, data_dir, "test_irfold_obj_fn")
        for struct_pred, _ in model_outputs
    ]

    for i in range(len(model_outputs) - 1):
        preceding_model_obj_fn = model_outputs[i][1]
        preceding_model_true_free_energy = true_free_energies[i]
        preceding_model_obj_fn_error = abs(
            preceding_model_obj_fn - preceding_model_true_free_energy
        )

        succeeding_model_obj_fn = model_outputs[i + 1][1]
        succeeding_model_true_free_energy = true_free_energies[i + 1]
        succeeding_model_obj_fn_error = abs(
            succeeding_model_obj_fn - succeeding_model_true_free_energy
        )

        assert preceding_model_obj_fn_error >= succeeding_model_obj_fn_error
