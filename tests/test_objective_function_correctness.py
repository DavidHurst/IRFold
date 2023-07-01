from irfold import (
    IRFoldBase,
    IRFoldVal1,
    IRFoldVal2,
    IRFoldCor2,
    IRFoldCor3,
)
from irfold.util import calc_free_energy


def test_objective_function_more_correct_with_increasing_model_version(
    sequence,
    sequence_length,
    sequence_name,
    data_dir,
    all_irs,
):
    model_outputs = [
        ir_fold_variant.fold(
            sequence, 2, sequence_length, sequence_name, out_dir=data_dir
        )
        for ir_fold_variant in [
            IRFoldBase,
            IRFoldVal1,
            IRFoldVal2,
            IRFoldCor2,
            IRFoldCor3,
        ]
    ]

    true_free_energies = [
        calc_free_energy(struct_pred, sequence, data_dir, sequence_name)
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
