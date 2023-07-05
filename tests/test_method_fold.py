import random

import pytest
from pathlib import Path

from irfold import (
    IRFoldBase,
    IRFoldVal1,
    IRFoldVal2,
    IRFoldCor2,
    IRFoldCor3,
    IRFoldCorX,
)


# Test that when a sub-child calls its own function, that is called and not a parent's


@pytest.mark.parametrize(
    "ir_fold_variant",
    [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2, IRFoldCor3, IRFoldCorX],
)
def test_output_not_none(
    ir_fold_variant,
    sequence,
    sequence_length,
    sequence_name,
    data_dir,
    all_irs,
):
    secondary_structure_pred, obj_fn_value = ir_fold_variant.fold(
        sequence,
        2,
        sequence_length,
        sequence_length - 1,
        sequence_name,
        out_dir=data_dir,
    )

    assert secondary_structure_pred is not None
    assert obj_fn_value is not None


@pytest.mark.parametrize(
    "ir_fold_variant",
    [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2, IRFoldCor3, IRFoldCorX],
)
def test_output_type_correct(
    ir_fold_variant,
    sequence,
    sequence_length,
    sequence_name,
    data_dir,
    all_irs,
):
    secondary_structure_pred, obj_fn_value = ir_fold_variant.fold(
        sequence,
        2,
        sequence_length,
        sequence_length - 1,
        sequence_name,
        out_dir=data_dir,
    )

    assert isinstance(secondary_structure_pred, str)
    assert isinstance(obj_fn_value, float)


@pytest.mark.parametrize(
    "ir_fold_variant",
    [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2, IRFoldCor3, IRFoldCorX],
)
def test_balanced_brackets(ir_fold_variant, data_dir):
    seq_len = 40
    for _ in range(10):
        seq = "".join(random.choice("ACGU") for _ in range(seq_len))
        secondary_structure_pred, _ = ir_fold_variant.fold(
            seq, 2, seq_len, seq_len - 1, "test_seq", out_dir=data_dir
        )

        assert secondary_structure_pred.count("(") == secondary_structure_pred.count(
            ")"
        )


@pytest.mark.parametrize(
    "ir_fold_variant",
    [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2, IRFoldCor3, IRFoldCorX],
)
def test_ir_fold_variant_performance_written_to_file(
    ir_fold_variant,
    sequence,
    sequence_length,
    sequence_name,
    data_dir,
    all_irs,
):
    _, _ = ir_fold_variant.fold(
        sequence,
        2,
        sequence_length,
        sequence_length - 1,
        sequence_name,
        save_performance=True,
        out_dir=data_dir,
    )

    # Performance file for IRFold variant created
    performance_file_path = (
        Path(data_dir) / f"{ir_fold_variant.__name__}_performance.csv"
    )
    assert performance_file_path.exists()

    # Should only have one entry written to the file (2 lines including column titles)
    with open(performance_file_path, "r") as file:
        lines = file.readlines()
    assert len(lines) == 2

    _, _ = ir_fold_variant.fold(
        sequence,
        2,
        sequence_length,
        sequence_length - 1,
        sequence_name,
        save_performance=True,
        out_dir=data_dir,
    )

    # Should now have two entries written to the file (3 lines including column titles)
    with open(performance_file_path, "r") as file:
        lines = file.readlines()

    assert len(lines) == 3
