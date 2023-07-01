from pathlib import Path

import pytest

from irfold import (
    IRFoldBase,
    IRFoldVal1,
    IRFoldVal2,
    IRFoldCor2,
    IRFoldCor3,
    IRFoldCorX,
)


# ToDo: Write test for not having iupacpal compiled


@pytest.mark.parametrize(
    "ir_fold_variant",
    [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2, IRFoldCor3, IRFoldCorX],
)
def test_not_none(ir_fold_variant, sequence, sequence_length, sequence_name, data_dir):
    assert (
        ir_fold_variant.find_irs(
            sequence,
            2,
            sequence_length,
            sequence_length - 1,
            seq_name=sequence_name,
            out_dir=data_dir,
        )
        is not None
    )


@pytest.mark.parametrize(
    "ir_fold_variant",
    [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2, IRFoldCor3, IRFoldCorX],
)
def test_number_of_irs_found(
    ir_fold_variant,
    sequence,
    sequence_length,
    sequence_name,
    data_dir,
    all_irs,
):
    found_irs = ir_fold_variant.find_irs(
        sequence,
        2,
        sequence_length,
        sequence_length - 1,
        seq_name=sequence_name,
        out_dir=data_dir,
    )

    assert len(found_irs) == len(all_irs)


@pytest.mark.parametrize(
    "ir_fold_variant",
    [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2, IRFoldCor3, IRFoldCorX],
)
def test_found_irs_output_file_created(
    ir_fold_variant, sequence, sequence_length, sequence_name, data_dir
):
    _ = ir_fold_variant.find_irs(
        sequence,
        2,
        sequence_length,
        sequence_length - 1,
        seq_name=sequence_name,
        out_dir=data_dir,
    )

    # Sequence file is created
    assert (Path(data_dir) / f"{sequence_name}.fasta").exists()

    # Found IRs file is created
    assert (Path(data_dir) / f"{sequence_name}_found_irs.txt").exists()

    # Correct sequence is written to file
    with open(
        str(Path(data_dir).resolve() / f"{sequence_name}.fasta"), "r"
    ) as seq_file:
        written_seq = seq_file.readlines()[1]
    assert written_seq == sequence
