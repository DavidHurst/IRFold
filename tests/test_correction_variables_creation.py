import re

import pytest

from irfold import IRFoldCor2, IRFoldCor3

irs = [
    ((0, 1), (28, 29)),
    ((0, 1), (23, 24)),
    ((0, 1), (19, 20)),
    ((0, 1), (6, 7)),
    ((2, 3), (9, 10)),
    ((3, 4), (28, 29)),
    ((3, 4), (23, 24)),
    ((3, 4), (19, 20)),
    ((3, 4), (6, 7)),
    ((6, 7), (10, 11)),
    ((7, 8), (13, 14)),
    ((8, 9), (13, 14)),
    ((10, 12), (27, 29)),
]

ir_pair_correction_variable_names = [
    "ir_0_ir_4_fe_corrector",
    "ir_0_ir_6_fe_corrector",
    "ir_0_ir_7_fe_corrector",
    "ir_0_ir_10_fe_corrector",
    "ir_0_ir_11_fe_corrector",
    "ir_1_ir_4_fe_corrector",
    "ir_1_ir_5_fe_corrector",
    "ir_1_ir_7_fe_corrector",
    "ir_1_ir_10_fe_corrector",
    "ir_1_ir_11_fe_corrector",
    "ir_1_ir_12_fe_corrector",
    "ir_2_ir_4_fe_corrector",
    "ir_2_ir_5_fe_corrector",
    "ir_2_ir_6_fe_corrector",
    "ir_2_ir_10_fe_corrector",
    "ir_2_ir_11_fe_corrector",
    "ir_2_ir_12_fe_corrector",
    "ir_5_ir_11_fe_corrector",
    "ir_6_ir_10_fe_corrector",
    "ir_6_ir_11_fe_corrector",
    "ir_6_ir_12_fe_corrector",
    "ir_7_ir_10_fe_corrector",
    "ir_7_ir_11_fe_corrector",
    "ir_7_ir_12_fe_corrector",
]
ir_triplet_correction_variable_names = [
    "ir_0_ir_6_ir_12_fe_corrector",
    "ir_0_ir_7_ir_12_fe_corrector",
    "ir_0_ir_10_ir_12_fe_corrector",
    "ir_0_ir_11_ir_12_fe_corrector",
    "ir_1_ir_6_ir_12_fe_corrector",
    "ir_1_ir_10_ir_11_fe_corrector",
    "ir_2_ir_3_ir_7_fe_corrector",
    "ir_2_ir_3_ir_10_fe_corrector",
    "ir_2_ir_3_ir_11_fe_corrector",
    "ir_3_ir_4_ir_6_fe_corrector",
    "ir_3_ir_4_ir_7_fe_corrector",
    "ir_3_ir_4_ir_12_fe_corrector",
    "ir_3_ir_5_ir_6_fe_corrector",
    "ir_3_ir_5_ir_7_fe_corrector",
]


def test_correction_pair_variables_created(rna_seq_30_bases_19_irs, data_dir):
    seq = rna_seq_30_bases_19_irs[0]
    seq_len = len(seq)
    _, variables = IRFoldCor2.get_solver(
        irs, seq_len, seq, data_dir, str(rna_seq_30_bases_19_irs)
    )

    correction_vars = [var for var in variables if "corrector" in var.Name()]

    assert len(correction_vars) > 0

    for v in correction_vars:
        ir_indices = re.findall(r"-?\d+\.?\d*", v.Name())
        assert (
            len(ir_indices) == 2
        )  # Correction variables have each IRs index in the name


def test_correction_triplet_variables_created(rna_seq_30_bases_19_irs, data_dir):
    seq = rna_seq_30_bases_19_irs[0]
    seq_len = len(seq)
    _, variables = IRFoldCor3.get_solver(
        irs, seq_len, seq, data_dir, str(rna_seq_30_bases_19_irs)
    )

    correction_vars = [var for var in variables if "corrector" in var.Name()]

    for v in correction_vars:
        print(f'"{v.Name()}",')

    assert len(correction_vars) > 0

    for v in correction_vars:
        ir_indices = re.findall(r"-?\d+\.?\d*", v.Name())
        assert (
            len(ir_indices) == 2 or len(ir_indices) == 3
        )  # Correction variables have each IRs index in the name


@pytest.mark.parametrize("irfold_cor", [IRFoldCor2, IRFoldCor3])
def test_expected_pair_correction_variables_created(
    rna_seq_30_bases_19_irs, data_dir, irfold_cor
):
    seq = rna_seq_30_bases_19_irs[0]
    seq_len = len(seq)
    _, variables = irfold_cor.get_solver(
        irs, seq_len, seq, data_dir, str(rna_seq_30_bases_19_irs)
    )

    ir_pair_corrector_variable_names = [
        v.Name() for v in variables if len(re.findall(r"-?\d+\.?\d*", v.Name())) == 2
    ]

    assert set(ir_pair_corrector_variable_names) == set(
        ir_pair_correction_variable_names
    )


def test_expected_triplet_correction_variables_created(
    rna_seq_30_bases_19_irs, data_dir
):
    seq = rna_seq_30_bases_19_irs[0]
    seq_len = len(seq)
    _, variables = IRFoldCor3.get_solver(
        irs, seq_len, seq, data_dir, str(rna_seq_30_bases_19_irs)
    )

    ir_triplet_corrector_variable_names = [
        v.Name() for v in variables if len(re.findall(r"-?\d+\.?\d*", v.Name())) == 3
    ]

    assert set(ir_triplet_corrector_variable_names) == set(
        ir_triplet_correction_variable_names
    )
