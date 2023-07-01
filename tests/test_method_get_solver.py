from ortools.sat.python.cp_model import CpModel, IntVar

import pytest

from irfold import IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2


@pytest.mark.parametrize("irfold", [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2])
def test_solver_generated_0_irs(data_dir, irfold, rna_seq_15_bases_0_irs):
    seq = rna_seq_15_bases_0_irs[0]
    seq_len = len(seq)

    solver, variables = irfold.get_solver(
        [],
        seq_len,
        seq,
        data_dir,
        "test_get_solver_seq",
    )

    assert solver is not None
    assert variables is not None

    assert isinstance(solver, CpModel)
    assert isinstance(variables, list)

    assert len(variables) == 0


@pytest.mark.parametrize("irfold", [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2])
def test_solver_generated_3_irs(data_dir, irfold, list_of_irs, rna_seq_15_bases_3_irs):
    seq = rna_seq_15_bases_3_irs[0]
    seq_len = len(seq)

    solver, variables = irfold.get_solver(
        list_of_irs,
        seq_len,
        seq,
        data_dir,
        "test_get_solver_seq",
    )

    assert solver is not None
    assert variables is not None

    assert isinstance(solver, CpModel)
    assert isinstance(variables, list)

    for var in variables:
        assert isinstance(var, IntVar)
