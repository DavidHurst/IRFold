from ortools.sat.python.cp_model import CpModel, IntVar

import pytest

from irfold import (
    IRFoldBase,
    IRfold,
)


@pytest.mark.parametrize(
    "ir_fold_variant",
    [
        IRFoldBase,
        IRfold,
    ],
)
def test_not_none(
    ir_fold_variant, all_irs, sequence, sequence_length, sequence_name, data_dir
):
    all_irs = list(all_irs)
    solver, variables = ir_fold_variant._get_cp_model(
        all_irs,
        sequence_length,
        sequence,
        data_dir,
        sequence_name,
    )

    assert solver is not None
    assert variables is not None

    assert isinstance(solver, CpModel)
    assert isinstance(variables, list)

    for var in variables:
        assert isinstance(var, IntVar)


@pytest.mark.parametrize(
    "ir_fold_variant, variable_names",
    [
        (IRFoldBase, pytest.lazy_fixture("all_irs_names")),
        (IRfold, pytest.lazy_fixture("ir_indicator_variables_names")),
    ],
)
def test_number_of_variables_generated(
    ir_fold_variant,
    variable_names,
    all_irs,
    sequence,
    sequence_length,
    sequence_name,
    data_dir,
):
    all_irs = list(all_irs)
    _, variables = ir_fold_variant._get_cp_model(
        all_irs,
        sequence_length,
        sequence,
        data_dir,
        sequence_name,
    )

    assert len(variables) == len(variable_names)


@pytest.mark.parametrize(
    "ir_fold_variant, variable_names",
    [
        (IRFoldBase, pytest.lazy_fixture("all_irs_names")),
        (IRfold, pytest.lazy_fixture("ir_indicator_variables_names")),
    ],
)
def test_variables_for_correct_variables_generated(
    ir_fold_variant,
    variable_names,
    all_irs,
    sequence,
    sequence_length,
    sequence_name,
    data_dir,
):
    print("\n", "=" * 60)
    all_irs = list(all_irs)
    _, variables = ir_fold_variant._get_cp_model(
        all_irs,
        sequence_length,
        sequence,
        data_dir,
        sequence_name,
    )

    assert [v.Name() for v in variables] == variable_names
