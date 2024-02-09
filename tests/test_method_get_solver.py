from ortools.linear_solver.pywraplp import Solver

import pytest

from irfold import IRfold
from irfold.util import ir_has_valid_gap_size


@pytest.mark.parametrize(
    "ir_fold_variant",
    [
        IRfold,
    ],
)
def test_not_none(
    ir_fold_variant, all_irs, sequence, sequence_length, sequence_name, data_dir
):
    all_irs = list(all_irs)
    solver, variables = ir_fold_variant._get_ilp_model(
        all_irs,
        sequence_length,
        sequence,
        data_dir,
        sequence_name,
    )

    assert solver is not None
    assert variables is not None

    assert isinstance(solver, Solver)
    assert isinstance(variables, list)


@pytest.mark.parametrize(
    "ir_fold_variant, variable_names",
    [
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
    _, variables = ir_fold_variant._get_ilp_model(
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
    _, variables = ir_fold_variant._get_ilp_model(
        all_irs,
        sequence_length,
        sequence,
        data_dir,
        sequence_name,
    )

    assert [v.name() for v in variables] == variable_names


def test_correct_number_of_constraints_generated(
    all_irs,
    sequence,
    sequence_length,
    sequence_name,
    data_dir,
    all_co_located_ir_pairs,
    all_partially_nested_ir_pairs,
    all_invalid_gap_size_irs,
):
    solver, variables = IRfold._get_ilp_model(
        list(all_irs),
        sequence_length,
        sequence,
        data_dir,
        sequence_name,
    )

    # There should be a constraint for all co-located IR pairs and all paritally nested IR pairs
    all_invalidly_relatively_positioned_ir_pairs = list(all_co_located_ir_pairs) + list(
        all_partially_nested_ir_pairs
    )

    # Remove pair with invalid gap size IRs as they will be removed before constraint generations
    invalidly_relatively_positioned_ir_pairs_with_valid_gap_size = []
    for (ir_a, ir_b) in all_invalidly_relatively_positioned_ir_pairs:
        if ir_has_valid_gap_size(ir_a) and ir_has_valid_gap_size(ir_b):
            invalidly_relatively_positioned_ir_pairs_with_valid_gap_size.append((ir_a, ir_b))

    num_constraints = len(solver.constraints())

    assert num_constraints == len(invalidly_relatively_positioned_ir_pairs_with_valid_gap_size)
