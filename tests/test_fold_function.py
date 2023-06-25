import pytest
from pathlib import Path

from irfold import IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2


# Test that when a sub-child calls its own function, that is called and not a parent's


@pytest.mark.parametrize("irfold", [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2])
def test_output_not_none(irfold, find_irs_params):
    pred, obj_fn_value = irfold.fold(**find_irs_params)

    assert pred is not None
    assert obj_fn_value is not None


@pytest.mark.parametrize("irfold", [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2])
def test_output_type_correct(irfold, find_irs_params):
    pred, obj_fn_value = irfold.fold(**find_irs_params)

    assert isinstance(pred, str)
    assert isinstance(obj_fn_value, float)


@pytest.mark.parametrize("irfold", [IRFoldBase, IRFoldVal1, IRFoldVal2, IRFoldCor2])
def test_output_written_to_file(irfold, find_irs_params, data_dir):
    find_irs_params["save_performance"] = True

    _, _ = irfold.fold(**find_irs_params)

    assert (
        Path(data_dir) / f"{irfold.__name__}_performance.csv"
    ).exists()  # Performance file for IRFold variant created
