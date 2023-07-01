from irfold.util import ir_has_valid_gap_size


def test_ir_with_invalid_gap_size(invalid_gap_size_ir):
    assert ir_has_valid_gap_size(invalid_gap_size_ir) == False


def test_ir_with_valid_gap_size(valid_gap_size_ir):
    assert ir_has_valid_gap_size(valid_gap_size_ir) == True
