from irfold.util import ir_pair_not_nested, irs_not_nested


def test_method_ir_pair_not_nested_with_not_nested_ir_pair(not_nested_ir_pair):
    ir_a = not_nested_ir_pair[0]
    ir_b = not_nested_ir_pair[1]

    assert ir_pair_not_nested(ir_a, ir_b) == True


def test_method_ir_pair_not_nested_with_nested_ir_pair(wholly_nested_ir_pair):
    ir_a = wholly_nested_ir_pair[0]
    ir_b = wholly_nested_ir_pair[1]

    assert ir_pair_not_nested(ir_a, ir_b) == False


def test_method_irs_not_nested_with_not_nested_ir_pair(not_nested_ir_pair):
    assert irs_not_nested(not_nested_ir_pair) == True


def test_method_irs_not_nested_with_nested_ir_pair(wholly_nested_ir_pair):
    assert irs_not_nested(wholly_nested_ir_pair) == False
