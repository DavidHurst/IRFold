from irfold.util import ir_pair_partially_nested


def test_method_ir_pair_partially_nested_with_partially_nested_ir_pair(
    partially_nested_ir_pair,
):
    ir_a = partially_nested_ir_pair[0]
    ir_b = partially_nested_ir_pair[1]

    assert ir_pair_partially_nested(ir_a, ir_b) == True


def test_method_ir_pair_partially_nested_with_wholly_nested_ir_pair(
    wholly_nested_ir_pair,
):
    ir_a = wholly_nested_ir_pair[0]
    ir_b = wholly_nested_ir_pair[1]

    assert ir_pair_partially_nested(ir_a, ir_b) == False


def test_method_ir_pair_partially_nested_with_not_nested_ir_pair(not_nested_ir_pair):
    ir_a = not_nested_ir_pair[0]
    ir_b = not_nested_ir_pair[1]

    assert ir_pair_partially_nested(ir_a, ir_b) == False
