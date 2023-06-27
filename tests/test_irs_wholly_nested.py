from irfold.util import ir_pair_wholly_nested


def test_irs_wholly_nested(wholly_nested_ir_pair):
    ir_a = wholly_nested_ir_pair[0]
    ir_b = wholly_nested_ir_pair[1]

    assert ir_pair_wholly_nested(ir_a, ir_b) == True


def test_irs_not_wholly_nested(entirely_disjoint_ir_pair):
    ir_a = entirely_disjoint_ir_pair[0]
    ir_b = entirely_disjoint_ir_pair[1]

    assert ir_pair_wholly_nested(ir_a, ir_b) == False
