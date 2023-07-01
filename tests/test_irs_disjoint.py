from irfold.util import ir_pair_disjoint, irs_disjoint


def test_method_ir_pair_disjoint_disjoint_ir_pair(disjoint_ir_pair):
    ir_a = disjoint_ir_pair[0]
    ir_b = disjoint_ir_pair[1]

    assert ir_pair_disjoint(ir_a, ir_b) == True


def test_method_ir_pair_disjoint_not_disjoint_ir_pair(wholly_nested_ir_pair):
    ir_a = wholly_nested_ir_pair[0]
    ir_b = wholly_nested_ir_pair[1]

    assert ir_pair_disjoint(ir_a, ir_b) == False


def test_method_irs_disjoint_disjoint_ir_pair(disjoint_ir_pair):
    assert irs_disjoint(disjoint_ir_pair) == True


def test_method_irs_disjoint_not_disjoint_ir_pair(wholly_nested_ir_pair):
    assert irs_disjoint(wholly_nested_ir_pair) == False