from irfold.util import ir_pair_disjoint


def test_irs_disjoint(entirely_disjoint_ir_pair):
    ir_a = entirely_disjoint_ir_pair[0]
    ir_b = entirely_disjoint_ir_pair[1]

    assert ir_pair_disjoint(ir_a, ir_b) == True
