from irfold.util import ir_pair_match_same_bases


def test_irs_matching_same_bases(matching_same_bases_ir_pair):
    ir_a = matching_same_bases_ir_pair[0]
    ir_b = matching_same_bases_ir_pair[1]

    assert ir_pair_match_same_bases(ir_a, ir_b) == True


def test_irs_not_matching_same_bases(not_matching_same_bases_ir_pair):
    ir_a = not_matching_same_bases_ir_pair[0]
    ir_b = not_matching_same_bases_ir_pair[1]

    assert ir_pair_match_same_bases(ir_a, ir_b) == False
