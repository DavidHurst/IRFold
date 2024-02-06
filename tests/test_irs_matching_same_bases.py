from irfold.util import ir_pair_co_located, irs_co_located


def test_method_ir_pair_match_same_bases_irs_matching_same_bases(
    matching_the_same_bases_ir_pair,
):
    ir_a = matching_the_same_bases_ir_pair[0]
    ir_b = matching_the_same_bases_ir_pair[1]

    assert ir_pair_co_located(ir_a, ir_b) == True


def test_method_ir_pair_match_same_bases_irs_not_matching_same_bases(
    not_matching_the_same_bases_ir_pair,
):
    ir_a = not_matching_the_same_bases_ir_pair[0]
    ir_b = not_matching_the_same_bases_ir_pair[1]

    assert ir_pair_co_located(ir_a, ir_b) == False


def test_method_irs_match_same_bases_irs_matching_same_bases(
    matching_the_same_bases_ir_pair,
):
    assert irs_co_located(matching_the_same_bases_ir_pair) == True


def test_method_irs_match_same_bases_irs_not_matching_same_bases(
    not_matching_the_same_bases_ir_pair,
):
    assert irs_co_located(not_matching_the_same_bases_ir_pair) == False
