import pytest

from irfold.util import ir_pair_disjoint


def test_irs_disjoint(
    list_of_ir_pairs_entirely_disjoint, list_of_ir_pairs_wholly_nested
):
    for ir_a, ir_b in list_of_ir_pairs_entirely_disjoint:
        assert ir_pair_disjoint(ir_a, ir_b) == True

    for ir_a, ir_b in list_of_ir_pairs_wholly_nested:
        assert ir_pair_disjoint(ir_a, ir_b) == False
