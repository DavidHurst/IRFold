from irfold.util import ir_pair_forms_valid_loop, irs_to_dot_bracket, calc_free_energy


# Note: Not checking if IRs have valid gap sizes, probably should check fixtures for this


def test_ir_pairs_forming_valid_loop_wholly_nested(wholly_nested_ir_pair):
    ir_a = wholly_nested_ir_pair[0]
    ir_b = wholly_nested_ir_pair[1]

    assert ir_pair_forms_valid_loop(ir_a, ir_b) == True


def test_ir_pairs_forming_valid_loop_entirely_disjoint(entirely_disjoint_ir_pair):
    ir_a = entirely_disjoint_ir_pair[0]
    ir_b = entirely_disjoint_ir_pair[1]

    assert ir_pair_forms_valid_loop(ir_a, ir_b) == True


def test_ir_pairs_forming_valid_loop_valid_num_bases_in_intersection(
    valid_num_bases_in_pair_intersection_ir_pair,
):
    ir_a = valid_num_bases_in_pair_intersection_ir_pair[0]
    ir_b = valid_num_bases_in_pair_intersection_ir_pair[1]

    assert ir_pair_forms_valid_loop(ir_a, ir_b) == True


def test_ir_pairs_forming_invalid_loop_invalid_num_bases_in_intersection(
    invalid_num_bases_in_pair_intersection_ir_pair,
):
    ir_a = invalid_num_bases_in_pair_intersection_ir_pair[0]
    ir_b = invalid_num_bases_in_pair_intersection_ir_pair[1]

    assert ir_pair_forms_valid_loop(ir_a, ir_b) == False

# ToDo: Add test calculating free energy of invalid loop forming IR pair, should be close to 100,000
