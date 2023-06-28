from irfold.util import ir_pair_forms_valid_loop, calc_free_energy, irs_to_dot_bracket


# Note: These tests don't check if IRs have valid gap sizes, probably should check fixtures for this


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


def test_ir_pairs_forming_invalid_loop_free_energy_proof(
    invalid_num_bases_in_pair_intersection_ir_pair, rna_seq_30_bases_19_irs, data_dir
):
    """An invalid loop forming will be detected by RNAlib's FE calculation.
    An invalid loop is assigned "infinity" free energy where infinity is 100k
    """
    seq = rna_seq_30_bases_19_irs[0]
    seq_len = len(seq)

    ir_a = invalid_num_bases_in_pair_intersection_ir_pair[0]
    ir_b = invalid_num_bases_in_pair_intersection_ir_pair[1]

    pairs_db_repr = irs_to_dot_bracket([ir_a, ir_b], seq_len)

    pairs_free_energy = calc_free_energy(pairs_db_repr, seq, data_dir)

    assert (
        pairs_free_energy >= 90_000
    )  # Test for 90k as valid loops will bring total FE down but not by 10k


def test_ir_pairs_forming_valid_loop_free_energy_proof(
    valid_num_bases_in_pair_intersection_ir_pair, rna_seq_30_bases_19_irs, data_dir
):
    """An invalid loop forming will be detected by RNAlib's FE calculation.
    An invalid loop is assigned "infinity" free energy where infinity is 100k
    """
    seq = rna_seq_30_bases_19_irs[0]
    seq_len = len(seq)

    ir_a = valid_num_bases_in_pair_intersection_ir_pair[0]
    ir_b = valid_num_bases_in_pair_intersection_ir_pair[1]

    pairs_db_repr = irs_to_dot_bracket([ir_a, ir_b], seq_len)

    pairs_free_energy = calc_free_energy(pairs_db_repr, seq, data_dir)

    assert (
        pairs_free_energy <= 90_000
    )  # Test for 90k as valid loops will bring total FE down but not by 10k
