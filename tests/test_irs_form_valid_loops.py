import pytest

from irfold.util import (
    ir_pair_forms_valid_loop,
    irs_form_valid_loop,
    calc_free_energy,
    irs_to_dot_bracket,
    irs_match_same_bases,
    ir_has_valid_gap_size,
)


@pytest.mark.parametrize(
    "valid_loop_forming_pair",
    [
        pytest.lazy_fixture("forms_valid_loop_ir_pair"),
        pytest.lazy_fixture("disjoint_ir_pair"),
        pytest.lazy_fixture("wholly_nested_ir_pair"),
        pytest.lazy_fixture("valid_number_of_bases_in_intersection_ir_pair"),
    ],
)
def test_method_ir_pair_forms_valid_loop_valid_loop_forming_pair(
    valid_loop_forming_pair,
):
    ir_a = valid_loop_forming_pair[0]
    ir_b = valid_loop_forming_pair[1]

    assert ir_pair_forms_valid_loop(ir_a, ir_b) == True


def test_method_ir_pair_forms_valid_loop_invalid_loop_forming_pair(
    forms_invalid_loop_ir_pair,
):
    ir_a = forms_invalid_loop_ir_pair[0]
    ir_b = forms_invalid_loop_ir_pair[1]

    assert ir_pair_forms_valid_loop(ir_a, ir_b) == False


@pytest.mark.parametrize(
    "valid_loop_forming_pair",
    [
        pytest.lazy_fixture("forms_valid_loop_ir_pair"),
        pytest.lazy_fixture("disjoint_ir_pair"),
        pytest.lazy_fixture("wholly_nested_ir_pair"),
        pytest.lazy_fixture("valid_number_of_bases_in_intersection_ir_pair"),
    ],
)
def test_method_irs_form_valid_loop_valid_loop_forming_pair(valid_loop_forming_pair):
    assert irs_form_valid_loop(valid_loop_forming_pair) == True


def test_method_irs_form_valid_loop_invalid_loop_forming_pair(
    forms_invalid_loop_ir_pair,
):
    assert irs_form_valid_loop(forms_invalid_loop_ir_pair) == False


def test_ir_pairs_forming_valid_loop_is_assigned_reasonable_free_energy_by_RNAlib(
    forms_valid_loop_ir_pair, sequence, sequence_length, data_dir
):
    """An invalid loop forming pair will be detected by RNAlib's FE calculation and is
    assigned "infinity" free energy where infinity is 100k
    """
    if irs_match_same_bases(forms_valid_loop_ir_pair) or any(
        [not ir_has_valid_gap_size(ir) for ir in forms_valid_loop_ir_pair]
    ):
        pytest.skip("Pair matches same bases")

    pair_db_repr = irs_to_dot_bracket(list(forms_valid_loop_ir_pair), sequence_length)
    pairs_free_energy = calc_free_energy(pair_db_repr, sequence, data_dir)

    # Test for 90k as valid loops will bring total FE down but not by 10k
    assert pairs_free_energy <= 90_000


def test_ir_pairs_forming_invalid_loop_is_assigned_near_100k_free_energy_by_RNAlib(
    forms_invalid_loop_ir_pair, sequence, sequence_length, data_dir
):
    """An invalid loop forming pair will be detected by RNAlib's FE calculation and is
    assigned "infinity" free energy where infinity is 100k
    """
    if irs_match_same_bases(forms_invalid_loop_ir_pair) or any(
        [not ir_has_valid_gap_size(ir) for ir in forms_invalid_loop_ir_pair]
    ):
        pytest.skip("Pair matches same bases or IR in pair has invalid gap size.")

    pair_db_repr = irs_to_dot_bracket(list(forms_invalid_loop_ir_pair), sequence_length)
    pairs_free_energy = calc_free_energy(pair_db_repr, sequence, data_dir)

    # Test for 90k as valid loops will bring total FE down but not by 10k
    assert pairs_free_energy >= 90_000
