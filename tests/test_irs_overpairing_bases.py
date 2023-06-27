import pytest
from irfold.util import ir_pair_match_same_bases


@pytest.fixture
def list_of_ir_pairs_with_expected_value_for_base_overpairing(
    list_of_ir_pairs_matching_same_bases, list_of_ir_pairs_not_matching_same_bases
):
    pairs_overmatching_with_expected = list(
        zip(
            list_of_ir_pairs_matching_same_bases,
            [True for _ in range(len(list_of_ir_pairs_matching_same_bases))],
        )
    )
    pairs_not_overmatching_with_expected = list(
        zip(
            list_of_ir_pairs_not_matching_same_bases,
            [False for _ in range(len(list_of_ir_pairs_not_matching_same_bases))],
        )
    )
    return pairs_not_overmatching_with_expected + pairs_overmatching_with_expected


def test_ir_pair_matches_same_bases(
    list_of_ir_pairs_with_expected_value_for_base_overpairing,
):
    for ir_pair, expected in list_of_ir_pairs_with_expected_value_for_base_overpairing:
        actual = ir_pair_match_same_bases(ir_pair[0], ir_pair[1])
        assert actual == expected
