import pytest

from irfold import IRFold0, IRFold1


@pytest.fixture(scope="module")
def entirely_disjoint_ir_pairs():
    return [
        (((2, 3), (7, 8)), ((10, 12), (16, 18))),  # Entirely disjoint
        (((2, 3), (7, 8)), ((9, 12), (16, 19))),  # Entirely disjoint
    ]


@pytest.fixture(scope="module")
def list_of_valid_and_invalid_ir_pairs(entirely_disjoint_ir_pairs):
    # Valid: pair does not match same bases or form invalid loop
    valid_pairs = entirely_disjoint_ir_pairs
    invalid_pairs = [
        (((2, 3), (7, 8)), ((2, 3), (16, 18))),  # Match same bases
        (((2, 3), (12, 13)), ((10, 11), (18, 19))),  # Invalid loop
    ]

    return zip(valid_pairs, [True for _ in range(len(valid_pairs))]) + zip(
        invalid_pairs, [False for _ in range(len(invalid_pairs))]
    )


def test_ir_gap_size_check():
    invalid_gap_size_irs = [
        ((0, 1), (2, 3)),
        ((0, 1), (3, 4)),
        ((0, 1), (4, 5)),
    ]

    valid_gap_size_irs = [
        ((0, 1), (5, 6)),
        ((2, 5), (10, 13)),
        ((3, 6), (12, 15)),
    ]

    for ir in invalid_gap_size_irs:
        assert IRFold1.ir_has_valid_gap_size(ir) == False

    for ir in valid_gap_size_irs:
        assert IRFold1.ir_has_valid_gap_size(ir) == True


def test_ir_pair_matches_same_bases(
    list_of_ir_pairs_with_expected_value_for_base_overpairing,
):
    for ir_pair, expected in list_of_ir_pairs_with_expected_value_for_base_overpairing:
        print(ir_pair, expected)
        actual = IRFold0.ir_pair_match_same_bases(ir_pair[0], ir_pair[1])
        assert actual == expected


# Test irs are disjoint

# Test
