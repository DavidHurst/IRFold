import pytest

from ir_fold import IRFold0, IRFold1


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


# Test irs pair same bases

# Test irs are disjoint

# Test
