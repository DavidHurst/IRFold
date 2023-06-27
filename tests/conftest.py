import sys
import pytest
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))


def setup_module():
    print("setup_module")


@pytest.fixture(scope="module")
def data_dir():
    return str(Path(__file__).parent / "tests_data")


@pytest.fixture(scope="module")
def rna_seq_15_bases_3_irs():
    return "CACCACCAUAAGGCU", 3


@pytest.fixture(scope="module")
def rna_seq_15_bases_1_irs():
    pass


@pytest.fixture(scope="module")
def rna_seq_15_bases_0_irs():
    return "AAAAAAAAAAAAAAA", 0


@pytest.fixture(scope="module")
def find_irs_params(rna_seq_15_bases_3_irs, data_dir):
    seq = rna_seq_15_bases_3_irs[0]
    seq_len = len(seq)
    return {
        "sequence": seq,
        "min_len": 2,
        "max_len": seq_len,
        "max_gap": seq_len - 1,
        "mismatches": 0,
        "seq_name": "test_irs_seq",
        "out_dir": data_dir,
    }


@pytest.fixture(scope="module")
def expected_dot_bracket_reprs(rna_seq_15_bases_3_irs):
    return ["..((.......))..", ".....((....))..", "..........((.))"]


@pytest.fixture(scope="module")
def sequence_lengths(expected_dot_bracket_reprs):
    return [len(db_repr) for db_repr in expected_dot_bracket_reprs]


# ====================  IR group definitions   ====================


@pytest.fixture(scope="module")
def list_of_irs(rna_seq_15_bases_3_irs):
    return [((2, 3), (11, 12)), ((5, 6), (11, 12)), ((10, 11), (13, 14))]


@pytest.fixture(scope="module")
def list_of_ir_pairs_not_matching_same_bases():
    return [
        (((0, 1), (5, 6)), ((10, 12), (15, 17))),
        (((3, 4), (16, 17)), ((9, 10), (19, 20))),
    ]


@pytest.fixture(scope="module")
def list_of_ir_pairs_matching_same_bases():
    return [
        (((0, 1), (5, 6)), ((5, 6), (10, 12))),
        (((3, 4), (16, 17)), ((14, 16), (22, 24))),
    ]


@pytest.fixture(scope="module")
def list_of_ir_pairs_invalid_gap_sizes():
    # Valid means IR has a >= 3 bases in hairpin
    return [
        ((0, 1), (2, 3)),
        ((0, 1), (3, 4)),
        ((0, 1), (4, 5)),
    ]


@pytest.fixture(scope="module")
def list_of_ir_pairs_valid_gap_sizes():
    # Valid means IR has a >= 3 bases in hairpin
    return [
        ((0, 1), (5, 6)),
        ((2, 5), (10, 13)),
        ((3, 6), (12, 15)),
    ]


@pytest.fixture(scope="module")
def list_of_ir_pairs_wholly_nested():
    return [
        (((7, 8), (10, 12)), ((2, 3), (16, 18))),  # First ir nested in second
        (((1, 3), (15, 18)), ((5, 7), (11, 13))),  # Second ir nested in first
    ]


@pytest.fixture(scope="module")
def list_of_ir_pairs_entirely_disjoint():
    return [
        (((2, 3), (7, 8)), ((10, 12), (16, 18))),  # Entirely disjoint
        (((2, 3), (7, 8)), ((9, 12), (16, 19))),  # Entirely disjoint
    ]


@pytest.fixture(scope="module")
def list_of_valid_and_invalid_ir_pairs(list_of_ir_pairs_entirely_disjoint):
    # Valid: pair does not match same bases or form invalid loop
    valid_pairs = list_of_ir_pairs_entirely_disjoint
    invalid_pairs = [
        (((2, 3), (7, 8)), ((2, 3), (16, 18))),  # Match same bases
        (((2, 3), (12, 13)), ((10, 11), (18, 19))),  # Invalid loop
    ]

    return zip(valid_pairs, [True for _ in range(len(valid_pairs))]) + zip(
        invalid_pairs, [False for _ in range(len(invalid_pairs))]
    )
