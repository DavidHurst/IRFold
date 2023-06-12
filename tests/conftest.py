import sys
import pytest
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))


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
def list_of_irs(rna_seq_15_bases_3_irs):
    return [((2, 3), (11, 12)), ((5, 6), (11, 12)), ((10, 11), (13, 14))]


@pytest.fixture
def list_of_ir_pairs_not_matching_same_bases():
    return [
        (((0, 1), (5, 6)), ((10, 12), (15, 17))),
        (((3, 4), (16, 17)), ((9, 10), (19, 20))),
    ]


@pytest.fixture
def list_of_ir_pairs_matching_same_bases():
    return [
        (((0, 1), (5, 6)), ((5, 6), (10, 12))),
        (((3, 4), (16, 17)), ((14, 16), (22, 24))),
    ]


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


@pytest.fixture(scope="module")
def expected_dot_bracket_reprs(rna_seq_15_bases_3_irs):
    return ["..((.......))..", ".....((....))..", "..........((.))"]


@pytest.fixture(scope="module")
def sequence_lengths(expected_dot_bracket_reprs):
    return [len(db_repr) for db_repr in expected_dot_bracket_reprs]
