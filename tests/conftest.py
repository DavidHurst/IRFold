import sys
import pytest
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

@pytest.fixture(scope="module")
def data_dir():
    return str(Path(__file__).parent / "tests_data")

# @pytest.fixture
# def rna_seq_15_bases_3_irs():
#     return "CACCACCAUAAGGCU", 3

@pytest.fixture
def rna_seq_15_bases_3_irs():
    return "CACCACCAUAAGGCU", 2


@pytest.fixture
def rna_seq_15_bases_0_irs():
    return "AAAAAAAAAAAAAAA", 0


@pytest.fixture
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
def list_of_irs():
    return [((3, 4), (12, 13)), ((6, 7), (12, 13)), ((11, 12), (14, 15))]


@pytest.fixture(scope="module")
def expected_dot_bracket_reprs():
    return ["..((.......))..", ".....((....))..", "..........((.))"]


@pytest.fixture(scope="module")
def sequence_lengths(expected_dot_bracket_reprs):
    return [len(db_repr) for db_repr in expected_dot_bracket_reprs]
