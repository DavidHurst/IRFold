"""File to store share fixtures between tests"""

import pytest
from pathlib import Path

# @pytest.fixture
# def rna_seq_15_bases_3_irs():
#     return "CACCACCAUAAGGCU", 3


@pytest.fixture()
def data_dir():
    return str(Path(__file__).parent / "tests_data")


@pytest.fixture
def ir_list():
    return ((3, 4), (12, 13))  # , ((6, 7), (12, 13))]


@pytest.fixture
def expected_dot_bracket_reprs():
    return "..((.......)).."  # , ".....((....))..", "..........((.))"]


@pytest.fixture
def sequence_lengths():
    return 15
