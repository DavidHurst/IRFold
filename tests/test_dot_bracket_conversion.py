import pytest
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from ir_fold import IRFold

DATA_DIR = str(Path(__file__).parent / 'tests_data')


@pytest.fixture
def rna_seq_15_bases_3_irs():
    return "CACCACCAUAAGGCU", 3


@pytest.fixture
def rna_seq_15_bases_0_irs():
    return "AAAAAAAAAAAAAAA", 0


@pytest.fixture
def find_irs_params(rna_seq_15_bases_3_irs):
    seq = rna_seq_15_bases_3_irs[0]
    seq_len = len(seq)
    return {
        "sequence": seq,
        "min_len": 2,
        "max_len": seq_len,
        "max_gap": seq_len - 1,
        "mismatches": 0,
        "out_dir": DATA_DIR,
    }


@pytest.fixture
def list_of_found_irs(find_irs_params):
    return IRFold.find_irs(**find_irs_params)


def test_ir_to_db_conversion(list_of_found_irs, find_irs_params):
    seq_len = len(find_irs_params["sequence"])
    expected_db_reprs = ["..((.......))..", ".....((....))..", "..........((.))"]

    for i, ir in enumerate(list_of_found_irs):
        generated_db_repr = IRFold.irs_to_dot_bracket([ir], seq_len)
        assert len(generated_db_repr) == seq_len
        assert generated_db_repr == expected_db_reprs[i]
