import pytest
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from src import IRFold

DATA_DIR = "./tests_data"

# ToDo: add conftest.py file to share these fixtures across tests
@pytest.fixture
def IRFold_obj():
    return IRFold(DATA_DIR)


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
    }


@pytest.fixture
def list_of_found_irs(IRFold_obj, find_irs_params):
    return IRFold_obj.find_irs(**find_irs_params)


# ToDo: Write test for not having iupacpal compiled
# ToDo: Paremetrise tests below to run multiple sequences with varying ir counts


def test_irs_found(list_of_found_irs):
    assert list_of_found_irs is not None


def test_irs_out_files_created(list_of_found_irs, find_irs_params):
    list_of_found_irs
    seq = find_irs_params["sequence"]

    assert (Path(DATA_DIR) / "seq.fasta").exists()  # Sequence file is created
    assert (Path(DATA_DIR) / "seq_found_irs.txt").exists()  # Found IRs file is created

    with open(str(Path(DATA_DIR) / "seq.fasta")) as seq_file:
        written_seq = seq_file.readlines()[1]
    assert written_seq == seq


def test_n_irs_found(list_of_found_irs, rna_seq_15_bases_3_irs):
    irs_found = list_of_found_irs

    assert len(irs_found) == rna_seq_15_bases_3_irs[1]
