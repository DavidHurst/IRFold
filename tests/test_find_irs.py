import pytest
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from ir_fold import IRFold0

DATA_DIR = str(Path(__file__).parent / "tests_data")



@pytest.fixture
def list_of_found_irs(find_irs_params):
    return IRFold0.find_irs(**find_irs_params)


# ToDo: Write test for not having iupacpal compiled
# ToDo: Parametrise tests below to run multiple sequences with varying ir counts


def test_irs_found(list_of_found_irs):
    assert list_of_found_irs is not None


def test_irs_out_files_created(list_of_irs, find_irs_params):
    seq = find_irs_params["sequence"]
    seq_name = find_irs_params["seq_name"]

    assert (Path(DATA_DIR) / f"{seq_name}.fasta").exists()  # Sequence file is created
    assert (
        Path(DATA_DIR) / f"{seq_name}_found_irs.txt"
    ).exists()  # Found IRs file is created

    with open(str(Path(DATA_DIR).resolve() / f"{seq_name}.fasta")) as seq_file:
        written_seq = seq_file.readlines()[1]
    assert written_seq == seq


def test_n_irs_found(list_of_irs, rna_seq_15_bases_3_irs):
    irs_found = list_of_irs

    assert len(irs_found) == rna_seq_15_bases_3_irs[1]


def test_found_irs_gap_over_3(list_of_irs):
    for ir in list_of_irs:
        gap_sz = ir[1][0] - ir[0][1] - 1
        assert gap_sz >= 3
