from pathlib import Path
from ir_fold import IRFold0

# ToDo: Write test for not having iupacpal compiled
# ToDo: Parametrise tests below to run multiple sequences with varying ir counts


def test_irs_found(find_irs_params):
    assert IRFold0.find_irs(**find_irs_params) is not None


def test_irs_out_files_created(list_of_irs, find_irs_params, data_dir):
    seq = find_irs_params["sequence"]
    seq_name = find_irs_params["seq_name"]

    assert (Path(data_dir) / f"{seq_name}.fasta").exists()  # Sequence file is created
    assert (
        Path(data_dir) / f"{seq_name}_found_irs.txt"
    ).exists()  # Found IRs file is created

    with open(str(Path(data_dir).resolve() / f"{seq_name}.fasta")) as seq_file:
        written_seq = seq_file.readlines()[1]
    assert written_seq == seq


def test_n_irs_found(list_of_irs, rna_seq_15_bases_3_irs):
    n_irs_in_seq = rna_seq_15_bases_3_irs[1]
    assert len(list_of_irs) == n_irs_in_seq
