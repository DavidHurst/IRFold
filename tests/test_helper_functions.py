from pathlib import Path

from irfold.util import irs_to_dot_bracket, calc_free_energy, write_performance_to_file


def test_db_conversion_lengths_match(
    list_of_irs, expected_dot_bracket_reprs, sequence_lengths
):
    for ir, expected, seq_len in zip(
        list_of_irs, expected_dot_bracket_reprs, sequence_lengths
    ):
        generated_db_repr = irs_to_dot_bracket([ir], seq_len)
        assert len(generated_db_repr) == len(expected)


def test_db_conversion_output_matches(
    list_of_irs, expected_dot_bracket_reprs, sequence_lengths
):
    for ir, expected_db_repr, seq_len in zip(
        list_of_irs, expected_dot_bracket_reprs, sequence_lengths
    ):
        generated_db_repr = irs_to_dot_bracket([ir], seq_len)
        assert generated_db_repr == expected_db_repr


def test_calc_free_energy(data_dir, rna_seq_15_bases_0_irs):
    seq = rna_seq_15_bases_0_irs[0]
    free_energy = calc_free_energy("...............", seq, data_dir, "test_fe_calc_seq")

    assert free_energy is not None
    assert free_energy == 0.0


def test_write_performance_to_file(data_dir):
    perf_file_name = "test_write_perf_to_file"
    perf_file = Path(data_dir) / f"{perf_file_name}_performance.csv"

    if perf_file.exists():
        perf_file.unlink()

    write_performance_to_file(
        dot_bracket_repr="first_sample",
        obj_fn_final_value=0.0,
        dot_bracket_repr_mfe=0.0,
        seq_len=0,
        out_dir=data_dir,
        ssp_model_name=perf_file_name,
        n_irs_found=0,
        solver_num_booleans=0,
        solver_solve_time=0.0,
        solver_num_branches_explored=0,
        solver_num_conflicts=0,
    )

    # Check only one line written
    with open(perf_file, "r") as file:
        lines = file.readlines()

    assert len(lines) == 2
    assert lines[1].split(",")[0] == "first_sample"

    # Check the next write is added to the same file
    write_performance_to_file(
        dot_bracket_repr="second_sample",
        obj_fn_final_value=0.0,
        dot_bracket_repr_mfe=0.0,
        seq_len=0,
        out_dir=data_dir,
        ssp_model_name=perf_file_name,
        n_irs_found=0,
        solver_num_booleans=0,
        solver_solve_time=0.0,
        solver_num_branches_explored=0,
        solver_num_conflicts=0,
    )

    with open(perf_file, "r") as file:
        lines = file.readlines()

    assert len(lines) == 3
    assert lines[2].split(",")[0] == "second_sample"
