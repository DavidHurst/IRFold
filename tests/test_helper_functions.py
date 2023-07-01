from pathlib import Path

from irfold.util import irs_to_dot_bracket, calc_free_energy, write_performance_to_file


def test_dot_bracket_conversion_lengths_match(all_irs, sequence_length):
    generated_db_repr = irs_to_dot_bracket(all_irs, sequence_length)
    assert len(generated_db_repr) == sequence_length


def test_db_conversion_output_matches_irs(
    all_irs, all_ir_dot_bracket_reprs, sequence_length
):
    # ToDo: This should use the ir and ir_dot_bracket_repr fixtures but weird combinatorial issue going on
    for ir, db_repr in zip(all_irs, all_ir_dot_bracket_reprs):
        generated_db_repr = irs_to_dot_bracket([ir], sequence_length)
        assert generated_db_repr == db_repr


# ToDo: Write tests for IR pair and triplet dot bracket conversions


def test_calc_free_energy(data_dir, sequence, sequence_name, sequence_length):
    free_energy = calc_free_energy(
        "." * sequence_length, sequence, data_dir, sequence_name
    )

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
