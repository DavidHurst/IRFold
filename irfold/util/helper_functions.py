import csv
import itertools
import subprocess
from typing import List, Tuple

from pathlib import Path

import RNA

IR = Tuple[Tuple[int, int], Tuple[int, int]]


def irs_to_dot_bracket(irs: List[IR], seq_len: int) -> str:
    """Does not support mismatches."""
    paired_bases: List[str] = ["." for _ in range(seq_len)]  # Initially none paired

    for ir in irs:
        n_base_pairs: int = ir[0][1] - ir[0][0] + 1  # Assumes no mismatches

        left_strand: Tuple[int, int]
        right_strand: Tuple[int, int]
        left_strand, right_strand = ir[0], ir[1]

        paired_bases[left_strand[0] : left_strand[1] + 1] = [
            "(" for _ in range(n_base_pairs)
        ]

        paired_bases[right_strand[0] : right_strand[1] + 1] = [
            ")" for _ in range(n_base_pairs)
        ]

    return "".join(paired_bases)


def get_valid_gap_sz_ir_n_tuples(
    n: int, num_irs: int, ir_list: List[IR], invalid_gap_sz_irs_idxs: List[int]
) -> Tuple[List[Tuple[IR, ...]], List[Tuple[int, ...]]]:
    """Returns all possible and valid (valid gap size) IR n-tuples i.e. tuples of size n that can be made from the
    provided IR list. E.g. all valid tuples of size 2 i.e. pairs that can be created from the IR list.
    Also returns the indices of each IR n-tuple."""
    all_unique_possible_idx_n_tuples: List[Tuple[int, ...]] = list(
        itertools.combinations([i for i in range(num_irs)], n)
    )
    valid_ir_idx_n_tuples: List[Tuple[int, ...]] = [
        ir_idx_n_tuple
        for ir_idx_n_tuple in all_unique_possible_idx_n_tuples
        if all([ir_idx_n_tuple[i] not in invalid_gap_sz_irs_idxs for i in range(n)])
    ]

    valid_ir_n_tuples: List[Tuple[IR, ...]] = [
        tuple(ir_list[ir_idx] for ir_idx in ir_idx_n_tuple)
        for ir_idx_n_tuple in valid_ir_idx_n_tuples
    ]

    return valid_ir_n_tuples, valid_ir_idx_n_tuples


def calc_free_energy(
    dot_brk_repr: str, sequence: str, out_dir: str, seq_name: str = "seq"
) -> float:
    out_dir_path: Path = Path(out_dir).resolve()
    if not out_dir_path.exists():
        out_dir_path = Path.cwd().resolve()

    out_file: str = str(out_dir_path / f"{seq_name}_calculated_ir_energies.txt")

    # RNAlib requires a file passed as parameter even if not writing to it
    with open(out_file, "a") as file:
        free_energy = RNA.eval_structure_simple(sequence, dot_brk_repr, -1, file)

    return free_energy


def create_seq_file(seq: str, seq_name: str, file_name: str) -> None:
    with open(file_name, "w") as file:
        file.write(f">{seq_name}\n")
        file.write(seq)


def run_cmd(cmd):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return proc.returncode, stdout, stderr


def write_performance_to_file(
    dot_bracket_repr: str,
    obj_fn_final_value: float,
    dot_bracket_repr_mfe: float,
    seq_len: int,
    out_dir: str,
    ssp_model_name: str,
    n_irs_found: int = None,
    solver_num_booleans: int = None,
    solver_solve_time: float = None,
    solver_num_branches_explored: int = None,
    solver_num_conflicts: int = None,
):
    out_dir_path: Path = Path(out_dir).resolve()
    if not out_dir_path.exists():
        out_dir_path = Path.cwd().resolve()

    performance_file_name: str = f"{ssp_model_name}_performance.csv"
    performance_file_path: Path = (Path(out_dir_path) / performance_file_name).resolve()

    if (
        not performance_file_path.exists()
    ):  # First performance sample being written, make new file
        performance_file_path = (out_dir_path / performance_file_name).resolve()
        column_names = [
            "dot_bracket_repr",
            "obj_fn_final_value",
            "dot_bracket_repr_mfe",
            "seq_len",
            "n_irs_found",
            "solver_num_booleans",
            "solver_solve_time",
            "solver_iterations",
            "solver_num_branches_explored",
            "solver_num_conflicts",
        ]

        with open(str(performance_file_path), "w") as perf_file:
            writer = csv.writer(perf_file)
            writer.writerow(column_names)
            writer.writerow(
                [
                    dot_bracket_repr,
                    obj_fn_final_value,
                    dot_bracket_repr_mfe,
                    seq_len,
                    n_irs_found,
                    solver_num_booleans,
                    solver_solve_time,
                    solver_num_branches_explored,
                    solver_num_conflicts,
                ]
            )
    else:  # Add new performance sample to existing file
        with open(str(performance_file_path), "a") as perf_file:
            writer = csv.writer(perf_file)
            writer.writerow(
                [
                    dot_bracket_repr,
                    obj_fn_final_value,
                    dot_bracket_repr_mfe,
                    seq_len,
                    n_irs_found,
                    solver_num_booleans,
                    solver_solve_time,
                    solver_num_branches_explored,
                    solver_num_conflicts,
                ]
            )
