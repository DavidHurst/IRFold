import csv
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


def calc_free_energy(
    dot_brk_repr: str, sequence: str, out_dir: str, seq_name: str = "seq"
) -> float:
    out_dir_path: Path = Path(out_dir).resolve()
    if not out_dir_path.exists():
        out_dir_path = Path.cwd().resolve()

    out_file: str = str(out_dir_path / f"{seq_name}_calculated_ir_energies.txt")

    with open(out_file, "a") as file:
        file.write(f"Evaluating IR:\n")
        for i in range(len(dot_brk_repr)):
            file.write(f"{i + 1:<3}")
        file.write("\n")
        for b in dot_brk_repr:
            file.write(f"{b:<3}")
        file.write("\n")

        free_energy = RNA.eval_structure_simple(sequence, dot_brk_repr, 1, file)
        file.write(f"\n\n")

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
    class_name: str,
    n_irs_found: int = None,
    solver_num_booleans: int = None,
    solver_solve_time: float = None,
    solver_num_branches_explored: int = None,
    solver_num_conflicts: int = None,
):
    out_dir_path: Path = Path(out_dir).resolve()
    if not out_dir_path.exists():
        out_dir_path = Path.cwd().resolve()

    performance_file_name: str = f"{class_name}_performance.csv"
    performance_file_path: Path = (Path(out_dir_path) / performance_file_name).resolve()

    if not performance_file_path.exists():
        performance_file_path = (out_dir_path / performance_file_name).resolve()
        column_names = [
            "dot_bracket_repr",
            "solution_mfe",
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
    else:
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
