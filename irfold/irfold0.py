__all__ = ["IRFold0"]

import csv
import subprocess
import RNA
import itertools
import re

from pathlib import Path
from typing import Tuple, List
from ortools.linear_solver import pywraplp

# ((left_strand_start, left_strand_end), (right_strand_start, right_strand_end))
IR = Tuple[Tuple[int, int], Tuple[int, int]]


class IRFold0:
    """Base IRFold Model: RNA secondary structure prediction based on extracting optimal
    inverted repeat configurations from the primary sequence"""

    @classmethod
    def fold(
        cls,
        sequence: str,
        min_len: int,
        max_len: int,
        max_gap: int,
        seq_name: str = "seq",
        mismatches: int = 0,
        solver_backend: str = "SCIP",
        out_dir: str = ".",
        *,
        save_performance: bool = False,
    ) -> Tuple[str, float]:

        # Find IRs in sequence
        found_irs: List[IR] = IRFold0.find_irs(
            sequence=sequence,
            seq_name=seq_name,
            min_len=min_len,
            max_len=max_len,
            max_gap=max_gap,
            mismatches=mismatches,
            out_dir=out_dir,
        )

        n_irs_found: int = len(found_irs)
        seq_len: int = len(sequence)
        if n_irs_found == 0:  # Return sequence if no IRs found
            db_repr, mfe = "".join(["." for _ in range(seq_len)]), 0
            if save_performance:
                cls.__write_performance_to_file(
                    db_repr,
                    mfe,
                    seq_len,
                    n_irs_found,
                    None,
                    None,
                    None,
                    None,
                    None,
                    out_dir,
                )
            return db_repr, mfe

        # Define ILP and solve
        solver = cls.get_lp_solver(
            found_irs, seq_len, sequence, out_dir, seq_name, solver_backend
        )
        status = solver.Solve()

        if status == pywraplp.Solver.OPTIMAL:
            # Return dot bracket repr and mfe of final solution
            active_ir_idxs: List[int] = [
                i
                for i, v in enumerate(solver.variables())
                if int(v.solution_value()) == 1
            ]
            db_repr: str = IRFold0.irs_to_dot_bracket(
                [found_irs[i] for i in active_ir_idxs], seq_len
            )
            solution_mfe: float = solver.Objective().Value()
            if save_performance:
                cls.__write_performance_to_file(
                    db_repr,
                    solution_mfe,
                    seq_len,
                    n_irs_found,
                    solver.NumVariables(),
                    solver.NumConstraints(),
                    solver.wall_time(),
                    solver.iterations(),
                    solver.nodes(),
                    out_dir,
                )
            return db_repr, solution_mfe
        else:
            # The optimisation problem does not have an optimal solution
            db_repr, mfe = "".join(["." for _ in range(seq_len)]), 0
            if save_performance:
                cls.__write_performance_to_file(
                    db_repr,
                    mfe,
                    seq_len,
                    n_irs_found,
                    None,
                    None,
                    None,
                    None,
                    None,
                    out_dir,
                )
            return db_repr, mfe

    @staticmethod
    def find_irs(
        sequence: str,
        min_len: int,
        max_len: int,
        max_gap: int,
        seq_name: str = "seq",
        mismatches: int = 0,  # not supported yet
        out_dir: str = ".",
    ) -> List[IR]:
        out_dir_path: Path = Path(out_dir).resolve()
        if not out_dir_path.exists():
            out_dir_path = Path.cwd().resolve()

        # Check IUPACpal has been compiled to this cwd
        iupacpal_exe: Path = Path(__file__).parent / "IUPACpal"
        if not iupacpal_exe.exists():
            raise FileNotFoundError("Could not find IUPACpal executable.")

        # Write sequence to file for IUPACpal
        seq_file: str = str(out_dir_path / f"{seq_name}.fasta")
        IRFold0.create_seq_file(sequence, seq_name, seq_file)
        irs_output_file: str = str(out_dir_path / f"{seq_name}_found_irs.txt")

        # ToDo: Refactor this to capture stdout of running IUPACpal instead of writing to file then extracting
        _, out, _ = IRFold0.__run_cmd(
            [
                str(iupacpal_exe),
                "-f",
                seq_file,
                "-s",
                seq_name,
                "-m",
                str(min_len),
                "-M",
                str(max_len),
                "-g",
                str(max_gap),
                "-x",
                str(mismatches),
                "-o",
                str(out_dir_path / f"{seq_name}_found_irs.txt"),
            ]
        )

        if "Error" not in str(out):
            # Extract IR indices from format IUPACpal outputs
            found_irs: List[IR] = []

            with open(irs_output_file) as f_in:
                lines: List[str] = list(
                    line for line in (l.strip() for l in f_in) if line
                )

            ir_lines: List[str] = lines[lines.index("Palindromes:") + 1 :]
            formatted_irs: List[List[str]] = [
                ir_lines[i : i + 3] for i in range(0, len(ir_lines), 3)
            ]

            for f_ir in formatted_irs:
                ir_idxs: List[int] = re.findall(r"-?\d+\.?\d*", "".join(f_ir))

                left_start, left_end = int(ir_idxs[0]) - 1, int(ir_idxs[1]) - 1
                right_start, right_end = int(ir_idxs[3]) - 1, int(ir_idxs[2]) - 1
                found_irs.append(((left_start, left_end), (right_start, right_end)))

            return found_irs
        else:
            print(str(out.decode("utf-8")))
            return []

    @staticmethod
    def create_seq_file(seq: str, seq_name: str, file_name: str) -> None:
        with open(file_name, "w") as file:
            file.write(f">{seq_name}\n")
            file.write(seq)

    @staticmethod
    def __run_cmd(cmd):
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        return proc.returncode, stdout, stderr

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def get_lp_solver(
        ir_list: List[IR],
        seq_len: int,
        sequence: str,
        out_dir: str,
        seq_name: str,
        backend: str = "SCIP",
    ) -> pywraplp.Solver:
        # Create solver
        solver: pywraplp.Solver = pywraplp.Solver.CreateSolver(backend)
        if solver is None:
            raise Exception("Failed to create solver.")
        n_irs: int = len(ir_list)

        # Evaluate free energy of each IR respectively
        db_reprs: List[str] = [
            IRFold0.irs_to_dot_bracket([ir_list[i]], seq_len) for i in range(n_irs)
        ]
        ir_free_energies: List[float] = [
            IRFold0.calc_free_energy(db_repr, sequence, out_dir, seq_name)
            for db_repr in db_reprs
        ]

        # Create binary indicator variables
        variables = [solver.IntVar(0, 1, f"ir_{i}") for i in range(n_irs)]

        # Add XOR between IRs that match the same bases
        unique_idx_pairs: List[Tuple[int, int]] = list(
            itertools.combinations([i for i in range(n_irs)], 2)
        )
        unique_ir_pairs: List[Tuple[IR, IR]] = [
            (ir_list[i], ir_list[j]) for i, j in unique_idx_pairs
        ]
        incompatible_ir_pair_idxs: List[Tuple[int, int]] = [
            idx_pair
            for ir_pair, idx_pair in zip(unique_ir_pairs, unique_idx_pairs)
            if IRFold0.ir_pair_incompatible(ir_pair[0], ir_pair[1])
        ]

        for inc_ir_a, inc_ir_b in incompatible_ir_pair_idxs:
            solver.Add(variables[inc_ir_a] + variables[inc_ir_b] <= 1)

        # Define objective function
        obj_fn = solver.Objective()
        for i in range(n_irs):
            obj_fn.SetCoefficient(variables[i], ir_free_energies[i])
        obj_fn.SetMinimization()

        return solver

    @staticmethod
    def ir_pair_match_same_bases(ir_a: IR, ir_b: IR) -> bool:
        # Check if IRs match the same bases
        ir_a_left_strand, ir_a_right_strand = ir_a[0], ir_a[1]
        paired_base_idxs_a = [
            idx for idx in range(ir_a_left_strand[0], ir_a_left_strand[1] + 1)
        ] + [idx for idx in range(ir_a_right_strand[0], ir_a_right_strand[1] + 1)]

        ir_b_left_strand, ir_b_right_strand = ir_b[0], ir_b[1]
        paired_base_idxs_b = [
            idx for idx in range(ir_b_left_strand[0], ir_b_left_strand[1] + 1)
        ] + [idx for idx in range(ir_b_right_strand[0], ir_b_right_strand[1] + 1)]

        any_shared_base_pairs: bool = any(
            [ir_b_bases in paired_base_idxs_a for ir_b_bases in paired_base_idxs_b]
        )
        return any_shared_base_pairs

    @staticmethod
    def ir_pair_incompatible(ir_a: IR, ir_b: IR) -> bool:
        return IRFold0.ir_pair_match_same_bases(ir_a, ir_b)

    @classmethod
    def __write_performance_to_file(
        cls,
        dot_bracket_repr: str,
        solution_mfe: float,
        seq_len: int,
        n_irs_found: int,
        solver_num_variables: int,
        solver_num_constraints: int,
        solver_solve_time: float,
        solver_iterations: int,
        solver_num_branch_bound_nodes: int,
        out_dir: str,
    ):
        out_dir_path: Path = Path(out_dir).resolve()
        if not out_dir_path.exists():
            out_dir_path = Path.cwd().resolve()

        performance_file_name: str = f"{cls.__name__}_performance.csv"
        performance_file_path: Path = (
            Path(out_dir_path) / performance_file_name
        ).resolve()

        if not performance_file_path.exists():
            performance_file_path = (out_dir_path / performance_file_name).resolve()
            column_names = [
                "dot_bracket_repr",
                "solution_mfe",
                "seq_len",
                "n_irs_found",
                "solver_num_variables",
                "solver_num_constraints",
                "solver_solve_time",
                "solver_iterations",
                "solver_num_branch_bound_nodes",
            ]

            with open(str(performance_file_path), "w") as perf_file:
                writer = csv.writer(perf_file)
                writer.writerow(column_names)
                writer.writerow(
                    [
                        dot_bracket_repr,
                        solution_mfe,
                        seq_len,
                        n_irs_found,
                        solver_num_variables,
                        solver_num_constraints,
                        solver_solve_time,
                        solver_iterations,
                        solver_num_branch_bound_nodes,
                    ]
                )
            #     print(f'\n >>>>>> {cls.__name__} <<<<<')
            # with open(str(performance_file_path), "r") as perf_file:
            #     print(f'File contents after create and write:')
            #     lines = perf_file.readlines()
            #     for l in lines:
            #         print(l)

        else:
            with open(str(performance_file_path), "a") as perf_file:
                writer = csv.writer(perf_file)
                writer.writerow(
                    [
                        dot_bracket_repr,
                        solution_mfe,
                        seq_len,
                        n_irs_found,
                        solver_num_variables,
                        solver_num_constraints,
                        solver_solve_time,
                        solver_iterations,
                        solver_num_branch_bound_nodes,
                    ]
                )
            #     print(f'\n >>>>>> {cls.__name__} <<<<<')
            #
            # with open(str(performance_file_path), "r") as perf_file:
            #     print(f'File contents after append:')
            #     lines = perf_file.readlines()
            #     for l in lines:
            #         print(l)
