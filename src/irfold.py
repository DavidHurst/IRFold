import typing
import subprocess
import RNA
import errno
import os
import itertools

from pathlib import Path
from typing import Tuple, List
from ortools.linear_solver import pywraplp


__all__ = ["IRFold"]

# ((left_strand_start, left_strand_end), (right_strand_start, right_strand_end))
IR = Tuple[Tuple[int, int], Tuple[int, int]]


class IRFold:
    """RNA secondary structure prediction based on extracting optimal
    inverted repeat configurations from the primary sequence"""

    def __init__(self, data_dir: str = None):
        if data_dir is None:
            self.data_dir = "."
        else:
            self.data_dir = data_dir

    def fold(
        self,
        sequence: str,
        seq_name: str,
        min_len: int,
        max_len: int,
        max_gap: int,
        input_file: str = None,
        mismatches: int = 0,
        irs_output_file: str = "found_irs.txt",
        ir_energies_output_file: str = None,
        solver_backend: str = "SCIP",
    ) -> Tuple[str, float]:
        seq_len: int = len(sequence)

        # Find IRs in sequence
        found_irs: List[IR] = self.find_irs(
            sequence=sequence,
            input_file=input_file,
            seq_name=seq_name,
            min_len=min_len,
            max_len=max_len,
            max_gap=max_gap,
            mismatches=mismatches,
            irs_output_file=irs_output_file,
        )
        n_irs_found: int = len(found_irs)
        if n_irs_found == 0:
            print(f"No IRs found, returning primary sequence.")
            return "".join(["." for _ in range(seq_len)]), 0

        # Evaluate free energy of each IR respectively
        db_reprs: List[str] = [
            self.__irs_to_dot_bracket([found_irs[i]], seq_len)
            for i in range(n_irs_found)
        ]
        ir_free_energies: List[float] = [
            self.__eval_free_energy(db_repr, sequence, ir_energies_output_file)
            for db_repr in db_reprs
        ]
        for i, (rpr, e) in enumerate(zip(db_reprs, ir_free_energies)):
            print(f"IR#{i:<2}:      ", rpr, e)

        # Define ILP and solve
        solver = self.__get_solver(
            n_irs_found, found_irs, ir_free_energies, solver_backend
        )
        status = solver.Solve()

        if status == pywraplp.Solver.OPTIMAL:
            # Return dot bracket repr and mfe of final solution
            active_ir_idxs: List[int] = [
                i
                for i, v in enumerate(solver.variables())
                if int(v.solution_value()) == 1
            ]
            dp_repr: str = self.__irs_to_dot_bracket(
                [found_irs[i] for i in active_ir_idxs], seq_len
            )
            return dp_repr, solver.Objective().Value()
        else:
            print("The problem does not have an optimal solution")
            return "".join(["." for _ in range(seq_len)]), 0

    def find_irs(
        self,
        sequence: str,
        input_file: str,
        seq_name: str,
        min_len: int,
        max_len: int,
        max_gap: int,
        mismatches: int = 0,  # not supported yet
        irs_output_file: str = "iupacpal_output.txt",
    ) -> List[IR]:
        # If no file specified create it for IUPACpal
        if input_file is None:
            input_file = Path(self.data_dir) / f"{seq_name}.fasta"
            with open(input_file, "w") as file:
                file.write(f">{seq_name}\n")
                file.write(sequence)
        else:  # Check file provided exists
            raise NotImplementedError

        _, out, err = self.__run(
            [
                "./IUPACpal",
                "-f",
                input_file,
                "-s",
                seq_name,
                "-m",
                str(min_len),
                "-M",
                str(max_len),
                "-g",
                str(max_gap),
                "-x",
                str(0),
                "-o",
                irs_output_file,
            ]
        )

        inverted_repeats = []

        valid_run = not "Error" in str(out)

        if valid_run:
            with open(irs_output_file) as f_in:
                lines = list(line for line in (l.strip() for l in f_in) if line)

            lines = lines[lines.index("Palindromes:") + 1 :]

            chunks = [lines[i : i + 3] for i in range(0, len(lines), 3)]

            for chunk in chunks:
                left_start, left_end = self.__extract_locs(chunk[0])
                right_end, right_start = self.__extract_locs(chunk[2])
                inverted_repeats.append(
                    ((left_start, left_end), (right_start, right_end))
                )

            return inverted_repeats
        else:
            print(str(out.decode("utf-8")))
            return None

    def __run(self, cmd):
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        return proc.returncode, stdout, stderr

    def __extract_locs(self, data):
        x = data.split(" ")
        a = int(x[0])
        b = int(x[-1])
        return a, b

    def __irs_to_dot_bracket(self, irs: List[IR], seq_len: int) -> str:
        """Does not support mismatches."""
        paired_bases: List[str] = ["." for _ in range(seq_len)]  # Initially none paired
        for ir in irs:
            n_base_pairs: int = ir[0][1] - ir[0][0] + 1

            lhs: Tuple[int, int]
            rhs: Tuple[int, int]
            lhs, rhs = ir[0], ir[1]

            # IUPACpal returns base pairings using 1-based indexing
            paired_bases[lhs[0] - 1 : lhs[1]] = ["(" for _ in range(n_base_pairs)]
            paired_bases[rhs[0] - 1 : rhs[1]] = [")" for _ in range(n_base_pairs)]

        return "".join(paired_bases)

    def __eval_free_energy(
        self, dot_brk_repr: str, sequence: str, out_file: str
    ) -> float:
        if out_file is None:
            out_file = f"{self.data_dir}/calculated_ir_energies.txt"

        with open(out_file, "a") as file:
            file.write(f"Evaluting IR:\n")
            for i in range(len(dot_brk_repr)):
                file.write(f"{i+1:<3}")
            file.write("\n")
            for b in dot_brk_repr:
                file.write(f"{b:<3}")
            file.write("\n")

            free_energy = RNA.eval_structure_simple(sequence, dot_brk_repr, 1, file)
            file.write(f"\n\n")

        return free_energy

    def __get_solver(
        self,
        n_irs: int,
        all_irs: List[IR],
        ir_free_energies: List[float],
        backend: str = "SCIP",
    ) -> pywraplp.Solver:
        # Create solver
        solver: pywraplp.Solver = pywraplp.Solver.CreateSolver(backend)
        if solver is None:
            raise Exception("Failed to create solver.")

        # Create binary indicator variables
        variables = [solver.IntVar(0, 1, f"ir_{i}") for i in range(n_irs)]

        # Add XOR between IRs that match the same bases
        unique_idx_pairs: List[Tuple[int, int]] = list(
            itertools.combinations([i for i in range(n_irs)], 2)
        )
        unique_ir_pairs: List[Tuple[IR, IR]] = [
            (all_irs[i], all_irs[j]) for i, j in unique_idx_pairs
        ]
        incompatible_ir_pair_idxs: List[Tuple[int, int]] = [
            idx_pair
            for ir_pair, idx_pair in zip(unique_ir_pairs, unique_idx_pairs)
            if self.irs_incompatible(ir_pair[0], ir_pair[1])
        ]

        for inc_ir_a, inc_ir_b in incompatible_ir_pair_idxs:
            solver.Add(variables[inc_ir_a] + variables[inc_ir_b] <= 1)

        # Define objective function
        obj_fn = solver.Objective()
        for i in range(n_irs):
            obj_fn.SetCoefficient(variables[i], ir_free_energies[i])
        obj_fn.SetMinimization()

        return solver

    def irs_incompatible(self, ir_a: IR, ir_b: IR) -> bool:
        paired_base_idxs_a = [idx for idx in range(ir_a[0][0], ir_a[0][1] + 1)] + [
            idx for idx in range(ir_a[1][0], ir_a[1][1] + 1)
        ]
        paired_base_idxs_b = [idx for idx in range(ir_b[0][0], ir_b[0][1] + 1)] + [
            idx for idx in range(ir_b[1][0], ir_b[1][1] + 1)
        ]

        return any(
            [ir_b_bases in paired_base_idxs_a for ir_b_bases in paired_base_idxs_b]
        )
