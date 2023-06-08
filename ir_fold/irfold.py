import subprocess
import RNA
import itertools
import re

from pathlib import Path
from typing import Tuple, List
from ortools.linear_solver import pywraplp

# ((left_strand_start, left_strand_end), (right_strand_start, right_strand_end))
IR = Tuple[Tuple[int, int], Tuple[int, int]]


class IRFold:
    """RNA secondary structure prediction based on extracting optimal
    inverted repeat configurations from the primary sequence"""

    @staticmethod
    def fold(
        sequence: str,
        min_len: int,
        max_len: int,
        max_gap: int,
        # seq_file: str = None,
        # seq_name: str = "seq",
        mismatches: int = 0,
        solver_backend: str = "SCIP",
        out_dir: str = ".",
    ) -> Tuple[str, float]:
        # Find IRs in sequence
        found_irs: List[IR] = IRFold.find_irs(
            sequence=sequence,
            # seq_file=seq_file,
            # seq_name=seq_name,
            min_len=min_len,
            max_len=max_len,
            max_gap=max_gap,
            mismatches=mismatches,
            out_dir=out_dir,
        )

        n_irs_found: int = len(found_irs)
        seq_len: int = len(sequence)
        if n_irs_found == 0:
            print(f"No IRs found, returning primary sequence.")
            return "".join(["." for _ in range(seq_len)]), 0

        # Evaluate free energy of each IR respectively
        db_reprs: List[str] = [
            IRFold.irs_to_dot_bracket([found_irs[i]], seq_len)
            for i in range(n_irs_found)
        ]
        ir_free_energies: List[float] = [
            IRFold.calc_free_energy(db_repr, sequence, out_dir) for db_repr in db_reprs
        ]

        # Define ILP and solve
        solver = IRFold.ir_ilp_solver(
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
            dp_repr: str = IRFold.irs_to_dot_bracket(
                [found_irs[i] for i in active_ir_idxs], seq_len
            )
            return dp_repr, solver.Objective().Value()
        else:
            print("The problem does not have an optimal solution")
            return "".join(["." for _ in range(seq_len)]), 0

    @staticmethod
    def find_irs(
        sequence: str,
        min_len: int,
        max_len: int,
        max_gap: int,
        mismatches: int = 0,  # not supported yet
        out_dir: str = ".",
    ) -> List[IR]:
        out_dir_path: Path = Path(out_dir).resolve()
        if not out_dir_path.exists():
            print(f"Output directory not found, defaulting to current directory.")
            out_dir_path = Path(".").resolve()

        # Check IUPACpal has been compiled to this cwd
        iupacpal_exe: Path = Path(__file__).parent / "IUPACpal"
        if not iupacpal_exe.exists():
            raise FileNotFoundError("Could not find IUPACpal executable.")

        # Write sequence to file for IUPACpal
        # ToDo: Parametrise these in/out files
        seq_name: str = "seq"
        seq_file: str = str(out_dir_path / f"{seq_name}.fasta")
        IRFold.create_seq_file(sequence, seq_name, seq_file)
        irs_output_file: str = str(out_dir_path / f"{seq_name}_found_irs.txt")

        # ToDo: Refactor this to capture stdout of running IUPACpal instead of writing to file then extracting
        _, out, _ = IRFold.__run_cmd(
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
                str(out_dir_path / "seq_found_irs.txt"),
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

                left_start, left_end = int(ir_idxs[0]), int(ir_idxs[1])
                right_start, right_end = int(ir_idxs[3]), int(ir_idxs[2])
                found_irs.append(((left_start, left_end), (right_start, right_end)))

            # Remove IRs with a gap < 3, these are sterically impossible
            valid_irs: List[IR] = [
                ir for ir in found_irs if ir[1][0] - ir[0][1] - 1 >= 3
            ]

            return valid_irs
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
            n_base_pairs: int = ir[0][1] - ir[0][0] + 1

            left_strand: Tuple[int, int]
            right_strand: Tuple[int, int]
            left_strand, right_strand = ir[0], ir[1]

            # IUPACpal returns base pairings using 1-based indexing
            paired_bases[left_strand[0] - 1 : left_strand[1]] = [
                "(" for _ in range(n_base_pairs)
            ]
            paired_bases[right_strand[0] - 1 : right_strand[1]] = [
                ")" for _ in range(n_base_pairs)
            ]

        return "".join(paired_bases)

    @staticmethod
    def calc_free_energy(dot_brk_repr: str, sequence: str, out_dir: str) -> float:
        out_dir_path: Path = Path(out_dir).resolve()
        if not out_dir_path.exists():
            print(f"Output directory not found, defaulting to current directory.")
            out_dir_path = Path(".").resolve()
        out_file: str = str(out_dir_path / "calculated_ir_energies.txt")

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
    def ir_ilp_solver(
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
            if IRFold.irs_share_base_pair(ir_pair[0], ir_pair[1]) and not IRFold.irs_form_valid_loop(ir_pair[0], ir_pair[1])
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
    def irs_share_base_pair(ir_a: IR, ir_b: IR) -> bool:
        # Check that IRs don't match the same bases
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
    def irs_disjoint(ir_a, ir_b):
        # ir_a comes entirely before ir_b
        ir_a_strictly_before_ir_b = ir_a[1][1] < ir_b[0][0]

        # ir_b comes entirely before ir_a
        ir_b_strictly_before_ir_a = ir_b[1][1] < ir_a[0][0]

        return ir_a_strictly_before_ir_b or ir_b_strictly_before_ir_a

    @staticmethod
    def irs_form_valid_loop(ir_a, ir_b):
        """This won't scale up, not entirely sure what kind of invalid loop this finds but
        matches up with loops that RNAlib assigns infinite free energy to so validated."""
        # Note: IUPACpal is 0-based indexing IRs
        latest_left_string_base_idx = (
            ir_a[0][1] if ir_a[0][1] > ir_b[0][1] else ir_b[0][1]
        ) - 1
        earliest_right_string_base_idx = (
            ir_a[1][0] if ir_a[1][0] < ir_b[1][0] else ir_b[1][0]
        ) - 1

        # Case 1: Disjoint
        if IRFold.irs_disjoint(ir_a, ir_b):
            return True

        # Case 2: Invalid loop
        bases_inbetween_parens = (
            earliest_right_string_base_idx - latest_left_string_base_idx - 1
        )
        if bases_inbetween_parens < 3:
            return False

        # Case 3: Valid loop
        return True
