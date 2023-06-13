__all__ = ["IRFold0"]

import csv
import subprocess
import itertools
import re

from pathlib import Path
from typing import Tuple, List
from ortools.sat.python.cp_model import (
    CpModel,
    CpSolver,
    IntVar,
    LinearExpr,
    OPTIMAL,
    FEASIBLE,
)

from irfold.util import (
    ir_pair_match_same_bases,
    irs_to_dot_bracket,
    calc_free_energy,
    IR,
    create_seq_file,
    write_performance_to_file,
    run_cmd,
)


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
        """Returns the MFE computed by simply adding the free energies of IRs which is incorrect but
        is illustrative of the IR free energy assumption additivity not holding consistently which is corrected
        for in class IRFold2."""

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
            db_repr, obj_fn_value = "".join(["." for _ in range(seq_len)]), 0
            if save_performance:
                write_performance_to_file(
                    db_repr, obj_fn_value, 0.0, seq_len, out_dir, cls.__name__
                )
            return db_repr, obj_fn_value

        # Define ILP and solve
        model, ir_variables = cls.get_solver(
            found_irs, seq_len, sequence, out_dir, seq_name, solver_backend
        )

        solver: CpSolver = CpSolver()
        status = solver.Solve(model)

        if status == OPTIMAL or status == FEASIBLE:
            # Return dot bracket repr and mfe of final solution
            active_ir_idxs: List[int] = [
                int(re.findall(r"-?\d+\.?\d*", v.Name())[0])
                for v in ir_variables
                if solver.Value(v) == 1
            ]
            db_repr: str = irs_to_dot_bracket(
                [found_irs[i] for i in active_ir_idxs], seq_len
            )

            obj_fn_value: float = solver.ObjectiveValue()
            dot_bracket_repr_mfe: float = calc_free_energy(
                db_repr, sequence, out_dir, seq_name
            )

            if save_performance:
                write_performance_to_file(
                    db_repr,
                    obj_fn_value,
                    dot_bracket_repr_mfe,
                    seq_len,
                    out_dir,
                    cls.__name__,
                )
            return db_repr, obj_fn_value
        else:
            # The optimisation problem does not have an optimal solution
            db_repr, obj_fn_value = "".join(["." for _ in range(seq_len)]), 0
            if save_performance:
                if save_performance:
                    write_performance_to_file(
                        db_repr, obj_fn_value, 0.0, seq_len, out_dir, cls.__name__
                    )
            return db_repr, obj_fn_value

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
        create_seq_file(sequence, seq_name, seq_file)
        irs_output_file: str = str(out_dir_path / f"{seq_name}_found_irs.txt")

        # ToDo: Refactor this to capture stdout of running IUPACpal instead of writing to file then extracting
        _, out, _ = run_cmd(
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
    def get_solver(
        ir_list: List[IR],
        seq_len: int,
        sequence: str,
        out_dir: str,
        seq_name: str,
        backend: str = "SCIP",
    ) -> Tuple[CpModel, List[IntVar]]:
        # Create solver
        model: CpModel = CpModel()

        n_irs: int = len(ir_list)

        # Evaluate free energy of each IR respectively
        db_reprs: List[str] = [
            irs_to_dot_bracket([ir_list[i]], seq_len) for i in range(n_irs)
        ]
        ir_free_energies: List[float] = [
            calc_free_energy(db_repr, sequence, out_dir, seq_name)
            for db_repr in db_reprs
        ]

        # Create binary indicator variables for IRs
        variables = [model.NewIntVar(0, 1, f"ir_{i}") for i in range(n_irs)]

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

        # All constraints and the objective must have integer coefficients for CP-SAT solver

        # Add XOR between IRs that match the same bases
        for ir_a_idx, ir_b_idx in incompatible_ir_pair_idxs:
            model.Add(variables[ir_a_idx] + variables[ir_b_idx] <= 1)

        # Define objective function
        variable_coefficients = [round(free_energy) for free_energy in ir_free_energies]
        obj_fn_expr = LinearExpr.WeightedSum(variables, variable_coefficients)

        model.Minimize(obj_fn_expr)

        return model, variables

    @staticmethod
    def ir_pair_incompatible(ir_a: IR, ir_b: IR) -> bool:
        return ir_pair_match_same_bases(ir_a, ir_b)
