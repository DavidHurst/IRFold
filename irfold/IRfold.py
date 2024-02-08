__all__ = ["IRfold"]

import re
from pathlib import Path
from typing import Tuple, List


from .util import (
    ir_has_valid_gap_size,
    IR,
    irs_to_dot_bracket,
    calc_free_energy,
    ir_pair_invalid_relative_pos,
    get_valid_gap_sz_ir_n_tuples,
    write_solver_performance_to_file,
    create_seq_file,
    run_cmd,
)
from ortools.sat.python.cp_model import (
    CpModel,
    CpSolver,
    IntVar,
    LinearExpr,
    OPTIMAL,
    FEASIBLE,
)
from tqdm import tqdm


class IRfold:
    """Extends IRFold model by validating found IRs and IRs in pairs before passing them to the solver."""

    @classmethod
    def fold(
        cls,
        sequence: str,
        out_dir: str = ".",
        *,
        seq_name: str = "seq",
        save_performance: bool = False,
    ) -> Tuple[str, float]:

        # Find IRs in sequence
        found_irs: List[IR] = cls._find_irs(
            sequence,
            out_dir,
            seq_name=seq_name,
        )

        n_irs_found: int = len(found_irs)
        seq_len: int = len(sequence)
        if n_irs_found == 0:  # Return sequence if no IRs found
            db_repr, obj_fn_value = "".join(["." for _ in range(seq_len)]), 0
            if save_performance:
                write_solver_performance_to_file(
                    db_repr, obj_fn_value, 0.0, seq_len, out_dir, cls.__name__
                )
            return db_repr, obj_fn_value

        # Define constraint programming problem and solve
        model, variables = cls._get_cp_model(
            found_irs, seq_len, sequence, out_dir, seq_name
        )

        solver: CpSolver = CpSolver()

        with tqdm(
            desc=f"{cls.__name__} - Running solver ({len(variables)} variables)"
        ) as _:
            status = solver.Solve(model)

        if status == OPTIMAL or status == FEASIBLE:
            # Return dot bracket repr and objective function's final value
            active_ir_idxs: List[int] = [
                int(re.findall(r"-?\d+\.?\d*", v.Name())[0])
                for v in variables
                if solver.Value(v) == 1 and "corrector" not in v.Name()
            ]
            db_repr: str = irs_to_dot_bracket(
                [found_irs[i] for i in active_ir_idxs], seq_len
            )
            obj_fn_value: float = solver.ObjectiveValue()

            if save_performance:
                dot_bracket_repr_mfe: float = calc_free_energy(
                    db_repr, sequence, out_dir, seq_name
                )
                write_solver_performance_to_file(
                    db_repr,
                    obj_fn_value,
                    dot_bracket_repr_mfe,
                    seq_len,
                    out_dir,
                    cls.__name__,
                    n_irs_found,
                    len(variables),
                    solver.WallTime(),
                    solver.NumBranches(),
                    solver.NumConflicts(),
                )
            return db_repr, obj_fn_value
        else:
            # The optimisation problem does not have an optimal solution
            db_repr, obj_fn_value = "".join(["." for _ in range(seq_len)]), 0
            if save_performance:
                write_solver_performance_to_file(
                    db_repr, obj_fn_value, 0.0, seq_len, out_dir, cls.__name__
                )
            return db_repr, obj_fn_value

    @staticmethod
    def _find_irs(
        sequence: str,
        out_dir: str = ".",
        *,
        seq_name: str = "seq",
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
                str(2),
                "-M",
                str(len(sequence)),
                "-g",
                str(len(sequence) - 1),
                "-x",
                str(0),
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
    def _get_cp_model(
        ir_list: List[IR],
        seq_len: int,
        sequence: str,
        out_dir: str,
        seq_name: str,
    ) -> Tuple[CpModel, List[IntVar]]:
        model: CpModel = CpModel()
        n_irs: int = len(ir_list)

        # Create binary indicator variables for IRs
        invalid_gap_sz_ir_idxs: List[int] = [
            i for i in range(n_irs) if not ir_has_valid_gap_size(ir_list[i])
        ]
        ir_indicator_variables = [
            model.NewIntVar(0, 1, f"ir_{i}")
            for i in range(n_irs)
            if i not in invalid_gap_sz_ir_idxs
        ]

        # If 1 or fewer variables, trivial or impossible optimisation problem, will be trivially handled by solver
        if len(ir_indicator_variables) <= 1:
            return model, ir_indicator_variables

        # Add XOR between IRs that are incompatible
        valid_ir_pairs, valid_idx_pairs = get_valid_gap_sz_ir_n_tuples(
            2, n_irs, ir_list, invalid_gap_sz_ir_idxs
        )
        incompatible_ir_pair_idxs: List[Tuple[int, int]] = [
            idx_pair
            for ir_pair, idx_pair in zip(valid_ir_pairs, valid_idx_pairs)
            if ir_pair_invalid_relative_pos(ir_pair[0], ir_pair[1])
        ]

        for ir_a_idx, ir_b_idx in incompatible_ir_pair_idxs:
            # Search required as some IR variables might not have had variables created as they were invalid so
            # list ordering of variables cannot be relied upon
            ir_a_var: IntVar = [
                var for var in ir_indicator_variables if str(ir_a_idx) in var.Name()
            ][0]
            ir_b_var: IntVar = [
                var for var in ir_indicator_variables if str(ir_b_idx) in var.Name()
            ][0]

            model.Add(ir_a_var + ir_b_var <= 1)

        # All constraints and the objective must have integer coefficients for CP-SAT solver

        # Obtain free energies of the IRs that are valid, they comprise the coefficients for ir vars
        variable_coefficients: List[float] = []
        for var in ir_indicator_variables:
            ir_idx: int = int(re.findall(r"-?\d+\.?\d*", var.Name())[0])
            ir_db_repr: str = irs_to_dot_bracket([ir_list[ir_idx]], seq_len)
            variable_coefficients.append(
                calc_free_energy(ir_db_repr, sequence, out_dir, seq_name)
            )

        # Define objective function
        obj_fn_expr = LinearExpr.WeightedSum(
            ir_indicator_variables, variable_coefficients
        )
        model.Minimize(obj_fn_expr)

        return model, ir_indicator_variables
