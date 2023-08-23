__all__ = ["IRFoldVal1"]

import re

from ortools.sat.python.cp_model import CpModel, IntVar, LinearExpr
from typing import Tuple, List
from irfold import IRFoldBase
from irfold.util import (
    ir_has_valid_gap_size,
    IR,
    irs_to_dot_bracket,
    calc_free_energy,
)


class IRFoldVal1(IRFoldBase):
    """Extends base IRFold model by validating found IRs in pairs before passing them to the solver."""

    @staticmethod
    def _get_cp_model(
        ir_list: List[IR],
        seq_len: int,
        sequence: str,
        out_dir: str,
        seq_name: str,
        max_n_tuple_sz_to_correct: int = 3,
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

        # All constraints and the objective must have integer coefficients for CP-SAT solver
        # Obtain free energies of the IRs that are valid, they comprise the coefficients for ir vars
        variable_coefficients: List[int] = []
        for var in ir_indicator_variables:
            ir_idx: int = int(re.findall(r"-?\d+\.?\d*", var.Name())[0])
            ir_db_repr: str = irs_to_dot_bracket([ir_list[ir_idx]], seq_len)
            ir_free_energy: float = calc_free_energy(
                ir_db_repr, sequence, out_dir, seq_name
            )
            variable_coefficients.append(round(ir_free_energy))

        # Define objective function
        obj_fn_expr = LinearExpr.WeightedSum(
            ir_indicator_variables, variable_coefficients
        )
        model.Minimize(obj_fn_expr)

        return model, ir_indicator_variables
