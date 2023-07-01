__all__ = ["IRFoldVal1"]

import itertools
import re

from ortools.sat.python.cp_model import CpModel, IntVar, LinearExpr
from typing import Tuple, List
from irfold import IRFoldBase
from irfold.util import (
    ir_has_valid_gap_size,
    IR,
    irs_to_dot_bracket,
    calc_free_energy,
    ir_pair_match_same_bases,
    ir_pair_forms_valid_loop,
)


class IRFoldVal1(IRFoldBase):
    """Extends base IRFold model by validating found IRs in pairs before passing them to the solver."""

    @staticmethod
    def get_solver(
        ir_list: List[IR],
        seq_len: int,
        sequence: str,
        out_dir: str,
        seq_name: str,
        max_n_tuple_sz_to_correct: int = 2
    ) -> Tuple[CpModel, List[IntVar]]:
        model: CpModel = CpModel()
        n_irs: int = len(ir_list)

        # Create binary indicator variables for IRs
        invalid_gap_ir_idxs: List[int] = [
            i for i in range(n_irs) if not ir_has_valid_gap_size(ir_list[i])
        ]
        ir_indicator_variables = [
            model.NewIntVar(0, 1, f"ir_{i}")
            for i in range(n_irs)
            if i not in invalid_gap_ir_idxs
        ]

        # If 1 or fewer variables, trivial or impossible optimisation problem, will be trivially handled by solver
        if len(ir_indicator_variables) <= 1:
            return model, ir_indicator_variables

        unique_possible_idx_pairs: List[Tuple[int, int]] = list(
            itertools.combinations([i for i in range(n_irs)], 2)
        )
        valid_idx_pairs: List[Tuple[int, int]] = [
            pair
            for pair in unique_possible_idx_pairs
            if pair[0] not in invalid_gap_ir_idxs and pair[1] not in invalid_gap_ir_idxs
        ]
        valid_ir_pairs: List[Tuple[IR, IR]] = [
            (ir_list[i], ir_list[j]) for i, j in valid_idx_pairs
        ]
        incompatible_ir_pair_idxs: List[Tuple[int, int]] = [
            idx_pair
            for ir_pair, idx_pair in zip(valid_ir_pairs, valid_idx_pairs)
            if IRFoldBase.ir_pair_incompatible(ir_pair[0], ir_pair[1])
        ]

        # Add XOR between IRs that are incompatible
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

    # @staticmethod
    # def ir_pair_incompatible(ir_a: IR, ir_b: IR) -> bool:
    #     return ir_pair_match_same_bases(ir_a, ir_b) or not ir_pair_forms_valid_loop(
    #         ir_a, ir_b
    #     )
