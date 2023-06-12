__all__ = ["IRFold1"]

import itertools
import re

from irfold import IRFold0
from typing import Tuple, List
from ortools.linear_solver import pywraplp

IR = Tuple[Tuple[int, int], Tuple[int, int]]


class IRFold1(IRFold0):
    """Extends base IRFold model by validating found IRs in pairs before passing them to the solver."""

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

        # Create binary indicator variables for IRs
        invalid_gap_ir_idxs: List[int] = [
            i for i in range(n_irs) if not IRFold1.ir_has_valid_gap_size(ir_list[i])
        ]
        variables = [
            solver.IntVar(0, 1, f"ir_{i}")
            for i in range(n_irs)
            if i not in invalid_gap_ir_idxs
        ]

        # If 1 or fewer variables, trivial or impossible optimisation problem
        if len(variables) <= 1:
            return solver

        # Add XOR between IRs that match the same bases
        unique_possible_idx_pairs: List[Tuple[int, int]] = list(
            itertools.combinations([i for i in range(n_irs)], 2)
        )
        valid_idx_pairs: List[Tuple[int, int]] = [
            pair
            for pair in unique_possible_idx_pairs
            if pair[0] not in invalid_gap_ir_idxs and pair[1] not in invalid_gap_ir_idxs
        ]

        unique_ir_pairs: List[Tuple[IR, IR]] = [
            (ir_list[i], ir_list[j]) for i, j in valid_idx_pairs
        ]
        incompatible_ir_pair_idxs: List[Tuple[int, int]] = [
            idx_pair
            for ir_pair, idx_pair in zip(unique_ir_pairs, valid_idx_pairs)
            if IRFold1.ir_pair_incompatible(ir_pair[0], ir_pair[1])
        ]

        for ir_a_idx, ir_b_idx in incompatible_ir_pair_idxs:
            ir_a_var = solver.LookupVariable(f"ir_{ir_a_idx}")
            ir_b_var = solver.LookupVariable(f"ir_{ir_b_idx}")
            solver.Add(
                ir_a_var + ir_b_var <= 1,
                f"ir_{ir_a_idx}_XOR_ir_{ir_b_idx}",
            )

        # Define objective function
        obj_fn = solver.Objective()
        for var in variables:
            ir_idx: int = int(re.findall(r"-?\d+\.?\d*", var.name())[0])
            ir_db_repr: str = IRFold0.irs_to_dot_bracket([ir_list[ir_idx]], seq_len)
            ir_free_energy: float = IRFold0.calc_free_energy(
                ir_db_repr, sequence, out_dir, seq_name
            )

            obj_fn.SetCoefficient(var, ir_free_energy)
        obj_fn.SetMinimization()

        return solver

    @staticmethod
    def ir_has_valid_gap_size(ir):
        left_strand_end_idx: int = ir[0][1]
        right_strand_start_idx: int = ir[1][0]

        return right_strand_start_idx - left_strand_end_idx - 1 >= 3

    @staticmethod
    def ir_pair_incompatible(ir_a: IR, ir_b: IR) -> bool:
        return IRFold0.ir_pair_match_same_bases(
            ir_a, ir_b
        ) or not IRFold1.ir_pair_forms_valid_loop(ir_a, ir_b)

    @staticmethod
    def ir_pair_forms_valid_loop(ir_a, ir_b):
        if IRFold1.ir_pair_wholly_nested(ir_a, ir_b) or IRFold1.ir_pair_wholly_nested(
            ir_b, ir_a
        ):
            return True

        ir_a_left_strand: Tuple[int, int]
        ir_b_left_strand: Tuple[int, int]

        ir_a_left_strand, ir_b_left_strand = ir_a[0], ir_b[0]

        latest_left_string_base_idx = (
            ir_a_left_strand[1]
            if ir_a_left_strand[1] > ir_b_left_strand[1]
            else ir_b_left_strand[1]
        )
        earliest_right_string_base_idx = (
            ir_a[1][0] if ir_a[1][0] < ir_b[1][0] else ir_b[1][0]
        )

        if IRFold1.ir_pair_disjoint(ir_a, ir_b):
            return True

        # Invalid loop
        bases_inbetween_parens = (
            earliest_right_string_base_idx - latest_left_string_base_idx - 1
        )
        if bases_inbetween_parens < 3:
            return False

        return True

    @staticmethod
    def ir_pair_disjoint(ir_a, ir_b):
        # ir_a comes entirely before ir_b
        ir_a_strictly_before_ir_b = ir_a[1][1] < ir_b[0][0]

        # ir_b comes entirely before ir_a
        ir_b_strictly_before_ir_a = ir_b[1][1] < ir_a[0][0]

        return ir_a_strictly_before_ir_b or ir_b_strictly_before_ir_a

    @staticmethod
    def ir_pair_wholly_nested(outer_ir, nested_ir):
        return outer_ir[0][0] < nested_ir[0][0] and nested_ir[1][1] < outer_ir[1][0]
