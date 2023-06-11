__all__ = ["IRFold1"]

import itertools

from irfold import IRFold0
from typing import Tuple, List
from ortools.linear_solver import pywraplp

IR = Tuple[Tuple[int, int], Tuple[int, int]]


class IRFold1(IRFold0):
    """Extends base IRFold model by validating found IRs in pairs before passing them to the solver."""

    @staticmethod
    def get_lp_solver(
        ir_list: List[IR],
        ir_free_energies: List[float],
        backend: str = "SCIP",
    ) -> pywraplp.Solver:
        # Create solver
        solver: pywraplp.Solver = pywraplp.Solver.CreateSolver(backend)
        if solver is None:
            raise Exception("Failed to create solver.")

        # Remove IRs with a gap < 3, these are sterically impossible
        valid_irs: List[IR] = [ir for ir in ir_list if ir[1][0] - ir[0][1] - 1 >= 3]
        n_irs: int = len(valid_irs)

        # Create binary indicator variables
        variables = [solver.IntVar(0, 1, f"ir_{i}") for i in range(n_irs)]

        # Add XOR constraint between IR pairs that are incompatible
        unique_idx_pairs: List[Tuple[int, int]] = list(
            itertools.combinations([i for i in range(n_irs)], 2)
        )
        unique_ir_pairs: List[Tuple[IR, IR]] = [
            (valid_irs[i], valid_irs[j]) for i, j in unique_idx_pairs
        ]
        incompatible_ir_pair_idxs: List[Tuple[int, int]] = [
            idx_pair
            for ir_pair, idx_pair in zip(unique_ir_pairs, unique_idx_pairs)
            if IRFold1.ir_pair_incompatible(ir_pair[0], ir_pair[1])
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

        # Case 1: Disjoint
        if IRFold1.ir_pair_disjoint(ir_a, ir_b):
            return True

        # Case 2: Invalid hairpin
        bases_inbetween_parens = (
            earliest_right_string_base_idx - latest_left_string_base_idx - 1
        )
        if bases_inbetween_parens < 3:
            return False

        # Case 3: Valid hairpin / loop?
        # hp_region_dp_repr = [" " for _ in range(seq_len)]
        # hp_region_dp_repr[latest_left_string_base_idx] = "("
        # hp_region_dp_repr[earliest_right_string_base_idx] = ")"
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
