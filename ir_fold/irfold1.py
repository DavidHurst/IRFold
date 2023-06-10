__all__ = ["IRFold1"]

import itertools

from ir_fold import IRFold0
from typing import Tuple, List
from ortools.linear_solver import pywraplp

IR = Tuple[Tuple[int, int], Tuple[int, int]]


class IRFold1(IRFold0):
    """Extends base IRFold model by validating found IRs before passing them to the solver
    and checks all IR pairs form valid loops when combined.
    """

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
        if IRFold1.irs_disjoint(ir_a, ir_b):
            return True

        # Case 2: Invalid loop
        bases_inbetween_parens = (
            earliest_right_string_base_idx - latest_left_string_base_idx - 1
        )
        if bases_inbetween_parens < 3:
            return False

        # Case 3: Valid loop
        return True

    @staticmethod
    # @override # Upgrade to python 3.12 for this, import from typing
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

        # Remove IRs with a gap < 3, these are sterically impossible
        valid_irs: List[IR] = [ir for ir in all_irs if ir[1][0] - ir[0][1] - 1 >= 3]

        # Add XOR between IRs that match the same bases
        unique_idx_pairs: List[Tuple[int, int]] = list(
            itertools.combinations([i for i in range(n_irs)], 2)
        )
        unique_ir_pairs: List[Tuple[IR, IR]] = [
            (valid_irs[i], valid_irs[j]) for i, j in unique_idx_pairs
        ]
        incompatible_ir_pair_idxs: List[Tuple[int, int]] = [
            idx_pair
            for ir_pair, idx_pair in zip(unique_ir_pairs, unique_idx_pairs)
            if IRFold1.irs_share_base_pair(ir_pair[0], ir_pair[1])
            and not IRFold1.irs_form_valid_loop(ir_pair[0], ir_pair[1])
        ]

        for inc_ir_a, inc_ir_b in incompatible_ir_pair_idxs:
            solver.Add(variables[inc_ir_a] + variables[inc_ir_b] <= 1)

        # Define objective function
        obj_fn = solver.Objective()
        for i in range(n_irs):
            obj_fn.SetCoefficient(variables[i], ir_free_energies[i])
        obj_fn.SetMinimization()

        return solver
