__all__ = ["IRFold2"]

import itertools

from irfold import IRFold1
from typing import Tuple, List
from ortools.linear_solver import pywraplp

IR = Tuple[Tuple[int, int], Tuple[int, int]]


class IRFold2(IRFold1):
    """Extends IRFold1 by adding solver variables to correct for the additivity of free energy IRs not holding
    consistently"""

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

        # Remove IRs with a gap < 3, these are sterically impossible
        valid_irs: List[IR] = [ir for ir in ir_list if ir[1][0] - ir[0][1] - 1 >= 3]
        n_irs: int = len(valid_irs)

        # Generate all pairs of valid IRs
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

        # Create variables representing IR inclusion in final structure and variables to correct free energy additivity
        ir_indicator_vars, _ = IRFold2.__get_solver_variables(
            valid_irs,
            solver,
            sequence,
            out_dir,
            seq_name,
            unique_idx_pairs,
            unique_ir_pairs,
        )

        # Add XOR constraint between IR pairs that are incompatible
        for inc_ir_a, inc_ir_b in incompatible_ir_pair_idxs:
            solver.Add(ir_indicator_vars[inc_ir_a] + ir_indicator_vars[inc_ir_b] <= 1)

        # Evaluate free energy of each IR respectively
        db_reprs: List[str] = [
            IRFold2.irs_to_dot_bracket([valid_irs[i]], seq_len) for i in range(n_irs)
        ]
        ir_free_energies: List[float] = [
            IRFold2.calc_free_energy(db_repr, sequence, out_dir, seq_name)
            for db_repr in db_reprs
        ]

        # Define objective function
        obj_fn = solver.Objective()
        for i in range(n_irs):
            obj_fn.SetCoefficient(ir_indicator_vars[i], ir_free_energies[i])
        obj_fn.SetMinimization()

        return solver

    @staticmethod
    def __get_solver_variables(
        ir_list: List[IR],
        solver: pywraplp.Solver,
        sequence: str,
        out_dir: str,
        seq_name: str,
        unique_idx_pairs: List[Tuple[int, int]],
        unique_ir_pairs: List[Tuple[IR, IR]],
    ):
        n_irs: int = len(ir_list)

        # Create binary indicator variables
        ir_binary_indicators: List[pywraplp.Solver.IntVar] = [
            solver.IntVar(0, 1, f"ir_{i}_indicator") for i in range(n_irs)
        ]

        # Create binary variables representing the correction
        # constant needed to remedy each ir pair's additive free energy

        # ToDo: Need a correction variable for ir pairs that assumption doesn't hold for, check experiment 2 results analysis
        # for list
        ir_pair_fe_correction_binary_indicators: List[pywraplp.Solver.IntVar] = [
            solver.IntVar(0, 1, f"ir_{i}_ir_{j}_fe_correction_indicator")
            for i, j in unique_idx_pairs
        ]

        return ir_binary_indicators, ir_pair_fe_correction_binary_indicators
