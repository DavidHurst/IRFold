__all__ = ["IRFold2"]

import itertools
import re

from irfold import IRFold1, IRFold0
from typing import Tuple, List
from ortools.linear_solver import pywraplp

IR = Tuple[Tuple[int, int], Tuple[int, int]]


class IRFold2(IRFold1):
    """Extends IRFold1 by adding solver variables to correct for the additivity of free energy IRs not holding
    consistently"""

    @staticmethod
    def get_solver(
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
        ir_indicator_variables = [
            solver.IntVar(0, 1, f"ir_{i}")
            for i in range(n_irs)
            if i not in invalid_gap_ir_idxs
        ]

        # If 1 or fewer variables, trivial or impossible optimisation problem, will be handled by ortools
        if len(ir_indicator_variables) <= 1:
            return solver

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
            if IRFold1.ir_pair_incompatible(ir_pair[0], ir_pair[1])
        ]

        # Add ir indicator variables and their coefficients (free energies of IRs)
        obj_fn = solver.Objective()
        for var in ir_indicator_variables:
            ir_idx: int = int(re.findall(r"-?\d+\.?\d*", var.name())[0])
            ir_db_repr: str = IRFold0.irs_to_dot_bracket([ir_list[ir_idx]], seq_len)
            ir_free_energy: float = IRFold0.calc_free_energy(
                ir_db_repr, sequence, out_dir, seq_name
            )

            obj_fn.SetCoefficient(var, ir_free_energy)

        # Add variable for each IR pair that forms an invalid loop which corrects the pairs additive free energy
        ir_pair_fe_correction_indicator_vars = []
        for (ir_a_idx, ir_b_idx), (ir_a, ir_b) in zip(valid_idx_pairs, valid_ir_pairs):
            if IRFold1.ir_pair_incompatible(ir_a, ir_b):
                continue
            ir_a_db_repr: str = IRFold0.irs_to_dot_bracket([ir_list[ir_a_idx]], seq_len)
            ir_b_db_repr: str = IRFold0.irs_to_dot_bracket([ir_list[ir_b_idx]], seq_len)
            ir_pair_db_repr: str = IRFold0.irs_to_dot_bracket(
                [ir_list[ir_a_idx], ir_list[ir_b_idx]], seq_len
            )
            print(
                f"ir_{ir_a_idx} and ir_{ir_b_idx} form invalid loop, adding correction variable"
            )

            ir_pair_additive_free_energy: float = sum(
                [
                    IRFold0.calc_free_energy(ir_db_repr, sequence, out_dir, seq_name)
                    for ir_db_repr in [ir_a_db_repr, ir_b_db_repr]
                ]
            )
            ir_pair_db_repr_free_energy: float = IRFold0.calc_free_energy(
                ir_pair_db_repr, sequence, out_dir, seq_name
            )

            free_energy_difference: float = (
                ir_pair_additive_free_energy - ir_pair_db_repr_free_energy
            )

            if free_energy_difference > 0:
                correction_var = solver.IntVar(
                    0, 1, f"ir_{ir_a_idx}_ir_{ir_b_idx}_fe_corrector"
                )
                ir_pair_fe_correction_indicator_vars.append(correction_var)

                obj_fn.SetCoefficient(correction_var, -free_energy_difference)

        print(f"after correction var generation")

        # Add XOR constraint between incompatible IR pairs
        for ir_a_idx, ir_b_idx in incompatible_ir_pair_idxs:
            ir_a_var = solver.LookupVariable(f"ir_{ir_a_idx}")
            ir_b_var = solver.LookupVariable(f"ir_{ir_b_idx}")
            solver.Add(
                ir_a_var + ir_b_var <= 1,
                f"ir_{ir_a_idx}_XOR_ir_{ir_b_idx}",
            )

        # Add constraints that only activate correction variables when the IR pair they represent is active
        for correction_var in ir_pair_fe_correction_indicator_vars:
            print(f"Adding activate {correction_var.name()}")
            ir_idxs: List[str] = re.findall(r"-?\d+\.?\d*", correction_var.name())
            ir_a_idx: int = int(ir_idxs[0])
            ir_b_idx: int = int(ir_idxs[1])

            ir_a_var = solver.LookupVariable(f"ir_{ir_a_idx}")
            ir_b_var = solver.LookupVariable(f"ir_{ir_b_idx}")

            # ir_pair_is_active_var == (ir_a + ir_b >= 2)
            ir_pair_is_active_var = solver.BoolVar(
                f"ir_{ir_a_idx}_ir_{ir_b_idx}_both_1"
            )
            solver.Add(ir_a_var + ir_b_var >= 2).OnlyEnforceIf(ir_pair_is_active_var)
            solver.Add(ir_a_var + ir_b_var < 2).OnlyEnforceIf(
                ir_pair_is_active_var.Not()
            )

            # ir_pair_is_active_var implies correction_var == 1
            solver.Add(correction_var == 1).OnlyEnforceIf(ir_pair_is_active_var)

            # not(ir_pair_is_active_var) implies correction_var == 0
            solver.Add(correction_var == 0).OnlyEnforceIf(ir_pair_is_active_var.Not())

        obj_fn.SetMinimization()

        return solver
