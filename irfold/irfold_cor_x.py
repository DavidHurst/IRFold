__all__ = ["IRFoldCorX"]

import itertools
import re

from irfold import IRFoldVal2
from typing import Tuple, List
from ortools.sat.python.cp_model import CpModel, IntVar, LinearExpr
from irfold.util import (
    ir_has_valid_gap_size,
    IR,
    irs_to_dot_bracket,
    calc_free_energy,
    get_valid_ir_n_tuples,
    irs_incompatible,
)


class IRFoldCorX(IRFoldVal2):

    # ToDo: Refactor parent classes to enable this
    # def __init__(self, max_n_tuple_sz_to_correct: int = 2):
    #     if max_n_tuple_sz_to_correct < 2:
    #         raise ValueError('Cannot corrct IR groups smaller than 2.')
    #     self.max_n_tuple_sz_to_correct = max_n_tuple_sz_to_correct
    #     self.__name__ = f'IRFoldCor{max_n_tuple_sz_to_correct}'

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

        if max_n_tuple_sz_to_correct < 2 or max_n_tuple_sz_to_correct > len(ir_list):
            raise ValueError('Cannot corrct IR groups smaller than 2 or larger than the number of IRs.')

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

        # Obtain free energies of the IRs that are valid, these will the coefficients for IR indicator vars
        ir_indicator_var_coeffs: List[int] = []
        for var in ir_indicator_variables:
            ir_idx: int = int(re.findall(r"-?\d+\.?\d*", var.Name())[0])
            ir_db_repr: str = irs_to_dot_bracket([ir_list[ir_idx]], seq_len)
            ir_free_energy: float = calc_free_energy(
                ir_db_repr, sequence, out_dir, seq_name
            )
            ir_indicator_var_coeffs.append(round(ir_free_energy))

        # Correct IR tuples of increasing size until max specified size to correct
        obj_fn_expressions: List[LinearExpr.WeightedSum] = []
        obj_fn_vars: List[IntVar] = []
        for tuple_sz in range(2, max_n_tuple_sz_to_correct + 1):
            valid_ir_n_tuples, valid_ir_idx_n_tuples = get_valid_ir_n_tuples(
                tuple_sz, n_irs, ir_list, invalid_gap_sz_ir_idxs
            )

            # Add XOR between IRs that are incompatible
            incompatible_ir_idx_n_tuples: List[Tuple[int, ...]] = [
                ir_idx_n_tuple
                for ir_n_tuple, ir_idx_n_tuple in zip(
                    valid_ir_n_tuples, valid_ir_idx_n_tuples
                )
                if irs_incompatible(ir_n_tuple)
            ]

            # Get ir indicator variables corresponding to incompatible IR indices found above
            for incomp_ir_idx_n_tuple in incompatible_ir_idx_n_tuples:
                # Search for indicator variables by name required as some IR validation process likely removed some IRs
                # so cannot rely on creation order to get IRs
                incomp_ir_variables: List[IntVar] = [
                    var
                    for var in ir_indicator_variables
                    if any(
                        [
                            True
                            for ir_idx in incomp_ir_idx_n_tuple
                            if str(ir_idx) in var.Name()
                        ]
                    )
                ]

                model.Add(
                    LinearExpr.Sum(incomp_ir_variables) <= tuple_sz
                )  # Prevents IRs being active simultaneously

            # Generate correction variable and coefficient for each IR n-tuple whose free energy additivity is erroneous
            (
                correction_vars,
                correction_var_coeffs,
            ) = IRFoldCorX.generate_ir_n_tuple_correction_variables_w_coeffs(
                model,
                seq_len,
                sequence,
                out_dir,
                seq_name,
                valid_ir_idx_n_tuples,
                valid_ir_n_tuples,
            )
            obj_fn_vars += correction_vars

            # Add constraints only activating correction variables if all variable in n-tuple are active
            IRFoldCorX.add_activation_constraints_for_correction_vars(
                correction_vars, ir_indicator_variables, model, tuple_sz
            )

            obj_fn_expressions.append(
                LinearExpr.WeightedSum(
                    correction_vars,
                    correction_var_coeffs,
                )
            )

        # Define objective function
        obj_fn_expressions.append(
            LinearExpr.WeightedSum(ir_indicator_variables, ir_indicator_var_coeffs)
        )
        model.Minimize(LinearExpr.Sum(obj_fn_expressions))

        return model, ir_indicator_variables + obj_fn_vars

    @staticmethod
    def generate_ir_n_tuple_correction_variables_w_coeffs(
        model: CpModel,
        seq_len: int,
        sequence: str,
        out_dir: str,
        seq_name: str,
        valid_ir_idx_n_tuples: List[Tuple[int, ...]],
        valid_ir_n_tuples: List[Tuple[IR, ...]],
    ) -> Tuple[List[IntVar], List[int]]:
        correction_vars: List[IntVar] = []
        correction_var_coeffs: List[int] = []

        # Add correction variable for each n-tuple that requires one
        for ir_idx_n_tuple, ir_n_tuple in zip(valid_ir_idx_n_tuples, valid_ir_n_tuples):
            if irs_incompatible(list(ir_n_tuple)):
                # ToDo: Pass list of compatible IRs to this function, already computed incompatible,
                #       set difference all with that to get compatible

                # This n-tuple will never be in final solution so seems efficient to not create correction var for it
                # but this might hurt search accuracy
                continue

            # Check if the addition of n-tuples' IR's free energies correctly represents the free energy of the n-tuple
            individual_ir_db_reprs: List[str] = [
                irs_to_dot_bracket([ir], seq_len) for ir in ir_n_tuple
            ]
            ir_n_tuple_db_repr: str = irs_to_dot_bracket(ir_n_tuple, seq_len)

            ir_n_tuple_additive_free_energy: float = round(
                sum(
                    [
                        calc_free_energy(ir_db_repr, sequence, out_dir, seq_name)
                        for ir_db_repr in individual_ir_db_reprs
                    ]
                ),
                4,
            )
            ir_n_tuple_true_free_energy: float = round(
                calc_free_energy(ir_n_tuple_db_repr, sequence, out_dir, seq_name), 4
            )
            free_energy_difference: float = round(
                ir_n_tuple_additive_free_energy - ir_n_tuple_true_free_energy, 4
            )

            # If additivity assumption does not hold i.e. diff > 0, add correction variable for n-tuple
            if abs(free_energy_difference) > 0:
                print(
                    f"IRFoldCorX - {ir_n_tuple} needs correcting, diff = {free_energy_difference}"
                )
                irs_idxs_corrected_repr = "".join(
                    [f"ir_{ir_idx}_" for ir_idx in ir_idx_n_tuple]
                )
                correction_var: IntVar = model.NewIntVar(
                    0, 1, f"{irs_idxs_corrected_repr}fe_corrector"
                )
                correction_vars.append(correction_var)

                # Add value required to make additive equal true free energy. Int rounding necessary for CpModel
                if ir_n_tuple_additive_free_energy > ir_n_tuple_true_free_energy:
                    correction_var_coeffs.append(-round(abs(free_energy_difference)))
                else:
                    correction_var_coeffs.append(round(abs(free_energy_difference)))

        return correction_vars, correction_var_coeffs

    @staticmethod
    def add_activation_constraints_for_correction_vars(
        correction_vars: List[IntVar],
        ir_indicator_variables: List[IntVar],
        model: CpModel,
        current_tuple_size: int,
    ) -> None:
        # Add constraints that only activate correction variables when the IR n-tuple the correction variable represents is active
        for cor_var in correction_vars:
            ir_idxs: List[str] = [
                ir_idx for ir_idx in re.findall(r"-?\d+\.?\d*", cor_var.Name())
            ]

            relevant_ir_indicator_vars: List[IntVar] = [
                var
                for var in ir_indicator_variables
                if any([True for ir_idx in ir_idxs if ir_idx in cor_var.Name()])
            ]

            # Create variable which is 1 only when all IRs are 1
            irs_idxs_corrected_repr = "".join([f"ir_{ir_idx}_" for ir_idx in ir_idxs])
            ir_n_tuple_is_active_var = model.NewBoolVar(
                f"{irs_idxs_corrected_repr}all_1"
            )
            model.Add(
                LinearExpr.Sum(relevant_ir_indicator_vars) == current_tuple_size
            ).OnlyEnforceIf(ir_n_tuple_is_active_var)
            model.Add(
                LinearExpr.Sum(relevant_ir_indicator_vars) != current_tuple_size
            ).OnlyEnforceIf(ir_n_tuple_is_active_var.Not())

            # Set IR n-tuple's correction variable to 1 only if all IR variables are active, 0 otherwise
            model.Add(cor_var == 1).OnlyEnforceIf(ir_n_tuple_is_active_var)
            model.Add(cor_var == 0).OnlyEnforceIf(ir_n_tuple_is_active_var.Not())
