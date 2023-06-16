import itertools
import random
import re

import pytest

from irfold import IRFoldBase, IRFoldVal1


@pytest.mark.parametrize("irfold_variant", [IRFoldBase, IRFoldVal1])
def test_correct_ir_xor_constraints_generated(data_dir, irfold_variant):
    print(f"\n  testing {irfold_variant.__name__}")
    seq_len = 15
    seq = "AUGUAACAACCCGAC"
    irs = [
        ((1, 3), (10, 12)),
        ((2, 3), (13, 14)),
        ((2, 3), (8, 9)),
    ]
    ir_pairs_matching_same_bases = [
        (0, 1),
        (0, 2),
        (1, 2),
    ]

    solver = irfold_variant.get_solver(
        irs,
        seq_len,
        seq,
        data_dir,
        "test_irfold0_xor_constraints_seq",
    )

    for constraint, incompatible_ir_pair_idx in zip(
        solver.constraints(), ir_pairs_matching_same_bases
    ):
        ir_a_idx = incompatible_ir_pair_idx[0]
        ir_b_idx = incompatible_ir_pair_idx[1]

        xor_constraint_ir_idxs = re.findall(r"-?\d+\.?\d*", constraint.name())
        constraint_ir_a = int(xor_constraint_ir_idxs[0])
        constraint_ir_b = int(xor_constraint_ir_idxs[1])

        assert ir_a_idx == constraint_ir_a
        assert ir_b_idx == constraint_ir_b

    assert solver is not None
