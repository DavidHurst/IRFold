from irfold.util import ir_pair_invalid_relative_pos


def test_ir_pair_invalid_relative_pos(
    all_co_located_ir_pairs, all_partially_nested_ir_pairs, all_ir_pairs
):
    all_ir_pairs = list(all_ir_pairs)
    all_invalid_relative_pos_ir_pairs = list(all_partially_nested_ir_pairs) + list(
        all_co_located_ir_pairs
    )

    invalid_relative_pos_irs = []
    for pair in all_ir_pairs:
        if ir_pair_invalid_relative_pos(pair[0], pair[1]):
            invalid_relative_pos_irs.append(pair)

    assert len(invalid_relative_pos_irs) == len(all_invalid_relative_pos_ir_pairs)
