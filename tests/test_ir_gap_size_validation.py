from irfold.util import ir_has_valid_gap_size


def test_ir_gap_size_validation_invalid_irs(list_of_ir_pairs_with_invalid_gap_sizes):
    for ir in list_of_ir_pairs_with_invalid_gap_sizes:
        assert ir_has_valid_gap_size(ir) == False

def test_ir_gap_size_validation_valid_irs(list_of_ir_pairs_with_valid_gap_sizes):
    for ir in list_of_ir_pairs_with_valid_gap_sizes:
        assert ir_has_valid_gap_size(ir) == True
