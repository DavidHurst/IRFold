from irfold import IRFoldBase


def test_conversion_lengths_match(
    list_of_irs, expected_dot_bracket_reprs, sequence_lengths
):
    for ir, expected, seq_len in zip(
        list_of_irs, expected_dot_bracket_reprs, sequence_lengths
    ):
        generated_db_repr = IRFoldBase.irs_to_dot_bracket([ir], seq_len)
        assert len(generated_db_repr) == len(expected)


def test_conversion_output_matches(
    list_of_irs, expected_dot_bracket_reprs, sequence_lengths
):
    for ir, expected_db_repr, seq_len in zip(
        list_of_irs, expected_dot_bracket_reprs, sequence_lengths
    ):
        generated_db_repr = IRFoldBase.irs_to_dot_bracket([ir], seq_len)
        assert generated_db_repr == expected_db_repr
