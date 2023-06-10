import pytest
from pytest_lazyfixture import lazy_fixture
from ir_fold import IRFold0


@pytest.mark.parametrize(
    "irs, expected_db_reprs, seq_lens",
    [
        (
            lazy_fixture("list_of_irs"),
            lazy_fixture("expected_dot_bracket_reprs"),
            lazy_fixture("sequence_lengths"),
        )
    ],
)
def test_conversion_lengths_match(irs, expected_db_reprs, seq_lens):
    for ir, expected, s_len in zip(irs, expected_db_reprs, seq_lens):
        generated_db_repr = IRFold0.irs_to_dot_bracket([ir], s_len)
        assert len(generated_db_repr) == len(expected)


@pytest.mark.parametrize(
    "irs, expected_db_reprs, seq_lens",
    [
        (
            lazy_fixture("list_of_irs"),
            lazy_fixture("expected_dot_bracket_reprs"),
            lazy_fixture("sequence_lengths"),
        )
    ],
)
def test_conversion_output_matches(irs, expected_db_reprs, seq_lens):
    for ir, expected, s_len in zip(irs, expected_db_reprs, seq_lens):
        generated_db_repr = IRFold0.irs_to_dot_bracket([ir], s_len)
        assert generated_db_repr == expected
