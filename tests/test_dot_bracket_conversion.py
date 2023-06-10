import pytest
from pytest_lazyfixture import lazy_fixture
from ir_fold import IRFold0


@pytest.mark.parametrize(
    "irs, expected_db_reprs, seq_lens",
    [
        (
            lazy_fixture("ir_list"),
            lazy_fixture("expected_dot_bracket_reprs"),
            lazy_fixture("sequence_lengths"),
        )
    ],
)
def test_ir_to_db_repr(irs, expected_db_reprs, seq_lens):
    generated_db_repr = IRFold0.irs_to_dot_bracket([irs], seq_lens)
    assert generated_db_repr == expected_db_reprs
