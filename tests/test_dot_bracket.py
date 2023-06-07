import pytest
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from src import IRFold


@pytest.fixture
def db_repr_two_bps():
    return "..((..)).."


@pytest.fixture
def db_repr_two_bps_at_end():
    return "((....))"


@pytest.fixture
def db_repr_three_bps_one_nested():
    return ".((..(..).)).."


def test_always_passes():
    assert True


def test_always_fails():
    assert False
