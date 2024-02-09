import sys
import pytest
from pathlib import Path

from sequence_fixtures import SequenceA

sys.path.append(str(Path(__file__).resolve().parents[1]))

DATA_DIR = str(Path(__file__).parent / "tests_data")


@pytest.fixture(scope="module")
def data_dir():
    return DATA_DIR


@pytest.fixture(scope="module", params=[SequenceA.all_irs.values()])
def all_irs(request):
    return request.param


@pytest.fixture(scope="module", params=[list(SequenceA.all_irs.keys())])
def all_irs_names(request):
    return request.param


@pytest.fixture(scope="module", params=[SequenceA.all_irs_dot_bracket.values()])
def all_ir_dot_bracket_reprs(request):
    return request.param


@pytest.fixture(scope="module", params=[SequenceA.sequence])
def sequence(request):
    return request.param


@pytest.fixture(scope="module", params=[SequenceA.sequence_name])
def sequence_name(request):
    return request.param


@pytest.fixture(scope="module", params=[SequenceA.sequence_length])
def sequence_length(request):
    return request.param


@pytest.fixture(scope="module", params=[ir for ir in list(SequenceA.all_irs.values())])
def ir(request):
    return request.param


@pytest.fixture(
    scope="module",
    params=[db_repr for db_repr in SequenceA.all_irs_dot_bracket.values()],
)
def ir_dot_bracket_repr(request):
    return request.param


@pytest.fixture(scope="module", params=list(SequenceA.all_ir_pairs.values()))
def ir_pair(request):
    return request.param


@pytest.fixture(scope="module", params=[SequenceA.all_ir_pairs.values()])
def all_ir_pairs(request):
    return request.param


@pytest.fixture(scope="module", params=[SequenceA.invalid_gap_size_irs.values()])
def all_invalid_gap_size_irs(request):
    return request.param


@pytest.fixture(scope="module", params=SequenceA.valid_gap_size_irs.values())
def valid_gap_size_ir(request):
    return request.param


@pytest.fixture(scope="module", params=SequenceA.invalid_gap_size_irs.values())
def invalid_gap_size_ir(request):
    return request.param


@pytest.fixture(scope="module", params=SequenceA.not_nested_ir_pairs.values())
def not_nested_ir_pair(request):
    return request.param


@pytest.fixture(scope="module", params=SequenceA.wholly_nested_ir_pairs.values())
def wholly_nested_ir_pair(request):
    return request.param


@pytest.fixture(scope="module", params=SequenceA.co_located_ir_pairs.values())
def co_located_ir_pair(request):
    return request.param


@pytest.fixture(scope="module", params=[SequenceA.co_located_ir_pairs.values()])
def all_co_located_ir_pairs(request):
    return request.param


@pytest.fixture(scope="module", params=SequenceA.non_co_located_ir_pairs.values())
def non_co_located_ir_pair(request):
    return request.param


@pytest.fixture(scope="module", params=SequenceA.partially_nested_ir_pairs.values())
def partially_nested_ir_pair(request):
    return request.param


@pytest.fixture(scope="module", params=[SequenceA.partially_nested_ir_pairs.values()])
def all_partially_nested_ir_pairs(request):
    return request.param


@pytest.fixture(scope="module", params=[SequenceA.ir_indicator_variables_names])
def ir_indicator_variables_names(request):
    return request.param


# Clear tests_data dir out for clean test environment
for file_name in Path(DATA_DIR).glob(("**/*")):
    if not str(file_name).endswith(".gitignore"):
        Path.unlink(Path(DATA_DIR / file_name))
