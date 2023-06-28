import sys
import pytest
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))


def setup_module():
    print("setup_module")


@pytest.fixture(scope="module")
def data_dir():
    return str(Path(__file__).parent / "tests_data")


@pytest.fixture(scope="module")
def rna_seq_15_bases_3_irs():
    return "CACCACCAUAAGGCU", 3


@pytest.fixture(scope="module")
def rna_seq_30_bases_19_irs():
    return "UGAUGACAAAUGCUUAACCCAAGCACGGCA", 19


@pytest.fixture(scope="module")
def rna_seq_15_bases_0_irs():
    return "AAAAAAAAAAAAAAA", 0


@pytest.fixture(scope="module")
def find_irs_params(rna_seq_15_bases_3_irs, data_dir):
    seq = rna_seq_15_bases_3_irs[0]
    seq_len = len(seq)
    return {
        "sequence": seq,
        "min_len": 2,
        "max_len": seq_len,
        "max_gap": seq_len - 1,
        "mismatches": 0,
        "seq_name": "test_irs_seq",
        "out_dir": data_dir,
    }


@pytest.fixture(scope="module")
def list_of_irs(rna_seq_15_bases_3_irs):
    return [((2, 3), (11, 12)), ((5, 6), (11, 12)), ((10, 11), (13, 14))]


@pytest.fixture(scope="module")
def expected_dot_bracket_reprs(rna_seq_15_bases_3_irs):
    return ["..((.......))..", ".....((....))..", "..........((.))"]


@pytest.fixture(scope="module")
def sequence_lengths(expected_dot_bracket_reprs):
    return [len(db_repr) for db_repr in expected_dot_bracket_reprs]


# ====================  IR sample definitions   ====================

# ToDo: Convert these fixtures to be classes where each class represents a sequence and stores IRs, IR pairs etc.
@pytest.fixture(
    scope="module",
    params=[
        ((0, 1), (28, 29)),  # Valid gap size = True
        ((0, 1), (23, 24)),  # Valid gap size = True
        ((0, 1), (19, 20)),  # Valid gap size = True
        ((0, 1), (6, 7)),  # Valid gap size = True
        ((2, 3), (9, 10)),  # Valid gap size = True
        ((3, 4), (28, 29)),  # Valid gap size = True
        ((3, 4), (23, 24)),  # Valid gap size = True
        ((3, 4), (19, 20)),  # Valid gap size = True
        ((3, 4), (6, 7)),  # Valid gap size = False
        ((6, 7), (10, 11)),  # Valid gap size = False
        ((7, 8), (13, 14)),  # Valid gap size = True
        ((8, 9), (13, 14)),  # Valid gap size = True
        ((10, 12), (27, 29)),  # Valid gap size = True
        ((10, 14), (20, 24)),  # Valid gap size = True
        ((10, 11), (19, 20)),  # Valid gap size = True
        ((13, 14), (15, 16)),  # Valid gap size = False
        ((17, 18), (26, 27)),  # Valid gap size = True
        ((18, 19), (26, 27)),  # Valid gap size = True
        ((22, 23), (27, 28)),  # Valid gap size = True
    ],
)
def ir(rna_seq_30_bases_19_irs, request):
    return request.param


@pytest.fixture(
    scope="module",
    params=[
        ((3, 4), (6, 7)),  # Valid gap size = False
        ((6, 7), (10, 11)),  # Valid gap size = False
        ((13, 14), (15, 16)),  # Valid gap size = False
    ],
)
def invalid_gap_size_ir(rna_seq_30_bases_19_irs, request):
    return request.param


@pytest.fixture(
    scope="module",
    params=[
        ((0, 1), (28, 29)),  # Valid gap size = True
        ((0, 1), (23, 24)),  # Valid gap size = True
        ((0, 1), (19, 20)),  # Valid gap size = True
        ((0, 1), (6, 7)),  # Valid gap size = True
        ((2, 3), (9, 10)),  # Valid gap size = True
        ((3, 4), (28, 29)),  # Valid gap size = True
        ((3, 4), (23, 24)),  # Valid gap size = True
        ((3, 4), (19, 20)),  # Valid gap size = True
        ((7, 8), (13, 14)),  # Valid gap size = True
        ((8, 9), (13, 14)),  # Valid gap size = True
        ((10, 12), (27, 29)),  # Valid gap size = True
        ((10, 14), (20, 24)),  # Valid gap size = True
        ((10, 11), (19, 20)),  # Valid gap size = True
        ((17, 18), (26, 27)),  # Valid gap size = True
        ((18, 19), (26, 27)),  # Valid gap size = True
        ((22, 23), (27, 28)),  # Valid gap size = True
    ],
)
def valid_gap_size_ir(rna_seq_30_bases_19_irs, request):
    return request.param


@pytest.fixture(
    scope="module",
    params=[
        (((0, 1), (28, 29)), ((2, 3), (9, 10))),
        (((0, 1), (28, 29)), ((3, 4), (23, 24))),
        (((0, 1), (28, 29)), ((3, 4), (19, 20))),
        (((0, 1), (28, 29)), ((7, 8), (13, 14))),
        (((0, 1), (28, 29)), ((8, 9), (13, 14))),
    ],
)
def not_matching_same_bases_ir_pair(rna_seq_30_bases_19_irs, request):
    return request.param


@pytest.fixture(
    scope="module",
    params=[
        (((0, 1), (28, 29)), ((0, 1), (23, 24))),  # Both match 0,1
        (((0, 1), (28, 29)), ((0, 1), (19, 20))),  # Both match 0,1
        (((0, 1), (28, 29)), ((0, 1), (6, 7))),  # Both match  0,1
        (((0, 1), (28, 29)), ((3, 4), (28, 29))),  # Both match 28,29
        (((0, 1), (28, 29)), ((10, 12), (27, 29))),  # Both match 28,29
    ],
)
def matching_same_bases_ir_pair(rna_seq_30_bases_19_irs, request):
    return request.param


@pytest.fixture(
    scope="module",
    params=[
        (((3, 4), (28, 29)), ((7, 8), (13, 14))),  # Second inside first
        (((3, 4), (28, 29)), ((8, 9), (13, 14))),  # Second inside first
        (((3, 4), (28, 29)), ((10, 11), (19, 20))),  # Second inside first
        (((0, 1), (19, 20)), ((8, 9), (13, 14))),  # Second inside first
        (((0, 1), (19, 20)), ((2, 3), (9, 10))),  # Second inside first
    ],
)
def wholly_nested_ir_pair(rna_seq_30_bases_19_irs, request):
    return request.param


@pytest.fixture(
    scope="module",
    params=[
        (((8, 9), (13, 14)), ((17, 18), (26, 27))),  # First comes before second
        (((7, 8), (13, 14)), ((18, 19), (26, 27))),  # First comes before second
        (((3, 4), (19, 20)), ((22, 23), (27, 28))),  # First comes before second
        (((2, 3), (9, 10)), ((22, 23), (27, 28))),  # First comes before second
        (((2, 3), (9, 10)), ((18, 19), (26, 27))),  # First comes before second
    ],
)
def entirely_disjoint_ir_pair(rna_seq_30_bases_19_irs, request):
    return request.param


@pytest.fixture(
    scope="module",
    params=[
        (
            ((3, 4), (23, 24)),
            ((10, 12), (27, 29)),
        ),  # Num bases between IR1 right and IR2 left = 10
        (
            ((3, 4), (23, 24)),
            ((17, 18), (26, 27)),
        ),  # Num bases between IR1 right and IR2 left = 4
        (
            ((3, 4), (23, 24)),
            ((18, 19), (26, 27)),
        ),  # Num bases between IR1 right and IR2 left = 3
        (
            ((3, 4), (19, 20)),
            ((10, 12), (27, 29)),
        ),  # Num bases between IR1 right and IR2 left = 6
    ],
)
def valid_num_bases_in_pair_intersection_ir_pair(rna_seq_30_bases_19_irs, request):
    """Invalid loop empirically determined to be when the
    number of bases enclosed by left_a & right_b or right_a & left_b
    is less than 3. Similarly to gap size validation.
    """
    return request.param


@pytest.fixture(
    scope="module",
    params=[
        (
            ((0, 1), (19, 20)),
            ((17, 18), (26, 27)),
        ),  # Num bases between IR1 right and IR2 left = 0
        (
            ((0, 1), (6, 7)),
            ((2, 3), (9, 10)),
        ),  # Num bases between IR1 right and IR2 left = 2
        (
            ((0, 1), (6, 7)),
            ((3, 4), (28, 29)),
        ),  # Num bases between IR1 right and IR2 left = 1
        (
            ((0, 1), (6, 7)),
            ((3, 4), (23, 24)),
        ),  # Num bases between IR1 right and IR2 left = 1
        (
            ((0, 1), (6, 7)),
            ((3, 4), (19, 20)),
        ),  # Num bases between IR1 right and IR2 left = 1
    ],
)
def invalid_num_bases_in_pair_intersection_ir_pair(rna_seq_30_bases_19_irs, request):
    """Invalid loop empirically determined to be when the
    number of bases enclosed by left_a & right_b or right_a & left_b
    is less than 3. Similarly to gap size validation.
    """
    return request.param
