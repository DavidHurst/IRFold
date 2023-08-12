import random
import sys
import pandas as pd
from math import comb

from pathlib import Path

from helper_functions import (
    get_all_ir_pairs_not_matching_same_bases_valid_gap_sz,
    eval_ir_pair_structure_and_mfe,
)

sys.path.append(str(Path(__file__).resolve().parents[1]))

from irfold import IRFoldBase
from irfold.util import (
    ir_pair_wholly_nested,
    ir_pair_disjoint,
    ir_has_valid_gap_size,
)

DATA_DIR = (Path(__file__).parent.parent / "data").resolve()
EXPERIMENT_1_DATA_DIR = (DATA_DIR / "experiment_1").resolve()

if __name__ == '__main__':
    seq = 'UGAUGACA'
    irs = IRFoldBase.__find_irs(seq, out_dir=str(EXPERIMENT_1_DATA_DIR))
    for ir in irs:
        print(ir)