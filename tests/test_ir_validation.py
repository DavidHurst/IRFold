import pytest
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from ir_fold import IRFold0

DATA_DIR = "./tests_data"
