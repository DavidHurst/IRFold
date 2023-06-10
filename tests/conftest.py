import sys

from fixtures import *
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

print(f"\n >>>> Ran conftest file")
