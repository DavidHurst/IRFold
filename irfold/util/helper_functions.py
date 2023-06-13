from typing import List, Tuple

from pathlib import Path

import RNA

IR = Tuple[Tuple[int, int], Tuple[int, int]]


def irs_to_dot_bracket(irs: List[IR], seq_len: int) -> str:
    """Does not support mismatches."""
    paired_bases: List[str] = ["." for _ in range(seq_len)]  # Initially none paired

    for ir in irs:
        n_base_pairs: int = ir[0][1] - ir[0][0] + 1  # Assumes no mismatches

        left_strand: Tuple[int, int]
        right_strand: Tuple[int, int]
        left_strand, right_strand = ir[0], ir[1]

        paired_bases[left_strand[0] : left_strand[1] + 1] = [
            "(" for _ in range(n_base_pairs)
        ]

        paired_bases[right_strand[0] : right_strand[1] + 1] = [
            ")" for _ in range(n_base_pairs)
        ]

    return "".join(paired_bases)


def calc_free_energy(
    dot_brk_repr: str, sequence: str, out_dir: str, seq_name: str = "seq"
) -> float:
    out_dir_path: Path = Path(out_dir).resolve()
    if not out_dir_path.exists():
        out_dir_path = Path.cwd().resolve()

    out_file: str = str(out_dir_path / f"{seq_name}_calculated_ir_energies.txt")

    with open(out_file, "a") as file:
        file.write(f"Evaluating IR:\n")
        for i in range(len(dot_brk_repr)):
            file.write(f"{i + 1:<3}")
        file.write("\n")
        for b in dot_brk_repr:
            file.write(f"{b:<3}")
        file.write("\n")

        free_energy = RNA.eval_structure_simple(sequence, dot_brk_repr, 1, file)
        file.write(f"\n\n")

    return free_energy


def create_seq_file(seq: str, seq_name: str, file_name: str) -> None:
    with open(file_name, "w") as file:
        file.write(f">{seq_name}\n")
        file.write(seq)
