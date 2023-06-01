from IUPACpal import find_inverted_repeats, config
import re
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi.graph.bulge_graph as fgb
import forgi
import itertools
import RNA


def shape_state_to_dot_bracket(included_irs, all_irs, seq_len):
    paired_bases = ["." for _ in range(seq_len)]  # Initially none paired

    for ir_idx in included_irs:
        # print(f'    Adding IR#{ir_idx} -> {all_irs[ir_idx]}')
        ir = all_irs[ir_idx]
        lhs = ir[0]
        rhs = ir[1]

        for i in range(lhs[0] - 1, lhs[1]):
            paired_bases[i] = "("

        for j in range(rhs[0] - 1, rhs[1]):
            paired_bases[j] = ")"

    dot_b_repr = "".join(paired_bases)
    right_brackets = dot_b_repr.count(")")
    left_brackets = dot_b_repr.count("(")

    if right_brackets == left_brackets:
        return dot_b_repr, 1
    else:
        return "".join(["." for _ in range(seq_len)]), -1


inverted_repeats = find_inverted_repeats(
    input_file="test_data/paper.fasta",
    seq_name="seq1",
    min_len=2,
    max_len=20,
    max_gap=20,
    mismatches=0,
    output_file="all_found.txt",
)

seq = "CGMSTACRCAGTCMCCGRAGMGY"
seq_len = len(seq)

# for x in dir(RNA):
#     print(x)

# PRINT FOUND INVERTED REPEATS / ERROR MESSAGE
if isinstance(inverted_repeats, str):
    print(f"IRs are type string")
    # print(inverted_repeats)
else:
    print(f"IRs are type: {type(inverted_repeats)}")
    print(
        "FORMAT: (left_strand_start, left_strand_end), (right_strand_start, right_strand_end)"
    )
    print("FOUND INVERTED REPEATS:")
    for i, ir in enumerate(inverted_repeats):
        gap_size = (ir[1][0] - ir[0][1]) - 1
        print(f'IR #{i:02d}: {str(ir).ljust(25, " ")}')

    # Exhaustively search all possible combinations of IRs
    irs_included = [ir_idx for ir_idx in range(len(inverted_repeats))]
    tot = 1000
    count = 0
    n = len(inverted_repeats)
    out_file = open("energy_out.txt", "w+")
    energy_values = []

    print(f"Seq. = {seq}")
    print(f"Seq. len = {seq_len}")
    print("=" * 40)

    # Check all possible combinations of IRs
    for i in range(1 << n):
        s = bin(i)[2:]
        s = "0" * (n - len(s)) + s

        state = list(map(int, list(s)))[::-1]
        active_irs = [ir for ir in irs_included if state[ir] == 1]
        dot_bracket_repr, valid = shape_state_to_dot_bracket(
            active_irs, inverted_repeats, seq_len
        )

        print(f"State -> {active_irs}")
        if valid == 1:
            free_energy = RNA.eval_structure_simple(seq, dot_bracket_repr, 1, out_file)
            energy_values.append(free_energy)

            print(f"  Dot bracket = {dot_bracket_repr}")
            print(f"  Free energy = {free_energy}")
        else:
            print("  Invalid state.")

        count += 1
        if count >= tot:
            break

    out_file.close()
    energy_values.sort()
    print("Free energies found")
    print(energy_values)

    # for ir_set in [ir_4, ir_14, irs_01_14]:
    #     print('IR structure evaluation'.center(70, '='))
    #     dot_bracket = irs_to_dot_bracket(ir_set, 23)
    #     print(f'Adding the pairings encoded in the following IRs: {ir_set}')
    #     print(dot_bracket)

    #     # Compute free energy of selected IR set
    #     out_file = open('energy_out.txt', 'w+')
    #     energy = RNA.eval_structure_simple(seq, dot_bracket, 0, out_file)
    #     out_file.close()
    #     print(f'Energy = {energy} kcal/mol')

    #     # Visualise selected IR set graphically
    #     bg = fgb.BulgeGraph.from_dotbracket(dot_bracket, seq)
    #     fvm.plot_rna(bg, text_kwargs={"fontweight":"black"}, lighten=0.7,
    #                      backbone_kwargs={"linewidth":3})
    #     plt.show()
