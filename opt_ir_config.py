from IUPACpal import find_inverted_repeats, config
import re
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi.graph.bulge_graph as fgb
import forgi
import itertools
import RNA

# import rnastructure
# import contrafold
from ortools.linear_solver import pywraplp
import random


def shape_state_to_dot_bracket(active_irs, all_irs, seq_len):
    """
    active_irs is list of indicies. A subset of all indices in all_irs, those
    indices that represent an IR being included in the shape e.g. [0,2,5]
    """
    paired_bases = ["." for _ in range(seq_len)]  # Initially none paired

    for ir_idx in active_irs:
        # print(f'    Adding IR#{ir_idx} -> {all_irs[ir_idx]}')
        ir = all_irs[ir_idx]
        lhs = ir[0]
        rhs = ir[1]

        for i in range(lhs[0] - 1, lhs[1]):
            paired_bases[i] = "("

        for j in range(rhs[0] - 1, rhs[1]):
            paired_bases[j] = ")"

    return "".join(paired_bases)


def eval_free_energy(dot_brk_repr, seq, out_fname):
    with open(f"{out_fname}.txt", "a") as out_file:
        out_file.write(f"Evaluating shape:\n")
        for i in range(len(dot_brk_repr)):
            out_file.write(f"{i+1:<3}")
        out_file.write("\n")
        for b in dot_brk_repr:
            out_file.write(f"{b:<3}")
        out_file.write("\n")

        free_energy = RNA.eval_structure_simple(seq, dot_brk_repr, 1, out_file)
        out_file.write(f"\n\n")


    return free_energy


def irs_incompatible(ir_a, ir_b):
    # print(f'  ir_a = {ir_a}')
    # print(f'  ir_b = {ir_b}')

    paired_base_idxs_a = [idx for idx in range(ir_a[0][0], ir_a[0][1] + 1)] + [
        idx for idx in range(ir_a[1][0], ir_a[1][1] + 1)
    ]
    paired_base_idxs_b = [idx for idx in range(ir_b[0][0], ir_b[0][1] + 1)] + [
        idx for idx in range(ir_b[1][0], ir_b[1][1] + 1)
    ]
    # print(f'  Paired base idxs in ir_a = {paired_base_idxs_a}')
    # print(f'  Paired base idxs in ir_b = {paired_base_idxs_b}')
    any_overlapping_pairings = any(
        [ir_b_bases in paired_base_idxs_a for ir_b_bases in paired_base_idxs_b]
    )

    return any_overlapping_pairings


def main():
    # print(help(RNA))
    # print(help(rnastructure))
    # print(help(contrafold))
    # for x in dir(RNA):
    #     print(x)
    # exit()

    # Generate random RNA sequence and write it to file for IUPACpal to use
    RNA.init_rand()
    seq_len = 20
    seq = RNA.random_string(seq_len, 'ACGU')
    seq_fname = "rand_rna"
    seq_energies_fname = f"{seq_fname}_energies"

    with open(f"{seq_fname}.fasta", "w") as file:
        file.write(">rand_rna\n")
        file.write(seq)

    print(f'Seq.: {seq}')
    print(f'Seq. len.: {seq_len}')
    print(f'Seq. file name: {seq_fname}')

    # Prep free energy values file
    with open(f"{seq_energies_fname}.txt", "w") as file:
        file.write(f'Energies for seq: {seq}\n\n'.center(50, "="))

    # Find inverted repeats in given sequence
    inverted_repeats = find_inverted_repeats(
        input_file="rand_rna.fasta",
        seq_name="rand_rna",
        min_len=2,
        max_len=20,
        max_gap=19,
        mismatches=0,
        output_file=f"{seq_fname}_irs.txt",
    )

    n_irs = len(inverted_repeats)

    print(f"Found IRs ({n_irs})".center(50, "="))
    print(inverted_repeats)

    # Compute free energy of each IR
    free_energies = []
    for i in range(n_irs):
        dot_b = shape_state_to_dot_bracket([i], inverted_repeats, seq_len)
        fe = eval_free_energy(dot_b, seq, seq_energies_fname)
        free_energies.append(fe)

    print("Free energies of IRs".center(50, "="))
    print(*zip([_ for _ in range(n_irs)], free_energies))

    # Instantiate MIP solver with SCIP backend
    solver = pywraplp.Solver.CreateSolver("SCIP")
    if not solver:
        print("Failed to init solver.")
        return

    # Define variables: binary indicators for each IR
    # (hand crafted for above IUPACpal settings and input for now)
    print("Variables".center(50, "="))
    variables = []
    for i in range(n_irs):
        variables.append(solver.IntVar(0, 1, f"ir_{i}"))

    # Add varibale which represents how many irs are included, n_irs - # included. Minimse this
    # included_irs_var = solver.IntVar(0, infinity, 'irs_included')

    print(f"Num variables: {solver.NumVariables()}")
    # for i, v in enumerate(variables):
    #     print(f'Var #{i}: {v.name()}')

    # Define constraints: XOR between those IRs that match the same bases
    print("Constraints".center(50, "="))

    # Find all pairs of IRs that are incompatible with each other i.e.
    # those that match the same bases
    unique_ir_pairs = itertools.combinations([_ for _ in range(n_irs)], 2)
    incompatible_ir_pairs = []
    n_unique_ir_pairs = 0
    for pair in unique_ir_pairs:
        n_unique_ir_pairs += 1
        ir_a = inverted_repeats[pair[0]]
        ir_b = inverted_repeats[pair[1]]
        incompatible = irs_incompatible(ir_a, ir_b)
        if incompatible:
            incompatible_ir_pairs.append(pair)

    # Add an XOR for all incompatiable IR pairs i.e. ir_i + ir_j <= 1
    for inc_pair in incompatible_ir_pairs:
        ir_a_idx = inc_pair[0]
        ir_b_idx = inc_pair[1]

        # constraint = solver.Constraint(0, 1, f"ir_{ir_a_idx} XOR ir_{ir_b_idx}")
        # constraint.SetCoefficient(variables[ir_a_idx], 1)
        # constraint.SetCoefficient(variables[ir_b_idx], 1)
        solver.Add(variables[ir_a_idx] + variables[ir_b_idx] <= 1)

    # Add constraint that at least one IR is included in any solution
    one_or_more_irs_constr = solver.Add(solver.Sum(variables) >= 1)

    # print(f'Num possible IR combinations = {n_unique_ir_pairs}, num valid = {n_unique_ir_pairs - len(list(incompatible_ir_pairs))}')
    print(f"Num constraints: {solver.NumConstraints()}")

    # Objective: free energy of current configuration of binary indicators
    # print("Objective".center(50, "="))
    obj_fn = solver.Objective()
    for i in range(n_irs):
        obj_fn.SetCoefficient(variables[i], free_energies[i])
    obj_fn.SetMinimization()

    print("Solution".center(50, "="))
    status = solver.Solve()
    if status == pywraplp.Solver.OPTIMAL:
        print("Objective value =", solver.Objective().Value())
        for v in variables:
            print(f'Var. "{v.name()}" opt. val. = {v.solution_value()}')
        print()
        print("Problem solved in %f milliseconds" % solver.wall_time())
        print("Problem solved in %d iterations" % solver.iterations())
        print("Problem solved in %d branch-and-bound nodes" % solver.nodes())
    else:
        print("The problem does not have an optimal solution.")

    print(f"Done.")


if __name__ == "__main__":
    main()
