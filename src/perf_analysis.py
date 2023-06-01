import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("performance.csv")

avg_df = df.groupby("seq_length").mean()


# Sequence length vs. num vars/solver iters/solver time
plt.plot(df.seq_length.unique(), avg_df.num_irs, label='IRs/Vars.')
plt.plot(df.seq_length.unique(), avg_df.num_solver_iters, label='Solver Iters.')
plt.plot(df.seq_length.unique(), avg_df.solution_time_msecs, label='Solver Time (Mili. Secs)')
plt.plot(df.seq_length.unique(), avg_df.num_b_n_b_nodes, label='Solver # B&B Nodes')
plt.legend()
plt.grid()
plt.xlabel('Sequence Length')
plt.savefig('seq_len_vs.png')
plt.show()

# Sequence length vs. num constraints
plt.plot(df.seq_length.unique(), avg_df.num_constraints, label='Constraints')
plt.legend()
plt.grid()
plt.xlabel('Sequence Length')
plt.ylabel('# Constraints')
plt.savefig('seq_len_vs_constraints.png')
plt.show()

# MFEs of solutions
our_mfes_anom_rm = [mfe if mfe < 10000 else 40 for mfe in df.our_preds_mfes]
mfe_data = [df.rnalib_preds_mfes, our_mfes_anom_rm]
fig, ax = plt.subplots()

ax.boxplot(mfe_data, labels=['RNAlib MFEs', 'Our MFEs'])

ax.title.set_text('MFE of Final Solution')
plt.savefig('us_vs_rnalib_mfe.png')
plt.show()
