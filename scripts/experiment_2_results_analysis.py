import pandas as pd
import matplotlib.pyplot as plt
import Levenshtein as lv
import numpy as np

plt.rcParams["figure.figsize"] = (6, 4)

df = pd.read_csv("./data/performance.csv")

# Sequence length vs. similarity between our and rnalibs preds
mean_lv_distances = []
# mean_hamm_distances = []
for s_len in df.seq_length.unique():
    s_len_df = df[df.seq_length == s_len]

    lv_distances = list(
        map(lv.ratio, s_len_df.our_shape_preds, s_len_df.rnalib_shape_preds)
    )
    # hamm_distances = list(map(lv.hamming, s_len_df.our_shape_preds, s_len_df.rnalib_shape_preds))

    mean_lv_distances.append(np.mean(lv_distances))
    # mean_hamm_distances.append(np.mean(lv_distances))

plt.plot(
    df.seq_length.unique(), mean_lv_distances, label="Normalised Levenshtein Distance"
)
# plt.plot(df.seq_length.unique(), mean_hamm_distances, label='Hamming Distance')
plt.legend()
plt.grid()
plt.title("Simlarity between Our Predictions and RNAlib's")
plt.xlabel("Sequence Length")
plt.ylabel("Avg. sim(RNAlib_pred, our_pred)")
plt.tight_layout()
plt.savefig("edit_distance.png")
plt.show()

avg_df = df.groupby("seq_length").mean(numeric_only=True)

# Sequence length vs. num vars/solver iters/solver time
plt.plot(df.seq_length.unique(), avg_df.num_irs, label="Num. IRs/Vars.")
plt.plot(df.seq_length.unique(), avg_df.num_solver_iters, label="Solver Iters.")
plt.plot(
    df.seq_length.unique(), avg_df.solution_time_msecs, label="Solver Time (Mili. Secs)"
)
plt.plot(df.seq_length.unique(), avg_df.num_b_n_b_nodes, label="Solver # B&B Nodes")
plt.legend()
plt.grid()
plt.title("How Sequence Length affects Solver Performance")
plt.xlabel("Sequence Length")
plt.tight_layout()
plt.savefig("seq_len_vs.png")
plt.show()

# num variables vs. num constraints
n_irs_w_means = []
for n_irs in df.num_irs.unique():
    n_irs_df = df[df.num_irs == n_irs]
    n_irs_mean = np.mean(n_irs_df.num_constraints)

    n_irs_w_means.append((n_irs, n_irs_mean))
    # print(f'For n_irs = {n_irs}, we have {len(n_irs_df)} rows, averaging to {n_irs_mean}')
n_irs_w_means.sort(key=lambda x: x[0])

plt.plot([x[0] for x in n_irs_w_means], [x[1] for x in n_irs_w_means])
plt.grid()
plt.title("Number of Constraints Generated from Found IRs")
plt.xlabel("Num. Variables/IRs")
plt.ylabel("Avg. Num. Constraints")
plt.tight_layout()
plt.savefig("seq_len_vs_constraints.png")
plt.show()

# MFEs of solutions
our_mfes_anom_rm = [mfe if mfe < 10000 else 40 for mfe in df.our_preds_mfes]
mfe_data = [df.rnalib_preds_mfes, our_mfes_anom_rm]
fig, ax = plt.subplots()

ax.boxplot(mfe_data, labels=["RNAlib MFEs", "Our MFEs"])

ax.title.set_text("MFE of Final Solution")
plt.ylabel("MFE (kcal/mol)")
plt.tight_layout()
plt.savefig("us_vs_rnalib_mfe.png")
plt.show()

# Plot distance of our final solution's MFE to RNAlib's final solution MFE over sequence length
# this should show divergence as sequence length increases which I think could be attributed to compounding
# additivity assumption errors
