#!/usr/bin/env python3
"""Immune infiltration via ssGSEA (Barbie et al. 2009) with MCPcounter gene sets.
Applies to human datasets (GSE36761 TOF; GSE57338 once probe-mapped).
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

OUT = "results"
os.makedirs(OUT, exist_ok=True)

def parse_mcpcounter(path="data/aux/MCPcounter_genes.txt"):
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip('"') for c in df.columns]
    sets = {}
    for _, row in df.iterrows():
        cell = row["Cell population"].strip('"')
        gene = row["HUGO symbols"].strip('"')
        sets.setdefault(cell, []).append(gene)
    return {k: list(set(v)) for k, v in sets.items()}

def ssgsea_score(expr, gene_set, alpha=0.25):
    """expr: Series of log2 expression indexed by gene. Returns ssGSEA enrichment score."""
    N = len(expr)
    ranks = expr.rank(ascending=True)  # 1..N
    set_genes = [g for g in gene_set if g in expr.index]
    if len(set_genes) == 0:
        return np.nan
    rank_set = ranks[set_genes].sort_values()
    weights = np.abs(rank_set.values) ** alpha
    sum_w = weights.sum()
    # KS-like walk
    step_up = weights / sum_w
    step_dn = 1.0 / (N - len(set_genes))
    walk = np.zeros(N + 1)
    pos = 0
    up = 0.0
    for i in range(1, N + 1):
        if i in rank_set.values:
            idx = int(np.searchsorted(rank_set.values, i))
            up += step_up[idx]
            walk[i] = up - (i - pos) * step_dn
            pos = i
        else:
            walk[i] = walk[i - 1]
    # adjust for trailing
    es = np.max(walk) - np.min(walk)
    return es

def run_ssgsea(mat, sets, label):
    scores = {}
    for cell, genes in sets.items():
        vals = mat.apply(lambda col: ssgsea_score(col, genes), axis=0)
        scores[cell] = vals
    sdf = pd.DataFrame(scores).T
    sdf.to_csv(f"{OUT}/ssgsea_{label}.tsv", sep="\t")
    print(f"[{label}] ssGSEA scores: {sdf.shape}")
    return sdf

# ---------------- GSE36761 human TOF ----------------
mat = pd.read_csv("data/processed/GSE36761_log2_matrix.tsv", sep="\t", index_col=0)
meta = pd.read_csv("data/processed/GSE36761_metadata.tsv", sep="\t").set_index("sample")
sets = parse_mcpcounter()
scores = run_ssgsea(mat, sets, "GSE36761_human")

# compare NH vs TOF per cell type
nh = meta.index[meta["group"]=="NH"].tolist(); tof = meta.index[meta["group"]=="TOF"].tolist()
rows = []
for cell in scores.index:
    x = scores.loc[cell, nh].values; y = scores.loc[cell, tof].values
    t, p = stats.ttest_ind(y, x, equal_var=False)
    fc = y.mean() - x.mean()
    rows.append({"cell_type": cell, "NH_mean": x.mean(), "TOF_mean": y.mean(),
                 "diff": fc, "pvalue": p})
imm = pd.DataFrame(rows)
imm["padj"] = multipletests(imm["pvalue"].values, method="fdr_bh")[1]
imm = imm.sort_values("pvalue")
imm.to_csv(f"{OUT}/immune_diff_GSE36761.tsv", sep="\t", index=False)
print("immune differences NH vs TOF:")
print(imm.to_string(index=False))

# correlate immune scores with C1QTNF9
c = mat.loc["C1QTNF9"]
corr_rows = []
for cell in scores.index:
    r_p, p_p = stats.pearsonr(scores.loc[cell], c)
    corr_rows.append({"cell_type": cell, "pearson_r": r_p, "pearson_p": p_p})
cim = pd.DataFrame(corr_rows).sort_values("pearson_p")
cim.to_csv(f"{OUT}/immune_ctrp9_corr_GSE36761.tsv", sep="\t", index=False)
print("\nimmune vs C1QTNF9 correlation:")
print(cim.to_string(index=False))
print("IMMUNE_DONE")
