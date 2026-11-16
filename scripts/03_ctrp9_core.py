#!/usr/bin/env python3
"""Core analysis 1: CTRP9 expression across datasets + correlations with
hypertrophy/fibrosis markers. Produces summary tables in results/.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

OUT = "results"
os.makedirs(OUT, exist_ok=True)

def log2fc_p(group_vals, ctrl_vals):
    """mean diff of log2 values (group - ctrl) + t-test."""
    y = np.asarray(group_vals, dtype=float)
    x = np.asarray(ctrl_vals, dtype=float)
    fc = y.mean() - x.mean()
    t, p = stats.ttest_ind(y, x, equal_var=False)
    return fc, p

def corr_table(df, ctrp9_gene, markers, label):
    """Correlate CTRP9 with marker genes across samples; save table."""
    rows = []
    if ctrp9_gene not in df.index:
        print(f"  [{label}] {ctrp9_gene} NOT in matrix")
        return None
    c = df.loc[ctrp9_gene]
    for m in markers:
        if m not in df.index:
            continue
        x = df.loc[m]
        # Pearson on log2 values
        r_p, p_p = stats.pearsonr(c, x)
        rho_s, p_s = stats.spearmanr(c, x)
        rows.append({"marker": m, "pearson_r": r_p, "pearson_p": p_p,
                     "spearman_rho": rho_s, "spearman_p": p_s})
    out = pd.DataFrame(rows).sort_values("pearson_p")
    out.to_csv(f"{OUT}/corr_{label}.tsv", sep="\t", index=False)
    print(f"  [{label}] saved corr_{label}.tsv ({len(out)} markers)")
    return out

# ---- hypertrophy / fibrosis marker panels ----
HYP_MARKERS = ["Nppb", "Myh7", "Myh6", "Acta1"]
FIB_MARKERS = ["Col1a1", "Col1a2", "Col3a1", "Postn", "Fn1", "Tgfb1"]

# ============ 1. GSE210180 rat PAB ============
print("### GSE210180 (rat PAB RV, n=3+3)")
mat = pd.read_csv("data/processed/GSE210180_log2fpkm_matrix.tsv", sep="\t", index_col=0)
meta = pd.read_csv("data/processed/GSE210180_metadata.tsv", sep="\t").set_index("sample")
sham = meta.index[meta["group"] == "sham"].tolist()
pab = meta.index[meta["group"] == "pab"].tolist()
g = "C1qtnf9"
print(f"  C1qtnf9: sham={mat.loc[g, sham].values.round(2)} pab={mat.loc[g, pab].values.round(2)}")
fc, p = log2fc_p(mat.loc[g, pab], mat.loc[g, sham])
print(f"  log2FC(PAB vs sham)={fc:+.3f}  t-test p={p:.4f}")
corr_table(mat, g, HYP_MARKERS + FIB_MARKERS, "GSE210180_rat")

# ============ 2. GSE36761 human TOF ============
print("### GSE36761 (human TOF RV, n=8+22)")
mat = pd.read_csv("data/processed/GSE36761_log2_matrix.tsv", sep="\t", index_col=0)
meta = pd.read_csv("data/processed/GSE36761_metadata.tsv", sep="\t").set_index("sample")
nh = meta.index[meta["group"] == "NH"].tolist()
tof = meta.index[meta["group"] == "TOF"].tolist()
# CTRP9 from transcript sum (separate file)
ctrp9 = pd.read_csv("data/processed/GSE36761_C1QTNF9_transcript_sum.tsv", sep="\t", index_col=0)
ctrp9_log = np.log2(ctrp9.iloc[0] + 1)  # Series indexed by sample
print(f"  C1QTNF9 (transcript sum): NH={ctrp9_log[nh].values.round(2)} TOF={ctrp9_log[tof].values.round(2)}")
fc, p = log2fc_p(ctrp9_log[tof].values, ctrp9_log[nh].values)
print(f"  log2FC(TOF vs NH)={fc:+.3f}  t-test p={p:.4f}")
u, p_u = stats.mannwhitneyu(ctrp9_log[tof].values, ctrp9_log[nh].values)
print(f"  Mann-Whitney p={p_u:.4f}")
# add CTRP9 into matrix for correlation
mat2 = mat.copy()
mat2.loc["C1QTNF9"] = ctrp9_log
corr_table(mat2, "C1QTNF9", ["NPPB", "MYH7", "MYH6", "ACTA1", "COL1A1", "COL1A2", "COL3A1", "POSTN", "FN1", "TGFB1"], "GSE36761_human")

# ============ 3. GSE119754 rat MCT ============
print("### GSE119754 (rat MCT-PAH RVF, n=6+5)")
cnt = pd.read_csv("data/processed/GSE119754_counts_symbol.tsv", sep="\t", index_col=0)
meta = pd.read_csv("data/processed/GSE119754_metadata.tsv", sep="\t").set_index("sample")
ctl = meta.index[meta["group"] == "control"].tolist()
mct = meta.index[meta["group"] == "MCT"].tolist()
# TMM-lite normalization: use library-size scaled logCPM
libsize = cnt.sum(axis=0)
cpm = cnt.div(libsize, axis=1) * 1e6
logcpm = np.log2(cpm + 1)
g = "C1qtnf9"
print(f"  C1qtnf9 logCPM: ctl={logcpm.loc[g, ctl].values.round(2)} mct={logcpm.loc[g, mct].values.round(2)}")
fc, p = log2fc_p(logcpm.loc[g, mct], logcpm.loc[g, ctl])
print(f"  log2FC(MCT vs CT)={fc:+.3f}  t-test p={p:.4f}")
corr_table(logcpm, g, HYP_MARKERS + FIB_MARKERS, "GSE119754_rat")

# ============ 4. GSE224428 pig PAB ============
print("### GSE224428 (pig PAB RV, n=4+4) - C1QTNF9 not annotated")
print("  (pig annotation lacks C1QTNF9; hypertrophy/fibrosis markers validated separately)")

# ============ summary table across datasets ============
rows = [
    {"dataset": "GSE210180", "species": "rat", "model": "neonatal PAB RV",
     "n": "3+3", "ctrp9_gene": "C1qtnf9", "direction": "PAB vs sham",
     "note": "thesis-matched model"},
    {"dataset": "GSE36761", "species": "human", "model": "TOF RV",
     "n": "8+22", "ctrp9_gene": "C1QTNF9", "direction": "TOF vs NH",
     "note": "transcript-level reconstruction"},
    {"dataset": "GSE119754", "species": "rat", "model": "MCT-PAH RVF",
     "n": "6+5", "ctrp9_gene": "C1qtnf9", "direction": "MCT vs CT",
     "note": ""},
    {"dataset": "GSE224428", "species": "pig", "model": "PAB RVF",
     "n": "4+4", "ctrp9_gene": "C1QTNF9", "direction": "PAB vs CT",
     "note": "not annotated in pig genome build"},
]
pd.DataFrame(rows).to_csv(f"{OUT}/dataset_summary.tsv", sep="\t", index=False)
print("\nsaved dataset_summary.tsv")
print("DONE")
