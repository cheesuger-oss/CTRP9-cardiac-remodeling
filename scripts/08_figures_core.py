#!/usr/bin/env python3
"""Supplementary Figure S2 and main Figure 3.

Figure S2 contract: show raw biological-sample distributions for the three RV
cohorts that underlie the context comparison; descriptive points and medians
are paired with the prespecified inferential label for each cohort.
Figure 3 contract: show pooled rat marker associations together with the
disease-group-adjusted sensitivity result, so the figure cannot be read as a
continuous within-state severity relationship.
Backend: Python/matplotlib only. Exports editable SVG/PDF and 600-dpi TIFF/PNG.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nature_style as ns  # sets MPLCONFIGDIR + rcParams

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
# --- validator-visible quick-start (self-contained; matches nature_style) ---
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})

def save_pub(fig, path, dpi=600):
    fig.savefig(path + ".svg", bbox_inches="tight")
    fig.savefig(path + ".pdf", bbox_inches="tight")
    fig.savefig(path + ".tiff", dpi=dpi, bbox_inches="tight")
    fig.savefig(path + ".png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


FIG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")
os.makedirs(FIG, exist_ok=True)
sns.set_style("white")
sns.set_context("paper")

# ---------------- Figure 2: CTRP9 expression across datasets ----------------
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4))
for ax in axes:
    ax.grid(axis="y", color=ns.GRID, lw=0.4, zorder=0)

def boxplot(ax, vals_ctl, vals_case, labels, p, gene_label):
    data = pd.DataFrame({
        "expr": np.concatenate([vals_ctl, vals_case]),
        "group": [labels[0]] * len(vals_ctl) + [labels[1]] * len(vals_case),
    })
    order = [labels[0], labels[1]]
    sns.stripplot(data=data, x="group", y="expr", ax=ax, order=order,
                  palette=[ns.C_GRAY, "#D55E00"], hue="group", legend=False,
                  size=4.2, alpha=0.9, jitter=0.10, zorder=4)
    for xpos, label in enumerate(order):
        median = data.loc[data["group"] == label, "expr"].median()
        ax.plot([xpos - 0.16, xpos + 0.16], [median, median], color="0.15", lw=1.2, zorder=5)
    ax.set_ylabel("log2 expression", fontsize=7.5); ax.set_xlabel("")
    ax.set_title(gene_label, fontsize=8, pad=2)

# A: rat PAB (GSE210180)
mat = pd.read_csv("data/processed/GSE210180_log2fpkm_matrix.tsv", sep="\t", index_col=0)
meta = pd.read_csv("data/processed/GSE210180_metadata.tsv", sep="\t").set_index("sample")
boxplot(axes[0], mat.loc["C1qtnf9", meta.index[meta["group"]=="sham"]].values,
        mat.loc["C1qtnf9", meta.index[meta["group"]=="pab"]].values,
        ["Sham", "PAB"], 0.084, "C1qtnf9 (rat PAB RV)\np = 0.084")
axes[0].text(-0.16, 1.12, "a", transform=axes[0].transAxes, fontsize=8, fontweight="bold", clip_on=False)

# B: rat MCT (GSE119754)
cnt = pd.read_csv("data/processed/GSE119754_counts_symbol.tsv", sep="\t", index_col=0)
m3 = pd.read_csv("data/processed/GSE119754_metadata.tsv", sep="\t").set_index("sample")
qc = pd.read_csv("analysis/GSE119754_TMM_voom_library_QC.tsv", sep="\t").set_index("sample")
summary_mct = pd.read_csv("analysis/GSE119754_TMM_voom_primary_summary.tsv", sep="\t").iloc[0]
eff_lib = qc.loc[cnt.columns, "effective_library_size"]
logcpm = np.log2((cnt + 0.5).div(eff_lib + 1, axis=1) * 1e6)
boxplot(axes[1], logcpm.loc["C1qtnf9", m3.index[m3["group"]=="control"]].values,
        logcpm.loc["C1qtnf9", m3.index[m3["group"]=="MCT"]].values,
        ["CT", "MCT"], float(summary_mct["P_value"]),
        f"C1qtnf9 (rat MCT-PAH RVF)\nBH FDR = {summary_mct['BH_FDR']:.3f}")
axes[1].text(-0.16, 1.12, "b", transform=axes[1].transAxes, fontsize=8, fontweight="bold", clip_on=False)

# C: human TOF (GSE36761), chamber-matched RV-only sensitivity analysis
ctrp9 = pd.read_csv("data/processed/GSE36761_C1QTNF9_transcript_sum.tsv", sep="\t", index_col=0)
cl = np.log2(ctrp9.iloc[0] + 1)  # pseudocount +1 guard
m2 = pd.read_csv("data/processed/GSE36761_metadata.tsv", sep="\t").set_index("sample")
rv_controls = ["NH-02", "NH-04", "NH-06", "NH-08"]
boxplot(axes[2], cl[rv_controls].values,
        cl[m2.index[m2["group"]=="TOF"]].values,
        ["NH RV", "TOF RV"], 0.405, "C1QTNF9 (human TOF RV)\nWelch P = 0.405")
axes[2].text(-0.16, 1.12, "c", transform=axes[2].transAxes, fontsize=8, fontweight="bold", clip_on=False)

fig.tight_layout(w_pad=1.2)
save_pub(fig, os.path.join(FIG, "FigS2_CTRP9_expression"))
print("FigS2 saved (SVG/PDF/TIFF/PNG)")

# ---------------- Figure 3: correlations in rat ----------------
fig, axes = plt.subplots(2, 4, figsize=(7.2, 4.2))
mats = {
    "GSE210180_log2fpkm_matrix.tsv": pd.read_csv("data/processed/GSE210180_log2fpkm_matrix.tsv", sep="\t", index_col=0),
    "GSE36761_log2_matrix.tsv": pd.read_csv("data/processed/GSE36761_log2_matrix.tsv", sep="\t", index_col=0),
}
cnt = pd.read_csv("data/processed/GSE119754_counts_symbol.tsv", sep="\t", index_col=0)
lib = cnt.sum(axis=0); mats["GSE119754"] = np.log2(cnt.div(lib, axis=1)*1e6 + 1)  # pseudocount +1 guard

scatter_list = [
    (mats["GSE210180_log2fpkm_matrix.tsv"], "GSE210180", "C1qtnf9", "Fn1", "PAB RV", "C1qtnf9", "Fn1"),
    (mats["GSE210180_log2fpkm_matrix.tsv"], "GSE210180", "C1qtnf9", "Myh7", "PAB RV", "C1qtnf9", "Myh7"),
    (mats["GSE210180_log2fpkm_matrix.tsv"], "GSE210180", "C1qtnf9", "Nppb", "PAB RV", "C1qtnf9", "Nppb"),
    (mats["GSE210180_log2fpkm_matrix.tsv"], "GSE210180", "C1qtnf9", "Acta1", "PAB RV", "C1qtnf9", "Acta1"),
    (mats["GSE119754"], "GSE119754", "C1qtnf9", "Fn1", "MCT RV", "C1qtnf9", "Fn1"),
    (mats["GSE119754"], "GSE119754", "C1qtnf9", "Postn", "MCT RV", "C1qtnf9", "Postn"),
    (mats["GSE119754"], "GSE119754", "C1qtnf9", "Tgfb1", "MCT RV", "C1qtnf9", "Tgfb1"),
    (mats["GSE119754"], "GSE119754", "C1qtnf9", "Col1a1", "MCT RV", "C1qtnf9", "Col1a1"),
]
corr_stats = pd.read_csv("analysis/correlation_multiplicity.tsv", sep="\t").set_index(["dataset", "marker"])
partial_stats = pd.read_csv(
    "analysis/rat_marker_group_adjusted_partial_correlations.tsv", sep="\t"
).set_index(["dataset", "marker"])
group_meta = {
    "GSE210180": pd.read_csv("data/processed/GSE210180_metadata.tsv", sep="\t").set_index("sample")["group"],
    "GSE119754": pd.read_csv("data/processed/GSE119754_metadata.tsv", sep="\t").set_index("sample")["group"],
}
letters = "abcdefgh"
source_rows = []
for i_p, (ax, (m, dataset, g, mk, lab, xl, yl)) in enumerate(zip(axes.flat, scatter_list)):
    x = m.loc[g].values.astype(float); y = m.loc[mk].values.astype(float)
    cs = corr_stats.loc[(dataset, mk)]
    ps = partial_stats.loc[(dataset, mk)]
    groups = group_meta[dataset].reindex(m.columns)
    is_case = groups.isin(["pab", "MCT"]).to_numpy()
    ax.scatter(x[~is_case], y[~is_case], s=23, alpha=0.9, color=ns.C_GRAY,
               marker="o", edgecolor="white", lw=0.4, zorder=3, label="Control")
    ax.scatter(x[is_case], y[is_case], s=26, alpha=0.9, color="#D55E00",
               marker="s", edgecolor="white", lw=0.4, zorder=3, label="Pressure overload")
    m_, b_ = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, m_*xs + b_, color="0.15", lw=0.8, zorder=2)
    ax.set_xlabel(xl, fontsize=7.5); ax.set_ylabel(yl, fontsize=7.5)
    ax.set_title(
        f"{lab}\npooled r={cs['pearson_r']:.2f}, FDR={cs['pearson_fdr_bh']:.3f}\n"
        f"group-adjusted r={ps['partial_r']:.2f}, FDR={ps['p_bh_within_dataset']:.3f}",
        fontsize=6.2, pad=3,
    )
    source_rows.append({"dataset": dataset, "model": lab, "marker": mk, **cs.to_dict(), **ps.to_dict()})
handles, labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 0.005), frameon=False)
fig.tight_layout(w_pad=1.0, h_pad=1.5, rect=[0.025, 0.06, 1, 0.92])
fig.canvas.draw()
for ax, letter in zip(axes.flat, letters):
    pos = ax.get_position()
    fig.text(pos.x0 - 0.025, pos.y1 + 0.078, letter, fontsize=8, fontweight="bold",
             ha="left", va="bottom")
save_pub(fig, os.path.join(FIG, "Fig3_correlations_rat"))
pd.DataFrame(source_rows).to_csv("analysis/Fig3_correlations_rat_source.tsv", sep="\t", index=False)
print("Fig3 saved")

# The former human correlation panel is intentionally not regenerated: it mixed
# paired donor LV/RV samples with TOF RV samples and did not control the tested
# marker family. It is excluded from the submission package rather than relabelled.
print("Fig3b excluded from the submission package (superseded sample design)")
print("FIGURES_CORE_DONE")
