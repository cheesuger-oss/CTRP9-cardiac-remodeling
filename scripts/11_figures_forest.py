#!/usr/bin/env python3
"""Main Figure 2: cross-dataset forest plot of CTRP9.

Figure contract
---------------
Core conclusion: myocardial C1QTNF9 effects differ by disease setting; the two
rat RV estimates are negative, TOF RV is imprecise and non-significant, and
age/sex-adjusted terminal LV estimates are positive.
Evidence role: one estimate and its model-matched 95% CI per prespecified
comparison; statistics are displayed without a pooled cross-dataset effect.
Archetype: quantitative comparison grid (forest plot).
Backend: Python/matplotlib only, including preview and export.
Export: 7.2-inch editable SVG/PDF plus 600-dpi TIFF/PNG.
"""
import csv
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nature_style as ns  # sets MPLCONFIGDIR + rcParams

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
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
    # Preserve the exact 7.2-inch (182.9-mm) journal width; layout is fixed above.
    fig.savefig(path + ".svg")
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".tiff", dpi=dpi)
    fig.savefig(path + ".png", dpi=dpi)
    plt.close(fig)


FIG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")
os.makedirs(FIG, exist_ok=True)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCE_OUT = os.path.join(ROOT, "analysis", "Fig2_CTRP9_context_forest_source.tsv")


def read_matrix(name):
    return pd.read_csv(os.path.join(ROOT, "data", "processed", name), sep="\t", index_col=0)


def read_meta(name):
    return pd.read_csv(os.path.join(ROOT, "data", "processed", name), sep="\t").set_index("sample")


def welch_summary(case, control):
    """Case-minus-control mean difference with a Welch 95% confidence interval."""
    case = np.asarray(case, dtype=float)
    control = np.asarray(control, dtype=float)
    effect = case.mean() - control.mean()
    v_case = case.var(ddof=1) / len(case)
    v_control = control.var(ddof=1) / len(control)
    se = np.sqrt(v_case + v_control)
    df = (v_case + v_control) ** 2 / (
        v_case ** 2 / (len(case) - 1) + v_control ** 2 / (len(control) - 1)
    )
    critical = stats.t.ppf(0.975, df)
    p_raw = stats.ttest_ind(case, control, equal_var=False).pvalue
    return effect, effect - critical * se, effect + critical * se, p_raw


def gse57338_groups():
    """Parse disease groups from frozen GEO sample titles and align by GSM accession."""
    series = os.path.join(ROOT, "data", "raw", "GSE57338_series_matrix.txt.gz")
    titles = accessions = None
    with gzip.open(series, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            elif line.startswith("!Sample_geo_accession"):
                accessions = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            if titles is not None and accessions is not None:
                break
    if titles is None or accessions is None or len(titles) != len(accessions):
        raise RuntimeError("Could not parse aligned GSE57338 titles and accessions")
    groups = {}
    for accession, title in zip(accessions, titles):
        if "Non-failing" in title:
            groups[accession] = "NF"
        elif "Ischemic" in title:
            groups[accession] = "ICM"
        elif "Idiopathic Dilated CMP" in title:
            groups[accession] = "DCM"
        else:
            raise ValueError(f"Unrecognized GSE57338 title: {title}")
    return groups


rows = []

# Rat neonatal PAB RV: log2 FPKM; Welch test.
matrix = read_matrix("GSE210180_log2fpkm_matrix.tsv")
meta = read_meta("GSE210180_metadata.tsv")
control = matrix.loc["C1qtnf9", meta.index[meta.group == "sham"]]
case = matrix.loc["C1qtnf9", meta.index[meta.group == "pab"]]
effect, low, high, p_raw = welch_summary(case, control)
rows.append(dict(dataset="GSE210180", label="Rat neonatal PAB RV", chamber="RV",
                 comparison="PAB vs sham", n_control=len(control), n_case=len(case),
                 effect=effect, ci_low=low, ci_high=high, p_raw=p_raw,
                 p_display=p_raw, p_label=f"Welch P={p_raw:.3f}",
                 inference="two-sided Welch t-test; no multiplicity correction"))

# Rat MCT-RVF: one filtered TMM+voom-limma path supplies effect, CI, P and FDR.
mct = pd.read_csv(
    os.path.join(ROOT, "analysis", "GSE119754_TMM_voom_primary_summary.tsv"), sep="\t"
).iloc[0]
effect = float(mct["log2FC"])
low = float(mct["CI_low"])
high = float(mct["CI_high"])
p_raw = float(mct["P_value"])
fdr = float(mct["BH_FDR"])
rows.append(dict(dataset="GSE119754", label="Rat MCT-PAH RVF", chamber="RV",
                 comparison="MCT vs control", n_control=int(mct["n_control"]),
                 n_case=int(mct["n_case"]),
                 effect=effect, ci_low=low, ci_high=high, p_raw=p_raw,
                 p_display=fdr, p_label=f"BH FDR={fdr:.3f}",
                 inference=("CPM filter, TMM effective libraries and voom-limma robust "
                            "empirical Bayes; Benjamini-Hochberg genome-wide FDR"),
                 ci_method="voom-limma 95% CI"))

# Human TOF RV: corrected RV-only comparison. The eight non-failing specimens
# are paired LV/RV samples from four donors, so only the four donor RV samples
# enter the chamber-matched comparison with the 22 TOF RV samples. Reuse the
# audited sensitivity-analysis summary, including its bootstrap 95% CI.
tof = pd.read_csv(
    os.path.join(ROOT, "analysis", "GSE36761_RV_only_C1QTNF9.tsv"), sep="\t"
).iloc[0]
effect = float(tof["mean_difference"])
low = float(tof["bootstrap_CI_low"])
high = float(tof["bootstrap_CI_high"])
p_raw = float(tof["welch_p"])
rows.append(dict(dataset="GSE36761", label="Human TOF RV", chamber="RV",
                 comparison="TOF vs non-failing RV",
                 n_control=int(tof["n_NH_RV"]), n_case=int(tof["n_TOF_RV"]),
                 effect=effect, ci_low=low, ci_high=high, p_raw=p_raw,
                 p_display=p_raw, p_label=f"Welch P={p_raw:.3f}",
                 inference=("two-sided Welch t-test; one prespecified gene-level "
                            "RV-only sensitivity comparison; no multiplicity correction"),
                 ci_method="bootstrap 95% CI"))

# Human end-stage LV failure: age/sex-adjusted HC3 estimates and intervals.
adjusted = pd.read_csv(
    os.path.join(ROOT, "analysis", "GSE57338_C1QTNF9_age_sex_adjusted.tsv"), sep="\t"
)
for label, contrast in [
    ("Human ischemic HF LV", "ICM vs NF"),
    ("Human dilated-CM HF LV", "DCM vs NF"),
]:
    row = adjusted.loc[(adjusted.gene == "C1QTNF9") & (adjusted.contrast == contrast)].iloc[0]
    effect = float(row["adjusted_difference"])
    low = float(row["ci_low"])
    high = float(row["ci_high"])
    p_raw = float(row["p_raw"])
    p_bonf = float(row["p_bonferroni"])
    rows.append(dict(dataset="GSE57338", label=label, chamber="LV",
                     comparison=contrast, n_control=int(row["n_NF"]),
                     n_case=int(row["n_ICM"] if contrast.startswith("ICM") else row["n_DCM"]),
                     effect=effect, ci_low=low, ci_high=high, p_raw=p_raw,
                     p_display=p_bonf, p_label=f"Bonf. P={p_bonf:.1e}",
                     inference=("ordinary least squares adjusted for age and sex with HC3 robust "
                                "covariance; Bonferroni across three contrasts"),
                     ci_method="age/sex-adjusted HC3 95% CI"))

source = pd.DataFrame(rows)
if "ci_method" not in source.columns:
    source.insert(9, "ci_method", "Welch 95% CI")
else:
    source["ci_method"] = source["ci_method"].fillna("Welch 95% CI")
source.to_csv(SOURCE_OUT, sep="\t", index=False)

# Color and marker shape redundantly encode chamber for color and grayscale reading.
rv_color = "#0072B2"
lv_color = "#D55E00"
styles = {"RV": (rv_color, "o"), "LV": (lv_color, "s")}
fig, (ax, ax_stats) = plt.subplots(
    1, 2, figsize=(7.2, 3.05), sharey=True,
    gridspec_kw={"width_ratios": [5.2, 1.35], "wspace": 0.03},
)
ypos = np.arange(len(source))[::-1]
for y, row in zip(ypos, source.itertuples(index=False)):
    color, marker = styles[row.chamber]
    ax.errorbar(row.effect, y,
                xerr=[[row.effect - row.ci_low], [row.ci_high - row.effect]],
                fmt=marker, color=color, markerfacecolor=color,
                markersize=6.5, capsize=3, linewidth=1.4, zorder=3)

ax.axvline(0, color="0.25", linestyle="--", linewidth=0.8, zorder=1)
ax.axhline(1.5, color="0.82", linewidth=0.8, zorder=1)
ax.set_yticks(ypos)
ax.set_yticklabels([
    f"{row.label}  (n={row.n_control}+{row.n_case})" for row in source.itertuples(index=False)
])
ax.set_xlim(-4.6, 2.25)
ax.set_ylim(-0.7, 4.7)
ax.set_xlabel("Difference in log2 expression (disease minus control)")
fig.suptitle("Context-dependent C1QTNF9 expression across cardiac remodelling datasets",
             fontsize=8, y=0.98)
ax_stats.set_xlim(0, 1)
ax_stats.set_ylim(-0.7, 4.7)
ax_stats.axis("off")
ax_stats.text(0.02, 4.55, "Statistical evidence", fontsize=6.5,
              fontweight="bold", ha="left", va="bottom")
for y, row in zip(ypos, source.itertuples(index=False)):
    ax_stats.text(0.02, y, row.p_label, va="center", ha="left", fontsize=6.5)

legend = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=rv_color,
           markeredgecolor=rv_color, markersize=6, label="Right ventricle"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor=lv_color,
           markeredgecolor=lv_color, markersize=6, label="Left ventricle"),
]
ax.legend(handles=legend, loc="lower left", ncol=2)
ax.text(-4.52, 1.63, "RV pressure-overload/remodeling", fontsize=6, color="0.4", va="bottom")
ax.text(-4.52, 1.37, "End-stage LV failure", fontsize=6, color="0.4", va="top")
fig.subplots_adjust(left=0.29, right=0.98, bottom=0.20, top=0.88)
save_pub(fig, os.path.join(FIG, "Fig2_CTRP9_context_forest"))

print(source[["label", "n_control", "n_case", "effect", "ci_low", "ci_high", "p_label"]].to_string(index=False))
print(f"Source data: {SOURCE_OUT}")
print("Main Figure 2 saved")
