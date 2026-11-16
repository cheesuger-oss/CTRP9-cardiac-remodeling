#!/usr/bin/env python3
"""Figures 1, 8, S1 (nature-figure compliant): workflow, mouse function, GTEx.
Backend: python (exclusive). Exports SVG/PDF/TIFF/PNG via nature_style.save_pub.
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
from matplotlib.patches import FancyBboxPatch
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

# ---------------- Figure 1: study workflow ----------------
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.set_xlim(0, 13); ax.set_ylim(0, 9); ax.axis("off")

def box(x, y, w, h, text, fc="#EAF2F8", ec="#2E86C1", fs=8, tc="k", weight="normal"):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=fc, ec=ec, lw=1.2)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, color=tc, fontweight=weight)

def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.0))

box(0.3, 7.6, 3.2, 1.1, "Public transcriptome\nrepositories", fc="#D6EAF8", ec="#2E86C1", fs=9, weight="bold")
box(4.2, 7.6, 4.4, 1.1, "Eight GEO series\nplus GTEx v8", fc="#D6EAF8", ec="#2E86C1", fs=8)
box(9.3, 7.6, 3.4, 1.1, "Predefined cardiac contrasts\nPAB / MCT / TOF / HF", fc="#D6EAF8", ec="#2E86C1", fs=8)
arrow(3.5, 8.1, 4.2, 8.1); arrow(8.6, 8.1, 9.3, 8.1)

box(0.3, 6.1, 3.2, 0.9, "Data parsing & QC\nID mapping (mygene/Ensembl)", fc="#FEF9E7", ec="#B7950B", fs=8)
box(4.2, 6.1, 4.4, 0.9, "Expression matrices\nFPKM / counts / RMA-normalized", fc="#FEF9E7", ec="#B7950B", fs=8)
box(9.3, 6.1, 3.4, 0.9, "Sample provenance\nand chamber matching", fc="#FEF9E7", ec="#B7950B", fs=8)
arrow(3.5, 6.6, 4.2, 6.6); arrow(8.6, 6.6, 9.3, 6.6)
for x0 in [1.9, 6.4, 11.0]:
    arrow(x0, 7.6, x0, 7.0)

box(0.3, 4.4, 3.2, 1.1, "CTRP9 expression\n+ hypertrophy/fibrosis\ncorrelations", fc="#E8F8F5", ec="#148F77", fs=8)
box(4.2, 4.4, 4.4, 1.1, "Differential expression\n(limma, 4 datasets)\n+ WGCNA co-expression", fc="#E8F8F5", ec="#148F77", fs=8)
box(9.3, 4.4, 3.4, 1.1, "Cell-type signatures\n+ predicted miRNA-mRNA\nnetwork (TargetScan)", fc="#E8F8F5", ec="#148F77", fs=8)
for x0 in [1.9, 6.4, 11.0]:
    arrow(x0, 6.1, x0, 5.5)

box(0.3, 2.6, 3.2, 1.1, "Mouse CTRP9 KO /\nAAV overexpression\n(GSE219101)", fc="#FDEDEC", ec="#C0392B", fs=8)
box(4.2, 2.6, 4.4, 1.1, "Enrichment\n(Enrichr online: GO/KEGG/Reactome)", fc="#FDEDEC", ec="#C0392B", fs=8)
box(9.3, 2.6, 3.4, 1.1, "GTEx heart baseline\nCTRP9 expression", fc="#FDEDEC", ec="#C0392B", fs=8)
for x0 in [1.9, 6.4, 11.0]:
    arrow(x0, 4.4, x0, 3.7)

box(0.3, 0.8, 12.4, 1.2, "Context-dependent C1QTNF9 expression and ECM-associated cardiac remodeling\nPublic transcriptome reanalysis | reproducible Python and R workflows",
    fc="#F4ECF7", ec="#7D3C98", fs=10, weight="bold")
for x0 in [1.9, 6.4, 11.0]:
    arrow(x0, 2.6, x0, 2.0)

fig.tight_layout()
save_pub(fig, os.path.join(FIG, "Fig1_workflow"))
print("Fig1 saved")

# ---------------- Figure 8: mouse CTRP9 functional validation ----------------
aav = pd.read_csv("data/processed/GSE219101_AAV_CTRP9.tsv", sep="\t")
ko = pd.read_csv("data/processed/GSE219101_WT_vs_KO.tsv", sep="\t", index_col=0)
markers = ["C1qtnf9", "Col1a1", "Col1a2", "Col3a1", "Fn1", "Postn", "Tgfb1", "Nppa", "Acta1", "Myh6"]
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))

def marker_bars(ax, df, symbol_col, fc_col, fdr_col, p_col, title):
    rows = []
    for g in markers:
        sub = df[df[symbol_col] == g] if symbol_col in df.columns else df.loc[[g]] if g in df.index else None
        if sub is None or len(sub) == 0:
            continue
        r = sub.iloc[0]
        rows.append({"gene": g, "logFC": r[fc_col], "FDR": r[fdr_col], "PValue": r[p_col]})
    d = pd.DataFrame(rows).set_index("gene")
    d = d.reindex([g for g in markers if g in d.index])
    y = np.arange(len(d))
    cols = ["#D55E00" if (not np.isnan(f) and f < 0.05) else ns.C_GRAY for f in d["FDR"]]
    ax.barh(y, d["logFC"], color=cols, edgecolor="white", lw=0.4)
    ax.set_yticks(y); ax.set_yticklabels(d.index, fontsize=7)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("log2FC", fontsize=7.5)
    ax.set_title(title, fontsize=8)
    for yi, (v, f, p) in enumerate(zip(d["logFC"], d["FDR"], d["PValue"])):
        if not np.isnan(f) or not np.isnan(p):
            lab = (f"FDR={f:.1e}" if f < 0.05 else f"FDR={f:.2f}") if not np.isnan(f) else f"P={p:.1e}; FDR NA"
            # Put statistics in the empty half of each diverging bar panel,
            # rather than on top of the bar or the zero line.
            x_frac = 0.98 if v < 0 else 0.02
            ha = "right" if v < 0 else "left"
            ax.text(x_frac, yi, lab, transform=ax.get_yaxis_transform(),
                    ha=ha, va="center", fontsize=5.5,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.2})

marker_bars(axes[0], aav, "mgi_symbol", "logFC", "FDR", "PValue", "CTRP9 overexpression (AAV vs Luc; n=4 vs 3)")
marker_bars(axes[1], ko, None, "log2FoldChange", "padj", "pvalue", "CTRP9 knockout (KO vs WT; n=4 vs 3)")
fig.tight_layout(w_pad=1.6, rect=[0.02, 0, 1, 1])
fig.canvas.draw()
for ax, letter in zip(axes, "ab"):
    pos = ax.get_position()
    fig.text(pos.x0 - 0.020, pos.y1 + 0.010, letter, fontsize=8, fontweight="bold",
             ha="left", va="bottom")
save_pub(fig, os.path.join(FIG, "Fig8_mouse_function"))
print("Fig8 saved")

# ---------------- Figure S1: GTEx C1QTNF9 tissue distribution ----------------
import gzip, io
txt = gzip.decompress(open("data/aux/GTEx_v8_gene_median_tpm.gct.gz", "rb").read()).decode()
lines = txt.splitlines()
gtex = pd.read_csv(io.StringIO("\n".join(lines[2:])), sep="\t", low_memory=False)
row = gtex[gtex[gtex.columns[1]] == "C1QTNF9"].iloc[0]
tpm = row.drop(gtex.columns[:2]).astype(float).sort_values(ascending=False)
fig, ax = plt.subplots(figsize=(7.2, 3.0))
tpm20 = tpm.head(20)
colors = ["#D55E00" if "Heart" in t else ns.C_GRAY for t in tpm20.index]
ax.barh(range(len(tpm20))[::-1], tpm20.values, color=colors, edgecolor="white", lw=0.4)
ax.set_yticks(range(len(tpm20))[::-1]); ax.set_yticklabels(tpm20.index, fontsize=6)
ax.set_xlabel("Median TPM (GTEx v8)", fontsize=7.5)
ax.set_title("C1QTNF9 (CTRP9) expression across human tissues — GTEx v8", fontsize=8)
fig.tight_layout()
save_pub(fig, os.path.join(FIG, "FigS1_GTEx_C1QTNF9_tissues"))
print("FigS1 saved")
print("ALL_EXTRA_FIGURES_DONE")
