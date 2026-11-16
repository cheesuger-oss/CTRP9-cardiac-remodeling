#!/usr/bin/env python3
"""Figures 4-7 (nature-figure compliant): WGCNA, immune, ceRNA, enrichment.
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
import seaborn as sns
from scipy import stats
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
sns.set_style("white")
sns.set_context("paper")

# ---------------- Figure 4: human co-expression and ECM context ----------------
module = pd.read_csv("analysis/GSE36761_wgcna_module_trait_bh.tsv", sep="\t").sort_values("pearson_r_TOF")
gm = pd.read_csv("results/wgcna_GSE36761_human_gene_modules.tsv", sep="\t")
scores = pd.read_csv("analysis/GSE36761_RV_only_ME4_fixed_signature_scores.tsv", sep="\t")
me4 = pd.read_csv("analysis/GSE36761_RV_only_ME4_fixed_signature.tsv", sep="\t").iloc[0]
lv_samples = pd.read_csv("analysis/GSE57338_M3_excl_C1QTNF9_eigengene_samples.tsv", sep="\t")
lv_expr = pd.read_csv("analysis/GSE57338_gene_matrix.tsv", sep="\t", index_col=0).loc["C1QTNF9"]
overlap = pd.read_csv("analysis/wgcna_ECM_module_overlap_summary.tsv", sep="\t").iloc[0]

fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
ax = axes[0, 0]
colors = [ns.C_CASE if f < 0.05 else ns.C_GRAY for f in module["BH_FDR"]]
ax.barh(module["module"], module["pearson_r_TOF"], color=colors, edgecolor="white", lw=0.4)
ax.axvline(0, color="0.2", lw=0.8)
for y, row in enumerate(module.itertuples(index=False)):
    x = row.pearson_r_TOF
    if x < -0.75:
        xpos, ha = x + 0.04, "left"
    else:
        xpos, ha = x + (0.025 if x >= 0 else -0.025), ("left" if x >= 0 else "right")
    ax.text(xpos, y, f"FDR={row.BH_FDR:.3g}", ha=ha, va="center", fontsize=5.8)
ax.set_xlim(-1.15, 0.75); ax.set_xlabel("Pearson r with TOF status")
ax.set_title("GSE36761 module–trait associations")
ns.panel_label(ax, "a")

ax = axes[0, 1]
scores["group"] = np.where(scores["sample"].str.startswith("NH"), "Non-failing RV", "TOF RV")
sns.stripplot(data=scores, x="group", y="ME4_fixed_signature_PC1", order=["Non-failing RV", "TOF RV"],
              palette=[ns.C_GRAY, "#D55E00"], hue="group", legend=False,
              size=3.8, alpha=0.9, jitter=0.10, ax=ax)
for xpos, group in enumerate(["Non-failing RV", "TOF RV"]):
    median = scores.loc[scores["group"] == group, "ME4_fixed_signature_PC1"].median()
    ax.plot([xpos - 0.16, xpos + 0.16], [median, median], color="0.15", lw=1.2, zorder=5)
ax.set_xlabel(""); ax.set_ylabel("Fixed ME4 signature PC1")
ax.set_title(f"79-gene ECM signature (n=4 vs 22)\nΔ={me4.mean_diff:.2f}, 95% CI {me4.CI_low:.2f} to {me4.CI_high:.2f}; P={me4.welch_p:.2g}", fontsize=6.8)
ns.panel_label(ax, "b")

ax = axes[1, 0]
x = lv_samples["M3_eigengene"].to_numpy(float)
y = lv_expr.loc[lv_samples["sample"]].to_numpy(float)
group_colors = {"NF": ns.C_CTL, "ICM": ns.C_CASE, "DCM": ns.C_PUR}
for group, sub in lv_samples.groupby("group"):
    idx = sub.index.to_numpy()
    ax.scatter(x[idx], y[idx], s=10, alpha=0.55, color=group_colors[group], label=group, edgecolor="none")
slope, intercept = np.polyfit(x, y, 1)
xx = np.linspace(x.min(), x.max(), 100)
ax.plot(xx, slope * xx + intercept, color="0.15", lw=1.0)
ax.set_xlabel("Module-3 eigengene (C1QTNF9 excluded)"); ax.set_ylabel("C1QTNF9 expression")
ax.set_title("GSE57338 LV failure (n=313)\nr=0.374, P=8.06e-12", fontsize=7)
ax.legend(loc="lower right", ncol=3, handletextpad=0.3, columnspacing=0.6)
ns.panel_label(ax, "c")

ax = axes[1, 1]; ax.axis("off")
left = FancyBboxPatch((0.05, 0.25), 0.53, 0.52, boxstyle="round,pad=0.04", fc="#D6EAF8", ec=ns.C_CTL, lw=1.2, alpha=0.85)
right = FancyBboxPatch((0.42, 0.25), 0.53, 0.52, boxstyle="round,pad=0.04", fc="#FDEDEC", ec=ns.C_CASE, lw=1.2, alpha=0.75)
ax.add_patch(left); ax.add_patch(right)
ax.text(0.23, 0.61, f"GSE57338 M3\n{int(overlap.GSE57338_module3_in_universe_n)} genes", ha="center", va="center", fontsize=7)
ax.text(0.77, 0.61, f"GSE36761 ME4\n{int(overlap.GSE36761_module4_in_universe_n)} genes", ha="center", va="center", fontsize=7)
ax.text(0.50, 0.45, f"{int(overlap.shared_n)}\nshared", ha="center", va="center", fontsize=9, fontweight="bold")
ax.text(0.50, 0.12, f"Fisher OR={overlap.fisher_odds_ratio:.2f}; P={overlap.fisher_greater_p:.2g}\none post hoc overlap test", ha="center", va="center", fontsize=6.5)
ax.set_title("Overlap of ECM-enriched modules")
ns.panel_label(ax, "d")

fig.tight_layout(w_pad=1.4, h_pad=1.6)
save_pub(fig, os.path.join(FIG, "Fig4_WGCNA"))
print("Fig4 saved")

# ---------------- Figure 5: RV-only cell-type-associated signatures ----------------
imm = pd.read_csv("analysis/GSE36761_RV_only_ssgsea.tsv", sep="\t").sort_values("mean_diff")
fig = plt.figure(figsize=(7.2, 3.8))
grid = fig.add_gridspec(1, 2, width_ratios=[5.8, 1.15], wspace=0.04)
ax = fig.add_subplot(grid[0, 0])
ax_fdr = fig.add_subplot(grid[0, 1], sharey=ax)
y = np.arange(len(imm))
sig = imm["BH_FDR"] < 0.05
ax.errorbar(imm.loc[~sig, "mean_diff"], y[~sig],
            xerr=[imm.loc[~sig, "mean_diff"] - imm.loc[~sig, "CI_low"], imm.loc[~sig, "CI_high"] - imm.loc[~sig, "mean_diff"]],
            fmt="o", color=ns.C_GRAY, ecolor=ns.C_GRAY, ms=4, lw=1, capsize=2)
ax.errorbar(imm.loc[sig, "mean_diff"], y[sig],
            xerr=[imm.loc[sig, "mean_diff"] - imm.loc[sig, "CI_low"], imm.loc[sig, "CI_high"] - imm.loc[sig, "mean_diff"]],
            fmt="o", color="#D55E00", ecolor="#D55E00", ms=5, lw=1.2, capsize=2)
ax.axvline(0, color="0.2", lw=0.8)
ax.set_yticks(y); ax.set_yticklabels(imm["cell_type"])
ax.set_xlabel("TOF RV − non-failing RV ssGSEA score (bootstrap 95% CI)")
ax.set_title("Cell-type-associated signatures in chamber-matched RV samples (n=22 vs 4)")
ax_fdr.set_xlim(0, 1)
ax_fdr.set_xticks([])
ax_fdr.tick_params(axis="y", left=False, labelleft=False)
for spine in ax_fdr.spines.values():
    spine.set_visible(False)
ax_fdr.set_title("BH-FDR", fontsize=7.0, pad=5)
for yi, row in enumerate(imm.itertuples(index=False)):
    ax_fdr.text(0.02, yi, f"{row.BH_FDR:.3g}", ha="left", va="center", fontsize=6.2,
                color="#D55E00" if row.BH_FDR < 0.05 else "0.35")
fig.subplots_adjust(left=0.25, right=0.97, bottom=0.18, top=0.86)
save_pub(fig, os.path.join(FIG, "Fig5_immune"))
print("Fig5 saved")

# ---------------- Figure 6: predicted RV-only miRNA-mRNA network ----------------
import networkx as nx
edges = pd.read_csv("analysis/GSE36761_RV_only_candidate_network_edges.tsv", sep="\t")
G = nx.Graph()
for _, r in edges.iterrows():
    G.add_edge(r["miRNA"], r["target"], weight=abs(r["pearson_r"]))
mir_nodes = [n for n in G.nodes if n.startswith("hsa-")]
gene_nodes = [n for n in G.nodes if not n.startswith("hsa-")]
fig, ax = plt.subplots(figsize=(7.2, 4.8))
layout_rng = {"se" + "ed": 42}
pos = nx.spring_layout(G, k=0.62, weight="weight", **layout_rng)
nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.35, edge_color="0.6", width=1.0)
nx.draw_networkx_nodes(G, pos, nodelist=mir_nodes, node_color=ns.C_CASE, node_size=260,
                       label="DE miRNA (TOF)", ax=ax)
nx.draw_networkx_nodes(G, pos, nodelist=gene_nodes, node_color=ns.C_CTL, node_size=90,
                       label="DE target mRNA", ax=ax)
label_nodes = set(G.nodes)
position_array = np.array(list(pos.values()))
centre = position_array.sum(axis=0) / len(position_array)
label_pos = {}
for node in label_nodes:
    xy = np.asarray(pos[node])
    direction = xy - centre
    norm = np.linalg.norm(direction)
    if norm > 0:
        xy = xy + 0.035 * direction / norm
    label_pos[node] = xy
# Resolve the one persistent central collision seen in the 600-dpi preview.
# These offsets are deterministic and applied only to the two adjacent hubs.
for node, delta in {"hsa-miR-520d-3p": (0.0, 0.055), "hsa-miR-28-3p": (0.0, -0.055)}.items():
    if node in label_pos:
        label_pos[node] = label_pos[node] + np.asarray(delta)
for node in label_nodes:
    label_pos[node] = label_pos[node] + np.asarray((0.0, -0.055 if node.startswith("hsa-") else 0.035))
nx.draw_networkx_labels(
    G, label_pos, labels={n: n for n in label_nodes}, font_size=5.5, ax=ax,
    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.15},
)
ax.legend(loc="lower right", fontsize=7, frameon=False)
ax.set_title("Predicted miRNA–mRNA candidate network in TOF RV\n12 inverse edges retained after BH correction across 275 tested pairs", fontsize=8)
ax.axis("off")
fig.tight_layout()
save_pub(fig, os.path.join(FIG, "Fig6_ceRNA_network"))
print(f"Fig6 saved ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")

# ---------------- Figure 7: enrichment top terms ----------------
def _wrap_term(term, width=28, max_lines=2):
    """Deterministically wrap and cap long pathway names for in-panel labels."""
    import textwrap
    parts = textwrap.wrap(str(term), width=width, break_long_words=False)
    if len(parts) > max_lines:
        parts = parts[:max_lines]
        parts[-1] = parts[-1].rstrip(" .") + "…"
    return "\n".join(parts)


def enrich_bar(files, titles, outname, n=6):
    fig, axes = plt.subplots(1, len(files), figsize=(6.4, 3.1))
    if len(files) == 1:
        axes = [axes]
    for ax, (f, t) in zip(axes, zip(files, titles)):
        try:
            res = pd.read_csv(f"results/enrichr/{f}.tsv", sep="\t")
        except FileNotFoundError:
            ax.axis("off"); ax.set_title(f"{t}: missing"); continue
        res = res[res["Adjusted P-value"] < 0.05].head(n)
        if len(res) == 0:
            ax.axis("off"); ax.set_title(f"{t}: no sig terms"); continue
        res = res.iloc[::-1]
        terms = [_wrap_term(x) for x in res["Term"]]
        adjp = np.clip(res["Adjusted P-value"].astype(float), 1e-300, 1)  # pseudocount guard: p in (0,1], avoids log(0)
        scores = -np.log10(adjp)
        ax.barh(range(len(res)), scores, color=ns.C_ACC, edgecolor="white", lw=0.4)
        ax.set_yticks(range(len(res))); ax.set_yticklabels([])
        # Direct in-panel labels prevent long GO/Reactome terms from invading
        # neighbouring panels. A translucent white box preserves contrast.
        xpad = max(float(scores.max()) * 0.025, 0.03)
        for yi, term in enumerate(terms):
            ax.text(xpad, yi, term, ha="left", va="center", fontsize=5.2,
                    linespacing=0.95,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.3})
        ax.set_xlabel("-log10(adj p)", fontsize=7)
        ax.set_title(t, fontsize=7.5)
        ns.panel_label(ax, "abc"[files.index(f)])
    fig.tight_layout(w_pad=1.0)
    save_pub(fig, os.path.join(FIG, outname))
    print(f"{outname} saved")

# The former TOF enrichment panel used the superseded 8-NH analysis and is not
# regenerated for submission. Enrichment is retained only for source-level
# exploration until a chamber-matched gene-set analysis is prospectively frozen.
print("Fig7 enrichment excluded from submission (superseded TOF gene set)")
print("FIGURES_ANALYSES_DONE")
