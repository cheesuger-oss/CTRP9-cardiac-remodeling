#!/usr/bin/env python3
"""Cross-cohort CTRP9 effect sizes and rat correlation sensitivity.

The analysis plan was locked in
``analysis/cross_cohort_effect_sensitivity_plan.md`` before these outputs were
inspected. All inputs are read-only. New tables are written to ``analysis/``
and the supplementary figure to ``figures/``.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image, ImageOps
from scipy import optimize, stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
FIGURES = ROOT / "figures"


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def metadata(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    unnamed = [column for column in table.columns if column.startswith("Unnamed:")]
    if unnamed:
        table = table.drop(columns=unnamed)
    return table.set_index("sample")


def gse57338_groups() -> dict[str, str]:
    titles = pd.read_csv(ROOT / "data/processed/GSE57338_titles.tsv", sep="\t")["sample"].tolist()
    accessions = pd.read_csv(
        ANALYSIS / "GSE57338_gene_matrix.tsv", sep="\t", nrows=1, index_col=0
    ).columns.tolist()
    if len(titles) != len(accessions):
        raise RuntimeError("GSE57338 titles are not aligned to expression columns")
    result: dict[str, str] = {}
    for accession, title in zip(accessions, titles):
        if "Non-failing" in title:
            result[accession] = "NF"
        elif "Ischemic" in title:
            result[accession] = "ICM"
        elif "Idiopathic Dilated CMP" in title:
            result[accession] = "DCM"
        else:
            raise ValueError(f"Unrecognised GSE57338 title: {title}")
    return result


def hedges_g(case: np.ndarray, control: np.ndarray) -> dict[str, float]:
    case = np.asarray(case, dtype=float)
    control = np.asarray(control, dtype=float)
    n1, n0 = len(case), len(control)
    df = n1 + n0 - 2
    pooled_var = ((n1 - 1) * np.var(case, ddof=1) + (n0 - 1) * np.var(control, ddof=1)) / df
    if pooled_var <= 0:
        raise ValueError("Pooled variance must be positive")
    d = (np.mean(case) - np.mean(control)) / np.sqrt(pooled_var)
    correction = 1 - 3 / (4 * df - 1)
    g = correction * d
    variance = (n1 + n0) / (n1 * n0) + g * g / (2 * df)
    se = np.sqrt(variance)
    return {
        "n_control": n0,
        "n_case": n1,
        "mean_control": float(np.mean(control)),
        "mean_case": float(np.mean(case)),
        "sd_control": float(np.std(control, ddof=1)),
        "sd_case": float(np.std(case, ddof=1)),
        "hedges_g": float(g),
        "sampling_variance": float(variance),
        "standard_error": float(se),
        "ci_low": float(g - stats.norm.ppf(0.975) * se),
        "ci_high": float(g + stats.norm.ppf(0.975) * se),
        "welch_p": float(stats.ttest_ind(case, control, equal_var=False).pvalue),
    }


def cohort_effects() -> pd.DataFrame:
    rows = []

    pab_matrix = read_matrix(ROOT / "data/processed/GSE210180_log2fpkm_matrix.tsv")
    pab_meta = metadata(ROOT / "data/processed/GSE210180_metadata.tsv")
    pab = pab_matrix.loc["C1qtnf9", pab_meta.index]
    result = hedges_g(pab[pab_meta.group.eq("pab")], pab[pab_meta.group.eq("sham")])
    rows.append({
        "dataset": "GSE210180", "context": "Rat neonatal PAB RV",
        "chamber": "RV", "comparison": "PAB vs sham",
        "expression_scale": "log2(FPKM+1)", **result,
    })

    mct_counts = read_matrix(ROOT / "data/processed/GSE119754_counts_symbol.tsv")
    mct_meta = metadata(ROOT / "data/processed/GSE119754_metadata.tsv")
    mct_logcpm = np.log2(mct_counts.div(mct_counts.sum(axis=0), axis=1) * 1e6 + 1)
    mct = mct_logcpm.loc["C1qtnf9", mct_meta.index]
    result = hedges_g(mct[mct_meta.group.eq("MCT")], mct[mct_meta.group.eq("control")])
    rows.append({
        "dataset": "GSE119754", "context": "Rat MCT RV failure",
        "chamber": "RV", "comparison": "MCT vs control",
        "expression_scale": "log2(CPM+1)", **result,
    })

    tof_matrix = read_matrix(ROOT / "data/processed/GSE36761_log2_matrix.tsv")
    provenance = pd.read_csv(
        ANALYSIS / "GSE36761_sample_provenance_corrected.tsv", sep="\t"
    ).set_index("sample")
    rv_samples = provenance.index[provenance.chamber.eq("RV")]
    if not set(rv_samples).issubset(tof_matrix.columns):
        raise RuntimeError("Corrected GSE36761 RV samples are not all in the expression matrix")
    tof = tof_matrix.loc["C1QTNF9", rv_samples]
    rv_groups = provenance.loc[rv_samples, "group"]
    result = hedges_g(tof[rv_groups.eq("TOF")], tof[rv_groups.eq("NH")])
    rows.append({
        "dataset": "GSE36761", "context": "Human TOF RV",
        "chamber": "RV", "comparison": "TOF vs non-failing RV",
        "expression_scale": "provider-normalized log2 expression", **result,
    })

    hf_matrix = read_matrix(ANALYSIS / "GSE57338_gene_matrix.tsv")
    hf_groups = pd.Series(gse57338_groups())
    if list(hf_matrix.columns) != list(hf_groups.index):
        raise RuntimeError("GSE57338 group mapping is not aligned to expression columns")
    hf = hf_matrix.loc["C1QTNF9"]
    result = hedges_g(hf[hf_groups.ne("NF")], hf[hf_groups.eq("NF")])
    rows.append({
        "dataset": "GSE57338", "context": "Human end-stage HF LV",
        "chamber": "LV", "comparison": "ICM/DCM vs non-failing",
        "expression_scale": "RMA log2 expression", **result,
    })

    table = pd.DataFrame(rows)
    table["welch_fdr_bh"] = multipletests(table["welch_p"], method="fdr_bh")[1]
    table["multiplicity_family"] = "four independent cohort-level CTRP9 contrasts"
    table["ci_method"] = "normal approximation from Hedges g sampling variance"
    return table


def random_effects(effect: np.ndarray, variance: np.ndarray, label: str) -> dict[str, float | str]:
    effect = np.asarray(effect, dtype=float)
    variance = np.asarray(variance, dtype=float)
    k = len(effect)
    if k < 2:
        raise ValueError("Random-effects synthesis needs at least two cohorts")

    def reml_objective(tau2: float) -> float:
        weights = 1 / (variance + tau2)
        mean = np.sum(weights * effect) / np.sum(weights)
        return 0.5 * (
            np.sum(np.log(variance + tau2))
            + np.log(np.sum(weights))
            + np.sum(weights * (effect - mean) ** 2)
        )

    upper = max(10.0, float(np.var(effect, ddof=1) * 20 + np.max(variance)))
    fit = optimize.minimize_scalar(reml_objective, bounds=(0, upper), method="bounded")
    tau2 = max(0.0, float(fit.x))
    weights = 1 / (variance + tau2)
    pooled = float(np.sum(weights * effect) / np.sum(weights))
    residual_q = float(np.sum(weights * (effect - pooled) ** 2))
    hk_scale = max(1.0, residual_q / (k - 1))
    se_hk = float(np.sqrt(hk_scale / np.sum(weights)))
    critical = float(stats.t.ppf(0.975, k - 1))
    ci_low, ci_high = pooled - critical * se_hk, pooled + critical * se_hk
    p_hk = float(2 * stats.t.sf(abs(pooled / se_hk), k - 1))

    fixed_weights = 1 / variance
    fixed_mean = float(np.sum(fixed_weights * effect) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * (effect - fixed_mean) ** 2))
    q_df = k - 1
    q_p = float(stats.chi2.sf(q, q_df))
    i2 = float(max(0.0, (q - q_df) / q) * 100) if q > 0 else 0.0
    return {
        "scope": label,
        "k": k,
        "pooled_hedges_g": pooled,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "p_modified_knapp_hartung": p_hk,
        "tau2_reml": tau2,
        "cochran_q": q,
        "q_df": q_df,
        "q_p": q_p,
        "i2_percent": i2,
        "interval_method": "REML with modified Knapp-Hartung 95% CI",
        "interpretation": "exploratory; contexts are biologically and technically heterogeneous",
    }


def meta_summaries(effects: pd.DataFrame) -> pd.DataFrame:
    rows = [
        random_effects(
            effects.hedges_g.to_numpy(), effects.sampling_variance.to_numpy(), "all four contexts"
        ),
        random_effects(
            effects.loc[effects.chamber.eq("RV"), "hedges_g"].to_numpy(),
            effects.loc[effects.chamber.eq("RV"), "sampling_variance"].to_numpy(),
            "RV cohorts only",
        ),
    ]
    return pd.DataFrame(rows)


def rat_matrices() -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    pab_matrix = read_matrix(ROOT / "data/processed/GSE210180_log2fpkm_matrix.tsv")
    pab_meta = metadata(ROOT / "data/processed/GSE210180_metadata.tsv")
    pab_matrix = pab_matrix.loc[:, pab_meta.index]

    mct_counts = read_matrix(ROOT / "data/processed/GSE119754_counts_symbol.tsv")
    mct_meta = metadata(ROOT / "data/processed/GSE119754_metadata.tsv")
    mct_matrix = np.log2(mct_counts.div(mct_counts.sum(axis=0), axis=1) * 1e6 + 1)
    mct_matrix = mct_matrix.loc[:, mct_meta.index]
    return {
        "GSE210180": (pab_matrix, pab_meta),
        "GSE119754": (mct_matrix, mct_meta),
    }


def correlation_sensitivity() -> tuple[pd.DataFrame, pd.DataFrame]:
    full = pd.read_csv(ANALYSIS / "correlation_multiplicity.tsv", sep="\t")
    partial = pd.read_csv(
        ANALYSIS / "rat_marker_group_adjusted_partial_correlations.tsv", sep="\t"
    )
    matrices = rat_matrices()
    loo_rows = []
    summary_rows = []
    for dataset, dataset_full in full.groupby("dataset", sort=False):
        matrix, meta = matrices[dataset]
        x = matrix.loc["C1qtnf9"].to_numpy(float)
        samples = matrix.columns.to_list()
        for row in dataset_full.itertuples(index=False):
            y = matrix.loc[row.marker].to_numpy(float)
            pearson_values, spearman_values = [], []
            for index, omitted in enumerate(samples):
                keep = np.arange(len(samples)) != index
                pearson_r, pearson_p = stats.pearsonr(x[keep], y[keep])
                spearman_r, spearman_p = stats.spearmanr(x[keep], y[keep])
                pearson_values.append(float(pearson_r))
                spearman_values.append(float(spearman_r))
                loo_rows.append({
                    "dataset": dataset,
                    "marker": row.marker,
                    "omitted_sample": omitted,
                    "n_after_deletion": int(keep.sum()),
                    "pearson_r": float(pearson_r),
                    "pearson_p_descriptive": float(pearson_p),
                    "spearman_rho": float(spearman_r),
                    "spearman_p_descriptive": float(spearman_p),
                    "inference_boundary": "descriptive leave-one-out diagnostic; P values not used for claims",
                })
            pearson_values = np.asarray(pearson_values)
            spearman_values = np.asarray(spearman_values)
            partial_row = partial.loc[
                partial.dataset.eq(dataset) & partial.marker.eq(row.marker)
            ].iloc[0]
            summary_rows.append({
                "dataset": dataset,
                "marker": row.marker,
                "n_full": int(row.n),
                "pearson_r_full": float(row.pearson_r),
                "pearson_fdr_bh": float(row.pearson_fdr_bh),
                "spearman_rho_full": float(row.spearman_rho),
                "spearman_fdr_bh": float(row.spearman_fdr_bh),
                "loo_pearson_min": float(pearson_values.min()),
                "loo_pearson_median": float(np.median(pearson_values)),
                "loo_pearson_max": float(pearson_values.max()),
                "loo_pearson_negative_fraction": float(np.mean(pearson_values < 0)),
                "loo_pearson_max_abs_shift": float(np.max(np.abs(pearson_values - row.pearson_r))),
                "loo_spearman_min": float(spearman_values.min()),
                "loo_spearman_median": float(np.median(spearman_values)),
                "loo_spearman_max": float(spearman_values.max()),
                "loo_spearman_negative_fraction": float(np.mean(spearman_values < 0)),
                "loo_spearman_max_abs_shift": float(np.max(np.abs(spearman_values - row.spearman_rho))),
                "pearson_direction_stable_all_deletions": bool(np.all(pearson_values < 0)),
                "spearman_direction_stable_all_deletions": bool(np.all(spearman_values < 0)),
                "group_adjusted_partial_r": float(partial_row.partial_r),
                "group_adjusted_partial_ci_low": float(partial_row.ci_low),
                "group_adjusted_partial_ci_high": float(partial_row.ci_high),
                "group_adjusted_partial_fdr_bh": float(partial_row.p_bh_within_dataset),
                "analysis_status": "post hoc sensitivity; complete prespecified marker family",
            })
    return pd.DataFrame(loo_rows), pd.DataFrame(summary_rows)


def plot_figure(effects: pd.DataFrame, meta: pd.DataFrame, sensitivity: pd.DataFrame) -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 6.5,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    })
    width_in = 183 / 25.4
    height_in = 145 / 25.4
    fig = plt.figure(figsize=(width_in, height_in), constrained_layout=True)
    grid = fig.add_gridspec(
        3, 2, height_ratios=[0.92, 0.30, 1.65], hspace=0.08, wspace=0.18
    )
    ax_a = fig.add_subplot(grid[0, :])
    ax_legend = fig.add_subplot(grid[1, :])
    ax_b = fig.add_subplot(grid[2, 0])
    ax_c = fig.add_subplot(grid[2, 1])
    ax_legend.axis("off")

    colors = {"RV": "#0072B2", "LV": "#D55E00"}
    markers = {"RV": "o", "LV": "s"}
    y = np.arange(len(effects))[::-1] + 2
    for ypos, row in zip(y, effects.itertuples(index=False)):
        ax_a.errorbar(
            row.hedges_g, ypos,
            xerr=[[row.hedges_g - row.ci_low], [row.ci_high - row.hedges_g]],
            fmt=markers[row.chamber], color=colors[row.chamber], markerfacecolor=colors[row.chamber],
            markersize=5.5, capsize=2.5, linewidth=1.1, zorder=3,
        )
    meta_y = [1, 0]
    for ypos, row in zip(meta_y, meta.itertuples(index=False)):
        x = row.pooled_hedges_g
        left, right = row.ci_low, row.ci_high
        diamond = np.array([[left, ypos], [x, ypos + 0.20], [right, ypos], [x, ypos - 0.20]])
        ax_a.fill(diamond[:, 0], diamond[:, 1], color="0.35", alpha=0.75, linewidth=0)
    labels = [
        f"{row.context}  (n={row.n_control}+{row.n_case})"
        for row in effects.itertuples(index=False)
    ] + [
        f"Exploratory summary: {row.scope}  (k={row.k})"
        for row in meta.itertuples(index=False)
    ]
    ax_a.set_yticks(list(y) + meta_y)
    ax_a.set_yticklabels(labels)
    ax_a.axvline(0, color="0.25", linestyle="--", linewidth=0.8)
    ax_a.axhline(1.5, color="0.82", linewidth=0.7)
    ax_a.set_xlabel("Hedges' g (disease minus reference), 95% CI")
    ax_a.set_title("Cross-cohort standardized C1QTNF9/C1qtnf9 effects", loc="left", fontweight="bold")
    ax_a.spines[["top", "right"]].set_visible(False)
    ax_a.margins(x=0.04, y=0.12)

    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["RV"],
               markeredgecolor=colors["RV"], label="RV cohort"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=colors["LV"],
               markeredgecolor=colors["LV"], label="LV cohort"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="0.35",
               markeredgecolor="0.35", label="Exploratory REML summary"),
    ]
    legend_top = ax_legend.legend(
        handles=legend, loc="upper center", bbox_to_anchor=(0.5, 1.05),
        frameon=False, ncol=3, handletextpad=0.4, columnspacing=1.2,
    )
    ax_legend.add_artist(legend_top)

    def sensitivity_panel(ax: plt.Axes, dataset: str, title: str) -> None:
        subset = sensitivity.loc[sensitivity.dataset.eq(dataset)].copy()
        subset = subset.sort_values("pearson_r_full")
        ypos = np.arange(len(subset))
        for pos, row in zip(ypos, subset.itertuples(index=False)):
            ax.hlines(pos, row.loo_pearson_min, row.loo_pearson_max, color="0.65", linewidth=2.0)
        ax.scatter(subset.pearson_r_full, ypos, s=22, color="#0072B2", label="Pearson", zorder=3)
        ax.scatter(subset.spearman_rho_full, ypos, s=20, color="#E69F00", marker="s",
                   label="Spearman", zorder=3)
        ax.scatter(subset.group_adjusted_partial_r, ypos, s=25, facecolors="white",
                   edgecolors="#6A3D9A", marker="^", linewidths=1.0,
                   label="Group-adjusted partial r", zorder=3)
        ax.axvline(0, color="0.25", linestyle="--", linewidth=0.7)
        ax.set_xlim(-1.05, 1.05)
        ax.set_yticks(ypos)
        ax.set_yticklabels(subset.marker)
        ax.set_xlabel("Correlation coefficient")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="0.92", linewidth=0.5)

    sensitivity_panel(ax_b, "GSE210180", "PAB RV (10 markers)")
    sensitivity_panel(ax_c, "GSE119754", "MCT RV failure (7 markers)")
    handles = [
        Line2D([0], [0], color="0.65", linewidth=2, label="Leave-one-out Pearson range"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#0072B2",
               markeredgecolor="#0072B2", label="Full Pearson"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#E69F00",
               markeredgecolor="#E69F00", label="Full Spearman"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor="white",
               markeredgecolor="#6A3D9A", label="Group-adjusted partial r"),
    ]
    ax_legend.legend(
        handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.05),
        ncol=4, frameon=False, handletextpad=0.4, columnspacing=1.1,
    )
    for label, ax in zip(["a", "b", "c"], [ax_a, ax_b, ax_c]):
        ax.text(-0.08, 1.05, label, transform=ax.transAxes, fontsize=8,
                fontweight="bold", va="top", ha="right")

    basename = FIGURES / "SuppFig_cross_cohort_effect_sensitivity"
    fig.savefig(basename.with_suffix(".svg"))
    fig.savefig(basename.with_suffix(".pdf"))
    fig.savefig(basename.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(basename.with_suffix(".png"), dpi=600)
    plt.close(fig)
    with Image.open(basename.with_suffix(".png")) as image:
        ImageOps.grayscale(image).save(
            basename.with_name(basename.name + "_grayscale").with_suffix(".png"),
            dpi=(600, 600),
        )


def main() -> None:
    ANALYSIS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    effects = cohort_effects()
    meta = meta_summaries(effects)
    loo, sensitivity = correlation_sensitivity()

    effects.to_csv(ANALYSIS / "cross_cohort_CTRP9_effect_sizes.tsv", sep="\t", index=False)
    meta.to_csv(ANALYSIS / "cross_cohort_CTRP9_meta_summary.tsv", sep="\t", index=False)
    loo.to_csv(ANALYSIS / "rat_marker_leave_one_out_correlations.tsv", sep="\t", index=False)
    sensitivity.to_csv(ANALYSIS / "rat_marker_sensitivity_summary.tsv", sep="\t", index=False)
    plot_figure(effects, meta, sensitivity)

    print("COHORT EFFECTS")
    print(effects[["dataset", "context", "n_control", "n_case", "hedges_g", "ci_low", "ci_high", "welch_p", "welch_fdr_bh"]].to_string(index=False))
    print("\nMETA-ANALYSIS")
    print(meta.to_string(index=False))
    print("\nCORRELATION SENSITIVITY SUMMARY")
    for dataset, table in sensitivity.groupby("dataset", sort=False):
        print(
            dataset,
            "Pearson direction-stable", int(table.pearson_direction_stable_all_deletions.sum()),
            "/", len(table),
            "Spearman direction-stable", int(table.spearman_direction_stable_all_deletions.sum()),
            "/", len(table),
            "partial FDR<0.05", int((table.group_adjusted_partial_fdr_bh < 0.05).sum()),
        )


if __name__ == "__main__":
    main()
