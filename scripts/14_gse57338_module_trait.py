#!/usr/bin/env python3
"""Test whether the GSE57338 C1QTNF9-containing WGCNA module differs by disease group.

Frozen inputs (read only):
  analysis/GSE57338_gene_matrix.tsv
  analysis/wgcna_GSE57338_gene_modules.tsv
  data/raw/GSE57338_series_matrix.txt.gz
  data/processed/GSE57338_chars.tsv

Outputs:
  analysis/GSE57338_M3_excl_C1QTNF9_eigengene_samples.tsv
  analysis/GSE57338_M3_excl_C1QTNF9_group_summary.tsv
  analysis/GSE57338_M3_excl_C1QTNF9_posthoc.tsv
  analysis/GSE57338_M3_excl_C1QTNF9_adjusted.tsv
  analysis/GSE57338_M3_excl_C1QTNF9_correlation.tsv

The eigengene is PC1 of sample-wise standardized expression for the 788 module-3
genes other than C1QTNF9. Its sign is oriented to correlate positively with
C1QTNF9, then standardized to
mean 0 and SD 1 across all 313 samples. Primary inference is one-way ANOVA plus
three pre-specified Welch pairwise contrasts with Bonferroni and BH correction.
An OLS model adjusted for age and sex is reported as a sensitivity analysis.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"


def sample_groups() -> dict[str, str]:
    titles = pd.read_csv(ROOT / "data/processed/GSE57338_titles.tsv", sep="\t")["sample"].tolist()
    accessions = pd.read_csv(
        ANALYSIS / "GSE57338_gene_matrix.tsv", sep="\t", nrows=1, index_col=0
    ).columns.tolist()
    if len(titles) != len(accessions):
        raise RuntimeError("GSE57338 titles are not aligned to expression columns")
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


def mean_diff_ci(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """Welch 95% CI for mean(a) - mean(b)."""
    diff = float(np.mean(a) - np.mean(b))
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    se = np.sqrt(va / len(a) + vb / len(b))
    df = (va / len(a) + vb / len(b)) ** 2 / (
        (va / len(a)) ** 2 / (len(a) - 1) + (vb / len(b)) ** 2 / (len(b) - 1)
    )
    crit = stats.t.ppf(0.975, df)
    return diff, diff - crit * se, diff + crit * se


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    """Bias-corrected standardized mean difference for a - b."""
    na, nb = len(a), len(b)
    pooled = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    d = (np.mean(a) - np.mean(b)) / pooled
    correction = 1 - 3 / (4 * (na + nb) - 9)
    return float(correction * d)


def main() -> None:
    matrix = pd.read_csv(ANALYSIS / "GSE57338_gene_matrix.tsv", sep="\t", index_col=0)
    modules = pd.read_csv(ANALYSIS / "wgcna_GSE57338_gene_modules.tsv", sep="\t")
    module_genes = modules.loc[modules["module"].eq(3), "gene"].tolist()
    available = [gene for gene in module_genes if gene in matrix.index]
    if "C1QTNF9" not in available or len(available) != len(module_genes):
        raise RuntimeError("Module-3 genes are not fully represented in the expression matrix")

    predictor_genes = [gene for gene in available if gene != "C1QTNF9"]
    x = matrix.loc[predictor_genes].T.astype(float)
    z = (x - x.mean(axis=0)) / x.std(axis=0, ddof=1)
    if z.isna().any().any():
        raise RuntimeError("Zero-variance or missing module genes encountered")
    pca = PCA(n_components=1)
    eigengene = pca.fit_transform(z).ravel()
    c1qtnf9 = matrix.loc["C1QTNF9", x.index].astype(float).values
    if stats.pearsonr(eigengene, c1qtnf9).statistic < 0:
        eigengene *= -1
    eigengene = (eigengene - eigengene.mean()) / eigengene.std(ddof=1)

    groups = sample_groups()
    if list(x.index) != list(groups):
        raise RuntimeError("Expression columns do not match GEO series sample order")
    sample_df = pd.DataFrame({"sample": x.index, "group": [groups[s] for s in x.index], "M3_eigengene": eigengene})

    chars = pd.read_csv(ROOT / "data/processed/GSE57338_chars.tsv", sep="\t")
    if len(chars) != len(sample_df):
        raise RuntimeError("GSE57338 characteristics are not aligned to the 313 samples")
    sample_df["sex"] = chars["3"].str.replace("gender: ", "", regex=False).str.lower().values
    sample_df["age"] = pd.to_numeric(chars["4"].str.replace("age: ", "", regex=False), errors="coerce").values
    sample_df.to_csv(ANALYSIS / "GSE57338_M3_excl_C1QTNF9_eigengene_samples.tsv", sep="\t", index=False)

    order = ["NF", "ICM", "DCM"]
    summary = sample_df.groupby("group")["M3_eigengene"].agg(n="size", mean="mean", sd="std", median="median").reindex(order)
    summary.to_csv(ANALYSIS / "GSE57338_M3_excl_C1QTNF9_group_summary.tsv", sep="\t")

    arrays = [sample_df.loc[sample_df.group.eq(group), "M3_eigengene"].values for group in order]
    f_stat, anova_p = stats.f_oneway(*arrays)
    comparisons = [("ICM", "NF"), ("DCM", "NF"), ("DCM", "ICM")]
    rows = []
    for case, control in comparisons:
        a = sample_df.loc[sample_df.group.eq(case), "M3_eigengene"].values
        b = sample_df.loc[sample_df.group.eq(control), "M3_eigengene"].values
        diff, low, high = mean_diff_ci(a, b)
        _, p = stats.ttest_ind(a, b, equal_var=False)
        rows.append({
            "contrast": f"{case} vs {control}", "n_case": len(a), "n_control": len(b),
            "mean_diff": diff, "ci_low": low, "ci_high": high,
            "hedges_g": hedges_g(a, b), "p_raw": p,
        })
    posthoc = pd.DataFrame(rows)
    posthoc["p_bonferroni"] = multipletests(posthoc.p_raw, method="bonferroni")[1]
    posthoc["p_bh"] = multipletests(posthoc.p_raw, method="fdr_bh")[1]
    posthoc.insert(0, "anova_F", f_stat)
    posthoc.insert(1, "anova_p", anova_p)
    posthoc.to_csv(ANALYSIS / "GSE57338_M3_excl_C1QTNF9_posthoc.tsv", sep="\t", index=False)

    adjusted_data = sample_df.dropna(subset=["age", "sex"]).copy()
    adjusted_data["group"] = pd.Categorical(adjusted_data["group"], categories=order)
    model = smf.ols("M3_eigengene ~ C(group, Treatment(reference='NF')) + age + C(sex)", data=adjusted_data).fit(cov_type="HC3")
    coefficients = model.summary2().tables[1].reset_index().rename(columns={"index": "term"})
    coefficients.insert(0, "n", int(model.nobs))
    coefficients.to_csv(ANALYSIS / "GSE57338_M3_excl_C1QTNF9_adjusted.tsv", sep="\t", index=False)

    corr = stats.pearsonr(sample_df.M3_eigengene, c1qtnf9)
    fisher_z = np.arctanh(corr.statistic)
    fisher_se = 1 / np.sqrt(len(sample_df) - 3)
    corr_low, corr_high = np.tanh(fisher_z + np.array([-1, 1]) * stats.norm.ppf(0.975) * fisher_se)
    pd.DataFrame([{
        "n": len(sample_df), "module_predictor_genes": len(predictor_genes),
        "pc1_variance_explained": pca.explained_variance_ratio_[0],
        "pearson_r": corr.statistic, "ci_low": corr_low, "ci_high": corr_high,
        "p_raw": corr.pvalue, "multiple_comparison_correction": "none (one pre-specified correlation)",
    }]).to_csv(ANALYSIS / "GSE57338_M3_excl_C1QTNF9_correlation.tsv", sep="\t", index=False)
    print(f"module predictors={len(predictor_genes)}; samples={len(sample_df)}; C1QTNF9 r={corr.statistic:.6f}, p={corr.pvalue:.6g}")
    print(f"ANOVA F={f_stat:.6f}, p={anova_p:.6g}")
    print(summary.to_string())
    print(posthoc.to_string(index=False))
    print("adjusted sensitivity:")
    print(coefficients.to_string(index=False))


if __name__ == "__main__":
    main()
