#!/usr/bin/env python3
"""Sensitivity analyses for the manuscript's central C1QTNF9 claims.

Read-only inputs
----------------
data/processed/GSE210180_log2fpkm_matrix.tsv
data/processed/GSE210180_metadata.tsv
data/processed/GSE119754_counts_symbol.tsv
data/processed/GSE119754_metadata.tsv
data/raw/GSE57338_series_matrix.txt.gz
data/processed/GSE57338_chars.tsv
analysis/GSE57338_gene_matrix.tsv
results/corr_GSE210180_rat.tsv
results/corr_GSE119754_rat.tsv

Outputs
-------
analysis/rat_marker_group_adjusted_partial_correlations.tsv
analysis/GSE57338_C1QTNF9_age_sex_adjusted.tsv
analysis/core_claim_sensitivity_summary.tsv

The rat analysis residualises both C1qtnf9 and each predefined marker on a
binary disease-group indicator, then calculates a partial Pearson correlation.
Benjamini-Hochberg correction is applied to the complete marker family within
each dataset. This evaluates whether pooled correlations persist beyond the
between-group contrast; it is exploratory and low-powered.

The GSE57338 analysis fits C1QTNF9 expression as a function of disease group,
age and sex with HC3 robust covariance. Three adjusted group contrasts are
reported with Bonferroni and BH correction.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"


def partial_correlation(x: np.ndarray, y: np.ndarray, group: np.ndarray) -> dict[str, float]:
    """Return the Pearson correlation between residuals after group adjustment."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    group = np.asarray(group, dtype=float)
    if not (len(x) == len(y) == len(group)):
        raise ValueError("x, y and group must have identical lengths")
    if len(x) < 5 or np.unique(group).size != 2:
        raise ValueError("partial correlation requires at least five observations in two groups")
    design = np.column_stack([np.ones(len(group)), group])
    x_resid = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    y_resid = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    if np.isclose(np.std(x_resid, ddof=1), 0) or np.isclose(np.std(y_resid, ddof=1), 0):
        raise ValueError("group-adjusted residual variance is zero")
    r = float(np.corrcoef(x_resid, y_resid)[0, 1])
    df = len(x) - 3
    if np.isclose(abs(r), 1):
        p_value = 0.0
    else:
        t_stat = r * np.sqrt(df / (1 - r * r))
        p_value = float(2 * stats.t.sf(abs(t_stat), df))
    fisher_se = 1 / np.sqrt(len(x) - 4)
    z = np.arctanh(np.clip(r, -0.999999999, 0.999999999))
    ci_low, ci_high = np.tanh(z + np.array([-1, 1]) * stats.norm.ppf(0.975) * fisher_se)
    return {
        "n": int(len(x)),
        "partial_r": r,
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "df": int(df),
        "p_raw": p_value,
    }


def parse_characteristics(chars: pd.DataFrame) -> pd.DataFrame:
    """Extract age and sex from the frozen GSE57338 characteristics table."""
    required = {"3", "4"}
    if not required.issubset(chars.columns):
        raise ValueError(f"characteristics table lacks columns {sorted(required - set(chars.columns))}")
    parsed = pd.DataFrame(index=chars.index)
    parsed["sex"] = chars["3"].str.replace("gender: ", "", regex=False).str.lower()
    parsed["age"] = pd.to_numeric(chars["4"].str.replace("age: ", "", regex=False), errors="coerce")
    return parsed


def fit_age_sex_adjusted_contrasts(data: pd.DataFrame, outcome: str) -> pd.DataFrame:
    """Fit an HC3 OLS model and return all three disease-group contrasts."""
    required = {"group", "age", "sex", outcome}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"adjusted model data lack columns {sorted(missing)}")
    model_data = data.dropna(subset=list(required)).copy()
    model_data["group"] = pd.Categorical(model_data["group"], categories=["NF", "ICM", "DCM"])
    model = smf.ols(
        f"{outcome} ~ C(group, Treatment(reference='NF')) + age + C(sex)",
        data=model_data,
    ).fit(cov_type="HC3")
    names = list(model.params.index)
    term_icm = "C(group, Treatment(reference='NF'))[T.ICM]"
    term_dcm = "C(group, Treatment(reference='NF'))[T.DCM]"
    if term_icm not in names or term_dcm not in names:
        raise RuntimeError("expected disease-group coefficients were not estimable")

    contrast_specs = [
        ("ICM vs NF", {term_icm: 1.0}),
        ("DCM vs NF", {term_dcm: 1.0}),
        ("DCM vs ICM", {term_dcm: 1.0, term_icm: -1.0}),
    ]
    rows = []
    counts = model_data["group"].value_counts()
    for label, weights in contrast_specs:
        vector = np.zeros(len(names))
        for term, weight in weights.items():
            vector[names.index(term)] = weight
        test = model.t_test(vector)
        ci = np.asarray(test.conf_int(alpha=0.05)).ravel()
        rows.append({
            "contrast": label,
            "model_n": int(model.nobs),
            "n_NF": int(counts.get("NF", 0)),
            "n_ICM": int(counts.get("ICM", 0)),
            "n_DCM": int(counts.get("DCM", 0)),
            "adjusted_difference": float(np.asarray(test.effect).ravel()[0]),
            "ci_low": float(ci[0]),
            "ci_high": float(ci[1]),
            "p_raw": float(np.asarray(test.pvalue).ravel()[0]),
            "covariance": "HC3",
            "covariates": "age + sex",
        })
    result = pd.DataFrame(rows)
    result["p_bonferroni"] = multipletests(result["p_raw"], method="bonferroni")[1]
    result["p_bh"] = multipletests(result["p_raw"], method="fdr_bh")[1]
    return result


def gse57338_groups() -> dict[str, str]:
    """Parse sample accession-to-group mapping from the frozen GEO series matrix."""
    titles = accessions = None
    series = ROOT / "data/raw/GSE57338_series_matrix.txt.gz"
    with gzip.open(series, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            elif line.startswith("!Sample_geo_accession"):
                accessions = next(csv.reader([line.rstrip("\n")], delimiter="\t"))[1:]
            if titles is not None and accessions is not None:
                break
    if titles is None or accessions is None or len(titles) != len(accessions):
        raise RuntimeError("could not parse aligned GSE57338 titles and accessions")
    mapping = {}
    for accession, title in zip(accessions, titles):
        if "Non-failing" in title:
            mapping[accession] = "NF"
        elif "Ischemic" in title:
            mapping[accession] = "ICM"
        elif "Idiopathic Dilated CMP" in title:
            mapping[accession] = "DCM"
        else:
            raise ValueError(f"unrecognised GSE57338 title: {title}")
    return mapping


def rat_partial_correlations() -> pd.DataFrame:
    """Compute group-adjusted partial correlations for both rat cohorts."""
    specs = {
        "GSE210180": {
            "matrix": ROOT / "data/processed/GSE210180_log2fpkm_matrix.tsv",
            "metadata": ROOT / "data/processed/GSE210180_metadata.tsv",
            "correlations": ROOT / "results/corr_GSE210180_rat.tsv",
            "case": "pab",
            "transform": "none",
        },
        "GSE119754": {
            "matrix": ROOT / "data/processed/GSE119754_counts_symbol.tsv",
            "metadata": ROOT / "data/processed/GSE119754_metadata.tsv",
            "correlations": ROOT / "results/corr_GSE119754_rat.tsv",
            "case": "MCT",
            "transform": "log2CPM+1",
        },
    }
    rows = []
    for dataset, spec in specs.items():
        matrix = pd.read_csv(spec["matrix"], sep="\t", index_col=0)
        metadata = pd.read_csv(spec["metadata"], sep="\t").set_index("sample")
        matrix = matrix.loc[:, metadata.index]
        if spec["transform"] == "log2CPM+1":
            matrix = np.log2(matrix.div(matrix.sum(axis=0), axis=1) * 1e6 + 1)
        markers = pd.read_csv(spec["correlations"], sep="\t")["marker"].tolist()
        group = metadata["group"].eq(spec["case"]).astype(float).to_numpy()
        x = matrix.loc["C1qtnf9"].to_numpy(float)
        dataset_rows = []
        for marker in markers:
            result = partial_correlation(x, matrix.loc[marker].to_numpy(float), group)
            dataset_rows.append({
                "dataset": dataset,
                "marker": marker,
                "adjustment": "binary disease group",
                "analysis_status": "post hoc sensitivity",
                **result,
            })
        table = pd.DataFrame(dataset_rows)
        table["p_bh_within_dataset"] = multipletests(table["p_raw"], method="fdr_bh")[1]
        rows.append(table)
    return pd.concat(rows, ignore_index=True)


def gse57338_adjusted() -> pd.DataFrame:
    """Fit the age/sex-adjusted GSE57338 C1QTNF9 model."""
    matrix = pd.read_csv(ANALYSIS / "GSE57338_gene_matrix.tsv", sep="\t", index_col=0)
    groups = gse57338_groups()
    if list(matrix.columns) != list(groups):
        raise RuntimeError("GSE57338 matrix columns do not match GEO sample order")
    chars = pd.read_csv(ROOT / "data/processed/GSE57338_chars.tsv", sep="\t")
    parsed = parse_characteristics(chars)
    if len(parsed) != matrix.shape[1]:
        raise RuntimeError("GSE57338 characteristics are not aligned to expression samples")
    model_data = pd.DataFrame({
        "sample": matrix.columns,
        "group": [groups[s] for s in matrix.columns],
        "expression": matrix.loc["C1QTNF9"].to_numpy(float),
        "age": parsed["age"].to_numpy(),
        "sex": parsed["sex"].to_numpy(),
    })
    result = fit_age_sex_adjusted_contrasts(model_data, "expression")
    result.insert(0, "gene", "C1QTNF9")
    result["analysis_status"] = "post hoc sensitivity"
    return result


def main() -> None:
    ANALYSIS.mkdir(exist_ok=True)
    rat = rat_partial_correlations()
    rat_path = ANALYSIS / "rat_marker_group_adjusted_partial_correlations.tsv"
    rat.to_csv(rat_path, sep="\t", index=False)

    adjusted = gse57338_adjusted()
    adjusted_path = ANALYSIS / "GSE57338_C1QTNF9_age_sex_adjusted.tsv"
    adjusted.to_csv(adjusted_path, sep="\t", index=False)

    summary_rows = []
    for dataset, table in rat.groupby("dataset", sort=False):
        summary_rows.append({
            "analysis": "group-adjusted rat marker partial correlations",
            "dataset": dataset,
            "tests": len(table),
            "fdr_lt_0_05": int((table["p_bh_within_dataset"] < 0.05).sum()),
            "interpretation_boundary": "post hoc; low power; evaluates pooled group separation",
        })
    summary_rows.append({
        "analysis": "age/sex-adjusted C1QTNF9 group contrasts",
        "dataset": "GSE57338",
        "tests": len(adjusted),
        "fdr_lt_0_05": int((adjusted["p_bh"] < 0.05).sum()),
        "interpretation_boundary": "post hoc; residual confounding remains possible",
    })
    pd.DataFrame(summary_rows).to_csv(
        ANALYSIS / "core_claim_sensitivity_summary.tsv", sep="\t", index=False
    )
    print(rat.to_string(index=False))
    print(adjusted.to_string(index=False))


if __name__ == "__main__":
    main()
