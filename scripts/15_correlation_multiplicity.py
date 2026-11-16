#!/usr/bin/env python3
"""Apply within-cohort BH correction to the prespecified CTRP9 correlations.

Frozen inputs
-------------
results/corr_GSE210180_rat.tsv
results/corr_GSE119754_rat.tsv

Output
------
analysis/correlation_multiplicity.tsv

Pearson and Spearman P values are corrected separately within each cohort
across every marker present in the original correlation result table.  The
script does not filter on nominal significance before correction.
"""

from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis" / "correlation_multiplicity.tsv"

COHORTS = {
    "GSE210180": {
        "path": ROOT / "results" / "corr_GSE210180_rat.tsv",
        "n": 6,
        "expected_tests": 10,
    },
    "GSE119754": {
        "path": ROOT / "results" / "corr_GSE119754_rat.tsv",
        "n": 11,
        "expected_tests": 7,
    },
}


def main() -> None:
    rows = []
    for dataset, spec in COHORTS.items():
        table = pd.read_csv(spec["path"], sep="\t")
        if len(table) != spec["expected_tests"]:
            raise ValueError(
                f"{dataset}: expected {spec['expected_tests']} marker tests, "
                f"found {len(table)}"
            )
        required = {
            "marker",
            "pearson_r",
            "pearson_p",
            "spearman_rho",
            "spearman_p",
        }
        missing = required.difference(table.columns)
        if missing:
            raise ValueError(f"{dataset}: missing columns {sorted(missing)}")

        table = table.copy()
        table["pearson_fdr_bh"] = multipletests(
            table["pearson_p"], method="fdr_bh"
        )[1]
        table["spearman_fdr_bh"] = multipletests(
            table["spearman_p"], method="fdr_bh"
        )[1]
        table.insert(0, "dataset", dataset)
        table.insert(1, "n", spec["n"])
        table["pearson_significant_fdr05"] = table["pearson_fdr_bh"] < 0.05
        table["spearman_significant_fdr05"] = table["spearman_fdr_bh"] < 0.05
        table["correction_family"] = (
            f"all {len(table)} prespecified available markers within {dataset}; "
            "Pearson and Spearman corrected separately"
        )
        rows.append(table)

    result = pd.concat(rows, ignore_index=True)
    result.to_csv(OUT, sep="\t", index=False)
    print(f"Saved {OUT.relative_to(ROOT)} ({len(result)} rows)")
    for dataset, group in result.groupby("dataset", sort=False):
        print(
            dataset,
            "Pearson FDR<0.05:",
            int(group["pearson_significant_fdr05"].sum()),
            "Spearman FDR<0.05:",
            int(group["spearman_significant_fdr05"].sum()),
        )


if __name__ == "__main__":
    main()
