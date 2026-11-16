#!/usr/bin/env python3
"""Apply BH correction to GSE36761 WGCNA module-trait tests.

Sources:
  results/wgcna_GSE36761_human_module_trait_cor.tsv
  results/wgcna_GSE36761_human_module_trait_p.tsv

Output:
  analysis/GSE36761_wgcna_module_trait_bh.tsv
"""

from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
COR_PATH = ROOT / "results" / "wgcna_GSE36761_human_module_trait_cor.tsv"
P_PATH = ROOT / "results" / "wgcna_GSE36761_human_module_trait_p.tsv"
OUT_PATH = ROOT / "analysis" / "GSE36761_wgcna_module_trait_bh.tsv"


def main() -> None:
    correlations = pd.read_csv(COR_PATH, sep="\t", index_col=0).iloc[:, 0]
    p_values = pd.read_csv(P_PATH, sep="\t", index_col=0).iloc[:, 0]
    if not correlations.index.equals(p_values.index):
        raise ValueError("Module order differs between correlation and P-value tables")

    result = pd.DataFrame(
        {
            "module": correlations.index,
            "pearson_r_TOF": correlations.to_numpy(),
            "PValue": p_values.to_numpy(),
        }
    )
    result["BH_FDR"] = multipletests(result["PValue"], method="fdr_bh")[1]
    result["n_module_tests"] = len(result)
    result.to_csv(OUT_PATH, sep="\t", index=False)
    print(result.to_string(index=False))
    print(f"Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()
