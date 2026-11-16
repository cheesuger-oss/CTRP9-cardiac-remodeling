#!/usr/bin/env python3
"""Cross-dataset overlap of the two ECM-enriched human WGCNA modules.

Sources
-------
analysis/wgcna_GSE57338_gene_modules.tsv (module 3; human LV failure)
results/wgcna_GSE36761_human_gene_modules.tsv (module 4; human TOF RV)

The test universe is the intersection of gene symbols present in both WGCNA
module tables. Fisher's exact test is one-sided (alternative='greater'); the
reported hypergeometric upper-tail probability is mathematically equivalent.
"""
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact, hypergeom


ROOT = Path(__file__).resolve().parents[1]
LEFT_PATH = ROOT / "analysis" / "wgcna_GSE57338_gene_modules.tsv"
RIGHT_PATH = ROOT / "results" / "wgcna_GSE36761_human_gene_modules.tsv"
OUT_SUMMARY = ROOT / "analysis" / "wgcna_ECM_module_overlap_summary.tsv"
OUT_GENES = ROOT / "analysis" / "wgcna_ECM_module_overlap_shared_genes.tsv"


def load_unique(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype={"gene": str, "module": int})
    frame = frame.dropna(subset=["gene", "module"]).copy()
    frame["gene"] = frame["gene"].str.strip()
    frame = frame[frame["gene"] != ""]
    conflicts = frame.groupby("gene")["module"].nunique()
    if (conflicts > 1).any():
        raise ValueError(f"Conflicting module assignments in {path}")
    return frame.drop_duplicates("gene").set_index("gene")


def main() -> None:
    left = load_unique(LEFT_PATH)
    right = load_unique(RIGHT_PATH)
    universe = sorted(set(left.index) & set(right.index))
    if not universe:
        raise ValueError("The two WGCNA tables have no shared gene-symbol universe")

    left_module_all = set(left.index[left["module"] == 3])
    right_module_all = set(right.index[right["module"] == 4])
    universe_set = set(universe)
    left_module = left_module_all & universe_set
    right_module = right_module_all & universe_set
    shared = sorted(left_module & right_module)

    a = len(shared)
    b = len(left_module - right_module)
    c = len(right_module - left_module)
    d = len(universe_set - (left_module | right_module))
    odds_ratio, fisher_p = fisher_exact([[a, b], [c, d]], alternative="greater")
    hypergeom_p = hypergeom.sf(a - 1, len(universe), len(left_module), len(right_module))
    expected = len(left_module) * len(right_module) / len(universe)
    fold_enrichment = a / expected if expected else float("nan")

    summary = pd.DataFrame([{
        "universe_intersection_n": len(universe),
        "GSE57338_module3_total_n": len(left_module_all),
        "GSE36761_module4_total_n": len(right_module_all),
        "GSE57338_module3_in_universe_n": len(left_module),
        "GSE36761_module4_in_universe_n": len(right_module),
        "shared_n": a,
        "expected_shared_n": expected,
        "fold_enrichment": fold_enrichment,
        "fisher_odds_ratio": odds_ratio,
        "fisher_greater_p": fisher_p,
        "hypergeom_upper_p": hypergeom_p,
        "fisher_table_a_shared": a,
        "fisher_table_b_left_only": b,
        "fisher_table_c_right_only": c,
        "fisher_table_d_neither": d,
        "multiple_comparison_correction": "none (single prespecified overlap test)",
    }])
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    pd.DataFrame({"gene": shared}).to_csv(OUT_GENES, sep="\t", index=False)
    print(summary.to_string(index=False))
    print(f"Shared genes: {', '.join(shared) if shared else '(none)'}")


if __name__ == "__main__":
    main()
