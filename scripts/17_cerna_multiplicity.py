#!/usr/bin/env python3
"""Rebuild the GSE36759/GSE36761 candidate miRNA-mRNA network with BH correction.

The original exploratory script filtered pairwise correlations at nominal
P < 0.05 and did not retain the full tested family. This audit deduplicates
TargetScan site-level predictions to unique miRNA-gene pairs, tests every
eligible pair across the 30 matched samples, and applies BH correction across
the complete correlation family.

Sources:
  analysis/TargetScanHuman_7.2_default_context_predictions.zip
  results/mirna_DEG_GSE36759.tsv
  results/limma_GSE36761_human.tsv
  data/processed/GSE36759_mirna_log2_matrix.tsv
  data/processed/GSE36761_log2_matrix.tsv

Outputs:
  analysis/ceRNA_all_candidate_correlations_bh.tsv
  analysis/ceRNA_edges_bh.tsv
  analysis/ceRNA_multiplicity_summary.tsv
"""

from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, ttest_ind
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
TARGETSCAN = ROOT / "analysis" / "TargetScanHuman_7.2_default_context_predictions.zip"
TARGETSCAN_MEMBER = "Predicted_Targets_Context_Scores.default_predictions.txt"


def main() -> None:
    mir_de = pd.read_csv(ROOT / "results/mirna_DEG_GSE36759.tsv", sep="\t")
    gene_de = pd.read_csv(ROOT / "results/limma_GSE36761_human.tsv", sep="\t")
    mir_expr = pd.read_csv(
        ROOT / "data/processed/GSE36759_mirna_log2_matrix.tsv", sep="\t", index_col=0
    )
    mir_meta = pd.read_csv(ROOT / "data/processed/GSE36759_metadata.tsv", sep="\t")
    gene_expr = pd.read_csv(
        ROOT / "data/processed/GSE36761_log2_matrix.tsv", sep="\t", index_col=0
    )

    nh_samples = mir_meta.loc[mir_meta["group"].eq("NH"), "sample"].tolist()
    tof_samples = mir_meta.loc[mir_meta["group"].eq("TOF"), "sample"].tolist()
    mir_recomputed_rows = []
    for mirna, values in mir_expr.iterrows():
        _, p_value = ttest_ind(
            values[tof_samples].to_numpy(dtype=float),
            values[nh_samples].to_numpy(dtype=float),
            equal_var=False,
        )
        mir_recomputed_rows.append(
            {
                "miRNA": mirna,
                "log2FC": values[tof_samples].mean() - values[nh_samples].mean(),
                "pvalue": p_value,
            }
        )
    mir_recomputed = pd.DataFrame(mir_recomputed_rows)
    mir_recomputed["padj"] = multipletests(mir_recomputed["pvalue"], method="fdr_bh")[1]
    mir_check = mir_recomputed.merge(mir_de, on="miRNA", suffixes=("_new", "_saved"))
    for column in ("log2FC", "pvalue", "padj"):
        if not np.allclose(
            mir_check[f"{column}_new"], mir_check[f"{column}_saved"], equal_nan=True
        ):
            raise ValueError(f"Saved miRNA differential-expression values differ for {column}")

    de_mir = mir_recomputed.loc[mir_recomputed["padj"] < 0.05].copy()
    de_gene = gene_de.loc[gene_de["adj.P.Val"] < 0.05].copy()
    mir_fc = de_mir.set_index("miRNA")["log2FC"]
    gene_fc = de_gene.set_index("gene")["logFC"]

    with zipfile.ZipFile(TARGETSCAN) as archive:
        with archive.open(TARGETSCAN_MEMBER) as handle:
            targetscan = pd.read_csv(handle, sep="\t")

    eligible = targetscan.loc[
        targetscan["miRNA"].isin(mir_fc.index)
        & targetscan["Gene Symbol"].isin(gene_fc.index)
    ].copy()
    # TargetScan is site/transcript level. A correlation is a unique
    # miRNA-gene hypothesis, so collapse duplicate sites before testing.
    pairs = (
        eligible.groupby(["miRNA", "Gene Symbol"], as_index=False)
        .agg(
            max_weighted_context_percentile=(
                "weighted context++ score percentile",
                "max",
            ),
            n_predicted_sites=("Site Type", "size"),
        )
        .rename(columns={"Gene Symbol": "target"})
    )

    common_samples = [sample for sample in gene_expr.columns if sample in mir_expr.columns]
    rows = []
    for pair in pairs.itertuples(index=False):
        if pair.miRNA not in mir_expr.index or pair.target not in gene_expr.index:
            continue
        x = mir_expr.loc[pair.miRNA, common_samples].to_numpy(dtype=float)
        y = gene_expr.loc[pair.target, common_samples].to_numpy(dtype=float)
        r, p_value = pearsonr(x, y)
        if not np.isfinite(r):
            continue
        mir_log2fc = float(mir_fc.loc[pair.miRNA])
        gene_log2fc = float(gene_fc.loc[pair.target])
        rows.append(
            {
                "miRNA": pair.miRNA,
                "target": pair.target,
                "mir_log2FC": mir_log2fc,
                "gene_log2FC": gene_log2fc,
                "pearson_r": r,
                "PValue": p_value,
                "opposite_group_direction": (mir_log2fc * gene_log2fc) < 0,
                "max_weighted_context_percentile": pair.max_weighted_context_percentile,
                "n_predicted_sites": pair.n_predicted_sites,
                "n_matched_samples": len(common_samples),
            }
        )

    tested = pd.DataFrame(rows).sort_values("PValue").reset_index(drop=True)
    tested["BH_FDR"] = multipletests(tested["PValue"], method="fdr_bh")[1]
    tested["negative_correlation"] = tested["pearson_r"] < 0
    tested["passes_network_filter"] = (
        tested["opposite_group_direction"]
        & tested["negative_correlation"]
        & (tested["BH_FDR"] < 0.05)
    )

    edges = tested.loc[tested["passes_network_filter"]].copy()
    c1qtnf9_predictions = targetscan.loc[
        targetscan["miRNA"].str.startswith("hsa-", na=False)
        & targetscan["Gene Symbol"].eq("C1QTNF9")
    ]
    summary = pd.DataFrame(
        [
            {
                "de_mirnas_bh_fdr_lt_0_05": len(de_mir),
                "mirna_welch_tests": len(mir_recomputed),
                "mirna_n_NH": len(nh_samples),
                "mirna_n_TOF": len(tof_samples),
                "de_mrnas_bh_fdr_lt_0_05": len(de_gene),
                "unique_targetscan_pairs_tested": len(tested),
                "nominal_negative_opposite_edges_p_lt_0_05": int(
                    (
                        tested["opposite_group_direction"]
                        & tested["negative_correlation"]
                        & (tested["PValue"] < 0.05)
                    ).sum()
                ),
                "bh_significant_negative_opposite_edges": len(edges),
                "targetscan_default_hsa_sites_for_C1QTNF9": len(c1qtnf9_predictions),
                "n_matched_samples": len(common_samples),
                "correlation_multiple_comparison_correction": (
                    "Benjamini-Hochberg across all unique eligible TargetScan miRNA-gene pairs"
                ),
            }
        ]
    )

    tested.to_csv(
        ROOT / "analysis/ceRNA_all_candidate_correlations_bh.tsv", sep="\t", index=False
    )
    edges.to_csv(ROOT / "analysis/ceRNA_edges_bh.tsv", sep="\t", index=False)
    summary.to_csv(ROOT / "analysis/ceRNA_multiplicity_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    if len(edges):
        print("\nBH-significant candidate edges:")
        print(edges.to_string(index=False))


if __name__ == "__main__":
    main()
