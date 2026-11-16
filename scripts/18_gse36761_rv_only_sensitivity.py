#!/usr/bin/env python3
"""GSE36761/GSE36759 sensitivity analysis restricted to right-ventricular controls.

This post-results sensitivity analysis corrects the public sample annotation:
NH-01/03/05/07 are LV and NH-02/04/06/08 are RV samples from four donors.
The frozen data directory is read only. Derived outputs are written to analysis/.

Sources:
  data/processed/GSE36761_log2_matrix.tsv
  data/processed/GSE36761_C1QTNF9_transcript_sum.tsv
  data/processed/GSE36759_mirna_log2_matrix.tsv
  results/ssgsea_GSE36761_human.tsv
  analysis/GSE36761_RV_only_limma.tsv
  analysis/TargetScanHuman_7.2_default_context_predictions.zip
"""

from pathlib import Path
import json
import zipfile

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis"
RV_CONTROLS = ["NH-02", "NH-04", "NH-06", "NH-08"]
LV_CONTROLS = ["NH-01", "NH-03", "NH-05", "NH-07"]
TOF = [f"TOF-{i:02d}" for i in range(1, 23)]
TARGETSCAN_MEMBER = "Predicted_Targets_Context_Scores.default_predictions.txt"
RNG = np.random.default_rng(20260823)


def hedges_g(case, control):
    case, control = np.asarray(case, float), np.asarray(control, float)
    n1, n0 = len(case), len(control)
    pooled = np.sqrt(((n1 - 1) * case.var(ddof=1) + (n0 - 1) * control.var(ddof=1)) / (n1 + n0 - 2))
    d = (case.mean() - control.mean()) / pooled
    return d * (1 - 3 / (4 * (n1 + n0) - 9))


def bootstrap_mean_diff(case, control, n_boot=20000):
    case, control = np.asarray(case, float), np.asarray(control, float)
    values = np.empty(n_boot)
    for i in range(n_boot):
        values[i] = RNG.choice(case, len(case), True).mean() - RNG.choice(control, len(control), True).mean()
    return np.quantile(values, [0.025, 0.975])


def main():
    gene = pd.read_csv(ROOT / "data/processed/GSE36761_log2_matrix.tsv", sep="\t", index_col=0)
    ctrp9_raw = pd.read_csv(ROOT / "data/processed/GSE36761_C1QTNF9_transcript_sum.tsv", sep="\t", index_col=0)
    gene.loc["C1QTNF9"] = np.log2(ctrp9_raw.loc["C1QTNF9"] + 1)
    mir = pd.read_csv(ROOT / "data/processed/GSE36759_mirna_log2_matrix.tsv", sep="\t", index_col=0)
    expected = set(RV_CONTROLS + LV_CONTROLS + TOF)
    if not expected.issubset(gene.columns) or not expected.issubset(mir.columns):
        raise ValueError("Expected GSE36761/GSE36759 samples are missing")

    provenance = pd.DataFrame([
        {"sample": s, "group": "NH", "chamber": "RV" if s in RV_CONTROLS else "LV", "donor_unit": f"NH_donor_{(int(s[-2:]) + 1)//2}"}
        for s in RV_CONTROLS + LV_CONTROLS
    ] + [{"sample": s, "group": "TOF", "chamber": "RV", "donor_unit": s} for s in TOF])
    provenance.to_csv(OUT / "GSE36761_sample_provenance_corrected.tsv", sep="\t", index=False)

    x = gene.loc["C1QTNF9", RV_CONTROLS].to_numpy(float)
    y = gene.loc["C1QTNF9", TOF].to_numpy(float)
    _, p = stats.ttest_ind(y, x, equal_var=False)
    lo, hi = bootstrap_mean_diff(y, x)
    expression = pd.DataFrame([{
        "contrast": "TOF_RV_vs_NH_RV", "n_NH_RV": len(x), "n_TOF_RV": len(y),
        "mean_NH_RV": x.mean(), "mean_TOF_RV": y.mean(), "mean_difference": y.mean()-x.mean(),
        "bootstrap_CI_low": lo, "bootstrap_CI_high": hi, "hedges_g": hedges_g(y, x),
        "welch_p": p, "multiple_comparison_correction": "none; one prespecified gene-level sensitivity test"
    }])
    expression.to_csv(OUT / "GSE36761_RV_only_C1QTNF9.tsv", sep="\t", index=False)

    # Cell-type-associated ssGSEA scores: BH across all ten prespecified populations.
    scores = pd.read_csv(ROOT / "results/ssgsea_GSE36761_human.tsv", sep="\t", index_col=0)
    rows = []
    for cell in scores.index:
        a = scores.loc[cell, RV_CONTROLS].to_numpy(float)
        b = scores.loc[cell, TOF].to_numpy(float)
        _, pv = stats.ttest_ind(b, a, equal_var=False)
        ci_lo, ci_hi = bootstrap_mean_diff(b, a)
        rows.append({"cell_type": cell, "n_NH_RV": len(a), "n_TOF_RV": len(b),
                     "NH_RV_mean": a.mean(), "TOF_RV_mean": b.mean(), "mean_diff": b.mean()-a.mean(),
                     "CI_low": ci_lo, "CI_high": ci_hi, "hedges_g": hedges_g(b, a), "PValue": pv})
    immune = pd.DataFrame(rows)
    immune["BH_FDR"] = multipletests(immune.PValue, method="fdr_bh")[1]
    immune.sort_values("PValue").to_csv(OUT / "GSE36761_RV_only_ssgsea.tsv", sep="\t", index=False)

    # Sensitivity of the previously discovered 79-gene ECM-enriched ME4 signature.
    # Module membership is held fixed; this is not a second WGCNA discovery run.
    modules = pd.read_csv(ROOT / "results/wgcna_GSE36761_human_gene_modules.tsv", sep="\t")
    me4_genes = modules.loc[modules.module.eq(4), "gene"].tolist()
    me4_genes = [g for g in me4_genes if g in gene.index]
    z = gene.loc[me4_genes, RV_CONTROLS + TOF].T
    z = (z - z.mean(axis=0)) / z.std(axis=0, ddof=1)
    u, s, _ = np.linalg.svd(z.to_numpy(float), full_matrices=False)
    pc1 = u[:, 0] * s[0]
    mean_signature = z.mean(axis=1).to_numpy(float)
    if np.corrcoef(pc1, mean_signature)[0, 1] < 0:
        pc1 = -pc1
    pc1 = pd.Series(pc1, index=z.index, name="ME4_fixed_signature_PC1")
    a, b = pc1.loc[RV_CONTROLS].to_numpy(float), pc1.loc[TOF].to_numpy(float)
    _, pv = stats.ttest_ind(b, a, equal_var=False)
    ci_lo, ci_hi = bootstrap_mean_diff(b, a)
    me4 = pd.DataFrame([{"signature": "original_GSE36761_module4_fixed_genes",
                         "genes_available": len(me4_genes), "n_NH_RV": len(a), "n_TOF_RV": len(b),
                         "variance_explained_PC1": (s[0] ** 2) / np.sum(s ** 2),
                         "NH_RV_mean": a.mean(), "TOF_RV_mean": b.mean(), "mean_diff": b.mean()-a.mean(),
                         "CI_low": ci_lo, "CI_high": ci_hi, "hedges_g": hedges_g(b, a), "welch_p": pv,
                         "multiple_comparison_correction": "none; one prespecified fixed-signature sensitivity test"}])
    me4.to_csv(OUT / "GSE36761_RV_only_ME4_fixed_signature.tsv", sep="\t", index=False)
    pc1.rename_axis("sample").reset_index().to_csv(OUT / "GSE36761_RV_only_ME4_fixed_signature_scores.tsv", sep="\t", index=False)

    # miRNA differential expression in the same RV-only contrast.
    mir_rows = []
    for name, values in mir.iterrows():
        a, b = values[RV_CONTROLS].to_numpy(float), values[TOF].to_numpy(float)
        _, pv = stats.ttest_ind(b, a, equal_var=False)
        mir_rows.append({"miRNA": name, "log2FC": b.mean()-a.mean(), "PValue": pv})
    mir_de = pd.DataFrame(mir_rows)
    mir_de["BH_FDR"] = multipletests(mir_de.PValue, method="fdr_bh")[1]
    mir_de.to_csv(OUT / "GSE36759_RV_only_mirna_DE.tsv", sep="\t", index=False)

    gene_de = pd.read_csv(OUT / "GSE36761_RV_only_limma.tsv", sep="\t")
    de_mir = mir_de[mir_de.BH_FDR < 0.05].set_index("miRNA")
    de_gene = gene_de[gene_de["adj.P.Val"] < 0.05].set_index("gene")
    with zipfile.ZipFile(OUT / "TargetScanHuman_7.2_default_context_predictions.zip") as archive:
        with archive.open(TARGETSCAN_MEMBER) as handle:
            ts = pd.read_csv(handle, sep="\t")
    eligible = ts[ts.miRNA.isin(de_mir.index) & ts["Gene Symbol"].isin(de_gene.index)]
    pairs = eligible.groupby(["miRNA", "Gene Symbol"], as_index=False).agg(
        max_weighted_context_percentile=("weighted context++ score percentile", "max"),
        n_predicted_sites=("Site Type", "size")).rename(columns={"Gene Symbol": "target"})
    common = RV_CONTROLS + TOF
    edge_rows = []
    for pair in pairs.itertuples(index=False):
        r, pv = stats.pearsonr(mir.loc[pair.miRNA, common], gene.loc[pair.target, common])
        edge_rows.append({"miRNA": pair.miRNA, "target": pair.target,
                          "mir_log2FC": de_mir.loc[pair.miRNA, "log2FC"],
                          "gene_log2FC": de_gene.loc[pair.target, "logFC"],
                          "pearson_r": r, "PValue": pv,
                          "opposite_group_direction": de_mir.loc[pair.miRNA, "log2FC"] * de_gene.loc[pair.target, "logFC"] < 0,
                          "n_matched_samples": len(common),
                          "max_weighted_context_percentile": pair.max_weighted_context_percentile,
                          "n_predicted_sites": pair.n_predicted_sites})
    tested = pd.DataFrame(edge_rows)
    if len(tested):
        tested["BH_FDR"] = multipletests(tested.PValue, method="fdr_bh")[1]
        tested["passes_network_filter"] = tested.opposite_group_direction & (tested.pearson_r < 0) & (tested.BH_FDR < 0.05)
        tested = tested.sort_values("PValue")
    tested.to_csv(OUT / "GSE36761_RV_only_candidate_network_all.tsv", sep="\t", index=False)
    edges = tested[tested.passes_network_filter].copy() if len(tested) else tested.copy()
    edges.to_csv(OUT / "GSE36761_RV_only_candidate_network_edges.tsv", sep="\t", index=False)

    summary = {
        "analysis_label": "post-results RV-only sensitivity analysis",
        "NH_RV_samples": RV_CONTROLS, "excluded_NH_LV_samples": LV_CONTROLS,
        "TOF_RV_n": len(TOF), "miRNA_tests": len(mir_de),
        "miRNA_DE_BH_FDR_lt_0.05": int((mir_de.BH_FDR < .05).sum()),
        "mRNA_DE_BH_FDR_lt_0.05": int((gene_de["adj.P.Val"] < .05).sum()),
        "candidate_pairs_tested": len(tested), "candidate_edges_retained": len(edges),
        "ME4_fixed_signature_genes": len(me4_genes),
        "random_seed": 20260823,
        "corrections": {"mRNA": "limma BH genome-wide", "miRNA": "BH across 626 miRNAs",
                        "ssGSEA": "BH across 10 cell types", "network": "BH across all eligible unique TargetScan pairs"}
    }
    (OUT / "GSE36761_RV_only_sensitivity_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("\nC1QTNF9\n", expression.to_string(index=False))
    print("\nssGSEA\n", immune.sort_values("PValue").to_string(index=False))
    print("\nME4 fixed signature\n", me4.to_string(index=False))


if __name__ == "__main__":
    main()
