#!/usr/bin/env python3
"""Donor-level C1QTNF9 localization and stage tests in GSE345646."""

from __future__ import annotations

import itertools
import importlib.metadata
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "analysis" / "external_validation" / "GSE345646"
RESULTS = PROJECT / "results"
COUNTS_PATH = RESULTS / "GSE345646_qc_filtered.h5ad"
CLUSTER_PATH = RESULTS / "GSE345646_clustered_hvg.h5ad"
ANNOTATION_PATH = PROJECT / "cluster_annotation_v1.tsv"
GROUP_ORDER = ["NF", "pRV", "RVF"]
STAGE = {"NF": 0, "pRV": 1, "RVF": 2}


def bh(values: list[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def exact_labelings(nf_n: int = 4, prv_n: int = 4, rvf_n: int = 3):
    indices = tuple(range(nf_n + prv_n + rvf_n))
    for nf in itertools.combinations(indices, nf_n):
        remaining = tuple(i for i in indices if i not in nf)
        for prv in itertools.combinations(remaining, prv_n):
            labels = np.full(len(indices), 2, dtype=int)
            labels[list(nf)] = 0
            labels[list(prv)] = 1
            yield labels


LABELINGS = list(exact_labelings())


def exact_tests(values: np.ndarray, observed_labels: np.ndarray) -> tuple[float, float, float, float, float]:
    groups = [values[observed_labels == i] for i in range(3)]
    h_obs, _ = stats.kruskal(*groups)
    rho_obs = stats.spearmanr(values, observed_labels).statistic
    h_null = np.empty(len(LABELINGS), dtype=float)
    rho_null = np.empty(len(LABELINGS), dtype=float)
    for i, labels in enumerate(LABELINGS):
        h_null[i] = stats.kruskal(*(values[labels == j] for j in range(3))).statistic
        rho_null[i] = stats.spearmanr(values, labels).statistic
    kw_p = float(np.mean(h_null >= h_obs - 1e-12))
    rho_p = float(np.mean(np.abs(rho_null) >= abs(rho_obs) - 1e-12))
    epsilon2 = max(0.0, float((h_obs - 3 + 1) / (len(values) - 3)))
    return float(h_obs), kw_p, float(rho_obs), rho_p, epsilon2


def main() -> None:
    annotations = pd.read_csv(ANNOTATION_PATH, sep="\t", dtype={"cluster": str})
    cluster_map = annotations.set_index("cluster")["cell_type"].to_dict()
    clustered = sc.read_h5ad(CLUSTER_PATH, backed="r")
    clusters = clustered.obs["leiden_0.5"].astype(str).copy()
    clustered.file.close()

    adata = sc.read_h5ad(COUNTS_PATH)
    clusters = clusters.reindex(adata.obs_names)
    if clusters.isna().any():
        raise ValueError("Cluster labels do not align with count matrix")
    cell_type = clusters.map(cluster_map)
    if cell_type.isna().any():
        raise ValueError("Unannotated primary cluster")
    if "C1QTNF9" not in adata.var_names:
        raise ValueError("C1QTNF9 absent after QC")
    target = np.asarray(adata[:, "C1QTNF9"].X.toarray()).ravel()
    total_counts = np.asarray(adata.X.sum(axis=1)).ravel()
    obs = adata.obs[["donor", "group"]].copy()
    obs["cell_type"] = cell_type.values
    obs["C1QTNF9_count"] = target
    obs["total_count"] = total_counts
    obs["C1QTNF9_positive"] = target > 0

    donor = obs.groupby(["donor", "group", "cell_type"], observed=True).agg(
        nuclei=("C1QTNF9_count", "size"),
        library_size=("total_count", "sum"),
        C1QTNF9_raw_count=("C1QTNF9_count", "sum"),
        C1QTNF9_positive_nuclei=("C1QTNF9_positive", "sum"),
    ).reset_index()
    donor["C1QTNF9_CPM"] = donor["C1QTNF9_raw_count"] / donor["library_size"] * 1e6
    donor["C1QTNF9_log2_CPM_plus_0.5"] = np.log2(donor["C1QTNF9_CPM"] + 0.5)
    donor["C1QTNF9_positive_fraction"] = donor["C1QTNF9_positive_nuclei"] / donor["nuclei"]
    donor["informative_nuclei_ge_20"] = donor["nuclei"] >= 20
    donor.to_csv(RESULTS / "C1QTNF9_pseudobulk_by_donor_celltype.tsv", sep="\t", index=False)

    localization = donor.groupby("cell_type", observed=True).agg(
        nuclei=("nuclei", "sum"),
        donors=("donor", "nunique"),
        C1QTNF9_raw_count=("C1QTNF9_raw_count", "sum"),
        donors_nonzero=("C1QTNF9_raw_count", lambda x: int((x > 0).sum())),
        C1QTNF9_positive_nuclei=("C1QTNF9_positive_nuclei", "sum"),
        library_size=("library_size", "sum"),
    ).reset_index()
    localization["C1QTNF9_CPM"] = localization["C1QTNF9_raw_count"] / localization["library_size"] * 1e6
    localization["C1QTNF9_positive_fraction"] = localization["C1QTNF9_positive_nuclei"] / localization["nuclei"]
    localization.to_csv(RESULTS / "C1QTNF9_celltype_localization.tsv", sep="\t", index=False)

    test_rows = []
    for cell, frame in donor.groupby("cell_type", observed=True):
        frame = frame[frame["informative_nuclei_ge_20"]].copy()
        group_n = frame.groupby("group", observed=True)["donor"].nunique().to_dict()
        total_target = int(frame["C1QTNF9_raw_count"].sum())
        eligible = total_target >= 20 and all(group_n.get(g, 0) >= 3 for g in GROUP_ORDER)
        row = {"cell_type": cell, "eligible": eligible, "C1QTNF9_raw_count": total_target,
               **{f"{g}_donor_n": int(group_n.get(g, 0)) for g in GROUP_ORDER}}
        for g in GROUP_ORDER:
            vals = frame.loc[frame["group"] == g, "C1QTNF9_log2_CPM_plus_0.5"]
            row[f"{g}_median_log2CPM_plus_0.5"] = float(vals.median()) if len(vals) else np.nan
        if eligible:
            # All eligible classes have the fixed 4/4/3 donor structure.
            frame["group"] = pd.Categorical(frame["group"], categories=GROUP_ORDER, ordered=True)
            frame = frame.sort_values(["group", "donor"])
            values = frame["C1QTNF9_log2_CPM_plus_0.5"].to_numpy()
            labels = frame["group"].map(STAGE).to_numpy(dtype=int)
            h, kw_p, rho, rho_p, epsilon2 = exact_tests(values, labels)
            row.update({"kruskal_H": h, "kruskal_exact_p": kw_p, "epsilon_squared": epsilon2,
                        "spearman_rho": rho, "spearman_exact_p": rho_p})
        test_rows.append(row)
    tests = pd.DataFrame(test_rows)
    eligible_mask = tests["eligible"].fillna(False)
    tests["kruskal_BH_FDR"] = np.nan
    tests["spearman_BH_FDR"] = np.nan
    if eligible_mask.any():
        tests.loc[eligible_mask, "kruskal_BH_FDR"] = bh(tests.loc[eligible_mask, "kruskal_exact_p"].tolist())
        tests.loc[eligible_mask, "spearman_BH_FDR"] = bh(tests.loc[eligible_mask, "spearman_exact_p"].tolist())
    tests.to_csv(RESULTS / "C1QTNF9_stage_tests.tsv", sep="\t", index=False)

    perturbation_candidates = localization[
        (localization["C1QTNF9_raw_count"] >= 50)
        & (localization["donors_nonzero"] >= np.ceil(0.25 * 11))
    ]["cell_type"].tolist()
    summary = {
        "retained_nuclei": int(len(obs)),
        "cell_types": int(localization.shape[0]),
        "eligible_stage_test_cell_types": int(eligible_mask.sum()),
        "comparison_correction": "BH separately across eligible cell types for Kruskal-Wallis and Spearman families",
        "virtual_perturbation_gate": "PASS" if perturbation_candidates else "FAIL_STOP",
        "virtual_perturbation_candidate_cell_types": perturbation_candidates,
        "software": {
            "python": platform.python_version(),
            **{name: importlib.metadata.version(name) for name in
               ("scanpy", "anndata", "scipy", "harmonypy", "umap-learn", "scrublet")},
        },
    }
    (RESULTS / "C1QTNF9_single_nucleus_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(localization.sort_values("C1QTNF9_CPM", ascending=False).to_string(index=False))
    print(tests.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
