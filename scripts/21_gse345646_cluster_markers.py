#!/usr/bin/env python3
"""Donor-aware cluster summaries for blinded GSE345646 cell annotation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "analysis" / "external_validation" / "GSE345646" / "results"
COUNTS_PATH = BASE / "GSE345646_qc_filtered.h5ad"
CLUSTER_PATH = BASE / "GSE345646_clustered_hvg.h5ad"
CLUSTER_KEY = "leiden_0.5"

MARKER_PANEL = [
    "TNNT2", "MYH6", "MYH7", "RYR2", "TTN", "NPPA", "NPPB",  # cardiomyocyte
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "POSTN",  # fibroblast
    "VWF", "PECAM1", "KDR", "EMCN", "FLT1", "CDH5", "PTPRB", "EGFL7", "ERG",
    "ESAM", "ENG", "RAMP2", "CA4", "ACKR1", "PLVAP", "EDNRB", "RGCC",  # endothelial
    "RGS5", "CSPG4", "MCAM", "PDGFRB", "ABCC9", "KCNJ8",  # pericyte
    "ACTA2", "TAGLN", "MYH11", "CNN1",  # smooth muscle
    "PTPRC", "LST1", "TYROBP", "CD68", "C1QA", "FCER1G",  # myeloid/immune
    "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "NKG7",  # lymphoid
    "KIT", "TPSAB1", "MS4A2",  # mast
    "PLP1", "S100B", "SOX10", "NRXN1", "MPZ", "SLC1A3", "GFAP",  # neural/Schwann
    "WT1", "TBX18", "ALDH1A2", "KRT19", "KRT8",  # mesothelial/epicardial
    "ADIPOQ", "PLIN1", "LPL", "FABP4",  # adipocyte
    "FBLN1", "FBLN2", "C7", "PI16", "DPT", "FAP", "THY1", "CILP", "COMP", "CFH",
    "CCL21", "PROX1", "FLT4", "PDPN",  # stromal and lymphatic distinction
]


def main() -> None:
    clustered = sc.read_h5ad(CLUSTER_PATH, backed="r")
    assignments = clustered.obs[[CLUSTER_KEY]].copy()
    clustered.file.close()

    adata = sc.read_h5ad(COUNTS_PATH)
    if not adata.obs_names.equals(assignments.index):
        assignments = assignments.reindex(adata.obs_names)
    if assignments[CLUSTER_KEY].isna().any():
        raise ValueError("Missing cluster assignment after index alignment")
    clusters = assignments[CLUSTER_KEY].astype(str)
    levels = sorted(clusters.unique(), key=int)
    codes = pd.Categorical(clusters, categories=levels).codes
    design = sparse.csr_matrix(
        (np.ones(adata.n_obs, dtype=np.float32), (np.arange(adata.n_obs), codes)),
        shape=(adata.n_obs, len(levels)),
    )
    aggregate = (design.T @ adata.X).tocsr()
    library = np.asarray(aggregate.sum(axis=1)).ravel()
    cpm = aggregate.multiply((1e6 / library)[:, None]).toarray()
    logcpm = np.log2(cpm + 1.0)

    # One-vs-rest cluster specificity; target gene is excluded from annotation.
    rest = (logcpm.sum(axis=0, keepdims=True) - logcpm) / (len(levels) - 1)
    specificity = logcpm - rest
    genes = np.asarray(adata.var_names)
    marker_rows = []
    for row, cluster in enumerate(levels):
        order = np.argsort(-specificity[row])
        rank = 0
        for idx in order:
            if genes[idx] == "C1QTNF9":
                continue
            marker_rows.append({
                "cluster": cluster,
                "rank": rank + 1,
                "gene": genes[idx],
                "log2CPM": logcpm[row, idx],
                "specificity_vs_mean_other_log2CPM": specificity[row, idx],
                "raw_count": int(aggregate[row, idx]),
            })
            rank += 1
            if rank == 40:
                break
    pd.DataFrame(marker_rows).to_csv(BASE / "cluster_top_markers_blinded.tsv", sep="\t", index=False)

    marker_rows = []
    gene_index = {gene: i for i, gene in enumerate(genes)}
    for cluster_i, cluster in enumerate(levels):
        for gene in MARKER_PANEL:
            if gene in gene_index:
                j = gene_index[gene]
                marker_rows.append({"cluster": cluster, "gene": gene,
                                    "raw_count": int(aggregate[cluster_i, j]),
                                    "CPM": float(cpm[cluster_i, j]),
                                    "log2CPM": float(logcpm[cluster_i, j])})
    pd.DataFrame(marker_rows).to_csv(BASE / "cluster_marker_panel_blinded.tsv", sep="\t", index=False)

    obs = adata.obs[["donor", "group"]].copy()
    obs["cluster"] = clusters.values
    donor_counts = obs.groupby(["cluster", "donor"], observed=True).size().rename("nuclei").reset_index()
    donor_counts["fraction_within_donor"] = donor_counts["nuclei"] / donor_counts.groupby("donor")["nuclei"].transform("sum")
    donor_counts.to_csv(BASE / "cluster_donor_representation.tsv", sep="\t", index=False)
    pd.DataFrame({"cluster": levels, "nuclei": np.bincount(codes, minlength=len(levels)),
                  "raw_library_size": library.astype(int)}).to_csv(
        BASE / "cluster_sizes.tsv", sep="\t", index=False
    )
    print(pd.DataFrame(marker_rows).pivot(index="cluster", columns="gene", values="log2CPM").round(2).to_string())


if __name__ == "__main__":
    main()
