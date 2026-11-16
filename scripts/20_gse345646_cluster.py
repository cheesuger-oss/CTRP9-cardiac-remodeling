#!/usr/bin/env python3
"""Memory-bounded clustering of QC-filtered GSE345646 RV nuclei."""

from __future__ import annotations

from pathlib import Path

import harmonypy as hm
import scanpy as sc
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "analysis" / "external_validation" / "GSE345646" / "results"
INPUT = BASE / "GSE345646_qc_filtered.h5ad"
OUTPUT = BASE / "GSE345646_clustered_hvg.h5ad"
SEED = 20261115

if not hasattr(sparse.csr_array, "getnnz"):
    sparse.csr_array.getnnz = sparse.csr_array._getnnz  # type: ignore[attr-defined]


def main() -> None:
    adata = sc.read_h5ad(INPUT)
    # The immutable input file remains the source of full-gene raw counts.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        flavor="cell_ranger",
        n_top_genes=3000,
        batch_key="donor",
    )
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.pca(adata, n_comps=50, svd_solver="arpack", random_state=SEED)

    harmony = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        "donor",
        random_state=SEED,
        max_iter_harmony=20,
        verbose=True,
    )
    corrected = harmony.Z_corr
    if corrected.shape != adata.obsm["X_pca"].shape:
        raise ValueError(
            f"Unexpected Harmony shape {corrected.shape}; expected {adata.obsm['X_pca'].shape}"
        )
    adata.obsm["X_pca_harmony"] = corrected
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    for resolution in (0.3, 0.5, 0.8):
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=f"leiden_{resolution}",
            flavor="igraph",
            directed=False,
            n_iterations=2,
            random_state=SEED,
        )
    adata.uns["analysis_contract"] = {
        "seed": SEED,
        "hvg_n": 3000,
        "hvg_flavor": "cell_ranger",
        "hvg_batch_key": "donor",
        "pca_n": 50,
        "neighbors_n": 15,
        "neighbors_pcs": 40,
        "integration": "Harmony on PCA; raw counts unchanged in separate QC file",
    }
    adata.write_h5ad(OUTPUT, compression="gzip")
    print(f"Wrote {OUTPUT}: {adata.n_obs} nuclei x {adata.n_vars} HVGs")
    for key in ("leiden_0.3", "leiden_0.5", "leiden_0.8"):
        print(key, adata.obs[key].value_counts().sort_index().to_dict())


if __name__ == "__main__":
    main()
