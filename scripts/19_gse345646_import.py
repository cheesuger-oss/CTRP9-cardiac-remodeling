#!/usr/bin/env python3
"""Import the 11 included GSE345646 CellBender matrices into AnnData."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "analysis" / "external_validation" / "GSE345646"
WORK = BASE / "work"
OUT = BASE / "results" / "GSE345646_included_raw.h5ad"

SAMPLES = [
    ("GSM10012320", "1343", "RVF", "GSM10012320_RV_1343_6K_100_2_6_cellbender_filtered.h5"),
    ("GSM10012321", "1392", "pRV", "GSM10012321_RV_1392_21K_100_2_6_cellbender_filtered.h5"),
    ("GSM10012322", "1467", "RVF", "GSM10012322_RV_1467_16K_100_2_6_cellbender_filtered.h5"),
    ("GSM10012323", "1561", "NF", "GSM10012323_RV_1561_20K_100_1_5_cellbender_filtered.h5"),
    ("GSM10012324", "1567", "pRV", "GSM10012324_RV_1567_25K_100_1_5_cellbender_filtered.h5"),
    ("GSM10012325", "1618", "pRV", "GSM10012325_RV_1618_17K_100_5_6_cellbender_filtered.h5"),
    ("GSM10012326", "1632", "RVF", "GSM10012326_RV_1632_25K_100_1_5_cellbender_filtered.h5"),
    ("GSM10012327", "1681", "NF", "GSM10012327_RV_1681_16K_100_2_6_cellbender_filtered.h5"),
    ("GSM10012328", "1691", "NF", "GSM10012328_RV_1691_17K_100_cellbender_5_6_filtered.h5"),
    ("GSM10012329", "1692", "pRV", "GSM10012329_RV_1692_17K_100_2_6_cellbender_filtered.h5"),
    ("GSM10012330", "1697", "NF", "GSM10012330_RV_1697_19K_100_2_6_cellbender_filtered.h5"),
]


def main() -> None:
    objects: dict[str, ad.AnnData] = {}
    metadata_rows = []
    for gsm, donor, group, filename in SAMPLES:
        obj = sc.read_10x_h5(WORK / filename, gex_only=True)
        obj.var_names_make_unique(join="__dup")
        obj.obs["gsm"] = gsm
        obj.obs["donor"] = donor
        obj.obs["group"] = group
        obj.obs["source_file"] = filename
        objects[gsm] = obj
        metadata_rows.append({"gsm": gsm, "donor": donor, "group": group, "source_file": filename,
                              "author_cellbender_nuclei": obj.n_obs})

    combined = ad.concat(
        objects,
        axis=0,
        join="inner",
        merge="same",
        label="library",
        index_unique="__",
        pairwise=False,
    )
    combined.obs_names_make_unique()
    combined.uns["provenance"] = {
        "accession": "GSE345646",
        "download_date": "2026-11-15",
        "archive_sha256": "02cd1dd99b4ac06132d4378a44d2b43b3d160c70b205b076b83c68d0d7217a61",
        "excluded_failed_qc_sample": "GSM10012319",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(OUT, compression="gzip")
    pd.DataFrame(metadata_rows).to_csv(BASE / "results" / "included_sample_metadata.tsv", sep="\t", index=False)
    print(f"Wrote {OUT}: {combined.n_obs} nuclei x {combined.n_vars} features")


if __name__ == "__main__":
    main()
