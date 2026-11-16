#!/usr/bin/env python3
"""Outcome-blind feasibility audit for GSE345646 CellBender matrices.

The audit reads only matrix dimensions, feature names, and cohort-wide raw
C1QTNF9 counts. It does not inspect condition-specific expression differences.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "analysis" / "external_validation" / "GSE345646"
WORK = BASE / "work"
RESULTS = BASE / "results"
ARCHIVE = BASE / "downloads" / "GSE345646_RAW.tar"

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


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in values])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inspect_matrix(path: Path) -> dict[str, object]:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(x) for x in matrix["shape"][:])
        genes = decode(matrix["features"]["name"][:])
        feature_ids = decode(matrix["features"]["id"][:])
        if len(set(feature_ids)) != len(feature_ids):
            raise ValueError(f"Non-unique feature IDs in {path.name}")
        matches = np.flatnonzero(genes == "C1QTNF9")
        if len(matches) != 1:
            raise ValueError(f"Expected one C1QTNF9 feature in {path.name}; found {len(matches)}")
        counts = sparse.csc_matrix(
            (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
            shape=shape,
        )
        target = counts[int(matches[0]), :]
        return {
            "feature_n": shape[0],
            "duplicated_gene_symbol_n": int(len(genes) - len(set(genes))),
            "nucleus_n": shape[1],
            "C1QTNF9_raw_umi": int(target.sum()),
            "C1QTNF9_positive_nuclei": int(target.getnnz()),
        }


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    for gsm, donor, group, filename in SAMPLES:
        path = WORK / filename
        if not path.exists():
            raise FileNotFoundError(path)
        rows.append({"gsm": gsm, "donor": donor, "group": group, "file": filename, **inspect_matrix(path)})

    table = pd.DataFrame(rows)
    table.to_csv(RESULTS / "feasibility_by_donor.tsv", sep="\t", index=False)
    group_n = table.groupby("group")["donor"].nunique().to_dict()
    total_umi = int(table["C1QTNF9_raw_umi"].sum())
    nonzero_donors = int((table["C1QTNF9_raw_umi"] > 0).sum())
    gates = {
        "eleven_readable_unique_feature_id_matrices": len(table) == 11,
        "at_least_three_donors_per_group": all(group_n.get(g, 0) >= 3 for g in ("NF", "pRV", "RVF")),
        "C1QTNF9_present_once_per_matrix": True,
        "C1QTNF9_at_least_20_umi_and_three_nonzero_donors": total_umi >= 20 and nonzero_donors >= 3,
        "at_least_500_nuclei_per_donor": bool((table["nucleus_n"] >= 500).all()),
    }
    summary = {
        "accession": "GSE345646",
        "download_date": "2026-11-15",
        "archive_sha256": sha256(ARCHIVE),
        "included_donors": len(table),
        "group_donor_n": group_n,
        "total_nuclei_author_cellbender_filtered": int(table["nucleus_n"].sum()),
        "C1QTNF9_total_raw_umi": total_umi,
        "C1QTNF9_nonzero_donors": nonzero_donors,
        "gates": gates,
        "overall_verdict": "PASS" if all(gates.values()) else "FAIL_STOP",
    }
    (RESULTS / "feasibility_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
