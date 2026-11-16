#!/usr/bin/env python3
"""Build a read-only provenance manifest for the final submission package."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis/submission_package_v2/reproducibility_manifest.json"

CORE_SCRIPTS = [
    "scripts/01_parse.py", "scripts/02_map_ids.py", "scripts/03_ctrp9_core.py",
    "scripts/05_immune_ssgsea.py", "scripts/06_ceRNA.py",
    "scripts/08_figures_core.py", "scripts/09_figures_analyses.py",
    "scripts/10_figures_extra.py", "scripts/11_figures_forest.py",
    "scripts/12_edger_lite.py", "scripts/13_module_overlap.py",
    "scripts/14_gse57338_module_trait.py", "scripts/15_correlation_multiplicity.py",
    "scripts/16_wgcna_trait_multiplicity.py", "scripts/17_cerna_multiplicity.py",
    "scripts/17_module_preservation.R", "scripts/18_gse36761_rv_only_limma.R",
    "scripts/18_gse36761_rv_only_sensitivity.py", "scripts/20_gse119754_tmm_voom_primary.R",
    "scripts/18_gse345646_feasibility.py", "scripts/19_gse345646_import.py",
    "scripts/20_gse345646_cluster.py", "scripts/21_gse345646_cluster_markers.py",
    "scripts/22_gse345646_ctrp9_pseudobulk.py", "scripts/23_gse345646_supp_figure.py",
    "scripts/24_cross_cohort_effect_sensitivity.py",
]

PACKAGE_SOURCES = [
    "analysis/submission_package_v2/manuscript_en.md",
    "analysis/submission_package_v2/manuscript_zh.md",
    "analysis/submission_package_v2/supplementary_information.md",
    "analysis/submission_package_v2/figure_legends.md",
    "analysis/submission_package_v2/claim_evidence_registry.md",
    "analysis/submission_package_v2/terminology_ledger.md",
    "analysis/submission_package_v2/build_docx_package.py",
    "analysis/submission_package_v2/build_jmcc_full_mock_preview.py",
    "analysis/submission_package_v2/audit_submission_package.py",
    "analysis/submission_package_v2/FINAL_CONCLUSION_LOCK.md",
    "analysis/submission_package_v2/INTERNAL_PEER_REVIEW.md",
    "analysis/submission_package_v2/REPRODUCIBILITY_README.md",
    "analysis/submission_package_v2/PROTEOMICS_PLACEMENT_DECISION.md",
    "analysis/submission_package_v2/SUBMISSION_RISK_REGISTER.md",
    "analysis/submission_package_v2/FINAL_PRESUBMISSION_REVIEW.md",
    "analysis/submission_package_v2/figure_legends_zh.md",
    "analysis/submission_package_v2/build_jmcc_bilingual_mock_preview.py",
    "scripts/33_build_submission_reproducibility_manifest.py",
]

KEY_OUTPUTS = [
    "analysis/GSE119754_TMM_voom_primary_summary.tsv",
    "analysis/GSE36761_RV_only_C1QTNF9.tsv",
    "analysis/GSE57338_C1QTNF9_age_sex_adjusted.tsv",
    "analysis/cross_cohort_CTRP9_effect_sizes.tsv",
    "analysis/cross_cohort_CTRP9_meta_summary.tsv",
    "analysis/rat_marker_sensitivity_summary.tsv",
    "analysis/wgcna_cross_dataset_module_preservation_targets.tsv",
    "analysis/external_validation/GSE345646/results/C1QTNF9_celltype_localization.tsv",
    "analysis/external_validation/GSE345646/results/C1QTNF9_stage_tests.tsv",
    "docs/data/md5_manifest.tsv",
    "analysis/2026-08-17_independent_reanalysis/core/R_sessionInfo.txt",
    "analysis/external_validation/GSE345646/requirements-scanpy.txt",
    "analysis/submission_package_v2/JMCC_bilingual_mock_publication_preview.pdf",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def record(rel: str) -> dict[str, object]:
    path = ROOT / rel
    return {
        "path": rel,
        "exists": path.is_file(),
        "bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256(path) if path.is_file() else None,
    }


def command_output(command: list[str]) -> str | None:
    try:
        return subprocess.run(
            command, cwd=ROOT, check=True, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def package_versions() -> dict[str, str | None]:
    names = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "seaborn", "PyYAML"]
    versions: dict[str, str | None] = {}
    for name in names:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def main() -> int:
    manifest = {
        "schema": "rovot.submission_reproducibility.v1",
        "generated_on": date.today().isoformat(),
        "project_root": str(ROOT),
        "frozen_data_policy": "No file under data/ may be modified; integrity is governed by docs/data/md5_manifest.tsv and scripts/verify_data.py --full.",
        "python": {"version": sys.version, "platform": platform.platform(), "packages": package_versions()},
        "r_version": command_output([
            str(ROOT / ".toolchain/miniconda3/bin/conda"), "run", "-n", "r_env",
            str(ROOT / ".toolchain/miniconda3/envs/r_env/bin/Rscript"), "--version"
        ]),
        "git_commit": command_output(["git", "rev-parse", "HEAD"]),
        "core_scripts": [record(x) for x in CORE_SCRIPTS],
        "package_sources": [record(x) for x in PACKAGE_SOURCES],
        "key_outputs": [record(x) for x in KEY_OUTPUTS],
        "excluded_internal_extension": {
            "paths": ["scripts/25-32", "analysis/proteomics_raw_reanalysis/"],
            "reason": "The prespecified PXD043768 direct-protein identification gate failed; no protein result enters the submitted evidence chain."
        },
    }
    missing = [
        item["path"] for group in ("core_scripts", "package_sources", "key_outputs")
        for item in manifest[group] if not item["exists"]
    ]
    warnings = []
    if manifest["git_commit"] is None:
        warnings.append("No Git commit was available in this restored workspace; a versioned public release is still required.")
    warnings.append("The Python package list is a runtime snapshot, not a fully resolved cross-platform lock file.")
    manifest["missing_required_files"] = missing
    manifest["warnings"] = warnings
    manifest["result"] = "PASS_WITH_WARNINGS" if not missing else "FAIL"
    OUT.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"output": str(OUT.relative_to(ROOT)), "result": manifest["result"], "missing": missing}, ensure_ascii=False))
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
