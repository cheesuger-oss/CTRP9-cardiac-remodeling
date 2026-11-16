# C1QTNF9/CTRP9 and context-dependent cardiac remodeling

This repository contains the analysis code, shareable derived data, source tables and vector figures supporting the public-data study of C1QTNF9/CTRP9 across right-ventricular pressure overload, tetralogy of Fallot and terminal left-ventricular failure.

The principal conclusion is deliberately bounded: CTRP9 is associated with cardiac remodeling, while myocardial C1QTNF9 expression varies by disease context. The repository does not claim a validated clinical biomarker, a universal mechanism or therapeutic efficacy.

## Contents

- `scripts/`: Python and R analysis and figure scripts.
- `data/processed/`: shareable frozen processed matrices and metadata used by the main analyses.
- `data/aux/`: GTEx v8 tissue summary and the predefined MCP-counter gene list.
- `analysis/`: statistical outputs and source tables used in the manuscript.
- `analysis/external_validation/GSE345646/`: donor-level single-nucleus summaries; large H5AD intermediates are not redistributed.
- `results/`: differential-expression, correlation, WGCNA, ssGSEA and ceRNA result tables.
- `figures/`: editable SVG versions of six main and five supplementary figures.
- `docs/`: data dictionary and figure-quality record.

## Public source datasets

GEO accessions: GSE210180, GSE119754, GSE224428, GSE36761, GSE36759, GSE217771, GSE219101, GSE57338 and GSE345646. GTEx v8 was used for tissue-level baseline expression. See `THIRD_PARTY_DATA.md` before redistributing any upstream file.

## Environment

Core analyses were run with Python 3.9 and R 4.5.3. Install the Python packages with:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

R scripts require `limma` and `WGCNA`. The independent single-nucleus workflow has a separate environment record at `analysis/external_validation/GSE345646/requirements-scanpy.txt`.

## Reproducing key corrected analyses

Run commands from the repository root:

```bash
.venv/bin/python scripts/13_module_overlap.py
.venv/bin/python scripts/14_gse57338_module_trait.py
.venv/bin/python scripts/15_correlation_multiplicity.py
.venv/bin/python scripts/16_wgcna_trait_multiplicity.py
.venv/bin/python scripts/24_cross_cohort_effect_sensitivity.py
```

The ceRNA and chamber-matched sensitivity workflows require the versioned TargetScan 7.2 input described in `THIRD_PARTY_DATA.md`. Large GSE345646 H5AD intermediates must be reconstructed from GEO with scripts 18–22; donor-level outputs required to inspect the reported statistics are included.

## Integrity

`release_manifest.json` records SHA-256 hashes and sizes for every archived file except the manifest itself. Verify it with:

```bash
.venv/bin/python scripts/build_release_manifest.py --verify
```

The manuscript project used a separate frozen-data gate; immediately before this release it returned `RESULT: 12 passed, 0 failed` for all 45 locked project inputs.

## Licensing

Code is released under the MIT License. Original derived tables and documentation are released under CC BY 4.0. Third-party source files remain governed by their source terms and are not relicensed; see `DATA_LICENSE.md` and `THIRD_PARTY_DATA.md`.

## AI-assisted preparation

OpenAI Codex assisted with language translation, computer programming, preparation of graphical illustrations and data visualizations, language editing and formatting. The authors reviewed and reran AI-assisted code using the specified source data, checked figures against the underlying results, and retain full responsibility. Generative AI was not used to fabricate scientific data, statistical results or experimental images.

## Archival status

This is release candidate `v0.1.0`. A Zenodo DOI should be minted from the final GitHub release only after the complete author/creator metadata have been confirmed.

