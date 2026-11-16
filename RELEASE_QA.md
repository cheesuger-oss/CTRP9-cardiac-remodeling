# Public release quality record

Date: 2026-11-16  
Release candidate: v0.1.0

## Verification

- Frozen project inputs: `RESULT: 12 passed, 0 failed` for 45 files.
- Public-release size before Git metadata: approximately 124 MB.
- Files larger than 95 MB: none.
- Python syntax compilation: PASS.
- R parsing: PASS for eight scripts.
- Private-path/contact/token scan: PASS; no author email, access token or personal absolute path is present.
- Key corrected analyses rerun from the isolated release directory: PASS.
- Exact output comparison against the frozen project passed for the module-overlap summary, correlation multiplicity table, WGCNA trait-multiplicity table and cross-cohort effect-size table.
- Internal submission documents, author-review files, large H5AD objects and the failed proteomics feasibility pilot are excluded.

## Nature-family review

**PASS.** The public bundle preserves the manuscript's evidence boundaries, including null findings, correction families and heterogeneous cross-context effects. It provides editable vector figures and source tables without adding stronger mechanistic or biomarker claims.

## Light workflow review

**PASS.** Scripts use relative repository paths, shareable derived inputs are separated from third-party sources, licenses and provenance are explicit, and every released file is covered by a SHA-256 manifest. Frozen `data/` in the source project was not modified.

The reviews agree; there is no conflict.

## Remaining archival action

GitHub provides a public version-controlled source location but not a DOI. A final GitHub release should be archived in Zenodo only after the complete creator list and affiliations are confirmed. The resulting version DOI must then replace the provisional repository wording in both journal packages.

