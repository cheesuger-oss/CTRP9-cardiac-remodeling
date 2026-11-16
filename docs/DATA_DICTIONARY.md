# DATA DICTIONARY — Frozen Data Package

Companion to `README.md`. Documents every file under `data/` (raw, aux, processed).
Checksums for all files: `docs/data/md5_manifest.tsv`.

---

## A. Processed matrices (`data/processed/`)

### A1. GSE210180 (rat neonatal PAB, RV)
- `GSE210180_log2fpkm_matrix.tsv` — rows = gene symbols (28,938), cols = 6 samples
  (`P21_Sham_1..3`, `P21_PAB_1..3`). Values = log2(FPKM+1). Gene annotation from the
  provider's `gene_name` column (duplicates collapsed by max mean).
- `GSE210180_metadata.tsv` — cols: `sample, group (sham|pab), species, model, tissue`.

### A2. GSE119754 (rat MCT-PAH RV failure)
- `GSE119754_counts_symbol.tsv` — rows = gene symbols or `MSTRG.*` ids (37,383), cols = 11 samples
  (`CT_RV_1..6`, `MCT_RV_1..4,6`; MCT_RV_5 absent in source). RAW COUNTS. Mixed IDs: ENSRNOG→symbol
  mapped via mygene; `MSTRG.*` rows kept with original ID (novel StringTie transcripts).
- `GSE119754_metadata.tsv` — cols: `sample, group (control|MCT), species, model, tissue`.
- Raw reference: `GSE119754_all_counts.csv.gz` (original mixed-ID counts),
  `GSE119754_condition_res_ct_rv-mct_rv.csv.gz` (author DE results; logFC/logCPM/PValue/FDR with
  mixed IDs).

### A3. GSE224428 (pig PAB RV failure)
- `GSE224428_counts_symbol.tsv` — rows = gene symbols or ENSSSCG ids (30,761), cols = 8 samples
  (`JBM220318_Pig_01..04` = control, `JBM220318_Pig_13..16` = PAB). RAW COUNTS.
  ENSSSCG→symbol via mygene; unmapped rows keep ENSSSCG id. **C1QTNF9 not annotated in pig.**
- `GSE224428_metadata.tsv` — cols: `sample, group (control|PAB), species, model, tissue`.
- Raw: `GSE224428_subread_counts_RV.txt.gz` (featureCounts output, 8 count columns + annotation
  cols Chr/Start/End/Strand/Length), `GSE224428_DEgenes_RV_control-band.txt.gz` (author DESeq2
  results; first column header is a single quoted string containing tabs — parse carefully).

### A4. GSE36761 (human TOF RV, gene-level)
- `GSE36761_log2_matrix.tsv` — rows = gene symbols (19,662 incl. added `C1QTNF9`), cols = 30 samples
  (`NH-01..08`, `TOF-01..22`). Values = log2(provider normalized value + 1).
  ENSG→symbol via mygene (19,914 mapped); unmapped ENSG rows dropped.
- `GSE36761_C1QTNF9_transcript_sum.tsv` — 1 row (`C1QTNF9`) × 30 samples = RAW SUM of
  ENST00000382071 + ENST00000332018 (transcript-level values). **This is the only source of CTRP9
  gene-level data in this dataset** (gene-level row was filtered out by the provider).
- `GSE36761_metadata.tsv` — cols: `sample, group (NH|TOF), species, model, tissue`.
- Raw: `GSE36761_gene_expression_levels_normalized.txt.gz` (26,521 ENSG × 30),
  `GSE36761_transcript_expression_levels_normalized.txt.gz` (50,593 ENST × 30),
  `GSE36761_gene_expression_levels.txt.gz` (non-normalized).

### A5. GSE36759 (human TOF RV, mature miRNA)
- `GSE36759_mirna_log2_matrix.tsv` — rows = mature miRNA names (626), cols = same 30 samples as
  A4. Values = log2(provider normalized + 1).
- `GSE36759_metadata.tsv` — cols: `sample, group (NH|TOF), species, model`.

### A6. GSE217771 (human TOF RV, miRNA-seq)
- `GSE217771_mirna_counts.tsv` — rows = miRNA names (1,668), cols = 10 samples
  (`control1..5`, `treat1..5`). RAW COUNTS.
- `GSE217771_metadata.tsv` — cols: `sample, group (control|treat), species, model`.
- Raw: `GSE217771_miRNA_count.txt.gz`.

### A7. GSE219101 (mouse CTRP9 function)
- `GSE219101_WT_vs_KO.tsv` — rows = gene symbols (1,230 unique), cols incl. per-sample counts,
  `log2FoldChange`, `pvalue`, `padj` (KO vs WT, HFD hearts). Source: author xlsx
  `GSE219101_results_diff_allsignificant_WT_vs_KO.xlsx`. Note: only 2,372 rows in source (pre-filtered);
  many genes have padj=NaN.
- `GSE219101_AAV_CTRP9.tsv` — rows = mouse gene symbols (13,219) aggregated from transcript-level
  xlsx `GSE219101_result_diff2_AAV_CTRP9vsLuc.xlsx`; cols: `mgi_symbol, logFC (mean), PValue (min),
  FDR (min), n_transcripts`. AAV-CTRP9 vs AAV-Luc, HFD hearts.
- Key QC: `C1qtnf9` log2FC = +1.78 (FDR=0.001) in AAV; −9.47 in KO.

### A8. GSE57338 (human heart failure, probe-level — INCOMPLETE)
- `GSE57338_probe_matrix.tsv` — 33,298 probes × 313 samples, RMA-normalized (log2) values.
  Sample titles in `GSE57338_titles.tsv`, characteristics in `GSE57338_chars.tsv`
  (fields: tissue, heart failure yes/no, disease status DCM/ICM/normal, gender, age).
- `mapping_gse57338_probe_symbol.tsv` — **TRUNCATED**: only ~16,214 of 54,675 probes mapped
  (from a truncated GPL570 platform table). Use with caution; verify individual genes.
- **Action needed before use**: obtain complete probe→symbol map (BioMart / other mirror).

### A9. ID mappings
- `mapping_human_ensg_symbol.tsv` — ENSG → symbol (19,914 entries; mygene, 2026-08).
- `mapping_rat_ensrnog_symbol.tsv` — ENSRNOG → symbol (17,527).
- `mapping_pig_enssscg_symbol.tsv` — ENSSSCG → symbol (18,879).
- `mapping_gse57338_probe_symbol.tsv` — see A8 (truncated).

---

## B. Auxiliary data (`data/aux/`)

| File | Source | Format | Notes |
|---|---|---|---|
| `MCPcounter_genes.txt` | github.com/ebecht/MCPcounter | TSV, quoted fields | cols: `HUGO symbols | Cell population | ENTREZID | ENSEMBL ID`; 10 cell populations |
| `GTEx_v8_gene_median_tpm.gct.gz` | storage.googleapis.com/adult-gtex (v8 RNA-seq) | GCT | header `#1.2`, dims `56200 54`; cols: `Name (ENSG) | Description (symbol) | 54 tissues`; values = median TPM |
| `miRDB_v6_prediction.txt.gz` | mirdb.org | TSV | cols: `miRNA | RefSeq target | score`; all species, ~4.3M rows — **targets are RefSeq IDs, not symbols** |

---

## C. Raw GEO files (`data/raw/`)

All downloaded from NCBI GEO FTP, unmodified (see md5 manifest). Naming convention
`<GSE>_<description>.<ext>`. Notes:
- Series matrix files (`*_series_matrix.txt.gz`) are **metadata-only stubs (~4 KB)** for these
  datasets (expression data lives in supplementary files) — do not expect expression values inside.
- `GSE119754_all_counts.csv.gz` — rownames in first column (quoted), mixed MSTRG/ENSRNOG IDs.
- `GSE224428_subread_counts_RV.txt.gz` — featureCounts format: line 1 starts with `# Program:...`
  (skiprows=1 to parse), then header `Geneid Chr Start End Strand Length <samples>`.
- `GSE57338_series_matrix.txt.gz` — 54 MB, contains the RMA expression matrix (probe rows).
- `GSE219101_*.xlsx` — Excel DEG result tables (openpyxl to read).

---

## D. Frozen integrity

- `docs/data/md5_manifest.tsv` — `<md5>\t<path>` for all 45 files under `data/`.
  Regenerate: `find data -type f | sort | while read f; do md5 -q "$f" | awk -v file="$f" '{print $1"\t"file}'; done > docs/data/md5_manifest.tsv`
- After verification, treat `data/` as immutable. Any new analysis writes to `results/` or a new
  `analysis/` directory, not into `data/`.
