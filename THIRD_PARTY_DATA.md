# Third-party data and retrieval record

The repository includes only shareable processed inputs and derived outputs needed for transparent review. Inclusion does not transfer ownership of upstream public datasets.

| Resource | Version/accession | Use | Retrieval |
|---|---|---|---|
| NCBI GEO | GSE210180 | Rat PAB RV transcriptome | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210180 |
| NCBI GEO | GSE119754 | Rat MCT-associated RV failure | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119754 |
| NCBI GEO | GSE224428 | Pig RV pressure overload | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224428 |
| NCBI GEO | GSE36761/GSE36759 | Human TOF cardiac mRNA/miRNA | https://www.ncbi.nlm.nih.gov/geo/ |
| NCBI GEO | GSE217771 | Supporting miRNA dataset | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217771 |
| NCBI GEO | GSE219101 | Public CTRP9 gain/loss summaries | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219101 |
| NCBI GEO | GSE57338 | Human terminal LV failure | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57338 |
| NCBI GEO | GSE345646 | Adult human RV single-nucleus data | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345646 |
| GTEx | v8 | Tissue median TPM | https://gtexportal.org/home/datasets |
| TargetScan Human | 7.2 default conserved predictions | Candidate miRNA-target pairs | https://www.targetscan.org/vert_72/ |

The TargetScan archive is not redistributed in this release. Place the official download at `analysis/TargetScanHuman_7.2_default_context_predictions.zip` before rerunning scripts 17 and 18. The frozen project copy had SHA-256 `24f9c6dc83a75da5b28839c26d3e77c37040f4af3ffac78d7d033d86b86faa80`; verify the downloaded file before use.

Large raw and H5AD files are intentionally omitted. Their derived donor-level and cohort-level outputs are included, and the reconstruction scripts preserve the public accession and processing logic.

