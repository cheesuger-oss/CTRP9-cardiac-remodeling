#!/usr/bin/env python3
"""Parse all GEO datasets into clean processed matrices + sample metadata.
Outputs to data/processed/.
"""
import os, gzip, json
import numpy as np
import pandas as pd

OUT = "data/processed"
os.makedirs(OUT, exist_ok=True)

def save(df, name):
    df.to_csv(f"{OUT}/{name}", sep="\t")
    print(f"saved {name}: {df.shape}")

# ============ 1. GSE210180: neonatal rat PAB RV (FPKM) ============
print("### GSE210180 (rat PAB RV)")
fpkm = pd.read_csv("data/raw/GSE210180_gene_fpkm.txt.gz", sep="\t", compression="gzip")
sample_cols = ["P21_Sham_1", "P21_Sham_2", "P21_Sham_3", "P21_PAB_1", "P21_PAB_2", "P21_PAB_3"]
mat = fpkm[sample_cols].copy()
mat["symbol"] = fpkm["gene_name"].astype(str)
# collapse duplicate symbols by max mean
mat = mat.groupby("symbol").mean(numeric_only=True).reset_index()
mat = mat.set_index("symbol").sort_index()
mat_log = np.log2(mat + 1)
save(mat_log, "GSE210180_log2fpkm_matrix.tsv")
meta = pd.DataFrame({
    "sample": sample_cols,
    "group": ["sham"] * 3 + ["pab"] * 3,
    "species": "rat", "model": "neonatal PAB RV", "tissue": "RV free wall"
})
save(meta, "GSE210180_metadata.tsv")

# ============ 2. GSE36761: human TOF RV (gene-level + transcript CTRP9) ============
print("### GSE36761 (human TOF RV)")
expr = pd.read_csv("data/raw/GSE36761_gene_expression_levels_normalized.txt.gz",
                   sep="\t", compression="gzip")
hm = pd.read_csv("data/processed/mapping_human_ensg_symbol.tsv", sep="\t", index_col=0)["symbol"].to_dict()
expr["symbol"] = expr["Ensembl Gene ID"].map(hm)
expr = expr.dropna(subset=["symbol"]).drop_duplicates(subset=["symbol"])
sample_cols = [c for c in expr.columns if c not in ("Ensembl Gene ID", "symbol")]
mat = expr.set_index("symbol")[sample_cols].sort_index()
# values appear to be raw-ish normalized counts (0..1350); log2 transform
mat_log = np.log2(mat + 1)
save(mat_log, "GSE36761_log2_matrix.tsv")
groups = ["NH" if c.startswith("NH") else "TOF" for c in sample_cols]
meta = pd.DataFrame({"sample": sample_cols, "group": groups,
                     "species": "human", "model": "TOF RV", "tissue": "RV myocardium"})
save(meta, "GSE36761_metadata.tsv")

# CTRP9 reconstructed from transcripts
tr = pd.read_csv("data/raw/GSE36761_transcript_expression_levels_normalized.txt.gz",
                 sep="\t", compression="gzip")
ctrp9_tr = tr[tr[tr.columns[0]].isin(["ENST00000382071", "ENST00000332018"])]
print("CTRP9 transcripts rows:", len(ctrp9_tr))
if len(ctrp9_tr) > 0:
    ctrp9_gene = ctrp9_tr[sample_cols].sum(axis=0)
    ctrp9_gene.name = "C1QTNF9"
    save(pd.DataFrame(ctrp9_gene).T, "GSE36761_C1QTNF9_transcript_sum.tsv")

# ============ 3. GSE36759: TOF mature miRNA ============
print("### GSE36759 (TOF miRNA)")
mir = pd.read_csv("data/raw/GSE36759_mature_miRNA_expression_levels_normalized.txt.gz",
                  sep="\t", compression="gzip")
mir = mir.set_index("miRNA").sort_index()
mir_log = np.log2(mir + 1)
save(mir_log, "GSE36759_mirna_log2_matrix.tsv")
mir_meta = pd.DataFrame({"sample": mir.columns,
                         "group": ["NH" if c.startswith("NH") else "TOF" for c in mir.columns],
                         "species": "human", "model": "TOF RV miRNA"})
save(mir_meta, "GSE36759_metadata.tsv")

# ============ 4. GSE217771: TOF miRNA-seq counts ============
print("### GSE217771 (TOF miRNA-seq)")
mir2 = pd.read_csv("data/raw/GSE217771_miRNA_count.txt.gz", sep="\t", compression="gzip")
mir2 = mir2.set_index("miRNA_name").sort_index()
save(mir2, "GSE217771_mirna_counts.tsv")
mir2_meta = pd.DataFrame({"sample": mir2.columns,
                          "group": ["control"] * 5 + ["treat"] * 5,
                          "species": "human", "model": "TOF miRNA-seq"})
save(mir2_meta, "GSE217771_metadata.tsv")

# ============ 5. GSE119754: rat MCT-PAH RV failure (counts, mixed IDs) ============
print("### GSE119754 (rat MCT-PAH RVF)")
cnt = pd.read_csv("data/raw/GSE119754_all_counts.csv.gz", compression="gzip", index_col=0)
rm = pd.read_csv("data/processed/mapping_rat_ensrnog_symbol.tsv", sep="\t", index_col=0)["symbol"].to_dict()
cnt["symbol"] = [rm.get(str(i), str(i)) for i in cnt.index]
# collapse by symbol (sum counts)
cnt = cnt.groupby("symbol").sum()
save(cnt, "GSE119754_counts_symbol.tsv")
mct_meta = pd.DataFrame({
    "sample": cnt.columns,
    "group": ["control" if c.startswith("CT") else "MCT" for c in cnt.columns],
    "species": "rat", "model": "MCT-PAH RV failure", "tissue": "RV"
})
save(mct_meta, "GSE119754_metadata.tsv")
print("C1qtnf9 in GSE119754:", "C1qtnf9" in cnt.index)

# ============ 6. GSE224428: pig PAB RV (counts) ============
print("### GSE224428 (pig PAB RV)")
pig = pd.read_csv("data/raw/GSE224428_subread_counts_RV.txt.gz", sep="\t",
                  compression="gzip", skiprows=1)
pm = pd.read_csv("data/processed/mapping_pig_enssscg_symbol.tsv", sep="\t", index_col=0)["symbol"].to_dict()
pig_cnt_cols = [c for c in pig.columns if c.startswith("JBM")]
pig_mat = pig[pig_cnt_cols].copy()
pig_mat["symbol"] = pig["Geneid"].map(pm).fillna(pig["Geneid"])
pig_mat = pig_mat.groupby("symbol").sum()
# sample order in count file: JBM220318_Pig_01..16; series matrix order: CDS30,31,32,33,15,18,19,20
# map by author DE file naming (control CDS30-33 -> Pig_01-04?); use positional mapping validated by DE file
pig_meta = pd.DataFrame({
    "sample": pig_cnt_cols,
    "group": ["control"] * 4 + ["PAB"] * 4,
    "species": "pig", "model": "PAB RV failure", "tissue": "RV free wall"
})
save(pig_mat, "GSE224428_counts_symbol.tsv")
save(pig_meta, "GSE224428_metadata.tsv")
print("C1QTNF9 in pig counts:", "C1QTNF9" in pig_mat.index)

# ============ 7. GSE57338: human DCM/ICM/control (RMA) ============
print("### GSE57338 (human heart failure)")
if os.path.exists("data/raw/GSE57338_series_matrix.txt.gz"):
    with gzip.open("data/raw/GSE57338_series_matrix.txt.gz", "rt") as f:
        lines = [l.rstrip("\n") for l in f]
    # parse sample characteristics
    chars = [l for l in lines if l.startswith("!Sample_characteristics_ch1")]
    titles = [l for l in lines if l.startswith("!Sample_title")]
    print("GSE57338 sample titles:", titles[0][:200] if titles else "none")
    print("GSE57338 char lines:", len(chars))
    for c in chars[:5]:
        print("  ", c[:250])
else:
    print("GSE57338 not yet downloaded")

print("\nALL PARSE DONE")
