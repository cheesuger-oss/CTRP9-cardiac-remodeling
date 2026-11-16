#!/usr/bin/env python3
"""ceRNA network construction:
DE miRNAs (GSE36759) -> TargetScan 7.2 predicted targets -> DE mRNAs (GSE36761)
with anti-correlation directionality + expression correlation filter.
"""
import os, gzip
import numpy as np
import pandas as pd
from scipy import stats

OUT = "results"
os.makedirs(OUT, exist_ok=True)

# ---------- load TargetScan predictions (hsa only) ----------
import zipfile
zf = zipfile.ZipFile("/tmp/targetscan.zip")
with zf.open("Predicted_Targets_Context_Scores.default_predictions.txt") as f:
    ts = pd.read_csv(f, sep="\t")
print("TargetScan rows:", len(ts))
ts_hsa = ts[ts["miRNA"].str.startswith("hsa-")].copy()
print("hsa rows:", len(ts_hsa))

# ---------- DE miRNAs from GSE36759 ----------
dm = pd.read_csv("results/mirna_DEG_GSE36759.tsv", sep="\t")
de_mir = dm[dm["padj"] < 0.05]
print("DE miRNAs:", len(de_mir))
de_mir_up = set(de_mir[de_mir["log2FC"] > 0]["miRNA"])
de_mir_dn = set(de_mir[de_mir["log2FC"] < 0]["miRNA"])

# ---------- DE mRNAs from GSE36761 (limma) ----------
deg = pd.read_csv("results/limma_GSE36761_human.tsv", sep="\t")
deg_sig = deg[deg["adj.P.Val"] < 0.05]
deg_sym = set(deg_sig["gene"])
print("DE mRNAs:", len(deg_sym))

# ---------- predicted targets of DE miRNAs that are DE mRNAs ----------
pred = ts_hsa[ts_hsa["miRNA"].isin(de_mir["miRNA"]) & ts_hsa["Gene Symbol"].isin(deg_sym)].copy()
pred["context_rank"] = pred["weighted context++ score percentile"]
print("predicted DE miRNA -> DE mRNA pairs:", len(pred))

# ---------- expression anti-correlation filter ----------
mat = pd.read_csv("data/processed/GSE36761_log2_matrix.tsv", sep="\t", index_col=0)
mir_mat = pd.read_csv("data/processed/GSE36759_mirna_log2_matrix.tsv", sep="\t", index_col=0)
# same 30 samples
common_samples = [c for c in mat.columns if c in mir_mat.columns]
print("common samples:", len(common_samples))

edges = []
for _, row in pred.iterrows():
    mir = row["miRNA"]; gene = row["Gene Symbol"]
    if mir not in mir_mat.index or gene not in mat.index:
        continue
    m_exp = mir_mat.loc[mir, common_samples].values
    g_exp = mat.loc[gene, common_samples].values
    r, p = stats.pearsonr(m_exp, g_exp)
    if np.isnan(r): continue
    # directionality: if miRNA up in TOF, target should be down (negative r), and vice versa
    mir_dir = 1 if mir in de_mir_up else -1
    gene_dir = 1 if deg_sig[deg_sig["gene"]==gene]["logFC"].iloc[0] > 0 else -1
    anti = (mir_dir * gene_dir) < 0  # opposite regulation
    edges.append({"miRNA": mir, "target": gene, "mir_log2FC": dm[dm["miRNA"]==mir]["log2FC"].iloc[0],
                  "gene_log2FC": deg_sig[deg_sig["gene"]==gene]["logFC"].iloc[0],
                  "pearson_r": r, "pearson_p": p, "anti_regulation": anti,
                  "context_score_pct": row["context_rank"]})

edge_df = pd.DataFrame(edges)
edge_df = edge_df[(edge_df["anti_regulation"]) & (edge_df["pearson_p"] < 0.05) & (edge_df["pearson_r"] < 0)]
edge_df = edge_df.sort_values("pearson_p")
print("filtered ceRNA edges (anti-correlated, p<0.05):", len(edge_df))
edge_df.to_csv(f"{OUT}/ceRNA_edges.tsv", sep="\t", index=False)

# keep top edges for network figure (top 60 by |r|)
top = edge_df.head(60)
top.to_csv(f"{OUT}/ceRNA_edges_top60.tsv", sep="\t", index=False)
print("\ntop edges:")
print(top[["miRNA","target","pearson_r","pearson_p","mir_log2FC","gene_log2FC"]].head(20).to_string(index=False))

# does any DE miRNA target C1QTNF9?
c1q = ts_hsa[ts_hsa["Gene Symbol"]=="C1QTNF9"]
print("\nTargetScan predictions for C1QTNF9:", len(c1q))
if len(c1q) > 0:
    print(c1q[["miRNA","context score","weighted context++ score percentile"]].sort_values("weighted context++ score percentile", ascending=False).head(10).to_string(index=False))
# and for fibrosis markers
for g in ["COL1A1","COL3A1","POSTN","FN1","TGFB1","NPPB"]:
    sub = ts_hsa[ts_hsa["Gene Symbol"]==g]
    n_de = sub[sub["miRNA"].isin(de_mir["miRNA"])]
    if len(n_de) > 0:
        print(f"  {g}: {len(n_de)} DE-miRNA predicted targets: {', '.join(n_de['miRNA'].unique()[:8])}")
print("CERNA_DONE")
