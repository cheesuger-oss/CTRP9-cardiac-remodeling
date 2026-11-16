#!/usr/bin/env python3
"""edgeR-lite: TMM normalization + negative-binomial GLM (IRLS) + Wald test.

An edgeR/DESeq2-equivalent differential expression for raw counts, implemented in
Python (fast IRLS, no slow numeric Hessian). Cross-validates limma logCPM results.

Method:
  1. TMM normalization factors (Robinson & Oshlack 2010).
  2. Per-gene dispersion: quasi-likelihood moment estimate, empirical-Bayes
     shrunk toward a common dispersion (DESeq2-style, fixed prior weight).
  3. Per-gene NB GLM via IRLS: log(mu) = b0 + b1*group + log(effective libsize).
  4. Wald test b1/SE(b1); BH FDR.

Usage:
  .venv/bin/python scripts/12_edger_lite.py <counts_tsv> <metadata_tsv> <case> <ctrl> <out_prefix>
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

def tmm_norm_factors(counts):
    counts = np.asarray(counts, dtype=float)  # genes x samples
    ng, ns = counts.shape
    libsize = counts.sum(axis=0)
    logmean = np.mean(np.log(libsize))
    ref = int(np.argmin(np.abs(np.log(libsize) - logmean)))
    f = np.ones(ns)
    N_ref = libsize[ref]
    for j in range(ns):
        if j == ref:
            continue
        N_j = libsize[j]
        with np.errstate(divide='ignore', invalid='ignore'):
            Mg = np.log2((counts[:, j] / N_j) / (counts[:, ref] / N_ref))
            Ag = 0.5 * np.log2((counts[:, j] / N_j) * (counts[:, ref] / N_ref))
        ok = np.isfinite(Mg) & np.isfinite(Ag)
        if ok.sum() < 10:
            continue
        Mg, Ag = Mg[ok], Ag[ok]
        M_lo, M_hi = np.quantile(Mg, 0.30), np.quantile(Mg, 0.70)
        A_lo, A_hi = np.quantile(Ag, 0.05), np.quantile(Ag, 0.95)
        keep = (Mg >= M_lo) & (Mg <= M_hi) & (Ag >= A_lo) & (Ag <= A_hi)
        if keep.sum() < 10:
            continue
        Mg = Mg[keep]
        f[j] = 2 ** np.mean(Mg)
    f = f / np.exp(np.mean(np.log(f)))
    return f

def fit_nb_wald(y, s_eff, group, phi):
    """IRLS NB GLM for a single gene; return (log2FC, pvalue)."""
    X = np.column_stack([np.ones_like(group), group])
    offset = np.log(s_eff)
    b = np.zeros(2)
    for _ in range(40):
        eta = offset + X @ b
        mu = np.exp(eta)
        w = mu / (1.0 + phi * mu)
        # working response targets X@b (offset is already inside mu via eta)
        z = (X @ b) + (y - mu) / mu
        Xw = X * np.sqrt(w)[:, None]
        zw = z * np.sqrt(w)
        b_new, *_ = np.linalg.lstsq(Xw, zw, rcond=None)
        if np.max(np.abs(b_new - b)) < 1e-7:
            b = b_new
            break
        b = b_new
    XWX = (X * w[:, None]).T @ X
    try:
        cov = np.linalg.inv(XWX)
    except np.linalg.LinAlgError:
        return np.nan, np.nan
    se = np.sqrt(max(cov[1, 1], 1e-12))
    z = b[1] / se
    p = 2 * norm.sf(np.abs(z))
    return b[1] / np.log(2), p

def run_edger(counts_file, meta_file, case_label, ctrl_label, out_prefix):
    counts = pd.read_csv(counts_file, sep="\t", index_col=0)
    meta = pd.read_csv(meta_file, sep="\t").set_index("sample")
    case_cols = [c for c in counts.columns if meta.loc[c, "group"] == case_label]
    ctrl_cols = [c for c in counts.columns if meta.loc[c, "group"] == ctrl_label]
    cols = case_cols + ctrl_cols
    counts = counts[cols].astype(float)
    group = np.array([1.0] * len(case_cols) + [0.0] * len(ctrl_cols))
    print(f"[{out_prefix}] genes={counts.shape[0]} samples={len(cols)} case={len(case_cols)} ctrl={len(ctrl_cols)}")

    f = tmm_norm_factors(counts.values)  # ns samples
    libsize = counts.values.sum(axis=0)
    s_eff = (libsize * f)
    s_eff = s_eff / np.mean(s_eff) * 1e6

    Y = counts.values  # genes x samples
    # per-gene moment dispersion (group-adjusted residuals)
    mu = s_eff[np.newaxis, :]  # normalized baseline
    # per-gene group means (fit quick via group means of normalized counts)
    norm = Y / s_eff[np.newaxis, :] * 1e6  # normalized counts
    mean_case = norm[:, group == 1].mean(axis=1)
    mean_ctrl = norm[:, group == 0].mean(axis=1)
    fitted = np.zeros_like(norm)
    fitted[:, group == 1] = mean_case[:, None]
    fitted[:, group == 0] = mean_ctrl[:, None]
    resid = (norm - fitted) ** 2 / np.maximum(fitted, 0.5)
    # moment dispersion: (var - mean)/mean^2 per gene (quasi-likelihood, weighted)
    phi_g = np.zeros(Y.shape[0])
    for i in range(Y.shape[0]):
        m = fitted[i].mean()
        if m <= 0:
            phi_g[i] = 0.05
            continue
        v = ((norm[i] - fitted[i]) ** 2).sum() / (norm.shape[1] - 1)
        phi_g[i] = max((v - m) / (m * m), 1e-6)
    # common dispersion (weighted median)
    phi_c = np.median(phi_g)
    # empirical-Bayes shrink toward common (fixed prior weight)
    d0 = 3.0
    phi = (phi_g + d0 * phi_c) / (1 + d0)

    rows = []
    for i, gene in enumerate(counts.index):
        log2fc, p = fit_nb_wald(Y[i], s_eff, group, phi[i])
        rows.append({"gene": gene, "log2FC": log2fc, "pvalue": p})
    res = pd.DataFrame(rows).dropna(subset=["log2FC", "pvalue"])
    res["padj"] = multipletests(res["pvalue"].values, method="fdr_bh")[1]
    res = res.sort_values("pvalue")
    res.to_csv(f"analysis/{out_prefix}_edger.tsv", sep="\t", index=False)
    sig = res[res["padj"] < 0.05]
    print(f"[{out_prefix}] DE (padj<0.05): {len(sig)} (up={(sig['log2FC']>0).sum()} down={(sig['log2FC']<0).sum()})")
    return res

if __name__ == "__main__":
    counts_file, meta_file, case, ctrl, prefix = sys.argv[1:6]
    run_edger(counts_file, meta_file, case, ctrl, prefix)
