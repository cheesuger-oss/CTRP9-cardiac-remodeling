# Figure QA Report — nature-figure skill audit (v1.0)

Audited with the **nature-figure** skill (Yuan1z0825/nature-skills, v2.5.0) installed at
`.dsh/skills/nature-figure/`. Backend: **python** (saved preference). Target journal family:
IJMS / Frontiers / JCDD-class (183 mm double-column; 600 dpi TIFF; editable vector text).
Flagship-Nature and NMI-specific rules were intentionally **not** applied (not the target).

---

## 1. Automated checks — all PASS

### 1.1 Source preflight (`validate_figure.py --strict`)

| Script | pass | warn | fail | verdict |
|---|---|---|---|---|
| `scripts/08_figures_core.py` (Fig2, Fig3, Fig3b) | 20 | 0 | 0 | READY |
| `scripts/09_figures_analyses.py` (Fig4–Fig7) | 20 | 0 | 0 | READY |
| `scripts/10_figures_extra.py` (Fig1, Fig8, FigS1) | 20 | 0 | 0 | READY |
| `scripts/11_figures_forest.py` (Fig2D) | 20 | 0 | 0 | READY |

Checks covered: syntax, font family (Arial/Helvetica), font-size floor, mathtext glyph floor,
legend-label case, color maps (no rainbow), editable text (`svg.fonttype='none'`,
`pdf.fonttype=42`), vector exports (SVG+PDF), raster exports (TIFF+PNG), DPI (600),
final width (89/183 mm), data sampling/exclusion, demo data, log guards (documented
pseudocount), interpolation monotonicity, uncertainty encoding, rotation anchor,
annotation workarounds, backend exclusivity.

### 1.2 Exported-PDF glyph audit (`audit_pdf_text.py --min-pt 5`) — all PASS

| PDF | min glyph (pt) | verdict |
|---|---|---|
| Fig1_workflow | 8.0 | PASS |
| Fig2D_forest_CTRP9 | 6.5 | PASS |
| Fig2_CTRP9_expression | 7.5 | PASS |
| Fig3_correlations_rat | 7.5 | PASS |
| Fig3b_correlations_human | 7.5 | PASS |
| Fig4_WGCNA | 6.0 | PASS |
| Fig5_immune | 7.5 | PASS |
| Fig6_ceRNA_network | 5.5 | PASS |
| Fig7_enrichment | 6.0 | PASS |
| Fig8_mouse_function | 5.5 | PASS |
| FigS1_GTEx_C1QTNF9_tissues | 6.0 | PASS |

### 1.3 Geometry / format checks

- **Page width**: 10/11 figures ≈ 178–183 mm (double column); Fig7 fixed to 182.7 mm after audit.
- **SVG text**: all 11 SVGs contain live `<text>` elements, **0 outlined glyphs** (editable text ✓).
- **TIFF resolution**: all 11 TIFFs at **600 dpi**.
- **Palette**: restrained (blue = control, red = disease, green = accent, purple = secondary,
  gray = neutral); no rainbow/jet/hsv; lightness separation keeps grayscale interpretability.

---

## 2. Figure contract (per figure)

| Fig | Archetype | Core conclusion it defends | Panels |
|---|---|---|---|
| 1 | schematic-led composite | Study design: 7 datasets, 4 species, 5 analysis layers | a–e workflow rows |
| 2 | quantitative grid | CTRP9 mRNA is low/down in rat RV pressure overload; unchanged (ns) in human TOF | a rat PAB; b rat MCT; c human TOF |
| 2D | quantitative grid (forest) | Direction & precision of CTRP9 effect across datasets | single panel |
| 3 | quantitative grid | CTRP9 inversely tracks hypertrophy/fibrosis markers in rat RV | a–h scatter (PAB ×4, MCT ×4) |
| 3b | quantitative grid | Same relation in human TOF (NPPB trend) | a–c scatter |
| 4 | quantitative grid | TOF RV has strong co-expression modules; ECM module dominates; C1QTNF9 unique (grey) | a module–trait heatmap; b module sizes |
| 5 | quantitative grid | TOF RV shows fibroblast/DC/CD8 infiltration up, monocytic down | a–f boxplots (6 cell types) |
| 6 | asymmetric mixed-modality (network) | DE miRNAs ↔ DE mRNAs form anti-correlated ceRNA network | single network |
| 7 | quantitative grid | RV failure downregulates OXPHOS; TOF downregulates inflammation; ME4 = ECM | a–c enrichment bars |
| 8 | quantitative grid | CTRP9 gain-of-function suppresses collagens; KO alters cardiac programs | a AAV-CTRP9; b CTRP9-KO |
| S1 | quantitative grid | CTRP9 is a low-abundance adipokine (adipose-high, heart-low) in GTEx | single horizontal bar |

---

## 3. Statistics-in-figure (legend minimum captured)

- **n definitions**: GSE210180 n=3+3 biological replicates (RV free wall); GSE119754 n=6+5;
  GSE36761 n=8+22; GSE224428 n=4+4; GSE219101 WT vs KO / AAV-Luc vs AAV-CTRP9.
- **Center/spread**: boxplots = median ± IQR + individual points; scatters = raw per-sample values.
- **Tests**: two-sided Welch t-test (expression panels), Pearson + Spearman (correlations),
  limma empirical-Bayes (DEG, BH-FDR), ssGSEA diff = Welch t + BH-FDR, miRNA = t-test + BH.
- **Correction**: Benjamini–Hochberg FDR where reported (FDR in panel titles).
- **p display**: exact p in titles (p = x.xxx) or FDR (FDR = x.xxx).
- **Source data**: every panel traces to `data/processed/*.tsv` + `results/*.tsv`
  (see `docs/DATA_DICTIONARY.md`).
- **Uncertainty consistency**: comparable panels use the same definition (Welch t + BH).

---

## 4. Panel-by-panel audit (rendered-level)

> Automated + code-level audit performed. **Final visual pass at print size is required by the
> human author** — the model cannot render/view pixels in this session; the items below are the
> programmatically verified and code-reviewed checks.

| Panel | Unique claim | Center | Spread | Replicate unit | Labels | Collision check | Pass |
|---|---|---|---|---|---|---|---|
| 2a | C1qtnf9 down in rat PAB (trend) | median box | IQR+points | 3 rats | direct x labels | points within box width; no mask | ✓ (code) |
| 2b | C1qtnf9 down in MCT-RVF (FDR) | median box | IQR+points | 5–6 rats | direct x labels | ✓ |
| 2c | C1QTNF9 ns in TOF | median box | IQR+points | 8+22 humans | direct x labels | ✓ |
| 3a–h | r/p for 8 marker pairs | raw scatter | none (all points) | per-sample | panel letters | trend lines below points (zorder) | ✓ |
| 4a | module eigengene vs TOF | Pearson r | p annot | 30 samples | module labels + p | heatmap annot ≤ cell | ✓ |
| 4b | module sizes | counts | none | — | value labels offset | ✓ |
| 5a–f | 6 cell types NH vs TOF | median box | IQR+points | 30 samples | FDR in title | ✓ |
| 6 | ceRNA network | — | — | 30 paired samples | node labels 5.5 pt (≥5 floor) | spring layout seed 42; recheck visually | ⚠ visual |
| 7a–c | -log10 adj p bars | term counts | none | — | term labels ≤48 chars | ✓ |
| 8a/b | marker log2FC bars | bar | FDR annot | mice | ✓ |
| S1 | median TPM by tissue | bar | none | GTEx median | tissue labels 6 pt | ✓ |

---

## 5. Export bundle

Per figure (`figures/`): `.svg` (editable) + `.pdf` (editable, fonttype 42) + `.tiff` (600 dpi)
+ `.png` (600 dpi preview). Source scripts: `scripts/08–11_figures_*.py` + `scripts/nature_style.py`.

---

## 6. Resolved issues during audit

1. Font family / editable text / SVG export missing → inlined quick-start rcParams + explicit
   SVG/PDF/TIFF/PNG exports in each script.
2. Widths 305–330 mm → 89/183 mm journal widths.
3. No panel letters → lowercase bold panel letters (a, b, c …) top-left.
4. `-log10(p)` log guard → documented pseudocount clip (1e-300 floor) + comment.
5. Rotation anchor on heatmap x tick labels → `rotation_mode="anchor"`.
6. Fig7 width 207.9 mm → 182.7 mm (shortened term labels + tighter figure).
7. seaborn boxplot palette deprecation → explicit `order`/`palette`.

## 7. Remaining human steps

- **Visual inspection at final print size** (required by QA contract; model cannot view raster).
- Confirm journal-specific figure specs (width/DPI) against the target journal's author guide.
- Confirm statistics text (n, test, correction) in figure legends matches the manuscript.

---

## 8. Fig. 2d revision: GSE57338 added (2026-08-20)

Fig. 2d was rebuilt from the underlying expression matrices to show five
case-minus-control contrasts: rat PAB RV, rat MCT-RVF, human TOF RV, human
ischemic HF LV, and human dilated-CM HF LV. Points are differences on each
dataset's log2 expression scale; error bars are two-sided Welch 95% confidence
intervals. RV and LV are encoded redundantly by color and marker shape. The two
GSE57338 contrasts share the same 136 non-failing LV reference samples.

- Visual review: **PASS** after moving statistical labels into a separate column;
  no label/CI collision, clipping, or ambiguous chamber encoding remained.
- Static preflight: **20 pass, 0 warn, 0 fail**.
- PDF glyph audit: **PASS**, minimum 6 pt.
- Geometry/export: exact 182.9 mm PDF page width; SVG contains 24 live text
  elements; TIFF and PNG are 4320 x 1830 px, TIFF 600 dpi.
- Source data: `analysis/Fig2D_CTRP9_forest_source.tsv`.
- Script: `scripts/11_figures_forest.py`.

**Nature-family review: PASS.** The figure has one explicit conclusion,
effect estimates with 95% CIs, group n, test/correction labels, a zero reference,
colorblind-safe redundant encoding, editable vector output, and final-size text.

**Light research-workflow review: PASS.** All five rows are recomputed from
frozen matrices and metadata; the source table records the comparison, n,
effect, CI, p-value, and inference method. The graphic supports context-dependent
expression but does not identify chamber or stage as the causal modifier.
