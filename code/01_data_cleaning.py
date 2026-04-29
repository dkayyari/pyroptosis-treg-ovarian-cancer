"""
=============================================================
Script 01: Data Download and Cleaning
Project: Pyroptosis-Treg Axis in Ovarian Cancer
Author:  Deeksha Kayyari
Date:    April 2026
=============================================================

Description:
    Downloads and cleans TCGA-OV gene expression, clinical,
    and survival data for downstream analysis.

Data Sources (download manually before running):
    1. TCGA-OV expression (log2 FPKM+1):
       https://xenabrowser.net/datapages/
       Dataset: TCGA-OV → gene expression RNAseq (IlluminaHiSeq)
       Save as: dataraw/tcga_expr.tsv.gz

    2. TCGA-OV survival:
       https://xenabrowser.net/datapages/
       Dataset: TCGA survival data
       Save as: dataraw/TCGA-OV.survival.tsv

    3. Gene probemap:
       https://xenabrowser.net/datapages/
       Save as: dataraw/gene_probemap.probemap

    4. GEO cohorts (run in terminal):
       wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE26nnn/GSE26193/matrix/GSE26193_series_matrix.txt.gz
       wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9891/matrix/GSE9891_series_matrix.txt.gz
       gunzip both files and place in dataraw/

    5. GDSC2 drug data:
       https://www.cancerrxgene.org/downloads/bulk_download
       Save as: dataraw/gdsc2.xlsx

    6. CCLE expression:
       https://depmap.org/portal/download/
       Save as: dataraw/ccle_expr.csv

Outputs:
    data_clean/TCGA_expr_clean.csv    — cleaned expression matrix
    data_clean/TCGA_survival_clean.csv — cleaned survival data
    data_clean/QC_plots.png           — quality control figures

Requirements:
    pip install pandas numpy matplotlib seaborn scikit-learn
=============================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import warnings
import os

warnings.filterwarnings("ignore")
np.random.seed(42)  # for reproducibility

# =============================================================
# CONFIGURATION — update paths for your system
# =============================================================
DATA_RAW   = r"C:\Users\deeksha\OneDrive - Indiana University\pyroptosis_treg_project\dataraw"
DATA_CLEAN = r"C:\Users\deeksha\OneDrive - Indiana University\pyroptosis_treg_project\data_clean"

# Create output directory if it does not exist
os.makedirs(DATA_CLEAN, exist_ok=True)

# =============================================================
# STEP 1 — Load TCGA expression data
# =============================================================
print("=" * 60)
print("SCRIPT 01 — DATA DOWNLOAD AND CLEANING")
print("=" * 60)

print("\n[1/4] Loading TCGA expression data...")

expr_path = os.path.join(DATA_RAW, "tcga_expr.tsv.gz")
if not os.path.exists(expr_path):
    raise FileNotFoundError(
        f"Expression file not found: {expr_path}\n"
        "Please download from UCSC Xena Browser (see script header)."
    )

tcga_raw = pd.read_csv(expr_path, sep="\t", index_col=0, compression="gzip")
print(f"  Raw expression matrix: {tcga_raw.shape}  (genes x samples)")
print(f"  Value range: {tcga_raw.values.min():.2f} — {tcga_raw.values.max():.2f}")

# =============================================================
# STEP 2 — Load survival data
# =============================================================
print("\n[2/4] Loading survival data...")

surv_path = os.path.join(DATA_RAW, "TCGA-OV.survival.tsv")
if not os.path.exists(surv_path):
    raise FileNotFoundError(f"Survival file not found: {surv_path}")

survival_raw = pd.read_csv(surv_path, sep="\t", index_col=0)
print(f"  Raw survival data: {survival_raw.shape}")
print(f"  Columns: {list(survival_raw.columns)}")

# Keep only OS and OS.time columns
# Rename to standard names
os_col      = [c for c in survival_raw.columns if "OS" in c and "time" not in c.lower()][0]
os_time_col = [c for c in survival_raw.columns if "time" in c.lower() and "OS" in c][0]

survival = survival_raw[[os_col, os_time_col]].copy()
survival.columns = ["OS", "OS_time"]
survival["OS"]      = pd.to_numeric(survival["OS"],      errors="coerce")
survival["OS_time"] = pd.to_numeric(survival["OS_time"], errors="coerce")
survival = survival.dropna()
survival = survival[survival["OS_time"] > 0]

print(f"  Cleaned survival: {survival.shape}")
print(f"  Events (OS=1): {survival['OS'].sum():.0f} / {len(survival)}")

# =============================================================
# STEP 3 — Quality control and cleaning
# =============================================================
print("\n[3/4] Quality control and cleaning...")

# 3a — Keep only samples with survival data
common_samples = list(set(tcga_raw.columns) & set(survival.index))
tcga = tcga_raw[common_samples].copy()
print(f"  Samples with expression AND survival: {len(common_samples)}")

# 3b — Remove zero-variance genes
gene_var = tcga.var(axis=1)
tcga     = tcga[gene_var > 0]
print(f"  Genes after variance filter: {tcga.shape[0]}")

# 3c — Remove outlier samples using PCA
pca  = PCA(n_components=2, random_state=42)
pca_coords = pca.fit_transform(tcga.T)

# Flag samples more than 3 standard deviations from mean on PC1 or PC2
pc1_mean, pc1_std = pca_coords[:, 0].mean(), pca_coords[:, 0].std()
pc2_mean, pc2_std = pca_coords[:, 1].mean(), pca_coords[:, 1].std()

outlier_mask = (
    (np.abs(pca_coords[:, 0] - pc1_mean) > 3 * pc1_std) |
    (np.abs(pca_coords[:, 1] - pc2_mean) > 3 * pc2_std)
)

outlier_samples = np.array(common_samples)[outlier_mask].tolist()
print(f"  Outlier samples removed (>3 SD on PCA): {len(outlier_samples)}")
print(f"  Outliers: {outlier_samples}")

good_samples = [s for s in common_samples if s not in outlier_samples]
tcga         = tcga[good_samples]
survival     = survival.loc[good_samples]

print(f"  Final expression matrix: {tcga.shape}  (genes x samples)")
print(f"  Final survival data:     {survival.shape}")

# 3d — Check all 33 pyroptosis-related genes are present
PRG_GENES = [
    "NLRP1","NLRP2","NLRP3","NLRP6","NLRP7","NLRC4","AIM2","PYCARD",
    "CASP1","CASP3","CASP4","CASP5","CASP6","CASP8","CASP9",
    "GSDMA","GSDMB","GSDMC","GSDMD","GSDME","PJVK",
    "IL1B","IL18","TNF","IL6","GPX4","NOD1","NOD2",
    "ELANE","PLCG1","TIRAP","PRKACA","SCAF11"
]
prgs_present = [g for g in PRG_GENES if g in tcga.index]
prgs_missing = [g for g in PRG_GENES if g not in tcga.index]
print(f"\n  PRGs present: {len(prgs_present)}/33")
if prgs_missing:
    print(f"  PRGs missing: {prgs_missing}")

# 3e — Check missing values
missing = tcga.isnull().sum().sum()
print(f"  Missing values: {missing}")

# =============================================================
# STEP 4 — Save cleaned data
# =============================================================
print("\n[4/4] Saving cleaned data...")

tcga.to_csv(os.path.join(DATA_CLEAN, "TCGA_expr_clean.csv"))
survival.to_csv(os.path.join(DATA_CLEAN, "TCGA_survival_clean.csv"))
print(f"  Saved: TCGA_expr_clean.csv   — {tcga.shape}")
print(f"  Saved: TCGA_survival_clean.csv — {survival.shape}")

# =============================================================
# QC PLOTS
# =============================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
fig.suptitle("Data Quality Control — TCGA-OV", fontsize=13, fontweight="bold")

# Plot 1 — Per-sample expression distribution
ax = axes[0]
sample_means = tcga.mean(axis=0)
ax.hist(sample_means, bins=40, color="#1D9E75", edgecolor="white", alpha=0.85)
ax.axvline(sample_means.mean(), color="red", linestyle="--",
           label=f"Mean={sample_means.mean():.2f}")
ax.set_title("Per-sample mean expression", fontsize=11)
ax.set_xlabel("Mean log2(FPKM+1)")
ax.set_ylabel("Number of samples")
ax.legend(fontsize=9)
sns.despine(ax=ax)

# Plot 2 — PRG expression levels
ax = axes[1]
prg_expr = tcga.loc[prgs_present].mean(axis=1).sort_values(ascending=False)
ax.barh(prg_expr.index[:10], prg_expr.values[:10], color="#534AB7", alpha=0.8)
ax.set_title("Top 10 PRG mean expression", fontsize=11)
ax.set_xlabel("Mean log2(FPKM+1)")
sns.despine(ax=ax)

# Plot 3 — Survival distribution
ax = axes[2]
ax.hist(survival["OS_time"][survival["OS"]==1]/365, bins=30,
        color="#D85A30", alpha=0.7, label=f"Deceased (n={int(survival['OS'].sum())})")
ax.hist(survival["OS_time"][survival["OS"]==0]/365, bins=30,
        color="#1D9E75", alpha=0.7, label=f"Censored (n={int((survival['OS']==0).sum())})")
ax.set_title("Survival time distribution", fontsize=11)
ax.set_xlabel("Time (years)")
ax.set_ylabel("Number of patients")
ax.legend(fontsize=9)
sns.despine(ax=ax)

plt.tight_layout()
plt.savefig(os.path.join(DATA_CLEAN, "QC_plots.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: QC_plots.png")

# =============================================================
# FINAL SUMMARY
# =============================================================
print("\n" + "=" * 60)
print("SCRIPT 01 COMPLETE — SUMMARY")
print("=" * 60)
print(f"  Final samples:      {tcga.shape[1]}")
print(f"  Final genes:        {tcga.shape[0]}")
print(f"  Survival events:    {int(survival['OS'].sum())} / {len(survival)}")
print(f"  Missing values:     {missing}")
print(f"  PRGs present:       {len(prgs_present)} / 33")
print(f"  Outliers removed:   {len(outlier_samples)}")
print(f"\n  Outputs: {DATA_CLEAN}")
print(f"  NEXT: Run 02_ssgsea_scoring.py")
print("=" * 60)
