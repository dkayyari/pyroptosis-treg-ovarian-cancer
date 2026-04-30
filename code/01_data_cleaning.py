"""
Script 01: Data Cleaning
Project: Pyroptosis-Treg Axis in Ovarian Cancer
Author:  Deeksha Kayyari | Indiana University | April 2026
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import warnings
import os

warnings.filterwarnings("ignore")
np.random.seed(42)

# =============================================================
# CHANGE THIS PATH to your project folder
# =============================================================
PROJECT_ROOT = r"C:\Users\deeksha\Downloads\pyroptosis_project_old"
# =============================================================

DATA_RAW   = os.path.join(PROJECT_ROOT, "dataraw")
DATA_CLEAN = os.path.join(PROJECT_ROOT, "data_clean")
os.makedirs(DATA_CLEAN, exist_ok=True)

print("=" * 60)
print("SCRIPT 01 — DATA CLEANING")
print("=" * 60)

print("\n[1/4] Loading TCGA expression data...")
expr_path = os.path.join(DATA_RAW, "tcga_expr.tsv.gz")
tcga_raw  = pd.read_csv(expr_path, sep="\t", index_col=0, compression="gzip")
print(f"  Raw expression: {tcga_raw.shape}")

print("\n[2/4] Loading survival data...")
surv_path    = os.path.join(DATA_RAW, "TCGA-OV.survival.tsv")
survival_raw = pd.read_csv(surv_path, sep="\t", index_col=0)
print(f"  Raw survival: {survival_raw.shape}")
print(f"  Columns: {list(survival_raw.columns)}")

os_col      = [c for c in survival_raw.columns if "OS" in c and "time" not in c.lower()][0]
os_time_col = [c for c in survival_raw.columns if "time" in c.lower() and "OS" in c][0]
survival    = survival_raw[[os_col, os_time_col]].copy()
survival.columns = ["OS", "OS_time"]
survival["OS"]      = pd.to_numeric(survival["OS"],      errors="coerce")
survival["OS_time"] = pd.to_numeric(survival["OS_time"], errors="coerce")
survival = survival.dropna()
survival = survival[survival["OS_time"] > 0]
print(f"  Cleaned survival: {survival.shape}")
print(f"  Events: {int(survival['OS'].sum())} / {len(survival)}")

print("\n[3/4] Quality control...")
common_samples = list(set(tcga_raw.columns) & set(survival.index))
tcga = tcga_raw[common_samples].copy()
gene_var = tcga.var(axis=1)
tcga = tcga[gene_var > 0]
print(f"  Samples with both expression and survival: {len(common_samples)}")
print(f"  Genes after variance filter: {tcga.shape[0]}")

# Remove outliers by PCA
pca = PCA(n_components=2, random_state=42)
pca_coords = pca.fit_transform(tcga.T)
pc1_mean, pc1_std = pca_coords[:,0].mean(), pca_coords[:,0].std()
pc2_mean, pc2_std = pca_coords[:,1].mean(), pca_coords[:,1].std()
outlier_mask = (
    (np.abs(pca_coords[:,0] - pc1_mean) > 3*pc1_std) |
    (np.abs(pca_coords[:,1] - pc2_mean) > 3*pc2_std)
)
outlier_samples = np.array(common_samples)[outlier_mask].tolist()
good_samples    = [s for s in common_samples if s not in outlier_samples]
tcga            = tcga[good_samples]
survival        = survival.loc[good_samples]
print(f"  Outliers removed: {len(outlier_samples)}")
print(f"  Final samples: {tcga.shape[1]}")

print("\n[4/4] Saving...")
tcga.to_csv(os.path.join(DATA_CLEAN, "TCGA_expr_clean.csv"))
survival.to_csv(os.path.join(DATA_CLEAN, "TCGA_survival_clean.csv"))
print("  Saved: TCGA_expr_clean.csv")
print("  Saved: TCGA_survival_clean.csv")

print("\n" + "="*60)
print("SCRIPT 01 COMPLETE")
print(f"  Samples: {tcga.shape[1]}  Genes: {tcga.shape[0]}")
print(f"  Events: {int(survival['OS'].sum())} / {len(survival)}")
print("="*60)
