"""
=============================================================
Script 02: ssGSEA Pathway Scoring
Project: Pyroptosis-Treg Axis in Ovarian Cancer
Author:  Deeksha Kayyari
Date:    April 2026
=============================================================

Description:
    Calculates single-sample Gene Set Enrichment Analysis (ssGSEA)
    scores for three gene sets per patient:
      A) 33 Pyroptosis-Related Genes (PRGs) → Pyroptosis score
      B) 13 Treg marker genes               → Treg score
      C) 10 CCL22-axis mechanism genes      → CCL22 score

    Also runs permutation test to confirm correlation is real.

Inputs:
    data_clean/TCGA_expr_clean.csv
    data_clean/TCGA_survival_clean.csv

Outputs:
    data_clean/ssgsea_scores.csv           — scores for all patients
    data_clean/Phase3_Figure1_scatter.png  — Figure 1 (main finding)
    data_clean/Phase3_distributions.png    — score distributions
    data_clean/Phase3_groups.png           — 4-group pie chart
    data_clean/Proof1_permutation_test.png — permutation test plot

Requirements:
    pip install pandas numpy matplotlib seaborn gseapy scipy
=============================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from scipy import stats
import warnings
import os

warnings.filterwarnings("ignore")
np.random.seed(42)  # for reproducibility

# =============================================================
# CONFIGURATION — update paths for your system
# =============================================================
DATA_CLEAN = r"C:\Users\deeksha\OneDrive - Indiana University\pyroptosis_treg_project\data_clean"

# =============================================================
# GENE SETS
# =============================================================
# Gene Set A — 33 Pyroptosis-Related Genes (PRGs)
# Source: Ye et al. (2021), Cao et al. (2022)
PYROPTOSIS_GENES = [
    "NLRP1","NLRP2","NLRP3","NLRP6","NLRP7","NLRC4","AIM2","PYCARD",
    "CASP1","CASP3","CASP4","CASP5","CASP6","CASP8","CASP9",
    "GSDMA","GSDMB","GSDMC","GSDMD","GSDME","PJVK",
    "IL1B","IL18","TNF","IL6","GPX4","NOD1","NOD2",
    "ELANE","PLCG1","TIRAP","PRKACA","SCAF11"
]

# Gene Set B — 13 Treg Marker Genes
# Key marker: FOXP3 (master transcription factor of Tregs)
TREG_GENES = [
    "FOXP3","IL2RA","CTLA4","IKZF2","TIGIT","CCR8",
    "TNFRSF18","LAYN","BATF","ENTPD1","CXCR5","TGFB1","IL10"
]

# Gene Set C — CCL22/CCR8 Axis Mechanism Genes
# Tests the mechanistic link: pyroptosis → CCL22 → Treg recruitment
CCL22_GENES = [
    "GSDMD","CASP1","IL1B","IL18","CCL22",
    "CCL17","CCL28","IL33","CCR8","IL1RL1"
]

# Confirm no overlap between gene sets (important for validity)
overlap = set(PYROPTOSIS_GENES) & set(TREG_GENES)
assert len(overlap) == 0, f"Gene sets overlap: {overlap}"
print(f"Gene set overlap check: PASSED (no overlap)")

# =============================================================
# STEP 1 — Load cleaned data
# =============================================================
print("=" * 60)
print("SCRIPT 02 — ssGSEA SCORING")
print("=" * 60)

print("\n[1/5] Loading cleaned data...")

tcga = pd.read_csv(
    os.path.join(DATA_CLEAN, "TCGA_expr_clean.csv"),
    index_col=0
)
survival = pd.read_csv(
    os.path.join(DATA_CLEAN, "TCGA_survival_clean.csv"),
    index_col=0
)
print(f"  Expression: {tcga.shape}  (genes x samples)")
print(f"  Survival:   {survival.shape}")

# =============================================================
# STEP 2 — Check gene coverage
# =============================================================
print("\n[2/5] Checking gene coverage...")
for name, genes in [("Pyroptosis", PYROPTOSIS_GENES),
                    ("Treg",       TREG_GENES),
                    ("CCL22_axis", CCL22_GENES)]:
    present = [g for g in genes if g in tcga.index]
    missing = [g for g in genes if g not in tcga.index]
    print(f"  {name}: {len(present)}/{len(genes)} genes found", end="")
    if missing:
        print(f"  (missing: {missing})")
    else:
        print(" ✓")

# =============================================================
# STEP 3 — Run ssGSEA separately for each gene set
# Note: run separately to avoid gseapy multi-set bug
# =============================================================
print("\n[3/5] Running ssGSEA (one gene set at a time)...")
print("  This takes approximately 20-30 minutes — please wait.")

def run_ssgsea(expr_df, gene_list, name):
    """
    Run ssGSEA for a single gene set.

    Parameters:
        expr_df   : DataFrame, genes x samples
        gene_list : list of gene symbols
        name      : str, name of the gene set

    Returns:
        pd.Series of NES scores per sample
    """
    # Keep only genes present in the expression matrix
    genes = [g for g in gene_list if g in expr_df.index]
    print(f"  {name}: {len(genes)} genes found — running ssGSEA...")

    result = gp.ssgsea(
        data               = expr_df,
        gene_sets          = {name: genes},
        outdir             = None,
        sample_norm_method = "rank",
        no_plot            = True,
        processes          = 1,
        min_size           = 1        # allow small gene sets (e.g., 10-13 genes)
    )

    scores = result.res2d[result.res2d["Term"] == name].set_index("Name")["NES"]
    scores.name = name
    print(f"  {name}: done — {len(scores)} patients scored")
    return scores

pyro_scores  = run_ssgsea(tcga, PYROPTOSIS_GENES, "Pyroptosis_score")
treg_scores  = run_ssgsea(tcga, TREG_GENES,       "Treg_score")
ccl22_scores = run_ssgsea(tcga, CCL22_GENES,      "CCL22_score")

# =============================================================
# STEP 4 — Combine scores and merge with survival
# =============================================================
print("\n[4/5] Combining scores and merging with survival data...")

scores_df = pd.concat([pyro_scores, treg_scores, ccl22_scores], axis=1)
print(f"  Combined scores: {scores_df.shape}")

# Match samples between scores and survival
common = list(set(scores_df.index) & set(survival.index))
master_df = pd.concat([scores_df.loc[common], survival.loc[common]], axis=1)
master_df.index.name = "sample"

# Convert to numeric — critical step to avoid type errors
for col in ["Pyroptosis_score", "Treg_score", "CCL22_score"]:
    master_df[col] = pd.to_numeric(master_df[col], errors="coerce")
master_df = master_df.dropna(subset=["Pyroptosis_score", "Treg_score"])

print(f"  Patients with scores and survival: {len(master_df)}")

# Score summary statistics
print(f"\n  Score summary:")
print(master_df[["Pyroptosis_score","Treg_score","CCL22_score"]].describe().round(3))

# Stratify patients into 4 groups using median split
pyro_median = master_df["Pyroptosis_score"].median()
treg_median = master_df["Treg_score"].median()

master_df["Pyro_group"] = master_df["Pyroptosis_score"].apply(
    lambda x: "High_Pyro" if x > pyro_median else "Low_Pyro"
)
master_df["Treg_group"] = master_df["Treg_score"].apply(
    lambda x: "High_Treg" if x > treg_median else "Low_Treg"
)
master_df["Combined_group"] = master_df["Pyro_group"] + "_" + master_df["Treg_group"]

print(f"\n  Group distribution (median split):")
print(master_df["Combined_group"].value_counts())

# KEY RESULT — Pearson correlation
r, p = stats.pearsonr(
    master_df["Pyroptosis_score"],
    master_df["Treg_score"]
)
print(f"\n{'='*50}")
print(f"  KEY RESULT — Pyroptosis vs Treg correlation")
print(f"  Pearson r = {r:.3f}")
print(f"  p-value   = {p:.4f}")
print(f"  Direction = {'POSITIVE' if r > 0 else 'NEGATIVE'}")
print(f"  Significant: {'YES ✓' if p < 0.05 else 'NO ✗'}")
print(f"{'='*50}")

# Validate Treg score using FOXP3 single-gene expression
if "FOXP3" in tcga.index:
    foxp3 = pd.to_numeric(tcga.loc["FOXP3", master_df.index], errors="coerce")
    r_fox, p_fox = stats.pearsonr(foxp3, master_df["Treg_score"])
    print(f"\n  Validation: FOXP3 vs Treg score: r = {r_fox:.3f}  p = {p_fox:.4f}")
    print(f"  {'✓ Treg score confirmed valid' if r_fox > 0.3 else '⚠ Check Treg score'}")

# Individual gene correlations (mechanistic validation)
print(f"\n  Individual gene correlations (no ssGSEA):")
gene_pairs = [("GSDMD","FOXP3"), ("GSDMD","CCL22"),
              ("CASP1","FOXP3"), ("CCL22","FOXP3")]
for g1, g2 in gene_pairs:
    if g1 in tcga.index and g2 in tcga.index:
        v1 = tcga.loc[g1].astype(float)
        v2 = tcga.loc[g2].astype(float)
        ri, pi = stats.pearsonr(v1, v2)
        print(f"  {g1} vs {g2}: r = {ri:.3f}  p = {pi:.4f}")

# =============================================================
# PERMUTATION TEST — confirm correlation is not a statistical artifact
# =============================================================
print(f"\n  Running permutation test (n=1,000)...")

n_permutations = 1000
null_r_values  = []
np.random.seed(42)

for i in range(n_permutations):
    shuffled = master_df["Treg_score"].sample(frac=1).values
    r_null, _ = stats.pearsonr(master_df["Pyroptosis_score"], shuffled)
    null_r_values.append(r_null)

null_r_values = np.array(null_r_values)
permutation_p = np.mean(np.abs(null_r_values) >= np.abs(r))

print(f"  Real r:             {r:.3f}")
print(f"  Null max r:         {null_r_values.max():.3f}")
print(f"  Null mean r:        {null_r_values.mean():.3f}")
print(f"  Permutation p:      {permutation_p:.4f}")
if permutation_p == 0:
    print(f"  → CONFIRMED: 0/{n_permutations} permutations reached r = {r:.3f}")
    print(f"  → Correlation is REAL biology, not a statistical artifact")

# =============================================================
# STEP 5 — Save outputs and generate figures
# =============================================================
print("\n[5/5] Saving outputs and generating figures...")

# Save master scores file
master_df.to_csv(os.path.join(DATA_CLEAN, "ssgsea_scores.csv"))
print("  Saved: ssgsea_scores.csv")

# --- Figure 1: Pyroptosis vs Treg scatter plot ---
fig, ax = plt.subplots(figsize=(7, 6))

color_map = {
    "High_Pyro_High_Treg": "#1D9E75",
    "High_Pyro_Low_Treg" : "#D85A30",
    "Low_Pyro_High_Treg" : "#534AB7",
    "Low_Pyro_Low_Treg"  : "#888780"
}

for group, color in color_map.items():
    sub = master_df[master_df["Combined_group"] == group]
    ax.scatter(sub["Pyroptosis_score"], sub["Treg_score"],
               c=color, alpha=0.65, s=30,
               label=f"{group} (n={len(sub)})")

# Regression line
m, b = np.polyfit(master_df["Pyroptosis_score"].astype(float),
                  master_df["Treg_score"].astype(float), 1)
x_line = np.linspace(master_df["Pyroptosis_score"].min(),
                     master_df["Pyroptosis_score"].max(), 100)
ax.plot(x_line, m * x_line + b, "k--", linewidth=1.5, alpha=0.7)

# Median split lines
ax.axvline(pyro_median, color="gray", linestyle=":", alpha=0.4)
ax.axhline(treg_median, color="gray", linestyle=":", alpha=0.4)

ax.set_xlabel("Pyroptosis ssGSEA score", fontsize=12)
ax.set_ylabel("Treg ssGSEA score", fontsize=12)
ax.set_title(
    f"Pyroptosis vs Treg score — TCGA-OV\n"
    f"Pearson r = {r:.3f}   p = {p:.4f}",
    fontsize=12
)
ax.legend(fontsize=8, loc="upper left")
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(os.path.join(DATA_CLEAN, "Phase3_Figure1_scatter.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: Phase3_Figure1_scatter.png  ← Figure 1 (main finding)")

# --- Score distributions ---
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
fig.suptitle(f"ssGSEA Score Distributions — TCGA-OV (n={len(master_df)})",
             fontsize=13, fontweight="bold")

for ax, col, color, title in zip(
    axes,
    ["Pyroptosis_score", "Treg_score", "CCL22_score"],
    ["#1D9E75", "#534AB7", "#BA7517"],
    ["Pyroptosis score (33 PRGs)",
     "Treg score (13 markers)",
     "CCL22 axis score (10 genes)"]
):
    ax.hist(master_df[col].dropna(), bins=40, color=color,
            edgecolor="white", linewidth=0.5, alpha=0.85)
    ax.axvline(master_df[col].median(), color="red", linestyle="--",
               linewidth=1.5, label=f"Median={master_df[col].median():.2f}")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("NES score")
    ax.set_ylabel("Number of patients")
    ax.legend(fontsize=9)
    sns.despine(ax=ax)

plt.tight_layout()
plt.savefig(os.path.join(DATA_CLEAN, "Phase3_distributions.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: Phase3_distributions.png")

# --- Permutation test plot ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(null_r_values, bins=50, color="#B5D4F4",
        edgecolor="white", label=f"Null distribution (n={n_permutations} permutations)")
ax.axvline(r,  color="#D85A30", linewidth=2.5, label=f"Real r = {r:.3f}")
ax.axvline(-r, color="#D85A30", linewidth=2.5, linestyle="--", alpha=0.5)
ax.set_xlabel("Pearson r", fontsize=12)
ax.set_ylabel("Frequency", fontsize=12)
ax.set_title(
    f"Permutation test — Pyroptosis vs Treg\n"
    f"Permutation p = {permutation_p:.4f}  "
    f"({'REAL biology confirmed' if permutation_p == 0 else 'check results'})",
    fontsize=12
)
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(DATA_CLEAN, "Proof1_permutation_test.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: Proof1_permutation_test.png")

# =============================================================
# FINAL SUMMARY
# =============================================================
print("\n" + "=" * 60)
print("SCRIPT 02 COMPLETE — SUMMARY")
print("=" * 60)
print(f"  Patients scored:      {len(master_df)}")
print(f"  Pyroptosis range:     {master_df['Pyroptosis_score'].min():.3f} — {master_df['Pyroptosis_score'].max():.3f}")
print(f"  Treg range:           {master_df['Treg_score'].min():.3f} — {master_df['Treg_score'].max():.3f}")
print(f"\n  KEY RESULT:")
print(f"  Pyroptosis vs Treg:   r = {r:.3f}   p = {p:.4f}")
print(f"  Permutation p:        {permutation_p:.4f}")
print(f"\n  Group sizes:")
for g, n in master_df["Combined_group"].value_counts().items():
    print(f"    {g}: {n} patients")
print(f"\n  Outputs saved to: {DATA_CLEAN}")
print(f"  NEXT: Run 03_core_analysis.py")
print("=" * 60)
