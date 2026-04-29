"""
=============================================================
Script 04: Drug Sensitivity Analysis
Project: Pyroptosis-Treg Axis in Ovarian Cancer
Author:  Deeksha Kayyari
Date:    April 2026
=============================================================

Description:
    Compares IC50 drug sensitivity between OC cell lines
    grouped by pyroptosis activity score using CCLE and GDSC2.
    This is an exploratory analysis (n=7 cell lines).

Inputs:
    dataraw/ccle_expr.csv   — CCLE expression (DepMap portal)
    dataraw/gdsc2.xlsx      — GDSC2 drug IC50 data

Outputs:
    data_clean/Figure5_drug_sensitivity.png
    data_clean/drug_sensitivity_results.csv

Requirements:
    pip install pandas numpy matplotlib seaborn scipy openpyxl
=============================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import warnings
import os

warnings.filterwarnings("ignore")
np.random.seed(42)

# =============================================================
# CONFIGURATION — update paths for your system
# =============================================================
DATA_RAW   = r"C:\Users\deeksha\OneDrive - Indiana University\pyroptosis_treg_project\dataraw"
DATA_CLEAN = r"C:\Users\deeksha\OneDrive - Indiana University\pyroptosis_treg_project\data_clean"

OC_DRUGS = ["Cisplatin","Paclitaxel","Olaparib","Gemcitabine"]

# Curated OC cell line mapping: CCLE ACH ID -> common name
KNOWN_OV = {
    "ACH-000672": "SKOV3",
    "ACH-000227": "OVCAR3",
    "ACH-000018": "A2780",
    "ACH-001209": "Caov3",
    "ACH-001028": "OVCAR4",
    "ACH-001029": "OVCAR8",
    "ACH-000925": "IGROV1"
}

# GDSC2 uses different cell line naming conventions
GDSC_NAMES = {
    "SKOV3":"SK-OV-3","OVCAR3":"OVCAR-3","A2780":"A2780",
    "Caov3":"Caov-3","OVCAR4":"OVCAR-4","OVCAR8":"OVCAR-8","IGROV1":"IGROV-1"
}

# Pyroptosis proxy genes for cell line scoring
PYRO_PROXY_GENES = ["GSDMD","GSDME","CASP1","NLRP3","AIM2","IL1B","IL18"]

print("=" * 60)
print("SCRIPT 04 — DRUG SENSITIVITY ANALYSIS")
print("=" * 60)

# Step 1 — Load CCLE
print("\n[1/4] Loading CCLE...")
ccle_raw  = pd.read_csv(os.path.join(DATA_RAW, "ccle_expr.csv"), index_col=0)
gene_cols = [c for c in ccle_raw.columns if "(" in c and ")" in c]
ccle_expr = ccle_raw[gene_cols].copy()
ccle_expr.columns = [c.split(" (")[0].strip() for c in gene_cols]
ccle_expr.index   = ccle_raw["ModelID"].values
print(f"  CCLE: {ccle_expr.shape}")

# Step 2 — Filter to OC cell lines
print("\n[2/4] Filtering to OC cell lines...")
available = {k:v for k,v in KNOWN_OV.items() if k in ccle_expr.index}
ccle_ov   = ccle_expr.loc[list(available.keys())].copy()
ccle_ov.index = [available[k] for k in ccle_ov.index]
ccle_ov   = ccle_ov[~ccle_ov.index.duplicated(keep="first")]
print(f"  After deduplication: {ccle_ov.index.tolist()}")

# Step 3 — Pyroptosis scores
print("\n[3/4] Calculating pyroptosis scores...")
pyro_genes = [g for g in PYRO_PROXY_GENES if g in ccle_ov.columns]
ccle_ov["Pyro_score"] = ccle_ov[pyro_genes].astype(float).mean(axis=1)
median_pyro = ccle_ov["Pyro_score"].median()
ccle_ov["Group"] = ccle_ov["Pyro_score"].apply(
    lambda x: "High_Pyroptosis" if x > median_pyro else "Low_Pyroptosis"
)
print(ccle_ov[["Pyro_score","Group"]].sort_values("Pyro_score", ascending=False).round(3))

# Step 4 — GDSC2 IC50
print("\n[4/4] Matching GDSC2...")
gdsc    = pd.read_excel(os.path.join(DATA_RAW, "gdsc2.xlsx"))
results = []
for cell_name, row_data in ccle_ov.iterrows():
    gdsc_name = GDSC_NAMES.get(cell_name)
    if not gdsc_name:
        continue
    cell_data = gdsc[(gdsc["CELL_LINE_NAME"]==gdsc_name) & (gdsc["DRUG_NAME"].isin(OC_DRUGS))]
    for _, dr in cell_data.iterrows():
        results.append({"drug":dr["DRUG_NAME"],"cell":cell_name,"group":row_data["Group"],"LN_IC50":float(dr["LN_IC50"])})

res_df = pd.DataFrame(results)
print(f"  Matches: {len(res_df)}")

# Figure 5
drugs_found = res_df["drug"].unique() if len(res_df) > 0 else OC_DRUGS
n_drugs     = len(drugs_found)
fig, axes   = plt.subplots(1, n_drugs, figsize=(4*n_drugs, 6), squeeze=False)
fig.suptitle("Drug sensitivity by pyroptosis activity — OC cell lines\nGDSC2 IC50 data", fontsize=13, fontweight="bold")

print(f"\n  Results:")
for idx, drug in enumerate(drugs_found):
    ax      = axes[0][idx]
    drug_df = res_df[res_df["drug"]==drug]
    high    = drug_df[drug_df["group"]=="High_Pyroptosis"]["LN_IC50"].values
    low     = drug_df[drug_df["group"]=="Low_Pyroptosis"]["LN_IC50"].values

    plot_data, plot_labels, plot_colors = [], [], []
    if len(high)>0: plot_data.append(high);  plot_labels.append(f"High Pyro\n(n={len(high)})"); plot_colors.append("#D85A30")
    if len(low)>0:  plot_data.append(low);   plot_labels.append(f"Low Pyro\n(n={len(low)})");  plot_colors.append("#1D9E75")

    if plot_data:
        bp = ax.boxplot(plot_data, patch_artist=True, medianprops=dict(color="black",linewidth=2.5), widths=0.5)
        for patch, color in zip(bp["boxes"], plot_colors):
            patch.set_facecolor(color); patch.set_alpha(0.75)
        for j, (data, color) in enumerate(zip(plot_data, plot_colors)):
            ax.scatter(np.random.normal(j+1, 0.05, len(data)), data, color=color, alpha=0.8, s=50, zorder=5)

    if len(high)>0 and len(low)>0:
        _, p = mannwhitneyu(high, low, alternative="two-sided")
        sig  = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
        y_max = max(np.max(high), np.max(low)) + 0.3
        ax.plot([1,2],[y_max,y_max], color="black", linewidth=1)
        ax.text(1.5, y_max+0.05, sig, ha="center", fontsize=12, fontweight="bold", color="red" if p<0.05 else "gray")
        ax.set_title(f"{drug}\np = {p:.4f}", fontsize=11)
        print(f"  {drug}: High={np.median(high):.2f}  Low={np.median(low):.2f}  p={p:.4f} {sig}")

    ax.set_xticks(range(1, len(plot_labels)+1))
    ax.set_xticklabels(plot_labels, fontsize=9)
    ax.set_ylabel("LN(IC50)\nlower = more sensitive", fontsize=9)
    sns.despine(ax=ax)

plt.tight_layout()
plt.savefig(os.path.join(DATA_CLEAN, "Figure5_drug_sensitivity.png"), dpi=150, bbox_inches="tight")
plt.close()
if len(res_df)>0:
    res_df.to_csv(os.path.join(DATA_CLEAN, "drug_sensitivity_results.csv"), index=False)

print("\n  Saved: Figure5_drug_sensitivity.png")
print("\n  Note: exploratory analysis — only 7 cell lines available.")
print(f"  NEXT: Run 05_extended_validation.py")
