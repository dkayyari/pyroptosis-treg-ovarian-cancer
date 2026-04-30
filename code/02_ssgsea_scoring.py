"""
Script 02: ssGSEA Scoring
Project: Pyroptosis-Treg Axis in Ovarian Cancer
Author:  Deeksha Kayyari | Indiana University | April 2026
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from scipy import stats
import requests
import warnings
import os

warnings.filterwarnings("ignore")
np.random.seed(42)

# =============================================================
# CHANGE THIS PATH to your project folder
# =============================================================
PROJECT_ROOT = r"C:\Users\deeksha\Downloads\pyroptosis_project_old"
# =============================================================

DATA_CLEAN = os.path.join(PROJECT_ROOT, "data_clean")
os.makedirs(DATA_CLEAN, exist_ok=True)

PYROPTOSIS_GENES = [
    "NLRP1","NLRP2","NLRP3","NLRP6","NLRP7","NLRC4","AIM2","PYCARD",
    "CASP1","CASP3","CASP4","CASP5","CASP6","CASP8","CASP9",
    "GSDMA","GSDMB","GSDMC","GSDMD","GSDME","PJVK",
    "IL1B","IL18","TNF","IL6","GPX4","NOD1","NOD2",
    "ELANE","PLCG1","TIRAP","PRKACA","SCAF11"
]
TREG_GENES = [
    "FOXP3","IL2RA","CTLA4","IKZF2","TIGIT","CCR8",
    "TNFRSF18","LAYN","BATF","ENTPD1","CXCR5","TGFB1","IL10"
]
CCL22_GENES = [
    "GSDMD","CASP1","IL1B","IL18","CCL22",
    "CCL17","CCL28","IL33","CCR8","IL1RL1"
]

print("=" * 60)
print("SCRIPT 02 — ssGSEA SCORING")
print("=" * 60)

# STEP 1 — Load data
print("\n[1/6] Loading clean data...")
tcga     = pd.read_csv(os.path.join(DATA_CLEAN, "TCGA_expr_clean.csv"), index_col=0)
survival = pd.read_csv(os.path.join(DATA_CLEAN, "TCGA_survival_clean.csv"), index_col=0)
print(f"  Expression: {tcga.shape}")
print(f"  Index sample: {tcga.index[:2].tolist()}")

# STEP 2 — Convert Ensembl IDs to gene symbols if needed
if tcga.index[0].startswith("ENSG"):
    print("\n[2/6] Converting Ensembl IDs to gene symbols (takes ~10 mins)...")
    tcga.index = tcga.index.str.split(".").str[0]  # remove version: ENSG000.15 -> ENSG000
    ensembl_ids = tcga.index.tolist()
    gene_map    = {}

    for i in range(0, len(ensembl_ids), 500):
        batch = ensembl_ids[i:i+500]
        try:
            resp = requests.post(
                "https://mygene.info/v3/query",
                data={"q":",".join(batch),"scopes":"ensembl.gene","fields":"symbol","species":"human"},
                timeout=30
            )
            if resp.status_code == 200:
                for hit in resp.json():
                    if "symbol" in hit and "query" in hit:
                        gene_map[hit["query"]] = hit["symbol"]
        except Exception as e:
            print(f"  Batch {i//500+1} error: {e}")
        if i % 5000 == 0:
            print(f"  Batch {i//500+1}: {len(gene_map)} mapped")

    tcga["gene_symbol"] = tcga.index.map(gene_map)
    tcga = tcga.dropna(subset=["gene_symbol"])
    tcga = tcga.groupby("gene_symbol").mean()
    tcga.to_csv(os.path.join(DATA_CLEAN, "TCGA_expr_clean_genesymbols.csv"))
    print(f"  After conversion: {tcga.shape}")
    print("  Saved: TCGA_expr_clean_genesymbols.csv")
else:
    print("\n[2/6] Gene symbols already present — skipping conversion")

found_p = [g for g in PYROPTOSIS_GENES if g in tcga.index]
found_t = [g for g in TREG_GENES if g in tcga.index]
print(f"  Pyroptosis genes: {len(found_p)}/33")
print(f"  Treg genes:       {len(found_t)}/13")

# STEP 3 — Run ssGSEA
print("\n[3/6] Running ssGSEA (takes 20-30 minutes)...")

def run_ssgsea(expr_df, gene_list, name):
    genes = [g for g in gene_list if g in expr_df.index]
    print(f"  {name}: {len(genes)} genes found — running...")
    result = gp.ssgsea(
        data=expr_df, gene_sets={name: genes},
        outdir=None, sample_norm_method="rank",
        no_plot=True, processes=1, min_size=1
    )
    scores = result.res2d[result.res2d["Term"]==name].set_index("Name")["NES"]
    scores.name = name
    print(f"  {name}: done")
    return scores

pyro_scores  = run_ssgsea(tcga, PYROPTOSIS_GENES, "Pyroptosis_score")
treg_scores  = run_ssgsea(tcga, TREG_GENES,       "Treg_score")
ccl22_scores = run_ssgsea(tcga, CCL22_GENES,      "CCL22_score")

# STEP 4 — Combine
print("\n[4/6] Combining scores...")
scores_df = pd.concat([pyro_scores, treg_scores, ccl22_scores], axis=1)

# STEP 5 — Merge with survival
print("\n[5/6] Merging with survival data...")
common    = list(set(scores_df.index) & set(survival.index))
master_df = pd.concat([scores_df.loc[common], survival.loc[common]], axis=1)
master_df.index.name = "sample"

for col in ["Pyroptosis_score","Treg_score","CCL22_score"]:
    master_df[col] = pd.to_numeric(master_df[col], errors="coerce")
master_df = master_df.dropna(subset=["Pyroptosis_score","Treg_score"])
print(f"  Patients: {len(master_df)}")

pyro_median = master_df["Pyroptosis_score"].median()
treg_median = master_df["Treg_score"].median()
master_df["Pyro_group"]     = master_df["Pyroptosis_score"].apply(lambda x:"High_Pyro" if x>pyro_median else "Low_Pyro")
master_df["Treg_group"]     = master_df["Treg_score"].apply(lambda x:"High_Treg" if x>treg_median else "Low_Treg")
master_df["Combined_group"] = master_df["Pyro_group"] + "_" + master_df["Treg_group"]
print(master_df["Combined_group"].value_counts())

r, p = stats.pearsonr(master_df["Pyroptosis_score"], master_df["Treg_score"])
print(f"\n  KEY RESULT: r = {r:.3f}  p = {p:.4f} {'✓' if p<0.05 else '✗'}")

null_r = []
np.random.seed(42)
for i in range(1000):
    shuffled  = master_df["Treg_score"].sample(frac=1).values
    r_null, _ = stats.pearsonr(master_df["Pyroptosis_score"], shuffled)
    null_r.append(r_null)
perm_p = np.mean(np.abs(np.array(null_r)) >= np.abs(r))
print(f"  Permutation p = {perm_p:.4f}")

# STEP 6 — Save and figures
print("\n[6/6] Saving outputs...")
master_df.to_csv(os.path.join(DATA_CLEAN, "ssgsea_scores.csv"))
print("  Saved: ssgsea_scores.csv")

fig, ax = plt.subplots(figsize=(7,6))
color_map = {"High_Pyro_High_Treg":"#1D9E75","High_Pyro_Low_Treg":"#D85A30",
             "Low_Pyro_High_Treg":"#534AB7","Low_Pyro_Low_Treg":"#888780"}
for group, color in color_map.items():
    sub = master_df[master_df["Combined_group"]==group]
    ax.scatter(sub["Pyroptosis_score"], sub["Treg_score"],
               c=color, alpha=0.65, s=30, label=f"{group} (n={len(sub)})")
m, b = np.polyfit(master_df["Pyroptosis_score"].astype(float),
                  master_df["Treg_score"].astype(float), 1)
xl = np.linspace(master_df["Pyroptosis_score"].min(), master_df["Pyroptosis_score"].max(), 100)
ax.plot(xl, m*xl+b, "k--", linewidth=1.5, alpha=0.7)
ax.axvline(pyro_median, color="gray", linestyle=":", alpha=0.4)
ax.axhline(treg_median, color="gray", linestyle=":", alpha=0.4)
ax.set_xlabel("Pyroptosis ssGSEA score", fontsize=12)
ax.set_ylabel("Treg ssGSEA score", fontsize=12)
ax.set_title(f"Pyroptosis vs Treg score — TCGA-OV\nPearson r = {r:.3f}   p = {p:.4f}", fontsize=12)
ax.legend(fontsize=8, loc="upper left")
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(os.path.join(DATA_CLEAN, "Phase3_Figure1_scatter.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: Phase3_Figure1_scatter.png")

print("\n" + "="*60)
print("SCRIPT 02 COMPLETE")
print(f"  r = {r:.3f}  p = {p:.4f}  n = {len(master_df)}")
print(f"  Outputs: {DATA_CLEAN}")
print("="*60)
