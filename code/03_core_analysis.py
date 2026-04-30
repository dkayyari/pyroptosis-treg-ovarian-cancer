"""
Script 03: Core Analysis
Project: Pyroptosis-Treg Axis in Ovarian Cancer
Author:  Deeksha Kayyari | Indiana University | April 2026
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.stats import kruskal
import warnings
import os

warnings.filterwarnings("ignore")
np.random.seed(42)

try:
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test, multivariate_logrank_test
except ImportError:
    os.system("pip install lifelines")
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test, multivariate_logrank_test

# =============================================================
# CHANGE THIS PATH to your project folder
# =============================================================
PROJECT_ROOT = r"C:\Users\deeksha\Downloads\pyroptosis_project_old"
# =============================================================

DATA_CLEAN = os.path.join(PROJECT_ROOT, "data_clean")

GROUP_COLORS = {
    "High_Pyro_High_Treg": "#1D9E75",
    "High_Pyro_Low_Treg" : "#D85A30",
    "Low_Pyro_High_Treg" : "#534AB7",
    "Low_Pyro_Low_Treg"  : "#888780"
}
GROUP_LABELS = {
    "High_Pyro_High_Treg": "High Pyroptosis + High Treg",
    "High_Pyro_Low_Treg" : "High Pyroptosis + Low Treg",
    "Low_Pyro_High_Treg" : "Low Pyroptosis + High Treg",
    "Low_Pyro_Low_Treg"  : "Low Pyroptosis + Low Treg"
}

print("=" * 60)
print("SCRIPT 03 — CORE ANALYSIS")
print("=" * 60)

# Load data
print("\n[Loading data]")

master_df = pd.read_csv(os.path.join(DATA_CLEAN, "ssgsea_scores.csv"), index_col=0)

# Use gene symbol version if available (created by script 02)
gene_sym_path = os.path.join(DATA_CLEAN, "TCGA_expr_clean_genesymbols.csv")
expr_path     = os.path.join(DATA_CLEAN, "TCGA_expr_clean.csv")
tcga = pd.read_csv(gene_sym_path if os.path.exists(gene_sym_path) else expr_path, index_col=0)
print(f"  Expression: {tcga.shape}")

for col in ["Pyroptosis_score","Treg_score","CCL22_score","OS","OS_time"]:
    if col in master_df.columns:
        master_df[col] = pd.to_numeric(master_df[col], errors="coerce")

master_df = master_df.dropna(subset=["OS","OS_time","Pyroptosis_score","Treg_score"])
master_df = master_df[master_df["OS_time"] > 0]

# Keep only samples present in both files
common_samples = list(set(tcga.columns) & set(master_df.index))
master_df      = master_df.loc[master_df.index.isin(common_samples)]
print(f"  Common samples: {len(common_samples)}")

if "Combined_group" not in master_df.columns:
    pm = master_df["Pyroptosis_score"].median()
    tm = master_df["Treg_score"].median()
    master_df["Pyro_group"]     = master_df["Pyroptosis_score"].apply(lambda x:"High_Pyro" if x>pm else "Low_Pyro")
    master_df["Treg_group"]     = master_df["Treg_score"].apply(lambda x:"High_Treg" if x>tm else "Low_Treg")
    master_df["Combined_group"] = master_df["Pyro_group"] + "_" + master_df["Treg_group"]

print(f"  Patients: {len(master_df)},  Deaths: {int(master_df['OS'].sum())}")
print(f"  Groups:\n{master_df['Combined_group'].value_counts()}")

# FIGURE 2 — KM Survival
print("\n[FIGURE 2] Kaplan-Meier Survival Curves")
fig, ax = plt.subplots(figsize=(10, 7))
kmf = KaplanMeierFitter()
for group, color in GROUP_COLORS.items():
    subset = master_df[master_df["Combined_group"] == group]
    if len(subset) == 0: continue
    kmf.fit(durations=subset["OS_time"]/365, event_observed=subset["OS"],
            label=f"{GROUP_LABELS[group]} (n={len(subset)})")
    kmf.plot_survival_function(ax=ax, color=color, linewidth=2, ci_show=True, ci_alpha=0.08)

results   = multivariate_logrank_test(master_df["OS_time"], master_df["Combined_group"], master_df["OS"])
p_logrank = results.p_value
print(f"  Log-rank p = {p_logrank:.4f}")

ax.set_xlabel("Time (years)", fontsize=13)
ax.set_ylabel("Overall survival probability", fontsize=13)
ax.set_title(f"Overall survival by pyroptosis-Treg group — TCGA-OV (n={len(master_df)})\nLog-rank p = {p_logrank:.4f}", fontsize=13)
ax.legend(loc="upper right", fontsize=9)
ax.set_ylim(0, 1.05); ax.set_xlim(left=0)
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(os.path.join(DATA_CLEAN, "Figure2_KM_survival.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: Figure2_KM_survival.png")

# FIGURE 3 — CCL22 Heatmap
print("\n[FIGURE 3] CCL22 Mechanism Heatmap")
EXECUTOR_GENES  = [g for g in ["GSDMD","CASP1","IL1B","IL18","NLRP3","AIM2"] if g in tcga.index]
RECRUITER_GENES = [g for g in ["CCL22","CCL17","CCR8","IL33","TGFB1","FOXP3"] if g in tcga.index]

mech_expr   = tcga.loc[EXECUTOR_GENES + RECRUITER_GENES, common_samples].T.astype(float)
corr_matrix = pd.DataFrame(index=EXECUTOR_GENES, columns=RECRUITER_GENES, dtype=float)
pval_matrix = pd.DataFrame(index=EXECUTOR_GENES, columns=RECRUITER_GENES, dtype=float)
annot       = pd.DataFrame(index=EXECUTOR_GENES, columns=RECRUITER_GENES, dtype=str)

for eg in EXECUTOR_GENES:
    for rg in RECRUITER_GENES:
        r, p = stats.pearsonr(mech_expr[eg], mech_expr[rg])
        corr_matrix.loc[eg, rg] = r
        pval_matrix.loc[eg, rg] = p
        stars = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else ""
        annot.loc[eg, rg] = f"{r:.2f}{stars}"

fig, ax = plt.subplots(figsize=(9, 5))
sns.heatmap(corr_matrix.astype(float), annot=annot, fmt="",
            cmap="RdYlGn", center=0, vmin=-0.5, vmax=0.5,
            linewidths=0.5, ax=ax, cbar_kws={"label":"Pearson r","shrink":0.8})
ax.set_title(f"Pyroptosis executors vs Treg recruitment genes\nPearson correlation — TCGA-OV (n={len(common_samples)})", fontsize=12)
ax.set_xlabel("Treg recruitment genes →", fontsize=11)
ax.set_ylabel("Pyroptosis executor genes →", fontsize=11)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
fig.text(0.99, 0.01, "* p<0.05  ** p<0.01  *** p<0.001", ha="right", va="bottom", fontsize=9, color="gray")
plt.tight_layout()
plt.savefig(os.path.join(DATA_CLEAN, "Figure3_CCL22_heatmap.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: Figure3_CCL22_heatmap.png")

# FIGURE 4 — Checkpoints
print("\n[FIGURE 4] Immune Checkpoint Boxplot")
CHECKPOINT_GENES = {"CD274":"PD-L1","PDCD1":"PD-1","CTLA4":"CTLA-4","HAVCR2":"TIM-3","LAG3":"LAG-3","TIGIT":"TIGIT"}
checkpoint_genes = [g for g in CHECKPOINT_GENES if g in tcga.index]
checkpoint_expr  = tcga.loc[checkpoint_genes, common_samples].T.astype(float)

kw_results  = {}
group_order = list(GROUP_COLORS.keys())
print("  Kruskal-Wallis test:")
for gene in checkpoint_genes:
    groups_data = [
        checkpoint_expr.loc[master_df[master_df["Combined_group"]==g].index, gene].dropna().values
        for g in group_order if len(master_df[master_df["Combined_group"]==g]) > 0
    ]
    groups_data = [g for g in groups_data if len(g) > 0]
    if len(groups_data) >= 2:
        _, p_kw = kruskal(*groups_data)
        kw_results[gene] = p_kw
        print(f"    {CHECKPOINT_GENES[gene]:8s}: p = {p_kw:.4f} {'✓' if p_kw<0.05 else '✗'}")

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle("Immune checkpoint expression by pyroptosis-Treg group\nTCGA-OV", fontsize=13, fontweight="bold")
axes = axes.flatten()
for idx, gene in enumerate(checkpoint_genes):
    ax   = axes[idx]
    p_kw = kw_results.get(gene, 1.0)
    plot_data = [
        checkpoint_expr.loc[checkpoint_expr.index.isin(master_df[master_df["Combined_group"]==g].index), gene].dropna().values
        for g in group_order
    ]
    bp = ax.boxplot(plot_data, patch_artist=True, medianprops=dict(color="black", linewidth=2))
    for patch, group in zip(bp["boxes"], group_order):
        patch.set_facecolor(GROUP_COLORS[group]); patch.set_alpha(0.75)
    ax.set_title(f"{CHECKPOINT_GENES[gene]}\np = {p_kw:.4f}{'*' if p_kw<0.05 else ''}", fontsize=11)
    ax.set_xticks(range(1, len(group_order)+1))
    ax.set_xticklabels([g.replace("_","\n") for g in group_order], fontsize=7)
    ax.set_ylabel("log2(FPKM+1)", fontsize=9)
    sns.despine(ax=ax)

legend_patches = [mpatches.Patch(color=c, label=GROUP_LABELS[g], alpha=0.75) for g,c in GROUP_COLORS.items()]
fig.legend(handles=legend_patches, loc="lower center", ncol=2, fontsize=9, bbox_to_anchor=(0.5,-0.02))
plt.tight_layout(rect=[0,0.06,1,1])
plt.savefig(os.path.join(DATA_CLEAN, "Figure4_checkpoints.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: Figure4_checkpoints.png")

print("\n" + "="*60)
print("SCRIPT 03 COMPLETE")
print(f"  KM log-rank p = {p_logrank:.4f}")
sig_ck = [CHECKPOINT_GENES.get(g,g) for g,p in kw_results.items() if p<0.05]
print(f"  Significant checkpoints: {sig_ck}")
print(f"  Outputs: {DATA_CLEAN}")
print("="*60)
