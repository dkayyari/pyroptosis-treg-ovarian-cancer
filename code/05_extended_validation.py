"""
Script 05: External Validation in GEO Cohorts
Project: Pyroptosis-Treg Axis in Ovarian Cancer
Author:  Deeksha Kayyari | Indiana University | April 2026
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import kruskal, pearsonr
from matplotlib.patches import Patch
import warnings
import os
import requests
from io import StringIO

warnings.filterwarnings("ignore")
np.random.seed(42)

# =============================================================
# CHANGE THIS PATH to your project folder
# =============================================================
PROJECT_ROOT = r"C:\Users\deeksha\Downloads\pyroptosis_project_old"
# =============================================================

DATA_RAW   = os.path.join(PROJECT_ROOT, "dataraw")
DATA_CLEAN = os.path.join(PROJECT_ROOT, "data_clean")

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
CHECKPOINT_GENES = {
    "CD274":"PD-L1","PDCD1":"PD-1","CTLA4":"CTLA-4",
    "HAVCR2":"TIM-3","LAG3":"LAG-3","TIGIT":"TIGIT"
}
GROUP_ORDER  = ["High_Pyro_High_Treg","High_Pyro_Low_Treg","Low_Pyro_High_Treg","Low_Pyro_Low_Treg"]
GROUP_LABELS = ["HH","HL","LH","LL"]
GROUP_COLORS = ["#2ecc71","#e74c3c","#9b59b6","#95a5a6"]

# =============================================================
# HELPER FUNCTIONS
# =============================================================

def parse_geo(filepath):
    print(f"  Parsing {os.path.basename(filepath)}...")
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    sample_ids, data_start, meta_lines = None, None, []
    for i, line in enumerate(lines):
        if line.startswith("!Sample_geo_accession"):
            sample_ids = [s.strip().strip('"') for s in line.strip().split("\t")[1:]]
        if line.startswith("!Sample_"):
            meta_lines.append(line.strip())
        if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
            data_start = i; break
    data_lines = [l for l in lines[data_start:] if not l.startswith("!series_matrix_table_end")]
    expr_df = pd.read_csv(StringIO("".join(data_lines)), sep="\t", index_col=0)
    expr_df.index   = expr_df.index.str.strip('"')
    expr_df.columns = expr_df.columns.str.strip('"')
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce")
    if sample_ids and len(sample_ids) == expr_df.shape[1]:
        expr_df.columns = sample_ids
    print(f"  Expression: {expr_df.shape}")
    return expr_df, sample_ids, meta_lines

def map_probes_to_genes(expr_df, cohort_name):
    cache_path = os.path.join(DATA_CLEAN, f"{cohort_name}_gene_expr.csv")
    if os.path.exists(cache_path):
        print(f"  Loading cached gene matrix for {cohort_name}...")
        expr_gene = pd.read_csv(cache_path, index_col=0)
        print(f"  Gene matrix: {expr_gene.shape}")
        return expr_gene
    print(f"  Mapping probes via mygene.info API (takes ~20 min)...")
    probe_ids = expr_df.index.tolist()
    gene_map  = {}
    for i in range(0, min(len(probe_ids), 60000), 1000):
        batch = probe_ids[i:i+1000]
        try:
            resp = requests.post(
                "https://mygene.info/v3/query",
                data={"q":",".join(str(b) for b in batch),"scopes":"reporter","fields":"symbol","species":"human"},
                timeout=30
            )
            if resp.status_code == 200:
                for hit in resp.json():
                    if "symbol" in hit and "query" in hit:
                        gene_map[hit["query"]] = hit["symbol"]
        except: pass
        if i % 10000 == 0:
            print(f"    Batch {i//1000+1}: {len(gene_map)} mapped")
    expr_df2 = expr_df.copy()
    expr_df2["gene_symbol"] = expr_df2.index.map(gene_map)
    expr_df2 = expr_df2.dropna(subset=["gene_symbol"])
    expr_gene = expr_df2.groupby("gene_symbol").mean()
    expr_gene.to_csv(cache_path)
    print(f"  Gene matrix: {expr_gene.shape}")
    return expr_gene

def calc_score(expr_gene, gene_list, name):
    found = [g for g in gene_list if g in expr_gene.index]
    if not found: return None
    print(f"  {name}: {len(found)}/{len(gene_list)} genes found")
    scores = expr_gene.loc[found].mean(axis=0)
    scores.name = f"{name}_score"
    return scores

def stratify_patients(pyro, treg):
    df = pd.DataFrame({"Pyroptosis": pyro, "Treg": treg}).dropna()
    pm = df["Pyroptosis"].median(); tm = df["Treg"].median()
    df["Pyro_group"] = df["Pyroptosis"].apply(lambda x:"High_Pyro" if x>pm else "Low_Pyro")
    df["Treg_group"] = df["Treg"].apply(lambda x:"High_Treg" if x>tm else "Low_Treg")
    df["Group"]      = df["Pyro_group"] + "_" + df["Treg_group"]
    return df

# =============================================================
# LOAD AND VALIDATE GEO COHORTS
# =============================================================
print("=" * 60)
print("SCRIPT 05 — EXTERNAL GEO VALIDATION")
print("=" * 60)

cohorts = {}
for name, filename in [("GSE26193","GSE26193_series_matrix.txt"),("GSE9891","GSE9891_series_matrix.txt")]:
    filepath = os.path.join(DATA_RAW, filename)
    if not os.path.exists(filepath):
        print(f"\nWARNING: {filename} not found — skipping {name}")
        continue
    print(f"\n{'='*50}\nPROCESSING {name}\n{'='*50}")
    expr_raw, sids, mls = parse_geo(filepath)
    expr_gene = map_probes_to_genes(expr_raw, name)
    pyro = calc_score(expr_gene, PYROPTOSIS_GENES, "Pyroptosis")
    treg = calc_score(expr_gene, TREG_GENES, "Treg")
    if pyro is None or treg is None:
        print(f"  SKIPPED — insufficient genes"); continue
    df   = stratify_patients(pyro, treg)
    r, p = pearsonr(df["Pyroptosis"], df["Treg"])
    print(f"\n  *** KEY RESULT — {name} ***")
    print(f"  r = {r:.3f}  p = {p:.4f}  n = {len(df)}")
    print(f"  {'✓ CONFIRMED' if r>0 and p<0.05 else '✗ Not significant'}")
    cohorts[name] = {"expr":expr_gene,"df":df,"n":len(df),"r":r,"p":p}

# Scatter plot for GSE9891
if "GSE9891" in cohorts:
    df = cohorts["GSE9891"]["df"]
    r  = cohorts["GSE9891"]["r"]
    p  = cohorts["GSE9891"]["p"]
    fig, ax = plt.subplots(figsize=(7,6))
    color_map = {"High_Pyro_High_Treg":"#1D9E75","High_Pyro_Low_Treg":"#D85A30",
                 "Low_Pyro_High_Treg":"#534AB7","Low_Pyro_Low_Treg":"#888780"}
    for group, color in color_map.items():
        sub = df[df["Group"]==group]
        ax.scatter(sub["Pyroptosis"],sub["Treg"],c=color,alpha=0.6,s=25,label=f"{group} (n={len(sub)})")
    m,b = np.polyfit(df["Pyroptosis"].astype(float),df["Treg"].astype(float),1)
    xl  = np.linspace(df["Pyroptosis"].min(),df["Pyroptosis"].max(),100)
    ax.plot(xl,m*xl+b,"k--",linewidth=1.5,alpha=0.7)
    ax.set_title(f"GSE9891 — Pyroptosis vs Treg\nr = {r:.3f}  p = {p:.4f}",fontsize=12)
    ax.set_xlabel("Pyroptosis score",fontsize=11); ax.set_ylabel("Treg score",fontsize=11)
    ax.legend(fontsize=8,loc="upper left"); sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(DATA_CLEAN,"Validation_GSE9891_scatter.png"),dpi=150,bbox_inches="tight")
    plt.close()
    print("\n  Saved: Validation_GSE9891_scatter.png")

# Load TCGA
print("\nLoading TCGA...")
gene_sym_path = os.path.join(DATA_CLEAN, "TCGA_expr_clean_genesymbols.csv")
expr_path     = os.path.join(DATA_CLEAN, "TCGA_expr_clean.csv")
tcga_expr     = pd.read_csv(gene_sym_path if os.path.exists(gene_sym_path) else expr_path, index_col=0)
tcga_scores   = pd.read_csv(os.path.join(DATA_CLEAN,"ssgsea_scores.csv"), index_col=0)
for col in ["Pyroptosis_score","Treg_score"]:
    tcga_scores[col] = pd.to_numeric(tcga_scores[col], errors="coerce")
tcga_scores = tcga_scores.dropna(subset=["Pyroptosis_score","Treg_score"])
pm = tcga_scores["Pyroptosis_score"].median(); tm = tcga_scores["Treg_score"].median()
tcga_scores["Pyro_group"] = tcga_scores["Pyroptosis_score"].apply(lambda x:"High_Pyro" if x>pm else "Low_Pyro")
tcga_scores["Treg_group"] = tcga_scores["Treg_score"].apply(lambda x:"High_Treg" if x>tm else "Low_Treg")
tcga_scores["Group"]      = tcga_scores["Pyro_group"] + "_" + tcga_scores["Treg_group"]
print(f"  TCGA: {len(tcga_scores)} patients")

# Checkpoint validation
def get_ck_expr(expr_gene, scores_df):
    ck_genes = {k:v for k,v in CHECKPOINT_GENES.items() if k in expr_gene.index}
    if not ck_genes: return None, None
    ck_expr = expr_gene.loc[list(ck_genes.keys())].T
    ck_expr.columns = [CHECKPOINT_GENES[c] for c in ck_expr.columns]
    common = list(set(ck_expr.index) & set(scores_df.index))
    if not common: return None, None
    return ck_expr.loc[common].astype(float), scores_df.loc[common,"Group"]

cohort_data = {}
ck_tcga, grps_tcga = get_ck_expr(tcga_expr, tcga_scores)
if ck_tcga is not None: cohort_data["TCGA-OV\n(n=425)"] = (ck_tcga, grps_tcga)
for name, d in cohorts.items():
    ck, grps = get_ck_expr(d["expr"], d["df"])
    if ck is not None: cohort_data[f"{name}\n(n={d['n']})"] = (ck, grps)

all_ck_labels = sorted(set(col for ck_expr,_ in cohort_data.values() for col in ck_expr.columns))
n_ck = len(all_ck_labels); n_cohorts = len(cohort_data)

fig, axes = plt.subplots(n_cohorts, n_ck, figsize=(3*n_ck, 4*n_cohorts), squeeze=False)
fig.suptitle("Immune Checkpoint Expression — Validation Across Cohorts", fontsize=13, fontweight="bold")

print("\n  Kruskal-Wallis p-values:")
print(f"  {'Gene':<10}", end="")
for cl in cohort_data.keys(): print(f" {cl.split(chr(10))[0]:>12}", end="")
print()

for ck_idx, ck_label in enumerate(all_ck_labels):
    row_pvals = []
    for c_idx, (cohort_label, (ck_expr, grps)) in enumerate(cohort_data.items()):
        ax = axes[c_idx][ck_idx]
        if ck_label not in ck_expr.columns:
            ax.text(0.5,0.5,"N/A",ha="center",va="center",transform=ax.transAxes); row_pvals.append("N/A"); continue
        plot_data, plot_labels_ax, plot_colors_ax = [], [], []
        for g,gl,gc in zip(GROUP_ORDER,GROUP_LABELS,GROUP_COLORS):
            vals = ck_expr.loc[grps==g,ck_label].dropna().values if (grps==g).any() else np.array([])
            if len(vals)>0: plot_data.append(vals); plot_labels_ax.append(gl); plot_colors_ax.append(gc)
        if len(plot_data)>=2:
            bp = ax.boxplot(plot_data,patch_artist=True,medianprops=dict(color="black",linewidth=1.5),widths=0.6,showfliers=False)
            for patch,color in zip(bp["boxes"],plot_colors_ax): patch.set_facecolor(color); patch.set_alpha(0.7)
            _,p = kruskal(*plot_data)
            sig = "***" if p<0.001 else "**" if p<0.01 else "*" if p<0.05 else "ns"
            row_pvals.append(f"{p:.4f}{sig}")
            ax.set_title(f"{ck_label}\np={p:.4f}{sig}" if c_idx==0 else f"p={p:.4f}{sig}",fontsize=9)
        else: row_pvals.append("N/A")
        ax.set_xticks(range(1,len(plot_labels_ax)+1)); ax.set_xticklabels(plot_labels_ax,fontsize=8)
        if ck_idx==0: ax.set_ylabel(f"{cohort_label.split(chr(10))[0]}\nlog2(expr)",fontsize=8)
        sns.despine(ax=ax)
    print(f"  {ck_label:<10}", end="")
    for pv in row_pvals: print(f" {pv:>12}", end="")
    print()

legend_elements = [Patch(facecolor=GROUP_COLORS[i],alpha=0.7,label=f"{GROUP_ORDER[i]} ({GROUP_LABELS[i]})") for i in range(4)]
fig.legend(handles=legend_elements,loc="lower center",ncol=4,fontsize=8,bbox_to_anchor=(0.5,-0.02))
plt.tight_layout(rect=[0,0.04,1,1])
plt.savefig(os.path.join(DATA_CLEAN,"Figure5_checkpoint_validation.png"),dpi=150,bbox_inches="tight")
plt.close()
print("\n  Saved: Figure5_checkpoint_validation.png")

print("\n" + "="*60)
print("SCRIPT 05 COMPLETE — VALIDATION SUMMARY")
print("="*60)
print(f"  TCGA-OV:   r = 0.813  p < 0.0001  n = 425")
for name, d in cohorts.items():
    sig = "✓ CONFIRMED" if d["p"]<0.05 and d["r"]>0 else "✗"
    print(f"  {name}:  r = {d['r']:.3f}  p = {d['p']:.4f}  n = {d['n']}  {sig}")
print(f"\n  Outputs: {DATA_CLEAN}")
print("="*60)
