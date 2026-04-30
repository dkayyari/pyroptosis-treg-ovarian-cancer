# Pyroptosis-Treg Axis in Ovarian Cancer

**Author:** Deeksha Kayyari | Indiana University | April 2026

---

## Overview

Bioinformatics analysis investigating whether pyroptosis activity explains the paradoxical enrichment of regulatory T cells (Tregs) in better-survival ovarian cancer patients.

**Key finding:** Strong positive correlation between pyroptosis activity and Treg abundance across 817 OC patients in three independent cohorts (r = 0.813, 0.517, 0.713; all p < 0.0001).

---

## Repository Structure

```
pyroptosis-treg-ovarian-cancer/
├── README.md
├── requirements.txt
└── code/
    ├── 01_data_cleaning.py
    ├── 02_ssgsea_scoring.py
    ├── 03_core_analysis.py
    ├── 04_drug_sensitivity_final.py
    └── 05_extended_validation.py
```

---

## How to Reproduce

### Step 1 — Clone the repository
```
git clone https://github.com/dkayyari/pyroptosis-treg-ovarian-cancer.git
cd pyroptosis-treg-ovarian-cancer
```

### Step 2 — Install packages
```
pip install -r requirements.txt
```

### Step 3 — Download raw data

Create a `dataraw/` folder inside the project and download these files:

| File | Source | Size |
|------|--------|------|
| `tcga_expr.tsv.gz` | [UCSC Xena](https://xenabrowser.net) → TCGA-OV → gene expression RNAseq (IlluminaHiSeq) | ~119 MB |
| `TCGA-OV.survival.tsv` | [UCSC Xena](https://xenabrowser.net) → TCGA-OV → survival data | small |
| `GSE26193_series_matrix.txt` | [GEO GSE26193](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26193) | ~45 MB |
| `GSE9891_series_matrix.txt` | [GEO GSE9891](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9891) | ~118 MB |
| `gdsc2.xlsx` | [GDSC2](https://www.cancerrxgene.org/downloads/bulk_download) → GDSC2 drug sensitivity | ~20 MB |
| `ccle_expr.csv` | [DepMap](https://depmap.org/portal/download/) → CCLE expression | ~291 MB |

For GSE9891 run in terminal:
```
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9891/matrix/GSE9891_series_matrix.txt.gz
gunzip GSE9891_series_matrix.txt.gz
```

### Step 4 — Set up folder structure

```
pyroptosis_project/
├── code/          ← scripts from this repo
├── dataraw/       ← downloaded files from Step 3
└── data_clean/    ← created automatically when scripts run
```

### Step 5 — Update path in each script

Open each script and change this one line at the top:

```python
PROJECT_ROOT = r"C:\Users\yourname\Downloads\pyroptosis_project"
```

Replace with the actual path to your `pyroptosis_project` folder.

### Step 6 — Run scripts in order

Open each script in Jupyter Notebook and click **Cell → Run All**

```
01_data_cleaning.py          (~5 min)   — cleans TCGA data
02_ssgsea_scoring.py         (~30 min)  — calculates ssGSEA scores
03_core_analysis.py          (~2 min)   — KM survival, heatmap, checkpoints
04_drug_sensitivity_final.py (~2 min)   — drug sensitivity analysis
05_extended_validation.py    (~20 min)  — GEO validation
```

All figures and cleaned data are saved to `data_clean/`.

---

## Key Results

| Cohort | Pearson r | p-value | n |
|--------|-----------|---------|---|
| TCGA-OV | 0.813 | < 0.0001 | 425 |
| GSE26193 | 0.517 | < 0.0001 | 107 |
| GSE9891 | 0.713 | < 0.0001 | 285 |

---

## Software

```
Python 3.11 | pandas 2.0 | numpy 1.24 | scipy 1.11
matplotlib 3.7 | seaborn 0.12 | gseapy 1.0
lifelines 0.27 | requests 2.31 | openpyxl 3.1
```
