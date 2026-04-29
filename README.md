# Pyroptosis-Treg Axis in Ovarian Cancer

**Author:** Deeksha Kayyari  
**Institution:** Indiana University  
**Course:** High-Throughput Bioinformatics  
**Date:** April 2026  

---

## Overview

This repository contains the complete bioinformatics analysis pipeline for the study:

> *"Pyroptosis Activity Drives Regulatory T Cell Recruitment in Ovarian Cancer via the CCL22/CCR8 Axis"*

We investigated whether pyroptosis activity explains the paradoxical enrichment of regulatory T cells (Tregs) in better-survival ovarian cancer patients — a finding reported but unexplained by Ye et al. (2021) and Cao et al. (2022).

**Key finding:** Strong positive correlation between pyroptosis activity and Treg abundance across 817 OC patients in three independent cohorts (TCGA r=0.813, GSE26193 r=0.517, GSE9891 r=0.713; all p<0.0001).

---

## Repository Structure

```
pyroptosis-treg-ovarian-cancer/
│
├── README.md                        ← this file
├── requirements.txt                 ← Python package versions
│
├── code/
│   ├── 01_data_cleaning.py          ← Phase 1+2: data download and QC
│   ├── 02_ssgsea_scoring.py         ← Phase 3: ssGSEA pathway scoring
│   ├── 03_core_analysis.py          ← Phase 4: KM survival, CCL22 heatmap, checkpoints
│   ├── 04_drug_sensitivity_final.py ← Phase 5A: drug sensitivity analysis
│   └── 05_extended_validation.py    ← Phase 5B: external GEO validation
│
├── figures/
│   ├── Figure1_scatter.png          ← Pyroptosis vs Treg scatter (TCGA)
│   ├── Figure2_checkpoints.png      ← Immune checkpoint boxplots
│   ├── Figure3_validation_GSE9891.png ← Validation scatter (Australia)
│   └── Figure4_checkpoint_validation.png ← Checkpoint validation (3 cohorts)
│
└── report/
    └── pyroptosis_report_FINAL.docx ← Final written report
```

---

## How to Reproduce

### Step 1 — Install dependencies

```bash
pip install -r requirements.txt
```

### Step 2 — Download raw data

Create the following folder structure:
```
pyroptosis_treg_project/
├── dataraw/     ← place downloaded files here
└── data_clean/  ← outputs will be saved here (created automatically)
```

Download these files and place in `dataraw/`:

| File | Source | Description |
|------|--------|-------------|
| `tcga_expr.tsv.gz` | [UCSC Xena](https://xenabrowser.net) | TCGA-OV RNA-seq expression |
| `TCGA-OV.survival.tsv` | [UCSC Xena](https://xenabrowser.net) | TCGA-OV survival data |
| `GSE26193_series_matrix.txt` | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26193) | GEO validation cohort 1 |
| `GSE9891_series_matrix.txt` | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9891) | GEO validation cohort 2 |
| `gdsc2.xlsx` | [GDSC](https://www.cancerrxgene.org/downloads/bulk_download) | Drug sensitivity data |
| `ccle_expr.csv` | [DepMap](https://depmap.org/portal/download/) | Cell line expression |

For GSE9891 (command line):
```bash
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9891/matrix/GSE9891_series_matrix.txt.gz
gunzip GSE9891_series_matrix.txt.gz
```

### Step 3 — Update file paths

Open each script and update the `DATA_RAW` and `DATA_CLEAN` variables at the top to match your system:

```python
DATA_RAW   = r"path/to/your/dataraw"
DATA_CLEAN = r"path/to/your/data_clean"
```

### Step 4 — Run scripts in order

```bash
python code/01_data_cleaning.py
python code/02_ssgsea_scoring.py       # takes ~30 minutes
python code/03_core_analysis.py
python code/04_drug_sensitivity_final.py
python code/05_extended_validation.py  # takes ~20 minutes first run (API call)
```

---

## Key Results

| Result | Value | Significance |
|--------|-------|-------------|
| Pyroptosis vs Treg (TCGA) | r = 0.813 | p < 0.0001 |
| Validation GSE26193 | r = 0.517 | p < 0.0001 |
| Validation GSE9891 | r = 0.713 | p < 0.0001 |
| Permutation test | p = 0.0000 | 0/1000 permutations exceeded real r |
| Immune checkpoints | all p < 0.0001 | 5/6 validated in GEO cohorts |

---

## Software

| Package | Version | Purpose |
|---------|---------|---------|
| Python | 3.11 | Core language |
| pandas | 2.0.3 | Data manipulation |
| numpy | 1.24.4 | Numerical computing |
| scipy | 1.11.4 | Statistical tests |
| gseapy | 1.0.6 | ssGSEA scoring |
| lifelines | 0.27.8 | Survival analysis |
| matplotlib | 3.7.3 | Plotting |
| seaborn | 0.12.2 | Statistical visualization |
| requests | 2.31.0 | mygene.info API (probe mapping) |

---

## References

1. Ye Z. et al. (2021). *Cell Death Discovery*, 7, 117.
2. Cao Y. et al. (2022). *Frontiers in Oncology*, 12, 857919.
3. Subramanian A. et al. (2005). *PNAS*, 102(43), 15545–15550.
4. Goldman M.J. et al. (2020). *Nature Biotechnology*, 38, 675–678.
