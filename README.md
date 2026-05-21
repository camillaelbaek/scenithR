# Scenith Flow Cytometry Analysis

A reproducible analysis pipeline for **Scenith** metabolic assays using flow cytometry. Includes both an interactive Shiny app and a self-contained R Markdown report.

---

## Overview

[Scenith](https://www.scenith.com/) (Single Cell ENergetic meTabolism by profilIng tHe protein synthesis) profiles cellular metabolism by measuring **puromycin incorporation** (translation) as a proxy for ATP availability, under targeted metabolic inhibition:

| Condition | Inhibitor | What it blocks |
|-----------|-----------|----------------|
| **Co** | — | Control (DMSO) |
| **DG** | 2-Deoxy-D-glucose | Glycolysis |
| **O** | Oligomycin A | OXPHOS (ATP synthase) |
| **DGO** | DG + Oligomycin | Both pathways |

Puromycin signal (APC-A) reflects translation, which is ATP-dependent. Four metabolic parameters are derived per genotype:

1. **Glucose dependence** = 100 × (Co − DG) / (Co − DGO)
2. **FAO/AAO capacity** = 100 − glucose dependence
3. **Mitochondrial dependence** = 100 × (Co − O) / (Co − DGO)
4. **Glycolytic capacity** = 100 − mitochondrial dependence

---

## Repository structure

```
.
├── ScenithApp/
│   ├── app.R                  # Shiny app (main analysis tool)
│   ├── run_app.R              # Launcher: installs packages + starts app
│   └── Mac_run_app.command    # Double-click launcher for macOS
├── analysis_v2.Rmd            # R Markdown report (scripted / offline use)
├── fcs-input/                 # FCS files — git-ignored
│   ├── exp1/
│   └── exp2/
├── .gitignore
└── README.md
```

> FCS files (`.fcs`), Excel files (`.xlsx`), Word documents, and FlowJo workspaces (`.wsp`) are listed in `.gitignore` and are not tracked by git.

---

## Quick start

### Option A — Shiny app (recommended)

**macOS:** double-click `ScenithApp/Mac_run_app.command`

**Windows / Linux:** open a terminal in the `ScenithApp/` folder and run:
```
Rscript run_app.R
```

The app will open in your browser. Then follow these steps:

1. **Upload FCS files** — select all `.fcs` files from your experiment folder
2. **Upload plate metadata** — a CSV or XLSX file describing your plate layout (see [Metadata format](#metadata-format) below)
3. Adjust gating thresholds if needed (see [Gating parameters](#gating-parameters))
4. Click **Run / Recompute analysis**

All plots have a **Download PNG** button.

### Option B — R Markdown report

Open `analysis_v2.Rmd` in RStudio and click **Knit**. Before knitting, update the two paths near the top of the document:

- `fcs.dir` — path to the folder containing your FCS files
- `read_xlsx(...)` — path to your metadata XLSX file

---

## Metadata format

The metadata file tells the app which wells belong to which experimental group. It must be a **CSV** or **XLSX** file with at least these three columns:

| Column | Description | Example values |
|--------|-------------|----------------|
| `well_code` | Well position: row letter + zero-padded column number | `B01`, `C12`, `G04` |
| `genotype` | Cell genotype or cell line label | `WT`, `KO`, `WT_HSV` |
| `treatment` | Metabolic condition applied to the well | `Co`, `DG`, `O`, `DGO` |

Additional columns are allowed (e.g. `replicate`, `passage`, `experimenter`) and will be carried through all tables.

**Download a pre-filled template** directly from the app sidebar ("Download CSV template"). It contains all rows B–H × columns 1–12 with empty `genotype` and `treatment` columns — fill those in and save.

### Well code format

- Use the format `B01` (uppercase letter + two-digit column number).
- Single-digit formats like `B1` are also accepted and normalized automatically.
- The app extracts well codes from FCS filenames using the pattern `[A-H][0-9]{1,2}` (e.g. `Specimen_004_B1_B01.fcs` → `B01`).

### Example metadata CSV

```csv
well_code,genotype,treatment
B01,WT,Co
B02,WT,Co
B03,WT,Co
B04,WT_HSV,Co
B05,WT_HSV,Co
C01,WT,DG
C02,WT,DG
D01,WT,O
E01,WT,DGO
G01,WT,UNST
```

### Treatment names for Scenith parameters

Scenith parameters are computed only when treatments spelled exactly **`Co`**, **`DG`**, **`O`**, and **`DGO`** are present in the metadata. You may have additional treatments — they will appear in QC and summary plots but not in the Scenith parameter calculations.

---

## Gating parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| FITC-A threshold | 4 000 | Live/Dead gate: keep cells **below** this value (low FITC-A = live) |
| APC-A threshold | 80 | Puromycin gate: keep cells **at or above** this value (puro+) |

Review the QC plots after running and adjust these in the sidebar if needed, then click **Run** again. The "Representative sample index" fields control which individual sample is shown in the singlet and live-gate diagnostic plots.

---

## App tabs

| Tab | What it shows |
|-----|---------------|
| **Overview** | Upload status checklist and experiment description |
| **Metadata** | Preview of uploaded metadata + per-FCS match status (green = matched, red = unmatched) |
| **Sample mapping** | Full table of FCS files joined to their metadata, with unmatched-well warnings |
| **Gating QC** | Singlet gate, live gate, and puromycin gate diagnostic plots |
| **Cell counts** | Bar plot and table of cell counts at each gating step |
| **Puromycin summary** | Mean APC-A per sample (all live cells and puro+ cells) |
| **Scenith parameters** | Derived metabolic parameters, density overlay plots, and per-genotype bar plots |

---

## Required R packages

`run_app.R` automatically installs any missing packages when you launch the app.

**CRAN packages**
```
shiny, dplyr, tidyr, ggplot2, DT, stringr, purrr, scales,
readr, readxl, ggridges, ggpubr, viridis, ggbeeswarm, sp
```

**Bioconductor packages** (installed via `BiocManager`)
```
flowCore, flowViz, ggcyto, openCyto
```

If Bioconductor packages fail to install automatically, run this in R:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("flowCore", "flowViz", "ggcyto", "openCyto"))
```

Requires **R ≥ 4.2**.

---

## Troubleshooting

**"Missing required column(s): genotype, treatment"**
Your metadata file is missing one or more required columns. Download the CSV template from the app sidebar, fill in the `genotype` and `treatment` columns, and re-upload.

**"X sample(s) could not be matched to the metadata" (orange warning)**
The well code extracted from those FCS filenames does not appear in the uploaded metadata. Check that the filename contains a standard well code (`[A-H][0-9]{1,2}`) and that the metadata covers all expected wells.

**Scenith plots show "treatments required" message**
The treatments `Co`, `DG`, `O`, and `DGO` must be present in the metadata and spelled exactly that way for the Scenith parameter calculations to run.

**App is slow on the all-samples FSC plot**
This plot loads raw events from every FCS file simultaneously. Keep the checkbox disabled until you specifically need it.

**Bioconductor packages fail to install**
Ensure R ≥ 4.2. Try `BiocManager::install(...)` manually (see above).

---

## Authors

- **Alba Perez Arribas** — experimental design and assay  
- **celbaek** — analysis pipeline and Shiny app
