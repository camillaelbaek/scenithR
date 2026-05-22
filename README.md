# scenithR — Scenith Flow Cytometry Analysis

A reproducible analysis pipeline for **Scenith** metabolic assays using flow cytometry. Includes an interactive Shiny app and a self-contained R Markdown report.

---

## Quick start

### Easiest — run directly from R (no download needed)

If you have R installed, paste this into your R console:

```r
shiny::runGitHub("scenithR", "camillaelbaek", subdir = "ScenithApp")
```

R fetches the app from GitHub and opens it in your browser. The first run installs missing packages automatically — this may take a few minutes.

---

### Alternative — download and run locally

1. Download the repository:
   - Click the green **Code** button on GitHub → **Download ZIP**
   - Unzip it anywhere on your computer

2. Open the `ScenithApp/` folder

3. Launch the app for your operating system:
   - **macOS:** double-click `Mac_run_app.command`
     *(first time: right-click → Open to bypass Gatekeeper)*
   - **Windows:** double-click `Win_run_app.bat`
   - **Any OS:** open a terminal in `ScenithApp/` and run `Rscript run_app.R`

The launcher installs any missing R packages automatically.

---

## Overview

[Scenith](https://www.scenith.com/) (Single Cell ENergetic meTabolism by profilIng tHe protein synthesis) profiles cellular metabolism by measuring **puromycin incorporation** (translation) as a proxy for ATP availability, under targeted metabolic inhibition:

| Condition | Inhibitor | What it blocks |
|-----------|-----------|----------------|
| **Co** | — | Control (DMSO) |
| **DG** | 2-Deoxy-D-glucose | Glycolysis |
| **O** | Oligomycin A | OXPHOS (ATP synthase) |
| **DGO** | DG + Oligomycin | Both pathways |

Four metabolic parameters are derived per genotype:

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
│   ├── Mac_run_app.command    # Double-click launcher for macOS
│   └── Win_run_app.bat        # Double-click launcher for Windows
├── analysis_v2.Rmd            # R Markdown report (scripted / offline use)
├── fcs-input/                 # FCS files — git-ignored
│   ├── exp1/
│   └── exp2/
├── .gitignore
└── README.md
```

> FCS files (`.fcs`), Excel files (`.xlsx`), Word documents, and FlowJo workspaces (`.wsp`) are listed in `.gitignore` and are not tracked by git.

---

## Metadata format

The metadata file tells the app which wells belong to which experimental group. It must be a **CSV** or **XLSX** file with at least these three columns:

| Column | Description | Example values |
|--------|-------------|----------------|
| `well_code` | Well position: row letter + zero-padded column number | `B01`, `C12`, `G04` |
| `genotype` | Cell genotype or cell line label | `WT`, `KO`, `WT_HSV` |
| `treatment` | Metabolic condition applied to the well | `Co`, `DG`, `O`, `DGO` |

Additional columns are allowed (e.g. `replicate`, `passage`) and will be carried through all tables.

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
C01,WT,DG
C02,WT,DG
D01,WT,O
E01,WT,DGO
G01,WT,UNST
```

### Treatment names for Scenith parameters

Scenith parameters are computed only when treatments spelled exactly **`Co`**, **`DG`**, **`O`**, and **`DGO`** are present. Additional treatments appear in QC and summary plots but not in the parameter calculations.

---

## Gating parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| FITC-A threshold | 4 000 | Live/Dead gate: keep cells **below** this value |
| APC-A threshold | 80 | Puromycin gate: keep cells **at or above** this value |

Adjust in the sidebar after reviewing the QC plots, then click **Run** again.

---

## App tabs

| Tab | What it shows |
|-----|---------------|
| **Overview** | Upload status and experiment description |
| **Metadata** | Uploaded metadata preview + per-FCS match status |
| **Sample mapping** | FCS files joined to metadata, with unmatched-well warnings |
| **Gating QC** | Singlet, live, and puromycin gate diagnostic plots |
| **Cell counts** | Cells retained at each gating step |
| **Puromycin summary** | Mean APC-A per sample (all live and puro+ cells) |
| **Scenith parameters** | Derived metabolic parameters, density plots, bar plots |

---

## Required R packages

`run_app.R` installs any missing packages automatically on launch.

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
Download the CSV template from the app sidebar, fill in the columns, and re-upload.

**"X sample(s) could not be matched"**
Check that FCS filenames contain a standard well code (`[A-H][0-9]{1,2}`) and that the metadata covers all expected wells.

**Scenith plots show "treatments required" message**
The treatments `Co`, `DG`, `O`, and `DGO` must be present and spelled exactly that way.

**App is slow on the all-samples FSC plot**
Keep that checkbox disabled until you specifically need it.

**Bioconductor packages fail to install**
Ensure R ≥ 4.2, then try `BiocManager::install(...)` manually (see above).

---

## Citation

If you use this pipeline, please cite the original Scenith method:

Argüello RJ, Combes AJ, Char R, Gigan JP, Baaziz AI, Bousiquot E, Camosseto V, Samad B, Tsui J, Yan P, Boissonneau S, Figarella-Branger D, Gatti E, Janssen E, Krummel MF, Pierre P.
**SCENITH: A Flow Cytometry-Based Method to Functionally Profile Energy Metabolism with Single-Cell Resolution.**
*Cell Metabolism.* 2020;32(6):1063–1075.e7.
DOI: [10.1016/j.cmet.2020.11.007](https://doi.org/10.1016/j.cmet.2020.11.007) · PMC: [PMC8407169](https://pmc.ncbi.nlm.nih.gov/articles/PMC8407169/)

---

## Authors

- **celbaek** — analysis pipeline and Shiny app
