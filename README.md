<p align="center">
  <img src="scenithr_logo.png" alt="scenithR" width="220"/>
</p>
<p align="center">
  <em>A reproducible analysis pipeline for Scenith metabolic assays using flow cytometry.</em>
</p>

---

## What you need to run the app

| Input | Format | Required | Notes |
|-------|--------|----------|-------|
| FCS files | `.fcs` | ✓ | All files from one experiment; select them all at once |
| Plate metadata | `.csv` or `.xlsx` | ✓ | Maps each well to genotype and perturbation (see below) |
| R ≥ 4.2 | — | ✓ | Packages installed automatically on first run |
| Panel/channel info | — | — | Choose a preset in the app; no file needed unless custom |

### Metadata columns

| Column | Required | Description | Example values |
|--------|----------|-------------|----------------|
| `well_code` | ✓ | Well position — letter + zero-padded number | `B01`, `C12`, `G04` |
| `genotype` | ✓ | Cell line or genotype label | `WT`, `KO`, `WT_HSV` |
| `perturbation` | ✓ | Metabolic probe condition | `Co`, `DG`, `O`, `DGO`, `UNST` |
| `treatment` | optional | Broader experimental condition | `fed`, `starved`, `drug_X` |
| `time` | optional | Timepoint | `0h`, `4h`, `24h` |

Extra columns (e.g. `replicate`, `passage`) are allowed and carried through all tables.

> **Download a pre-filled template** from the app sidebar — it covers all rows B–H × columns 1–12. Fill in `genotype` and `perturbation` and save as CSV.

### FCS file naming

The app extracts the well position from the FCS filename using the pattern `[A-H][0-9]{1,2}`. Standard instrument-exported names work automatically:

```
Specimen_004_B1_B01.fcs   →  well B01
Specimen_002_C12_C12.fcs  →  well C12
```

If your filenames do not contain a well code in this format, the files will still load but will not match the metadata.

### Perturbation labels for Scenith parameters

Four specific labels trigger the metabolic parameter calculations. Use these exact spellings in the `perturbation` column:

| Label | Condition |
|-------|-----------|
| `Co` | Control (DMSO) |
| `DG` | 2-Deoxy-D-glucose — glycolysis inhibitor |
| `O` | Oligomycin A — OXPHOS / ATP synthase inhibitor |
| `DGO` | DG + Oligomycin — combined inhibition |

Other perturbation labels (e.g. `UNST`, `DMSO_25uL`) are kept in QC plots but excluded from Scenith parameter calculations.

### Channel configuration

Configure which fluorescence channel plays which role in the **Channels** tab. Select a preset or choose Custom to assign channels from your FCS file.

| Role | Default | What it is |
|------|---------|------------|
| Scatter X | `FSC.A` | X-axis for singlet gate |
| Scatter Y | `FSC.H` | Y-axis for singlet gate |
| Live/Dead | `FITC.A` | Dead cell exclusion channel (set to None to skip) |
| Signal | `APC.A` | Puromycin readout channel |

Available presets: APC/FITC, APC/BV421, PE/FITC, and "no live/dead stain".

---

## Quick start

### Easiest — run directly from R (no download needed)

Paste this into your R console:

```r
shiny::runGitHub("scenithR", "camillaelbaek", subdir = "ScenithApp")
```

R fetches the app from GitHub and opens it in your browser. The first run installs missing packages automatically — this may take a few minutes.

---

### Alternative — download and run locally

1. Click the green **Code** button on GitHub → **Download ZIP**, then unzip it anywhere.
2. Open the `ScenithApp/` folder.
3. Launch for your OS:
   - **macOS:** double-click `Mac_run_app.command` *(first time: right-click → Open)*
   - **Windows:** double-click `Win_run_app.bat`
   - **Any OS:** run `Rscript run_app.R` in a terminal inside `ScenithApp/`

---

## Workflow in the app

1. Upload FCS files and metadata in the sidebar
2. **Channels tab** — select your panel preset (or configure custom channels)
3. **Gate 1: Singlets** — adjust the polygon gate; preview updates live
4. **Gate 2: Live/Dead** — set the threshold; skipped if no live/dead channel configured
5. **Gate 3: Signal** — set the minimum puromycin signal threshold
6. Click **Run analysis** → results appear in the remaining tabs

---

## Overview of the Scenith method

[Scenith](https://www.scenith.com/) profiles cellular metabolism by measuring **puromycin incorporation** (translation) as a proxy for ATP availability under targeted metabolic inhibition. Four parameters are derived per group:

1. **Glucose dependence** = 100 × (Co − DG) / (Co − DGO)
2. **FAO/AAO capacity** = 100 − glucose dependence
3. **Mitochondrial dependence** = 100 × (Co − O) / (Co − DGO)
4. **Glycolytic capacity** = 100 − mitochondrial dependence

When `treatment` and/or `time` columns are present in the metadata with more than one unique value, parameters are computed per **genotype × treatment × time** group.

---

## Repository structure

```
.
├── ScenithApp/
│   ├── app.R                  # Shiny app
│   ├── run_app.R              # Launcher: installs packages + starts app
│   ├── Mac_run_app.command    # macOS double-click launcher
│   └── Win_run_app.bat        # Windows double-click launcher
├── analysis_v2.Rmd            # R Markdown report (offline / scripted use)
├── fcs-input/                 # FCS files — git-ignored
├── .gitignore
└── README.md
```

> `.fcs`, `.xlsx`, `.wsp`, and `.docx` files are in `.gitignore` and not tracked by git.

---

## App tabs

| Tab | What it shows |
|-----|---------------|
| **Overview** | Status checklist and workflow guide |
| **Metadata** | Uploaded metadata preview + per-FCS match status (green/red) |
| **Channels** | Panel preset or custom channel assignment |
| **Gate 1: Singlets** | Polygon gate on scatter channels — 4 adjustable vertices |
| **Gate 2: Live/Dead** | Threshold gate on live/dead channel |
| **Gate 3: Signal** | Minimum threshold on signal channel |
| **Cell counts** | Cells retained at each gating step |
| **Summary** | Mean signal per sample (all live and signal-positive cells) |
| **Scenith parameters** | Derived parameters, density plots, and bar plots |

---

## Required R packages

Installed automatically by `run_app.R`.

**CRAN:** `shiny`, `dplyr`, `tidyr`, `ggplot2`, `DT`, `stringr`, `purrr`, `scales`, `readr`, `readxl`, `ggridges`, `ggpubr`, `viridis`, `ggbeeswarm`, `sp`

**Bioconductor** (via `BiocManager`): `flowCore`, `flowViz`, `ggcyto`, `openCyto`

If Bioconductor packages fail to install automatically:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("flowCore", "flowViz", "ggcyto", "openCyto"))
```

Requires **R ≥ 4.2**.

---

## Troubleshooting

**"Missing required column(s): genotype, perturbation"**
Download the CSV template from the sidebar, fill in the columns, and re-upload.

**"X sample(s) could not be matched"**
Check that FCS filenames contain a well code matching `[A-H][0-9]{1,2}` and that the metadata covers all those wells.

**Scenith plots show "perturbations required"**
The perturbations `Co`, `DG`, `O`, and `DGO` must all be present and spelled exactly that way.

**Gate preview is blank**
Upload FCS files and configure channels before using the gating tabs.

**Bioconductor packages fail to install**
Ensure R ≥ 4.2 and try `BiocManager::install(...)` manually (see above).

---

## Citation

### Citing this software

A `CITATION.cff` file is included in the root of this repository. GitHub surfaces it automatically as a **"Cite this repository"** button in the sidebar — click it to export a ready-made citation in APA, BibTeX, or other formats.

To cite manually:

> Elbaek C. (2025). *scenithR: a Shiny-based analysis pipeline for Scenith flow cytometry assays* (Version 1.0.0). GitHub. https://github.com/camillaelbaek/scenithR

### Citing the Scenith method

If you use the metabolic profiling approach, please also cite the original method paper:

Argüello RJ, Combes AJ, Char R, Gigan JP, Baaziz AI, Bousiquot E, Camosseto V, Samad B, Tsui J, Yan P, Boissonneau S, Figarella-Branger D, Gatti E, Janssen E, Krummel MF, Pierre P.
**SCENITH: A Flow Cytometry-Based Method to Functionally Profile Energy Metabolism with Single-Cell Resolution.**
*Cell Metabolism.* 2020;32(6):1063–1075.e7.
DOI: [10.1016/j.cmet.2020.11.007](https://doi.org/10.1016/j.cmet.2020.11.007) · PMC: [PMC8407169](https://pmc.ncbi.nlm.nih.gov/articles/PMC8407169/)

---

## Authors

- **celbaek** — analysis pipeline and Shiny app
