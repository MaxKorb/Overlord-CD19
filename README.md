# B-cell depletion analysis in the OVERLORD-MS trial

This repository contains the R code used for the analyses and figure generation for the study of peripheral B-cell depletion in the OVERLORD-MS laboratory sub-study.

The analyses compare B-cell depletion patterns in participants treated with rituximab or ocrelizumab and evaluate baseline characteristics associated with depletion.

## Analysis

The analysis script (`01_analysis.R`) includes code for:

- Definition of the analysis cohort
- Classification of sustained B-cell depletion
- Classification of early depletion at month 3
- Sensitivity analysis using an ≥80% depletion criterion
- Descriptive analyses of participant characteristics
- Assessment of B-cell depletion across follow-up visits
- Assessment of B-cell reappearance after initial depletion
- Statistical comparison of baseline characteristics between depletion groups
- Effect size estimation
- Correction for multiple testing using the false discovery rate (FDR) and Bonferroni methods
- Generation of Table 1 and Supplementary Tables 1–4

The figure script (`02_figures.R`) includes code for generation of the manuscript and supplementary figures.

B-cell depletion was defined using peripheral CD19+ B-cell counts, with a threshold of ≤5 cells/µL.

## Repository structure

```text
OVERLORD_B_cell_analysis/
├── scripts/
│   ├── 01_analysis.R
│   └── 02_figures.R
├── README.md
├── .gitignore
└── OVERLORD_B_cell_analysis.Rproj
```

## Data

Participant-level data are not included in this repository.

To reproduce the analyses, the analysis dataset should be stored locally as:

```text
data/OVERLORD_B_cell_analysis_data.csv
```

The expected local project structure is therefore:

```text
OVERLORD_B_cell_analysis/
├── data/
│   └── OVERLORD_B_cell_analysis_data.csv
├── scripts/
│   ├── 01_analysis.R
│   └── 02_figures.R
├── output/
├── README.md
├── .gitignore
└── OVERLORD_B_cell_analysis.Rproj
```

The `data/` directory and its contents should not be committed to the repository. Generated files in the `output/` directory are also excluded from version control.

## Requirements

The analyses and figures were generated in R using the packages required by the two scripts, including:

- `dplyr`
- `openxlsx`
- `rstatix`
- `rcompanion`

Additional packages used for figure generation are specified in `02_figures.R`.

Package and R version information for the final analyses is provided by `sessionInfo()` in the scripts.

## Running the analysis

Set the working directory to the root of the project.

Run the analysis and generate the manuscript tables with:

```r
source("scripts/01_analysis.R")
```

Generate the manuscript and supplementary figures with:

```r
source("scripts/02_figures.R")
```

Both scripts read the analysis dataset from the local `data/` directory and save generated files in the `output/` directory.

## Output

The analysis script generates:

- Table 1: Participant characteristics
- Supplementary Table 1
- Supplementary Table 2
- Supplementary Table 3
- Supplementary Table 4

Tables are exported in both `.csv` and `.xlsx` formats.

The figure script generates:

- Figure 1A: Study design
- Figure 2: CD19+ B-cell counts over time
- Figure 3: Proportion of participants with B-cell depletion over time
- Figure 4: CD19+ B-cell counts by treatment over time
- Supplementary Figure 1: Participants without sustained depletion
- Supplementary Figure 2: Participants without early depletion

Figures are exported in both `.pdf` and `.tiff` formats.

## Contributors

- Hannah Sofie Sjo – statistical analyses, generation of manuscript tables and figures, and preparation of the analysis and figure scripts.
- Max Korbmacher – statistical and methodological guidance, contribution to the analysis code, and preparation of the GitHub repository.

## Study

The analyses are based on data from the OVERLORD-MS trial (ClinicalTrials.gov identifier: NCT04578639).