<p align="center">
  <img src="man/figures/logo_pandora_seom.png" alt="Pandora SEOM" width="320"/>
</p>

# pantheia_model

**PANTHEIA-SIRI prognostic model — metastatic pancreatic cancer**
*PANTHEIA-SEOM Research Group*

---

## Overview

`pantheia_model` accompanies the external-validation study of the **PANTHEIA-SIRI**
prognostic model for overall survival (OS) in metastatic pancreatic ductal
adenocarcinoma. It bundles, in one place, everything needed to **use** the model and
to **reproduce** the paper:

- the fitted models (Weibull accelerated failure time for OS and PFS, logistic for
  objective response) using the 3-level tumour-burden coding (>5 cm / ≤5 cm /
  non-measurable);
- an interactive Shiny **calculator** — the function `pantheia()`;
- de-identified analysis datasets and three one-command scripts that regenerate the
  paper's **tables**, **figures** and **supplementary tables**.

## Installation

```r
# install.packages("remotes")
remotes::install_github("albertocarm/pantheia_model")
```

The reproduction scripts additionally use a few packages (declared in `Suggests`):

```r
install.packages(c("survminer", "cowplot", "dplyr", "Hmisc", "rms", "car"))
```

## 1. Run the prognostic calculator

`pantheia()` launches the Shiny app in your browser.

```r
library(pantheia_model)
pantheia()                      # opens the calculator locally

# non-interactive options are passed through to shiny::runApp(), e.g.
pantheia(launch.browser = FALSE, port = 8080)
```

## 2. Reproduce the paper

Each script can be run **two ways** — from an installed package, or from a clone of
this repository.

| Script | Reproduces |
|--------|------------|
| `reproduce_tables.R` | Table 2 (coefficients) + Table 3 (discrimination/calibration) |
| `reproduce_figures.R` | Figure 1 + Supplementary Figures S1–S2 |
| `reproduce_supplementary_tables.R` | Supplementary Tables S1–S2 |

**From the installed package** (any working directory):

```r
source(system.file("reproduce_tables.R",               package = "pantheia_model"))
source(system.file("reproduce_figures.R",              package = "pantheia_model"))
source(system.file("reproduce_supplementary_tables.R", package = "pantheia_model"))
```

**From a clone of the repo** (run from the package root):

```sh
Rscript inst/reproduce_tables.R
Rscript inst/reproduce_figures.R
Rscript inst/reproduce_supplementary_tables.R
```

### Expected output

```
# reproduce_tables.R
Table 2:  (Intercept) 3.805 ... Weibull sigma = 0.874
Table 3:  Derivation (n=593, pooled over 10 imputations)  C = 0.654 (0.627-0.681)
          Validation (n=62)                               C = 0.603 (0.518-0.687)

# reproduce_figures.R  -> writes to ./figures/
Figure_1.pdf         validation OS by risk group (KM, panels A-B) + calibration (C-D)
Supp_Figure_S1.pdf   probability density of death times (validation cohort)
Supp_Figure_S2.pdf   variable clustering of the model predictors (derivation)

# reproduce_supplementary_tables.R
Supp Table S1  GVIF and redundancy  (all GVIF^(1/2df) < 1.06; no redundant variables)
Supp Table S2  pairwise associations among cachexia, ECOG and SIRI
```

## Data

The de-identified datasets in `inst/extdata/` contain **only the model variables**
(`os_time, os_event, logsiri, diam3, regimen_cat, ecog_cat_3, CACS`); no patient
codes, centres, dates or other identifiers are included.

- `derivation_os_imputed.csv` — the 10 multiply-imputed datasets of the
  survival-analysis population (n = 593; column `imp` = 1–10). Derivation performance
  is therefore reproduced **over the imputations and pooled with Rubin's rules**, not
  on complete cases only.
- `validation_os.csv` — the independent external-validation cohort (n = 62).

`CACS` is the baseline symptom composite (anorexia, cachexia, asthenia, or
weight loss > 5%).

## Transparency

The source code, fitted model objects and de-identified data are released so that
every reported coefficient and performance metric can be independently reproduced,
in line with the STROBE/REMARK and TRIPOD reporting standards of the PANTHEIA-SEOM
study.
