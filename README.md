<p align="center">
  <img src="man/figures/logo_pandora_seom.png" alt="Pandora SEOM" width="320"/>
</p>

# pantheia_model

**PANTHEIA-SIRI prognostic model — metastatic pancreatic cancer**
*PANTHEIA-SEOM Research Group*

---

## Overview

`pantheia_model` accompanies the external-validation study of the PANTHEIA-SIRI
prognostic model for overall survival (OS) in metastatic pancreatic ductal
adenocarcinoma. It contains:

- the fitted models (Weibull accelerated failure time for OS and PFS, logistic
  for objective response) using the 3-level tumour-burden coding
  (>5 cm / <=5 cm / non-measurable);
- an interactive Shiny calculator, `pantheia()`;
- de-identified analysis datasets and a one-command script that reproduces
  **Table 2** (model coefficients) and **Table 3** (discrimination and
  calibration in the derivation and external-validation cohorts).

## Reproduce Table 2 and Table 3

```r
# from the package root
Rscript inst/reproduce_tables.R
# Expected (Table 3):
#   Derivation (n=593, pooled over 10 imputations)  C=0.654 (0.627-0.681)
#   Validation (n=62)                               C=0.603 (0.518-0.687)
```

## Reproduce the figures

```r
Rscript inst/reproduce_figures.R
# writes figures/Figure_1.pdf (validation KM + calibration),
#        figures/Supp_Figure_S1.pdf (death-time density),
#        figures/Supp_Figure_S2.pdf (variable clustering)
```

## Reproduce the supplementary tables

```r
Rscript inst/reproduce_supplementary_tables.R
# Supp Table S1 (collinearity: GVIF + redundancy) and
# Supp Table S2 (pairwise associations cachexia/ECOG/SIRI)
```

## Launch the calculator

```r
# remotes::install_github("albertocarm/pantheia_model")
pantheia_model::pantheia()
```

## Data

The de-identified datasets in `inst/extdata/` contain only the model variables
(`os_time, os_event, logsiri, diam3, regimen_cat, ecog_cat_3, CACS`); no patient
codes, centres, dates or other identifiers are included. `derivation_os_imputed.csv`
holds the 10 multiply-imputed datasets of the survival-analysis population (n=593,
`imp` column 1-10), so derivation performance is reproduced over the imputations and
pooled with Rubin's rules rather than on complete cases only; `validation_os.csv` is
the independent external cohort (n=62). "CACS" is the baseline
symptom composite (anorexia, cachexia, asthenia, or weight loss >5%).

## Transparency

The source code, fitted model objects and de-identified data are released so that
every reported coefficient and performance metric can be independently reproduced,
in line with the STROBE/REMARK and TRIPOD reporting standards of the
PANTHEIA-SEOM study.
