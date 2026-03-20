# Analysis of methacholine dose–response curves and rhythmic airway resistance profiles

![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Last Commit](https://img.shields.io/github/last-commit/robertmaidstone/asthma-shiftwork-ukb-female)

This repository contains an analysis pipeline for a respiratory study investigating how airway resistance varies across time, genotype, and allergen exposure. Airway resistance is measured using two different experimental techniques; DCP (measuring sRAW, sGAW and EF50) and Flexivent. Analysis focuses on comparison of methacholine dose–response curves and the temporal rhythmicity of derived resistance metrics. Analytical techniques used include basic statistical tests, non-linear modelling, harmonic regression, and production of publication ready figures.

---

## Overview

The study compares airway resistance in:
- Wild‑type (WT) vs Rev‑erbα knockout (KO) mice
- PBS vs House Dust Mite (HDM) exposure
- Multiple time points across the circadian cycle
Two major analytical components are included:
- Methacholine dose–response curve analysis
- Rhythmic modelling of summary resistance metrics over time
All analyses are performed in R using reproducible scripts contained in this repository.

## Methacholine Dose–Response Analysis
### Pairwise comparisons
Methacholine dose–response curves were compared pairwise across time points using a Mann–Whitney U test.
### Curve modelling
Individual mouse dose–response curves were modelled using the drc R package [1].
The slope parameter from each fitted curve was extracted and compared using a two‑way ANOVA (genotype × treatment) via the stats package [2].

## Summary Metrics and Rhythmic Modelling
Two summary statistics were derived from each methacholine curve:
- Area Under the Curve (AUC)
- Maximum resistance
These metrics were modelled over time as a sine wave with unknown:
- amplitude
- phase
- mean value
Fitting was performed using nonlinear least squares (nls, stats package [2]).

### Model comparison
To determine whether rhythmicity was present, each sine model was compared with a constant model using a likelihood ratio test (lmtest package [3]).
- Significant rhythms (p < 0.05) → plotted with solid lines
- Non‑significant rhythms → plotted with dashed line

### Joint modelling across conditions
To test for differences between WT/KO and PBS/HDM groups, data were modelled jointly with additional parameters for:
- genotype‑specific amplitude, phase, mean
- treatment‑specific amplitude, phase, mean
Significance of these added parameters (p < 0.05) indicates group‑level differences in rhythmic profiles

---

## Requirements

### **Software**
- **R ≥ 4.2**

### **Key Packages**
- `drc` – dose–response curve fitting
- `stats` – ANOVA, nonlinear least squares
- `lme4` - model fitting
- `lmtest` – likelihood ratio testing
- `tidyverse` – data manipulation and plotting
- `openxlsx` - reading in excel files
- `grid`/`patchwork`/`ggtext` - plotting
- `rstatix` -
- `purrr` -

## Running the Analysis
- Clone the repository
- Set the working directory
- Inputs from data/ folder
- Run an _analysis.R file (DCP and flexivent examples given) which will use functions from functions.R to run the analysis, some data manipulating may be needed.
- Outputs will be written to plots/
- 
## Citation
If using this code please cite the repository.

## Contributing
Contributions, suggestions, and extensions are welcome.
Please open an issue or submit a pull request.
## References

1.	Ritz, C., Baty, F., Streibig, J. C., Gerhard, D. (2015) Dose-Response Analysis Using R PLOS ONE,  10(12), e0146021
2.	R Core Team (2024). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. 
3.	Achim Zeileis, Torsten Hothorn (2002). Diagnostic Checking in Regression Relationships. R News 2(3), 7-10. URL https://CRAN.R-project.org/doc/Rnews/
