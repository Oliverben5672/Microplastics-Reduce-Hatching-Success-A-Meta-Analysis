# Microplastics Reduce Hatching Success in Oviparous Taxa: A Meta-Analysis

Data and R code accompanying:

> Warren & Hemmings (submitted for review). Microplastics Reduce Hatching Success in Oviparous Taxa: A Meta-Analysis.

**Version:** v2.0 (final analysis). 

The original analysis (29 studies, Hedges' g only) is archived at release [`v1.0-original-analysis`](https://github.com/Oliverben5672/Microplastics-Reduce-Hatching-Success-A-Meta-Analysis/releases/tag/v1.0-original-analysis) and is superseded by the code and data in this repository.

## Summary of the analysis

A multilevel meta-analysis of 39 studies (150 records), filtered using the CRED reliability framework. The primary dataset is pure microplastic exposures (102 records, 32 studies); records with a co-contaminant (secondary dataset) are analysed separately. The primary effect size is the log response ratio (lnRR, reported as % change); Hedges' g is a secondary scale. Models are multilevel (`rma.mv`) with study and record random effects, with species and phylogenetic structures as sensitivity analyses. Inference uses cluster-robust (CR2) variance estimation, reported only where df ≥ 4. Moderator, dose-response, sensitivity and publication-bias analyses follow.

## Repository structure

- `Data/` — dataset used in the meta-analysis [file name], and the PlotDigitiser validation data
- `Scripts/` — R scripts for all analyses, tables and figures
- `results_39studies_v5/` — saved outputs from the final run (results tables as CSV, `CRED_results.rds`, figures)

## Reproducing the analysis

[package list unchanged]

Run `CRED_39studies_v5.R` first: it reads from `Data/` and writes to `results_39studies_v5/`, and the other scripts read from that folder.

All analyses were conducted in R. Key packages:

- `metafor` — multilevel meta-analytic models
- `clubSandwich` — cluster-robust (CR2) inference
- `rotl` — Open Tree of Life species matching
- `ape` — phylogenetic tree manipulation
- `phytools` — phylogenetic signal estimation
- `dplyr`, `tidyverse` — data wrangling and plotting
- `kableExtra` — producing output tables

Run `CRED_39studies_v5.R` first, because it writes the results used by the other scripts.

Scripts:

- `CRED_39studies_v5.R` — main analysis: data preparation, primary and secondary models, moderators, dose-response, sensitivity analyses, plot-digitiser validation and figures
- `Publication_Bias_v5.R` — publication bias: multilevel Egger-type test, study-level regression, time-lag test and funnel plots
- `plotdigi.R` — Validation of values read using PlotDigitiser
- `Table_S2.R` — supplementary table of all model structures across all dataset configurations

## Contact

Oliver Warren — oliverben5672@gmail.com