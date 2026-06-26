# Microplastic effects on aquatic species hatching success — phylogenetic meta-analysis

Data and R code accompanying:

> Warren O. et al. (submitted for review). Microplastics Reduce Hatching Success in Oviparous Taxa: A Meta-Analysis.

## Repository structure

- `Data/` — dataset used in the meta-analysis
- `Scripts/` — R scripts for phylogenetic meta-analysis, sensitivity analyses, and figures
- `Outputs/` — saved outputs

## Reproducing the analysis

All analyses were conducted in R. Key packages:

- `metafor` — meta-analytic models
- `rotl` — Open Tree of Life species matching
- `ape` — phylogenetic tree manipulation
- `phytools` — phylogenetic signal estimation
- `dplyr` — data wrangling
- `tidyverse` — base data manipulation and plotting
- `ggalt` — dumbbell plots
- `ggthemes` — editing graph visuals
- `kableExtra` — producing output tables

Scripts include:

- `Simple Model.R` - Runs the simplest model and visualises all studies raw findings
- `Phylogenetic Model.R` - Runs the phylogenetic model by creating a correlation matrix
- `Taxonomic Model.R` - Runs the taxonomically structured model by including class as a moderator after filtering to classes with at least 3 observations
- `Publication Bias Assessment.R` - Assess publication bias via funnel plots and trim and fill
- `Test for Outliers.R` - Identifies outliers in each core model using Cook's distance and studentised residuals
- `Results Summary Table.R` - Produces a clean table with each core result
- `All Model Summary with Exclusions.R` - Produces a large table comparine all core models in full, with outliers removed, with co-contaminant studies removed, and with outliers and co-contaminant studies removed

## Contact

Oliver Warren — oliverben5672@gmail.com