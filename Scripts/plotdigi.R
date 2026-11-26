library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

PDval <- read_csv("Data/PDval.csv", show_col_types = FALSE)

tol_mhs <- 1     # a reading within ±1 percentage point of the paper counts as agreeing (adjust)
tol_sd  <- 10    # an SD within ±10% (relative) of the paper counts as agreeing (adjust)

# ---- 1. Controls: one row per distinct control per study (shared controls are repeated across treatment rows) ----
ctrl <- PDval %>%
  distinct(Study_ID, Reference, Paper_Control_MHS, Paper_Control_SD, PD_Control_MHS, PD_Control_SD)

# Studies with more than one distinct control row: either genuinely different controls or inconsistent digitising -> inspect
ctrl %>% count(Study_ID) %>% filter(n > 1)

# ---- 2. Long format: one row per comparison ----
long <- bind_rows(
  ctrl %>% transmute(Study_ID, Reference, group = "Control",   stat = "MHS", paper = Paper_Control_MHS, pd = PD_Control_MHS),
  ctrl %>% transmute(Study_ID, Reference, group = "Control",   stat = "SD",  paper = Paper_Control_SD,  pd = PD_Control_SD),
  PDval %>% transmute(Study_ID, Reference, group = "Treatment", stat = "MHS", paper = Paper_Treat_MHS,  pd = PD_Treat_MHS),
  PDval %>% transmute(Study_ID, Reference, group = "Treatment", stat = "SD",  paper = Paper_Treat_SD,   pd = PD_Treat_SD)
) %>%
  filter(!is.na(paper), !is.na(pd)) %>%                      # only pairs where both values exist
  mutate(diff = pd - paper,                                  # PlotDigitiser minus paper
         rel  = if_else(paper != 0, 100 * diff / paper, NA_real_),
         agree = if_else(stat == "MHS", abs(diff) <= tol_mhs, abs(rel) <= tol_sd))

# ---- 3. Agreement statistics ----
ccc <- function(x, y) 2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)   # Lin's concordance correlation

summ <- long %>% group_by(group, stat) %>%
  summarise(n = n(),
            bias            = mean(diff),                       # mean (PD - paper): systematic over- or under-reading
            sd_diff         = sd(diff),
            mean_abs_err    = mean(abs(diff)),
            median_abs_err  = median(abs(diff)),
            max_abs_err     = max(abs(diff)),
            median_rel_pct  = median(abs(rel), na.rm = TRUE),
            r               = cor(paper, pd),
            CCC             = ccc(paper, pd),
            pct_agree       = 100 * mean(agree),
            p_wilcoxon      = suppressWarnings(wilcox.test(pd, paper, paired = TRUE)$p.value),
            .groups = "drop")
print(summ, width = Inf)

# Largest discrepancies, to check by hand against the figures
long %>% arrange(desc(abs(if_else(stat == "MHS", diff, rel)))) %>%
  filter(!agree) %>% select(Study_ID, Reference, group, stat, paper, pd, diff, rel) %>% print(n = 30)

# ---- 4. Plots ----
p_scatter <- ggplot(long, aes(paper, pd, colour = group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.7) + facet_wrap(~ stat, scales = "free") +
  labs(x = "Value reported in paper", y = "PlotDigitiser reading", colour = NULL) + theme_classic()

p_ba <- long %>% group_by(stat) %>% mutate(m = (paper + pd) / 2, lo = mean(diff) - 1.96 * sd(diff), hi = mean(diff) + 1.96 * sd(diff)) %>%
  ggplot(aes(m, diff, colour = group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(aes(yintercept = lo), colour = "grey60") + geom_hline(aes(yintercept = hi), colour = "grey60") +
  geom_point(alpha = 0.7) + facet_wrap(~ stat, scales = "free") +
  labs(x = "Mean of paper and PlotDigitiser", y = "PlotDigitiser minus paper", colour = NULL) + theme_classic()

(p_scatter / p_ba) + plot_layout(guides = "collect")

# ---- 5. Does digitising change the effect size that goes into the meta-analysis? ----
fl <- 0.5   # same floor as the main analysis (percent)
eff <- PDval %>%
  filter(if_all(c(Paper_Control_MHS, Paper_Treat_MHS, PD_Control_MHS, PD_Treat_MHS), ~ !is.na(.x))) %>%
  mutate(lnRR_paper = log(pmax(Paper_Treat_MHS, fl) / pmax(Paper_Control_MHS, fl)),
         lnRR_pd    = log(pmax(PD_Treat_MHS,    fl) / pmax(PD_Control_MHS,    fl)),
         d_lnRR     = lnRR_pd - lnRR_paper)
c(mean_diff = mean(eff$d_lnRR), sd_diff = sd(eff$d_lnRR), r = cor(eff$lnRR_paper, eff$lnRR_pd), CCC = ccc(eff$lnRR_paper, eff$lnRR_pd))
eff %>% arrange(desc(abs(d_lnRR))) %>% select(Study_ID, Reference, lnRR_paper, lnRR_pd, d_lnRR) %>% head(15)
