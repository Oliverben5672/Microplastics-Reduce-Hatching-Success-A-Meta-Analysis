# =============================================================================
# CRED meta-analysis: microplastic exposure and hatching success  (39-study update, v5)
#   PRIMARY   = pure plastic only, CRED-filtered (R3 removed)
#   SECONDARY = co-contaminant rows included (contam tested as a moderator)
#   Effect size: lnRR (log response ratio, reported as % change in hatching) is PRIMARY;
#                Hedges' g (negative = reduced hatching) is the labelled SECONDARY scale.
#                g standardises by the SD between replicates, which differs by replicate unit
#                (single embryo vs tank), so it is not on a common scale across studies.
#
# [V4] marks a change made for v4: lnRR primary; sample-size-based sampling covariances (no fixed rho);
#      robust CIs/p only where Satterthwaite df >= 4 (model-based p always shown too); outliers defined on
#      each scale's own model; imputed SDs capped at sqrt(m (100 - m)); lnRR floor and zero-hatching checks.
# [V5] marks a change made for v5: the SECONDARY analysis (and the contaminant contrast) now uses a separate observation-level
#      variance for pure-plastic and co-contaminant records (fit_cred(resid = "auto")); the common-variance model is kept as a
#      labelled sensitivity. Table 2 reports by-type effects from that model instead of one pooled "all records" mean.
#      Also new: single-polymer PS vs PE test (18b), contrast under each lnRR floor x structure x residual variance.
# Sections 1-18 keep their numbers and headings.
# =============================================================================

library(tidyverse)
library(janitor)
library(metafor)
library(clubSandwich)   # CR2 small-sample robust inference
library(rotl)           # Open Tree of Life
library(ape)            # masks dplyr::where (section 1 uses tidyselect::where)
library(phytools)       # masks purrr::map (purrr::map_*() is always called with its prefix)
library(kableExtra)     # [V3] moved up from section 15
library(patchwork)      # [V3] moved up from section 15

out_dir <- "results_39studies_v5"    # [V5] new folder, so earlier outputs are not overwritten
df_min  <- 4                          # [V4] robust (CR2) CIs and p are reported only where Satterthwaite df >= df_min
secondary_resid <- "by_contam"        # [V5] datasets containing both contaminant types: "by_contam" = separate observation-level variance per type
#      (default); "common" = the earlier single variance (kept as a labelled sensitivity)
min_grp <- 5                          # [V5] fall back to the common variance if either contaminant type has fewer records than this
dir.create(out_dir, showWarnings = FALSE)

CRED_FINAL <- read_csv("Data/CRED_FINAL.csv")
# View(CRED_FINAL)

# 1. Clean names, drop empty Excel rows, tidy text -------------------------
dat <- CRED_FINAL %>%
  clean_names() %>%
  filter(!is.na(study_id)) %>%
  rename(contam    = contaminant_category,
         pathway   = exposure_pathway,
         lifestage = exposure_lifestage) %>%
  select(-any_of(c("x25", "x26", "x27"))) %>%
  mutate(across(tidyselect::where(is.character),
                ~ .x %>% str_replace_all("\u00a0", " ") %>%
                  str_replace_all("\u00c2", "") %>% str_squish()),
         es_id = row_number())

# 1b. Data validation --------------------------------------------------------
stopifnot(all(dat$cred_r %in% c("R1", "R2", "R3")))    # an R4 would otherwise pass filter(cred_r != "R3")

known <- list(
  taxa_class   = c("Actinopterygii","Amphibia","Crustacea","Mollusca","Nematoda","Oligochaeta"),
  polymer      = c("PS","LDPE","PE","CTR","Mix","HDPE","PA","AFH","PVC","PES","PET","PHB","PLA","PP","RFP"),
  concentration_unit = c("mg/L","ug/L","g/kg","MP/ml","mg MP/L","fibres/fish/day","mg/kg","% sed dw","MP/L","g/L"),
  pathway      = c("Acute_solution","Chronic_sediment","Chronic_solution","Ingestion",
                   "Ingestion + Chronic_sediment","Ingestion + Chronic_solution","Short_Term"),
  lifestage    = c("Combined","Embryonic","Larval","Maternal","Parental"),
  data_source  = c("PlotDigitiser","Author Correspondance","Primary Paper"),
  contam       = c("Pure Plastic","Plastic+Co-Contaminant","Plastic+Co-contaminant"),
  env_relevant = c("Yes","No"))
new_levels <- purrr::imap_dfr(known, \(lv, col) tibble(column = col, new_level = setdiff(unique(dat[[col]]), lv))) %>%
  filter(!is.na(new_level))
print(new_levels)
stopifnot(nrow(new_levels) == 0)     # [V3] restored: a new label stops the script until it has a decision

# impossible values, and hatching entered as a proportion (the lnRR floor of 0.5 assumes percent)
dat %>% filter(control_mhs < 0 | control_mhs > 100 | treat_mhs < 0 | treat_mhs > 100 |
                 control_sd < 0 | treat_sd < 0 | control_reps < 1 | treat_reps < 1) %>% print()
dat %>% group_by(study_id) %>% filter(max(c(control_mhs, treat_mhs)) <= 1) %>% distinct(study_id, reference) %>% print()

# 2. Exclusions and SD imputation -------------------------------------------
impute_ids <- c(3, 27, 38, 39)     # studies reporting no SDs (flagged IMPUTE in the spreadsheet)

# Studies whose replicates are single embryos (binary outcome per replicate): the SD follows from the mean,
# 100 * sqrt(p (1 - p) n / (n - 1)). Study 38 (Prata et al. 2022: n = 10 wells per treatment, 20 control wells) is the known case.
# [V4] If the spreadsheet has a replicate_unit column ("embryo" / "tank"), it decides; otherwise Study 38 is used.
binomial_ids <- if ("replicate_unit" %in% names(dat)) {
  dat %>% filter(str_to_lower(replicate_unit) %in% c("embryo", "individual")) %>% pull(study_id) %>% unique()
} else 38
bin_sd <- \(m, n) 100 * sqrt(pmin(pmax(m / 100, 0), 1) * (1 - pmin(pmax(m / 100, 0), 1)) * n / (n - 1))
# [V4] maximum SD of a variable bounded at 0-100 with mean m (Bhatia-Davis); imputed class-CV SDs are capped at this
sd_cap <- \(m) sqrt(pmax(m, 0) * (100 - pmin(m, 100)))

# arm-level rule: a record needs >= 2 replicates in both groups
arm_log <- dat %>%
  mutate(drop_arm = control_reps < 2 | treat_reps < 2) %>%
  group_by(study_id, reference) %>%
  summarise(arms = n(), arms_dropped = sum(drop_arm), .groups = "drop") %>%
  filter(arms_dropped > 0) %>% mutate(study_lost = arms_dropped == arms)
print(arm_log)

n_raw <- c(records = nrow(dat), studies = n_distinct(dat$study_id))
dat <- dat %>% filter(control_reps >= 2, treat_reps >= 2)
stopifnot(!any(dat$control_reps < 2 | dat$treat_reps < 2))
rbind(raw = n_raw, after_arm_exclusion = c(records = nrow(dat), studies = n_distinct(dat$study_id)))

# [V4] control groups (records sharing a control share a sampling covariance). Uses a control_group_id column when the
# spreadsheet has one; otherwise a group is inferred from the study and the control mean, SD and replicates.
has_cg <- "control_group_id" %in% names(dat)
if (has_cg) stopifnot(!anyNA(dat$control_group_id))
dat <- dat %>% mutate(control_id = if (has_cg) paste(study_id, control_group_id, sep = "_")
                      else paste(study_id, control_mhs, control_sd, control_reps, sep = "_"))
message("control groups taken from ", if (has_cg) "the control_group_id column" else "control mean/SD/replicates (inferred)")

# [V4] validation: an observed SD cannot exceed the maximum possible for a 0-100 variable (Bhatia-Davis bound, times
# sqrt(n/(n-1)) because these are sample SDs). Anything listed here is a possible data-entry error.
sd_cap_s <- \(m, n) sd_cap(m) * sqrt(n / (n - 1))
tab_sd_cap_check <- bind_rows(
  dat %>% filter(!is.na(control_sd), control_sd > sd_cap_s(control_mhs, control_reps) + 0.5) %>%
    transmute(study_id, reference, group = "control", mean = control_mhs, sd = control_sd, max_sd = sd_cap_s(control_mhs, control_reps)),
  dat %>% filter(!is.na(treat_sd), treat_sd > sd_cap_s(treat_mhs, treat_reps) + 0.5) %>%
    transmute(study_id, reference, group = "treatment", mean = treat_mhs, sd = treat_sd, max_sd = sd_cap_s(treat_mhs, treat_reps)))
print(tab_sd_cap_check, n = Inf, width = Inf)

# studies with missing SDs must be exactly those flagged for imputation
na_sd_ids <- dat %>% filter(is.na(control_sd) | is.na(treat_sd)) %>% pull(study_id) %>% unique()
print(na_sd_ids)
stopifnot(all(na_sd_ids %in% impute_ids))
if (length(setdiff(impute_ids, na_sd_ids)) > 0)
  message("impute_ids with no missing SD: ", paste(setdiff(impute_ids, na_sd_ids), collapse = ", "))

# [V3] SD imputation: class CV x the arm's own mean (replaces the constant class-mean SD)
sd_floor    <- 0.5     # percentage points: smallest SD an imputed arm can take
min_mean    <- 1       # arms with a mean below this (%) are left out of the CV calculation
dedupe_ctrl <- FALSE   # TRUE = count each control group once when averaging control CVs
min_n_cv    <- 5       # minimum observed arms per class before a class may be used for imputation

obs      <- dat %>% filter(!study_id %in% impute_ids)          # observed SDs only
ctrl_obs <- if (dedupe_ctrl) distinct(obs, study_id, taxa_class, control_mhs, control_sd, control_reps) else obs

class_cv <- full_join(
  ctrl_obs %>% filter(!is.na(control_sd), control_mhs >= min_mean) %>%
    group_by(taxa_class) %>%
    summarise(n_cv_ctrl = n(), cv_ctrl = mean(control_sd / control_mhs),
              cv_ctrl_median = median(control_sd / control_mhs), .groups = "drop"),
  obs %>% filter(!is.na(treat_sd), treat_mhs >= min_mean) %>%
    group_by(taxa_class) %>%
    summarise(n_cv_trt = n(), cv_trt = mean(treat_sd / treat_mhs),
              cv_trt_median = median(treat_sd / treat_mhs), .groups = "drop"),
  by = "taxa_class")
print(class_cv)

# the old constant class-mean SD is kept as a sensitivity alternative (section 13)
class_sd <- obs %>%
  select(taxa_class, control_sd, treat_sd) %>%
  pivot_longer(c(control_sd, treat_sd), values_to = "sd") %>%
  filter(!is.na(sd)) %>%
  group_by(taxa_class) %>%
  summarise(mean_sd = mean(sd), n_sd = n(), .groups = "drop")
print(class_sd)

# guard: every class that needs imputation must rest on enough observed arms
need <- dat %>% filter(study_id %in% impute_ids, !study_id %in% binomial_ids) %>% distinct(taxa_class) %>%
  left_join(class_cv, by = "taxa_class")
print(need)
stopifnot(isTRUE(all(need$n_cv_ctrl >= min_n_cv)), isTRUE(all(need$n_cv_trt >= min_n_cv)))

dat <- dat %>%
  left_join(select(class_cv, taxa_class, cv_ctrl, cv_trt), by = "taxa_class") %>%
  left_join(select(class_sd, taxa_class, mean_sd), by = "taxa_class") %>%
  mutate(sd_imputed = study_id %in% impute_ids & (is.na(control_sd) | is.na(treat_sd)),
         ctrl_na    = sd_imputed & is.na(control_sd),
         trt_na     = sd_imputed & is.na(treat_sd),
         binom      = study_id %in% binomial_ids,
         control_sd_cm = if_else(ctrl_na, mean_sd, control_sd),       # old rule (all imputed studies alike), kept for a sensitivity set
         treat_sd_cm   = if_else(trt_na,  mean_sd, treat_sd),
         # [V4] class-CV SDs are capped at the maximum SD possible for a 0-100 variable at that mean (arms at 0% or 100% fall to the floor)
         cap_hit    = sd_imputed & !binom & ((ctrl_na & cv_ctrl * control_mhs > sd_cap(control_mhs)) |
                                               (trt_na & cv_trt * treat_mhs > sd_cap(treat_mhs))),
         control_sd = if_else(ctrl_na, pmax(if_else(binom, bin_sd(control_mhs, control_reps), pmin(cv_ctrl * control_mhs, sd_cap(control_mhs))), sd_floor), control_sd),
         treat_sd   = if_else(trt_na,  pmax(if_else(binom, bin_sd(treat_mhs,  treat_reps),  pmin(cv_trt  * treat_mhs,  sd_cap(treat_mhs))),  sd_floor), treat_sd),
         floor_hit  = sd_imputed & (control_sd <= sd_floor | treat_sd <= sd_floor)) %>%
  select(-cv_ctrl, -cv_trt, -mean_sd, -ctrl_na, -trt_na)

stopifnot(!anyNA(dat$control_sd), !anyNA(dat$treat_sd))

# every imputed arm, with the SD it received (check these before going on)
imp_log <- dat %>% filter(sd_imputed) %>%
  select(study_id, reference, taxa_class, control_mhs, control_sd, control_sd_cm,
         treat_mhs, treat_sd, treat_sd_cm, cap_hit, floor_hit)
print(imp_log, n = Inf, width = Inf)
imp_log %>% group_by(study_id, reference, taxa_class) %>%
  summarise(arms = n(), floor_hits = sum(floor_hit), cap_hits = sum(cap_hit), min_treat_sd = min(treat_sd), max_treat_sd = max(treat_sd),
            .groups = "drop") %>% print(width = Inf)

# 3. Effect sizes (negative g = reduced hatching) ---------------------------
dat <- escalc(measure = "SMD",
              m1i = treat_mhs,   sd1i = treat_sd,   n1i = treat_reps,
              m2i = control_mhs, sd2i = control_sd, n2i = control_reps,
              data = dat)
stopifnot(!anyNA(dat$yi), !anyNA(dat$vi), all(dat$vi > 0))
dat <- tibble::as_tibble(dat)    # drop the escalc class: printing such objects with n = / width = fails ("invalid 'na.print'")
dat <- dat %>% mutate(es_scale = "g",                  # [V4] scale tag read by fit_cred() to build the right sampling covariances
                      zero_treat = treat_mhs <= 0)     # [V4] zero-hatching arms (lnRR needs a floor for these)

# 4. Concentration units -> mg/L (aqueous only) -----------------------------
unit_lookup <- tribble(
  ~concentration_unit, ~conc_class, ~to_mg_L,
  "mg/L",    "aqueous",  1,
  "mg MP/L", "aqueous",  1,      # confirm mass concentration (Pannetier et al. 2019)
  "ug/L",    "aqueous",  1e-3,
  "g/L",     "aqueous",  1e3,
  "g/kg",    "solid",    NA,
  "mg/kg",   "solid",    NA,
  "% sed dw","solid",    NA,
  "MP/ml",   "particle", NA,
  "MP/L",    "particle", NA,
  "fibres/fish/day", "other", NA
)

dat <- dat %>%
  left_join(unit_lookup, by = "concentration_unit") %>%
  mutate(conc_mg_L = concentration * to_mg_L) %>%
  select(-to_mg_L)
stopifnot(!anyNA(dat$conc_class))

# 5. Recode moderators ------------------------------------------------------
dat <- dat %>%
  mutate(
    contam = case_when(
      str_detect(str_to_lower(contam), "co-contaminant") ~ "Co-contaminant",
      contam == "Pure Plastic"                           ~ "Pure plastic"),
    polymer_grp = case_when(
      polymer == "PS"                    ~ "PS",
      polymer %in% c("PE","LDPE","HDPE") ~ "PE",
      polymer %in% c("Mix","CTR","AFH","RFP") ~ polymer,
      TRUE                               ~ "Other"),
    lifestage_grp = if_else(lifestage == "Maternal", "Parental", lifestage),
    route = case_when(
      str_detect(pathway, "^Ingestion \\+")               ~ "Combined",
      pathway == "Ingestion"                              ~ "Ingestion",
      pathway == "Chronic_sediment"                       ~ "Sediment",
      pathway %in% c("Acute_solution","Chronic_solution") ~ "Water",
      TRUE                                                ~ pathway),
    log_size = log10(size_um),
    size_bin = cut(size_um, breaks = c(0, 1, 100, Inf), right = FALSE,
                   labels = c("nano (<1 um)", "small micro (1-<100 um)", "large micro (>=100 um)"))
  )
stopifnot(!anyNA(dat$contam), !anyNA(dat$size_bin))

# Testable levels only (>= 3 studies in the 29-study data); the rest stay descriptive.
# run_mod() (section 9) applies its own >= 3-study rule within each dataset; this block is
# kept (comparable with the earlier run) and section 9 prints which held-out levels now qualify.
dat <- dat %>%
  mutate(
    polymer_test   = factor(if_else(polymer_grp   %in% c("PS","PE","Other"), polymer_grp, NA_character_),
                            levels = c("PS","PE","Other")),
    lifestage_test = factor(if_else(lifestage_grp %in% c("Embryonic","Parental","Combined"), lifestage_grp, NA_character_),
                            levels = c("Embryonic","Parental","Combined")),
    route_test     = factor(if_else(route %in% c("Water","Combined","Ingestion"), route, NA_character_),
                            levels = c("Water","Combined","Ingestion")),
    contam         = factor(contam, levels = c("Pure plastic","Co-contaminant")),
    size_bin       = relevel(size_bin, ref = "small micro (1-<100 um)"),
    env_relevant   = factor(env_relevant, levels = c("No","Yes"))
  )

# 6. Phylogeny and shared-control IDs ---------------------------------------
species_list <- sort(unique(dat$species))
taxa <- tnrs_match_names(species_list, context_name = "Animals")

failed <- is.na(taxa$ott_id)
if (any(failed)) {
  retry <- tnrs_match_names(taxa$search_string[failed])
  taxa$ott_id[failed]      <- retry$ott_id
  taxa$unique_name[failed] <- retry$unique_name
}
stopifnot(!anyNA(taxa$ott_id))

lookup <- tibble(species = species_list,
                 species_phylo = gsub("\\s*\\(.*\\)$", "", taxa$unique_name))

tree <- tol_induced_subtree(ott_ids = taxa$ott_id)
tree$tip.label <- tree$tip.label %>%
  gsub("_ott[0-9]+", "", .) %>% gsub("_", " ", .) %>% gsub("\\s*\\(.*\\)$", "", .)
stopifnot(setequal(lookup$species_phylo, tree$tip.label))

if (is.null(tree$edge.length)) tree <- compute.brlen(tree)
if (!is.ultrametric(tree))     tree <- force.ultrametric(tree, method = "extend")
phylo_cor <- vcv(tree, corr = TRUE)

dat <- dat %>% left_join(lookup, by = "species")     # [V4] control_id is now created in section 2

stopifnot(!anyNA(dat$species_phylo), all(dat$species_phylo %in% rownames(phylo_cor)))
dat %>% count(control_id) %>% count(n, name = "n_control_groups")

# [V4] control-group check: arms that share a control should share species, lifestage and exposure pathway. Flagged groups need
# confirming in the spreadsheet (a control_group_id column settles them); large groups are listed too.
tab_control_check <- dat %>% group_by(study_id, reference, control_id) %>%
  summarise(arms = n(), n_species = n_distinct(species), n_lifestage = n_distinct(lifestage),
            n_pathway = n_distinct(pathway), .groups = "drop") %>%
  mutate(flag = arms > 1 & (n_species > 1 | n_lifestage > 1 | n_pathway > 1))
tab_control_check %>% filter(flag) %>% print(n = Inf, width = Inf)
tab_control_check %>% filter(arms >= 5) %>% print(n = Inf, width = Inf)

# 7. Analysis datasets (always last, after every column exists) -------------
prep <- function(d) d %>% mutate(log_conc = log10(conc_mg_L))   # NA for non-aqueous units

dat_all       <- dat %>% prep()
dat_cred      <- dat_all %>% filter(cred_r != "R3")
dat_dose      <- dat_all %>% filter(conc_class == "aqueous")
dat_dose_cred <- dat_dose %>% filter(cred_r != "R3")   # aqueous, co-contaminants included (SECONDARY)

stopifnot(all(c("control_id", "species_phylo", "conc_mg_L", "log_conc") %in% names(dat_cred)))

# primary / secondary / sensitivity datasets on the g scale
sets_g <- list(
  "Primary: pure plastic, CRED"           = dat_cred %>% filter(contam == "Pure plastic"),
  "R3 retained (pure plastic)"            = dat_all  %>% filter(contam == "Pure plastic"),
  "PlotDigitiser excluded (pure plastic)" = dat_cred %>% filter(contam == "Pure plastic", data_source != "PlotDigitiser"),
  "Secondary: co-contaminants included"   = dat_cred
)
sapply(sets_g, \(d) c(k_es = nrow(d), k_studies = n_distinct(d$study_id)))

# log response ratio version of any dataset (zero hatching floored at floor_pct)
make_lnrr <- function(d, floor_pct = 0.5) {
  d %>% mutate(treat_floor = pmax(treat_mhs, floor_pct)) %>%
    escalc(measure = "ROM", m1i = treat_floor, sd1i = treat_sd, n1i = treat_reps,
           m2i = control_mhs, sd2i = control_sd, n2i = control_reps, data = .,
           var.names = c("yi_l", "vi_l")) %>%
    mutate(yi = yi_l, vi = vi_l, es_scale = "lnRR") %>%   # fit_cred() reads yi, vi and es_scale
    tibble::as_tibble()
}
sets_l <- purrr::map(sets_g, make_lnrr)

# [V3] the same data with the old constant class-mean SD for the imputed arms (sensitivity set)
swap_sd <- function(d) {
  d %>% select(-yi, -vi) %>%
    mutate(control_sd = control_sd_cm, treat_sd = treat_sd_cm) %>%
    escalc(measure = "SMD", m1i = treat_mhs, sd1i = treat_sd, n1i = treat_reps,
           m2i = control_mhs, sd2i = control_sd, n2i = control_reps, data = .) %>%
    tibble::as_tibble()
}

# shorthand used from here on
dat_pure      <- sets_g[["Primary: pure plastic, CRED"]]            # PRIMARY
dat_sec       <- sets_g[["Secondary: co-contaminants included"]]    # SECONDARY
dat_pure_dose <- dat_pure %>% filter(conc_class == "aqueous")       # PRIMARY, aqueous units only
# [V4] lnRR is the primary scale: these are the primary / secondary datasets for every headline analysis
dat_pure_l      <- sets_l[["Primary: pure plastic, CRED"]]
dat_sec_l       <- sets_l[["Secondary: co-contaminants included"]]
dat_pure_l_dose <- dat_pure_l %>% filter(conc_class == "aqueous")

# [V3] record flow for the Methods
flow <- tibble(
  step    = c("Records extracted", "After excluding arms with < 2 replicates", "After removing CRED R3 records",
              "Primary: pure plastic only", "Co-contaminant records (secondary analysis only)"),
  records = c(unname(n_raw["records"]), nrow(dat_all), nrow(dat_cred), nrow(dat_pure),
              sum(dat_cred$contam == "Co-contaminant")),
  studies = c(unname(n_raw["studies"]), n_distinct(dat_all$study_id), n_distinct(dat_cred$study_id),
              n_distinct(dat_pure$study_id), n_distinct(dat_cred$study_id[dat_cred$contam == "Co-contaminant"])))
print(flow)

# 8. Model function, helpers and the base CRED models -----------------------

# [V4] Sampling covariance matrix built from the data (replaces the fixed rho = 0.5).
# Records that share a control group are correlated through that control:
#   lnRR: Cov(i, j) = SD_c^2 / (n_c * m_c^2)                          the variance of the log control mean
#   g   : Cov(i, j) = 1/n_c + g_i g_j / (2 (n_c + n_ti + n_tj))        Gleser & Olkin (2009) approximation, several treatments vs one control
# Records with different controls are independent. If a block is not positive definite its off-diagonals are shrunk (with a warning).
build_V <- function(d) {
  sc <- unique(d$es_scale); stopifnot(length(sc) == 1, sc %in% c("lnRR", "g"))
  n  <- nrow(d); V <- diag(d$vi, nrow = n)
  for (cg in unique(d$control_id)) {
    i <- which(d$control_id == cg)
    if (length(i) < 2) next
    nc <- d$control_reps[i[1]]
    V[i, i] <- if (sc == "lnRR") {
      d$control_sd[i[1]]^2 / (nc * d$control_mhs[i[1]]^2)
    } else {
      1 / nc + outer(d$yi[i], d$yi[i]) / (2 * (nc + outer(d$treat_reps[i], d$treat_reps[i], "+")))
    }
    diag(V)[i] <- d$vi[i]
  }
  D <- diag(diag(V), nrow = n); off <- V - D; k <- 0
  while (inherits(try(chol(D + off), silent = TRUE), "try-error") && k < 25) { off <- off * 0.9; k <- k + 1 }
  if (k > 0) warning("build_V: off-diagonals shrunk by 0.9^", k, " to keep V positive definite")
  D + off
}

fit_cred <- function(data, struct = c("simple","phylo","phylo_sp","species"), mods = ~ 1, rho = 0.5,
                     method = "REML", vmethod = c("exact", "rho"), resid = c("auto", "common", "by_contam")) {
  # simple = study + observation | species = + species (no tree; the old paper's no-phylogeny model)
  # phylo = + phylogenetic species effect | phylo_sp = + species + phylogenetic species effect
  # vmethod "exact" = sample-size-based covariances (default); "rho" = fixed correlation rho (kept only for comparison)
  # [V5] resid "common" = one observation-level variance; "by_contam" = a separate one for pure-plastic and co-contaminant records
  #      (a diagonal structure, ~ contam | es_id); "auto" (default) = by_contam when the data contain both types (>= min_grp records
  #      each) and secondary_resid is "by_contam", otherwise common. The variance structure used is stored in attr(model, "resid").
  struct <- match.arg(struct); vmethod <- match.arg(vmethod); resid <- match.arg(resid)
  stopifnot("es_scale" %in% names(data))
  if (resid == "auto") {
    tt <- table(as.character(data$contam))
    resid <- if (secondary_resid == "by_contam" && length(tt) == 2 && all(tt >= min_grp)) "by_contam" else "common"
  }
  V <- if (vmethod == "exact") build_V(data) else vcalc(vi, cluster = control_id, obs = es_id, rho = rho, data = data)
  random <- list(~ 1 | study_id, ~ 1 | es_id)
  if (resid == "by_contam") {
    data$contam <- factor(as.character(data$contam), levels = c("Pure plastic", "Co-contaminant"))   # fixes the order of the two variances
    random[[2]] <- ~ contam | es_id
  }
  if (struct %in% c("phylo_sp", "species")) random <- c(random, list(~ 1 | species))
  if (struct %in% c("phylo", "phylo_sp"))   random <- c(random, list(~ 1 | species_phylo))
  R <- if (struct %in% c("phylo", "phylo_sp")) {
    sp <- unique(data$species_phylo)
    list(species_phylo = phylo_cor[sp, sp, drop = FALSE])
  }
  m <- if (resid == "by_contam") {
    rma.mv(yi, V, mods = mods, random = random, struct = "DIAG", R = R,
           data = data, method = method, test = "t")
  } else {
    rma.mv(yi, V, mods = mods, random = random, R = R,
           data = data, method = method, test = "t")
  }
  attr(m, "resid") <- resid
  m
}

# [V5] Variance helpers that work for both residual-variance structures.
resid_of    <- \(m) { r <- attr(m, "resid"); if (is.null(r)) "common" else r }
obs_by_type <- \(m) if (resid_of(m) == "by_contam") as.numeric(m$tau2) else c(NA_real_, NA_real_)   # (pure plastic, co-contaminant)
obs_var_vec <- function(m, d) {      # observation-level variance of each record
  if (resid_of(m) == "by_contam") as.numeric(m$tau2[match(as.character(d$contam), c("Pure plastic", "Co-contaminant"))])
  else rep(m$sigma2[match("es_id", m$s.names)], nrow(d))
}
extra_var <- \(m) sum(m$sigma2[m$s.names != "es_id"])                 # study (+ species, phylogeny) variance
s2_total  <- \(m, d) extra_var(m) + mean(obs_var_vec(m, d))            # all non-sampling variance, record-averaged

# helpers (each defined once)
add_wb <- function(d, x) {      # within/between-study split of a numeric moderator
  b <- paste0(x, "_b"); w <- paste0(x, "_w")
  d %>% group_by(study_id) %>%
    mutate(!!b := mean(.data[[x]], na.rm = TRUE),
           !!w := .data[[x]] - mean(.data[[x]], na.rm = TRUE)) %>% ungroup()
}
rob <- function(m, d) robust(m, cluster = d$study_id, clubSandwich = TRUE)
rob_safe <- function(m, cluster) {   # CR2, falling back to CR1 (e.g. models with a species or tree term)
  tryCatch({ r <- robust(m, cluster = cluster, clubSandwich = TRUE); attr(r, "infer") <- "CR2"; r },
           error = function(e) { r <- robust(m, cluster = cluster); attr(r, "infer") <- "CR1"; r })
}
sc_lab       <- \(d) if (d$es_scale[1] == "lnRR") "lnRR" else "Hedges' g"
scale_levels <- c("lnRR", "Hedges' g")                      # lnRR (primary) always first

# [V4] Robust inference, reported only where Satterthwaite df >= df_min; the model-based p and CI are always returned alongside.
inf_tbl <- function(m, rb, i = NULL) {
  nb <- length(rb$beta); i <- if (is.null(i)) seq_len(nb) else i
  dfv <- rb$dfs; if (is.null(dfv)) dfv <- rb$ddf; if (is.null(dfv)) dfv <- m$k - m$p
  df <- rep_len(as.numeric(dfv), nb)[i]
  ok <- !is.na(df) & df >= df_min
  tibble(term = rownames(rb$beta)[i], est = as.numeric(rb$beta)[i], df = df, reportable = ok,
         ci_lb = if_else(ok, rb$ci.lb[i], NA_real_), ci_ub = if_else(ok, rb$ci.ub[i], NA_real_),
         p_robust = if_else(ok, rb$pval[i], NA_real_),
         p_model = m$pval[i], ci_lb_model = m$ci.lb[i], ci_ub_model = m$ci.ub[i])
}

ov <- function(m, label, d) tibble(scale = sc_lab(d), model = label, resid = resid_of(m), k_es = m$k, k_studies = n_distinct(d$study_id),
                                   est = as.numeric(coef(m)), ci_lb = m$ci.lb, ci_ub = m$ci.ub, p_model = m$pval)   # model-based
vc <- function(m, label, d) {
  out <- tibble(scale = sc_lab(d), model = label, term = m$s.names, sigma2 = round(m$sigma2, 3))
  if (resid_of(m) == "by_contam")
    out <- bind_rows(out, tibble(scale = sc_lab(d), model = label, term = c("es_id: Pure plastic", "es_id: Co-contaminant"),
                                 sigma2 = round(m$tau2, 3)))
  mutate(out, resid = resid_of(m))
}
tidy_lrt <- function(a, label) tibble(comparison = label, LRT = as.numeric(a$LRT), p = a$pval)

red <- function(d, label, struct = "simple") {   # overall effect (robust, df-gated, + model-based p) and variance components
  m <- fit_cred(d, struct); rb <- rob_safe(m, d$study_id)
  bind_cols(tibble(scale = sc_lab(d), analysis = label, k_es = m$k, k_studies = n_distinct(d$study_id)),
            select(inf_tbl(m, rb, 1), -term),
            tibble(s2_study = m$sigma2[1], s2_obs = mean(obs_var_vec(m, d)),
                   s2_obs_pure = obs_by_type(m)[1], s2_obs_co = obs_by_type(m)[2],
                   resid = resid_of(m), infer = attr(rb, "infer")))
}

ov_dose <- function(d, struct = "simple") {      # overall effect + within-study concentration slope (aqueous records)
  m <- fit_cred(d, struct); rb <- rob_safe(m, d$study_id); o <- inf_tbl(m, rb, 1)
  w <- tryCatch({
    da  <- d %>% filter(conc_class == "aqueous", !is.na(log_conc)) %>% add_wb("log_conc")
    m2  <- fit_cred(da, struct, mods = ~ log_conc_b + log_conc_w); rb2 <- rob_safe(m2, da$study_id); i2 <- inf_tbl(m2, rb2, 3)
    tibble(k_dose = m2$k, within = i2$est, w_lb = i2$ci_lb, w_ub = i2$ci_ub, p_within = i2$p_robust,
           p_within_model = i2$p_model, df_within = i2$df)
  }, error = function(e) tibble(k_dose = NA_integer_, within = NA_real_, w_lb = NA_real_, w_ub = NA_real_,
                                p_within = NA_real_, p_within_model = NA_real_, df_within = NA_real_))
  bind_cols(tibble(scale = sc_lab(d), k_es = m$k, k_studies = n_distinct(d$study_id),
                   est = o$est, ci_lb = o$ci_lb, ci_ub = o$ci_ub, p_robust = o$p_robust, p_model = o$p_model,
                   df = o$df, resid = resid_of(m), infer = attr(rb, "infer")), w)
}

fmt_p  <- \(p) if_else(p < 0.001, "<0.001", sprintf("%.3f", p))
f      <- \(x, d) ifelse(is.na(x), NA_character_, sprintf("%.*f", d, x))
pct    <- \(x) 100 * (exp(x) - 1)
fmt_ci <- \(lo, hi, d) if_else(is.na(lo) | is.na(hi), "NE (df < 4)", paste0(f(lo, d), ", ", f(hi, d)))
fmt_pr <- \(p) if_else(is.na(p), "NE (df < 4)", fmt_p(p))

struct_labels <- c(simple = "simple", species = "species", phylo = "phylo", phylo_sp = "phylo + species")
# If a structure cannot be fitted with separate residual variances, it falls back to the common variance (the fit's "resid" attribute and
# the resid column of the tables record which was used) with a warning.
base_fits <- \(d, resid = "auto") purrr::map(purrr::set_names(names(struct_labels)), \(s)
                                             tryCatch(fit_cred(d, s, resid = resid), error = function(e) {
                                               warning("base_fits: ", s, " with resid = ", resid, " failed (", conditionMessage(e), "); using the common residual variance")
                                               fit_cred(d, s, resid = "common") }))

# datasets by scale (lnRR first everywhere)
d_P <- list(lnRR = dat_pure_l, g = dat_pure)     # PRIMARY
d_S <- list(lnRR = dat_sec_l,  g = dat_sec)      # SECONDARY

## ---- Class support (>= 3 STUDIES) ----
class_support <- dat_pure %>% group_by(taxa_class) %>%
  summarise(k_es = n(), k_studies = n_distinct(study_id), .groups = "drop")
print(class_support)
class_support_sec <- dat_sec %>% group_by(taxa_class) %>%
  summarise(k_es = n(), k_studies = n_distinct(study_id), .groups = "drop")
print(class_support_sec)

## ---- Base models: four random-effects structures x two datasets x two scales ----
fits_P <- purrr::map(d_P, base_fits)
fits_S <- purrr::map(d_S, base_fits)                                        # [V5] separate residual variance by contaminant type (default)
fits_S_common <- purrr::map(d_S, \(d) base_fits(d, resid = "common"))       # [V5] the earlier common-variance model, kept as a labelled sensitivity
fits_P$lnRR$phylo$s.names      # study_id, es_id, species_phylo
fits_P$lnRR$species$s.names    # study_id, es_id, species

tax_fit <- function(d) {       # taxa class as moderator, classes with >= 3 studies only
  cs <- d %>% group_by(taxa_class) %>% summarise(k = n_distinct(study_id), .groups = "drop")
  dt <- d %>% filter(taxa_class %in% cs$taxa_class[cs$k >= 3])
  list(data = dt, m0 = fit_cred(dt, "simple"), m = fit_cred(dt, "simple", mods = ~ 0 + taxa_class),
       omni = fit_cred(dt, "simple", mods = ~ taxa_class))
}
tax_P <- purrr::map(d_P, tax_fit)
tax_S <- purrr::map(d_S, tax_fit)

for (nm in names(d_P)) {
  cat("\n=====", nm, "=====\n")
  print(AIC(fits_P[[nm]]$simple, fits_P[[nm]]$species, fits_P[[nm]]$phylo, fits_P[[nm]]$phylo_sp))
  print(anova(fits_P[[nm]]$simple, fits_P[[nm]]$species))           # does a species term help?
  print(anova(fits_P[[nm]]$simple, fits_P[[nm]]$phylo))
  print(anova(fits_P[[nm]]$species, fits_P[[nm]]$phylo_sp))         # does phylogeny help once species is in?
  print(anova(fits_P[[nm]]$phylo, fits_P[[nm]]$phylo_sp))           # does the species term help once phylogeny is in?
  print(tax_P[[nm]]$omni)
}

## ---- Overall estimates and variance components (model-based p; robust inference is in Table 1) ----
overall_rows <- function(nm) bind_rows(
  purrr::imap_dfr(struct_labels, \(lab, s) ov(fits_P[[nm]][[s]], paste0("Primary: ", lab), d_P[[nm]])),
  ov(tax_P[[nm]]$m0, "Primary: taxa subset (no class)", tax_P[[nm]]$data),
  purrr::imap_dfr(struct_labels, \(lab, s) ov(fits_S[[nm]][[s]], paste0("Secondary: ", lab), d_S[[nm]])),
  ov(tax_S[[nm]]$m0, "Secondary: taxa subset (no class)", tax_S[[nm]]$data))
tab_overall <- purrr::map_dfr(names(d_P), overall_rows) %>%
  mutate(across(c(est, ci_lb, ci_ub), ~ round(.x, 3)), p_model = signif(p_model, 2))
print(tab_overall, n = Inf)

vc_rows <- function(nm) bind_rows(
  purrr::imap_dfr(struct_labels, \(lab, s) vc(fits_P[[nm]][[s]], paste0("Primary: ", lab), d_P[[nm]])),
  purrr::imap_dfr(struct_labels, \(lab, s) vc(fits_S[[nm]][[s]], paste0("Secondary: ", lab), d_S[[nm]])),
  purrr::imap_dfr(struct_labels, \(lab, s) vc(fits_S_common[[nm]][[s]], paste0("Secondary (common residual variance): ", lab), d_S[[nm]])))
tab_varcomp <- purrr::map_dfr(names(d_P), vc_rows)
print(tab_varcomp, n = Inf)

aic_rows <- function(fits, ds, d) purrr::imap_dfr(fits, \(m, s)
                                                  tibble(scale = sc_lab(d), dataset = ds, model = struct_labels[[s]], resid = resid_of(m), df = attr(logLik(m), "df"), AIC = AIC(m)))
tab_aic <- purrr::map_dfr(names(d_P), \(nm) bind_rows(aic_rows(fits_P[[nm]], "Primary", d_P[[nm]]),
                                                      aic_rows(fits_S[[nm]], "Secondary", d_S[[nm]])))
lrt_rows <- function(f, ds, d) bind_rows(
  tidy_lrt(anova(f$simple,  f$species),  "simple vs species"),
  tidy_lrt(anova(f$simple,  f$phylo),    "simple vs phylo"),
  tidy_lrt(anova(f$species, f$phylo_sp), "species vs phylo + species (phylogeny given species)"),
  tidy_lrt(anova(f$phylo,   f$phylo_sp), "phylo vs phylo + species (species given phylogeny)")) %>%
  mutate(dataset = ds, scale = sc_lab(d), .before = 1)
tab_lrt <- purrr::map_dfr(names(d_P), \(nm) bind_rows(lrt_rows(fits_P[[nm]], "Primary", d_P[[nm]]),
                                                      lrt_rows(fits_S[[nm]], "Secondary", d_S[[nm]])))
# [V5] secondary: separate vs common residual variance, every structure and scale (estimates and AIC)
sec_resid_rows <- function(nm) purrr::imap_dfr(struct_labels, \(lab, s) bind_rows(
  ov(fits_S[[nm]][[s]], lab, d_S[[nm]])        %>% mutate(AIC = as.numeric(AIC(fits_S[[nm]][[s]]))),
  ov(fits_S_common[[nm]][[s]], lab, d_S[[nm]]) %>% mutate(AIC = as.numeric(AIC(fits_S_common[[nm]][[s]])))))
tab_secondary_resid <- purrr::map_dfr(names(d_S), sec_resid_rows) %>%
  select(scale, model, resid, est, ci_lb, ci_ub, p_model, AIC)
print(tab_secondary_resid %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = Inf)

class_row <- \(t, ds) tibble(dataset = ds, scale = sc_lab(t$data), F = as.numeric(t$omni$QM),
                             df1 = t$omni$QMdf[1], df2 = t$omni$QMdf[2], p_model = t$omni$QMp)
tab_class <- purrr::map_dfr(names(d_P), \(nm) bind_rows(class_row(tax_P[[nm]], "Primary"), class_row(tax_S[[nm]], "Secondary")))
print(tab_lrt, n = Inf); print(tab_class)

## ---- Diagnostics (PRIMARY unless stated) ----
# 1. Are the variance components identifiable? (a clear peak = yes, flat = no)
profile(fits_P$lnRR$phylo)

# 2. Phylogenetic signal on species mean effect, as in the old paper
phy_sig <- function(d, ds) {
  set.seed(1)      # seeded so the K p-value is the same wherever it is reported
  me  <- d %>% group_by(species_phylo) %>% summarise(m = mean(yi), .groups = "drop") %>% deframe()
  tr_ <- keep.tip(tree, names(me)); me <- me[tr_$tip.label]
  l <- phylosig(tr_, me, method = "lambda", test = TRUE)
  k <- phylosig(tr_, me, method = "K", test = TRUE, nsim = 1000)
  tibble(dataset = ds, scale = sc_lab(d), lambda = l$lambda, lambda_p = l$P, K = k$K, K_p = k$P)
}
tab_phylosig <- bind_rows(phy_sig(dat_pure_l, "Primary"), phy_sig(dat_pure, "Primary"),
                          phy_sig(dat_sec_l, "Secondary"), phy_sig(dat_sec, "Secondary"))
print(tab_phylosig)

# 3. Sample-size-based covariances vs a fixed correlation rho (the fixed-rho fits are kept only as a comparison)
rho_one <- function(d, ds) {
  fits <- c(list(fit_cred(d, "simple")), purrr::map(c(0.2, 0.5, 0.8), \(r) fit_cred(d, "simple", rho = r, vmethod = "rho")))
  tibble(dataset = ds, scale = sc_lab(d), method = c("exact (sample-size based)", rep("fixed rho", 3)),
         rho = c(NA, 0.2, 0.5, 0.8),
         est = purrr::map_dbl(fits, \(m) as.numeric(coef(m))), ci_lb = purrr::map_dbl(fits, \(m) m$ci.lb),
         ci_ub = purrr::map_dbl(fits, \(m) m$ci.ub), s2_study = purrr::map_dbl(fits, \(m) m$sigma2[1]),
         s2_obs = purrr::map_dbl(fits, \(m) mean(obs_var_vec(m, d))), resid = purrr::map_chr(fits, resid_of))
}
tab_rho <- bind_rows(rho_one(dat_pure_l, "Primary"), rho_one(dat_pure, "Primary"),
                     rho_one(dat_sec_l, "Secondary"), rho_one(dat_sec, "Secondary"))
tab_rho_sec <- filter(tab_rho, dataset == "Secondary")
print(tab_rho, n = Inf)

# 4. Robust SEs, clustered on study (CR2)
rob(fits_P$lnRR$simple, dat_pure_l)

# 5. Influence: leave one study out (robust inference, df-gated)
tab_loo_overall <- purrr::map_dfr(d_P, \(d) purrr::map_dfr(unique(d$study_id), \(s)
                                                           red(filter(d, study_id != s), "drop") %>% mutate(dropped = s))) %>%
  left_join(distinct(dat_pure, study_id, reference), by = c("dropped" = "study_id")) %>%
  arrange(factor(scale, levels = scale_levels), est)
print(tab_loo_overall, n = 12)

# Descriptive summary for the classes that can't be tested (< 3 studies)
small_classes <- class_support$taxa_class[class_support$k_studies < 3]
tab_class_desc <- purrr::map_dfr(d_P, \(d) d %>% filter(taxa_class %in% small_classes) %>% group_by(scale = sc_lab(d), taxa_class) %>%
                                   summarise(k_es = n(), k_studies = n_distinct(study_id), mean_effect = mean(yi), min_effect = min(yi),
                                             max_effect = max(yi), .groups = "drop"))
print(tab_class_desc)

# 9. Single moderators ------------------------------------------------------
# For each moderator: intercept-only and moderator model on IDENTICAL rows, model-based and cluster-robust omnibus tests
# (CR2, CR1 fallback), robust level estimates, pseudo-R2. [V4] Robust F / CIs / p are reported only where the Satterthwaite
# df >= df_min (otherwise NA, printed "NE (df < 4)"); model-based p is always returned too; Holm is applied among the
# reportable tests.
run_mod <- function(var, d, continuous = FALSE, min_k = 3, struct = "simple") {
  d <- d %>% filter(!is.na(.data[[var]]))
  if (!continuous) {
    lv_all <- if (is.factor(d[[var]])) levels(d[[var]]) else sort(unique(as.character(d[[var]])))
    ks <- d %>% group_by(lvl = as.character(.data[[var]])) %>%
      summarise(ks = n_distinct(study_id), .groups = "drop")
    ok <- intersect(lv_all, ks$lvl[ks$ks >= min_k])
    if (length(ok) < 2) stop("not estimable: fewer than 2 levels with >= ", min_k, " studies")
    d <- d %>% filter(as.character(.data[[var]]) %in% ok)
    d[[var]] <- factor(as.character(d[[var]]), levels = ok)
  } else if (n_distinct(d$study_id) < min_k) stop("not estimable: fewer than ", min_k, " studies")
  
  m0  <- fit_cred(d, struct)                              # intercept only, identical rows
  m1  <- fit_cred(d, struct, mods = reformulate(var))     # omnibus
  rb1 <- rob_safe(m1, d$study_id)
  a   <- anova(rb1, btt = 2:length(coef(m1)))
  qd  <- as.numeric(a$QMdf)
  df2 <- if (length(qd) > 1) qd[2] else NA_real_
  okf <- !is.na(df2) & df2 >= df_min
  
  r2 <- \(x0, x1) if (x0 > 1e-8) max(0, (x0 - x1) / x0) else NA_real_
  fv <- var(as.numeric(model.matrix(m1) %*% coef(m1)))
  
  omni <- tibble(moderator = var, k_es = m1$k, k_studies = n_distinct(d$study_id),
                 n_levels = if (continuous) NA_integer_ else nlevels(d[[var]]),
                 F_model = as.numeric(m1$QM), p_model = m1$QMp,
                 F_robust = if (okf) as.numeric(a$QM) else NA_real_, df1 = qd[1], df2 = df2,
                 p_robust = if (okf) a$QMp else NA_real_, reliable = okf, infer = attr(rb1, "infer"),
                 R2_total = r2(s2_total(m0, d), s2_total(m1, d)),
                 R2_study = r2(m0$sigma2[1], m1$sigma2[1]),
                 R2_obs   = r2(mean(obs_var_vec(m0, d)), mean(obs_var_vec(m1, d))),
                 R2_marg  = fv / (fv + s2_total(m1, d)))
  
  if (continuous) {
    lev <- bind_cols(tibble(level = "slope per unit", k_es = m1$k, k_studies = n_distinct(d$study_id)),
                     select(inf_tbl(m1, rb1, 2), -term))
  } else {
    m2  <- fit_cred(d, struct, mods = reformulate(c("0", var)))      # level-specific means
    rb2 <- rob_safe(m2, d$study_id)
    n_tab <- d %>% group_by(level = as.character(.data[[var]])) %>%
      summarise(k_es = n(), k_studies = n_distinct(study_id), .groups = "drop")
    lev <- inf_tbl(m2, rb2) %>% mutate(level = sub(paste0("^", var), "", term)) %>% select(-term) %>%
      left_join(n_tab, by = "level") %>% relocate(level)
  }
  list(omnibus = omni, levels = lev)
}

specs <- tribble(
  ~var,             ~continuous, ~dose_only,
  "taxa_class",     FALSE, FALSE,
  "polymer_test",   FALSE, FALSE,
  "env_relevant",   FALSE, FALSE,
  "lifestage_test", FALSE, FALSE,
  "route_test",     FALSE, FALSE,
  "size_bin",       FALSE, FALSE,
  "log_size",       TRUE,  FALSE,
  "log_conc",       TRUE,  TRUE,       # aqueous units only, mg/L
  "contam",         FALSE, FALSE       # meaningful only in the secondary set
)

run_grid <- function(sets, specs, scale, struct = "simple") {
  om <- list(); lv <- list()
  for (sn in names(sets)) for (i in seq_len(nrow(specs))) {
    sp <- specs[i, ]; d <- sets[[sn]]
    if (sp$dose_only) d <- filter(d, conc_class == "aqueous")
    key <- paste(sn, sp$var, sep = " | ")
    r <- tryCatch(run_mod(sp$var, d, sp$continuous, struct = struct), error = function(e) conditionMessage(e))
    if (is.character(r)) {
      om[[key]] <- tibble(analysis = sn, scale = scale, struct = struct, moderator = sp$var, note = r)
    } else {
      om[[key]] <- bind_cols(tibble(analysis = sn, scale = scale, struct = struct), r$omnibus)
      lv[[key]] <- bind_cols(tibble(analysis = sn, scale = scale, struct = struct, moderator = sp$var), r$levels)
    }
  }
  list(omnibus = bind_rows(om), levels = bind_rows(lv))
}

res_l <- run_grid(sets_l, specs, "lnRR")           # [V4] lnRR first: all four datasets, simple structure
res_g <- run_grid(sets_g, specs, "Hedges' g")      # then g

drop_contam <- \(x) filter(x, !(moderator == "contam" & !grepl("^Secondary", analysis)))

omni <- bind_rows(res_l$omnibus, res_g$omnibus) %>% drop_contam() %>%
  group_by(analysis, scale, struct) %>% mutate(p_holm = p.adjust(p_robust, "holm")) %>% ungroup()   # Holm among reportable tests
if (!"note" %in% names(omni)) omni$note <- NA_character_      # every fit succeeded -> no note column

levs <- bind_rows(res_l$levels, res_g$levels) %>% drop_contam() %>%
  mutate(pct_change = if_else(scale == "lnRR", 100 * (exp(est) - 1), NA_real_),
         pct_lb     = if_else(scale == "lnRR", 100 * (exp(ci_lb) - 1), NA_real_),
         pct_ub     = if_else(scale == "lnRR", 100 * (exp(ci_ub) - 1), NA_real_))

omni %>% select(analysis, scale, moderator, k_es, k_studies, n_levels, F_robust, df1, df2,
                p_robust, p_holm, p_model, reliable, starts_with("R2"), any_of("note")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)
levs %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% print(n = Inf, width = Inf)

# held-out levels: which descriptive-only levels now have >= 3 studies?
level_support <- function(d, grp, test, moderator, dataset) d %>%
  filter(!is.na(.data[[grp]])) %>%
  group_by(level = as.character(.data[[grp]])) %>%
  summarise(k_es = n(), k_studies = n_distinct(study_id),
            in_test = any(!is.na(.data[[test]])), .groups = "drop") %>%
  mutate(moderator = moderator, dataset = dataset)
tab_level_support <- bind_rows(
  level_support(dat_pure, "polymer_grp",   "polymer_test",   "polymer",   "Primary"),
  level_support(dat_pure, "lifestage_grp", "lifestage_test", "lifestage", "Primary"),
  level_support(dat_pure, "route",         "route_test",     "route",     "Primary"),
  level_support(dat_sec,  "polymer_grp",   "polymer_test",   "polymer",   "Secondary"),
  level_support(dat_sec,  "lifestage_grp", "lifestage_test", "lifestage", "Secondary"),
  level_support(dat_sec,  "route",         "route_test",     "route",     "Secondary"))
tab_level_support %>% filter(!in_test, k_studies >= 3)    # held out of the tests but would now qualify

# Same moderators on the PRIMARY lnRR set under each random-effects structure (model-based and CR2/CR1 tests).
# CR1 rows (species / phylogenetic structures) are anti-conservative with this few clusters: compare estimates, not p.
res_struct <- purrr::map(purrr::set_names(c("simple", "species", "phylo", "phylo_sp")),
                         \(s) run_grid(sets_l["Primary: pure plastic, CRED"], specs, "lnRR", struct = s))
omni_struct <- purrr::map_dfr(res_struct, "omnibus") %>% drop_contam() %>%
  group_by(struct) %>% mutate(p_holm = p.adjust(p_robust, "holm")) %>% ungroup()
if (!"note" %in% names(omni_struct)) omni_struct$note <- NA_character_
omni_struct %>% arrange(moderator, struct) %>%
  select(moderator, struct, k_es, k_studies, F_robust, df2, p_robust, p_holm, p_model, infer, any_of("note")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# Leave-one-study-out on the pooled concentration slope (PRIMARY, both scales)
dat_dose_cred_l <- make_lnrr(dat_dose_cred)        # secondary aqueous records on the lnRR scale
tab_loo_dose <- purrr::map_dfr(list(dat_pure_l_dose, dat_pure_dose), \(d) purrr::map_dfr(unique(d$study_id), \(s) {
  dd <- filter(d, study_id != s); m <- fit_cred(dd, "simple", mods = ~ log_conc); ii <- inf_tbl(m, rob_safe(m, dd$study_id), 2)
  tibble(scale = sc_lab(dd), dropped = s, slope = ii$est, df = ii$df, p_robust = ii$p_robust, p_model = ii$p_model) })) %>%
  arrange(factor(scale, levels = scale_levels), desc(p_model))
print(tab_loo_dose, n = 10)

# Overlap checks
# SECONDARY: contaminant type vs relevance / lifestage
with(dat_sec, table(contam, env_relevant))
with(dat_sec, table(lifestage_test, contam, useNA = "ifany"))
dat_sec %>% group_by(contam) %>%
  summarise(k_studies = n_distinct(study_id), prop_env_yes = mean(env_relevant == "Yes"))
# PRIMARY aqueous: relevance vs concentration range
dat_pure_dose %>% group_by(env_relevant) %>%
  summarise(k_studies = n_distinct(study_id), med = median(conc_mg_L), min = min(conc_mg_L), max = max(conc_mg_L))

# B. Taxa-adjusted moderator models (class + moderator), with robust p (df-gated) and model-based p
adj <- function(var, data) {
  d <- data %>% filter(!is.na(.data[[var]])) %>% droplevels()
  m <- fit_cred(d, "simple", mods = as.formula(paste("~ taxa_class +", var)))
  idx <- grep(paste0("^", var), rownames(m$beta))
  rb <- rob_safe(m, d$study_id); a <- anova(rb, btt = idx)
  df2 <- as.numeric(a$QMdf)[2]; ok <- !is.na(df2) & df2 >= df_min
  tibble(scale = sc_lab(d), moderator = var, model = "taxa + moderator", k_es = m$k, k_studies = n_distinct(d$study_id),
         p_model = anova(m, btt = idx)$QMp, df2 = df2, p_robust = if (ok) a$QMp else NA_real_)
}
mods_adj <- c("polymer_test","env_relevant","lifestage_test","route_test","size_bin","log_size","log_conc")
adj_run <- \(tx, vars, ds) purrr::map_dfr(tx, \(t) purrr::map_dfr(vars, purrr::possibly(\(v) adj(v, t$data), otherwise = NULL))) %>%
  mutate(dataset = ds)
tab_adj <- bind_rows(adj_run(tax_P, mods_adj, "Primary"), adj_run(tax_S, c(mods_adj, "contam"), "Secondary")) %>%
  arrange(dataset, factor(scale, levels = scale_levels))
tab_adj %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf)

# 10. In-study variation, within/between models, non-linear dose ------------
within_var <- function(var, d) d %>% group_by(study_id) %>%
  summarise(n = n_distinct(.data[[var]], na.rm = TRUE), .groups = "drop") %>%
  summarise(moderator = var, studies_varying = sum(n >= 2), studies = n())

vars_pure <- c("polymer_test","env_relevant","lifestage_test","route_test","size_bin","log_size")
tab_within_var <- bind_rows(
  bind_rows(purrr::map_dfr(vars_pure, within_var, d = dat_pure),
            within_var("log_conc", dat_pure_dose)) %>% mutate(dataset = "Primary"),
  bind_rows(purrr::map_dfr(c(vars_pure, "contam"), within_var, d = dat_sec),
            within_var("log_conc", dat_dose_cred)) %>% mutate(dataset = "Secondary"))
print(tab_within_var, n = Inf)

# Within/between models: PRIMARY and SECONDARY, both scales (robust CI / p only where df >= 4; model-based p alongside)
wb_fit <- function(d, x, label) {
  d <- d %>% add_wb(x); m <- fit_cred(d, "simple", mods = reformulate(paste0(x, c("_b", "_w")))); rb <- rob(m, d)
  inf_tbl(m, rb) %>% mutate(model = label, scale = sc_lab(d), .before = 1)
}
env_n    <- \(d) d %>% mutate(env_n = as.numeric(env_relevant == "Yes"))
contam_n <- \(d) d %>% mutate(contam_n = as.numeric(contam == "Co-contaminant"))
tab_wb <- bind_rows(
  wb_fit(dat_pure_l_dose, "log_conc", "Primary: concentration"),   wb_fit(dat_pure_dose, "log_conc", "Primary: concentration"),
  wb_fit(env_n(dat_pure_l), "env_n",  "Primary: relevance"),       wb_fit(env_n(dat_pure), "env_n",  "Primary: relevance"),
  wb_fit(dat_dose_cred_l,  "log_conc", "Secondary: concentration"), wb_fit(dat_dose_cred,  "log_conc", "Secondary: concentration"),
  wb_fit(env_n(dat_sec_l),  "env_n",  "Secondary: relevance"),     wb_fit(env_n(dat_sec),  "env_n",  "Secondary: relevance"),
  wb_fit(contam_n(dat_sec_l), "contam_n", "Secondary: co-contaminant"), wb_fit(contam_n(dat_sec), "contam_n", "Secondary: co-contaminant"))
print(tab_wb %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = Inf, width = Inf)

# Contaminant-type contrast (secondary): pooled, within-study and between-study; used by Table 2 and the sensitivity tables
contam_contrast <- function(d, struct = "simple", resid = "auto") {
  d <- d %>% mutate(contam_n = as.numeric(contam == "Co-contaminant")); sc <- sc_lab(d)
  m  <- fit_cred(d, struct, mods = ~ contam_n, resid = resid); rb <- rob_safe(m, d$study_id); it <- inf_tbl(m, rb, 2)
  dw <- add_wb(d, "contam_n")
  mw <- fit_cred(dw, struct, mods = ~ contam_n_b + contam_n_w, resid = resid); rbw <- rob_safe(mw, dw$study_id); iw <- inf_tbl(mw, rbw)
  st_both  <- tapply(d$contam, d$study_id, \(x) n_distinct(x) == 2)
  n_both   <- sum(st_both)                                                    # studies that tested both types
  rec_both <- sum(as.character(d$study_id) %in% names(st_both)[st_both])      # records in those studies (the data informing the within-study contrast)
  bind_rows(
    mutate(it,      contrast = "Pooled: co-contaminant minus pure plastic", records = m$k,  studies = n_distinct(d$study_id), infer = attr(rb, "infer")),
    mutate(iw[3, ], contrast = paste0("Within-study (", n_both, " studies with both types)"), records = rec_both, studies = n_both, infer = attr(rbw, "infer")),
    mutate(iw[2, ], contrast = "Between-study", records = mw$k, studies = n_distinct(dw$study_id), infer = attr(rbw, "infer"))) %>%
    mutate(scale = sc, resid = resid_of(m), .before = 1) %>% select(-term)
}

# [V5] Effect of each contaminant type from ONE model (shared study variance, type-specific observation-level variance)
bytype <- function(d, struct = "simple", resid = "auto") {
  m <- fit_cred(d, struct, mods = ~ 0 + contam, resid = resid); rb <- rob_safe(m, d$study_id)
  se_m <- sqrt(diag(vcov(m))); lv <- c("Pure plastic", "Co-contaminant")
  inf_tbl(m, rb) %>%
    mutate(type = sub("^contam", "", term), scale = sc_lab(d), resid = resid_of(m), infer = attr(rb, "infer"),
           records = purrr::map_int(type, \(ty) sum(d$contam == ty)),
           studies = purrr::map_int(type, \(ty) n_distinct(d$study_id[d$contam == ty])),
           s2_study = m$sigma2[1],
           s2_obs = if (resid_of(m) == "by_contam") m$tau2[match(type, lv)] else m$sigma2[match("es_id", m$s.names)],
           tau = sqrt(s2_study + s2_obs),
           pi_lb = est - qt(0.975, m$k - m$p) * sqrt(se_m^2 + s2_study + s2_obs),
           pi_ub = est + qt(0.975, m$k - m$p) * sqrt(se_m^2 + s2_study + s2_obs)) %>%
    select(-term)
}

# Non-linear dose: linear within/between model plus a quadratic term, pooled and within-study;
# quadratic vs linear by likelihood-ratio test on ML fits, CR2 test on the term (df-gated)
run_dose <- function(d, struct = "simple") {
  d <- d %>% filter(conc_class == "aqueous", !is.na(log_conc)) %>%
    mutate(lc = log_conc - mean(log_conc)) %>% add_wb("log_conc")
  if (n_distinct(d$study_id) < 5) stop("fewer than 5 studies")
  
  m_wb <- fit_cred(d, struct, mods = ~ log_conc_b + log_conc_w); it <- inf_tbl(m_wb, rob(m_wb, d))
  linear <- tibble(scale = sc_lab(d), k_es = m_wb$k, k_studies = n_distinct(d$study_id),
                   within = it$est[3], w_lb = it$ci_lb[3], w_ub = it$ci_ub[3], p_within = it$p_robust[3],
                   p_within_model = it$p_model[3], df_within = it$df[3],
                   between = it$est[2], b_lb = it$ci_lb[2], b_ub = it$ci_ub[2], p_between = it$p_robust[2],
                   p_between_model = it$p_model[2], df_between = it$df[2])
  
  curve <- function(label, f_lin, f_quad) {
    l_ml <- fit_cred(d, struct, mods = f_lin,  method = "ML")
    q_ml <- fit_cred(d, struct, mods = f_quad, method = "ML")
    q    <- fit_cred(d, struct, mods = f_quad)
    j    <- grep("\\^2", names(coef(q))); iq <- inf_tbl(q, rob(q, d), j)
    tibble(scale = sc_lab(d), term = label, k_es = q$k, quad = iq$est, q_lb = iq$ci_lb, q_ub = iq$ci_ub,
           p_quad_robust = iq$p_robust, p_quad_model = iq$p_model, df = iq$df,
           LRT_p = anova(l_ml, q_ml)$pval, AIC_linear = as.numeric(AIC(l_ml)), AIC_quad = as.numeric(AIC(q_ml)))
  }
  list(linear = linear,
       quad = bind_rows(curve("pooled",       ~ lc, ~ lc + I(lc^2)),
                        curve("within-study", ~ log_conc_b + log_conc_w,
                              ~ log_conc_b + log_conc_w + I(log_conc_w^2))))
}

dose_in  <- list("Primary lnRR" = dat_pure_l, "Primary g" = dat_pure, "Secondary lnRR" = dat_sec_l, "Secondary g" = dat_sec)
dose_res <- purrr::imap(dose_in, \(d, nm) tryCatch(run_dose(d), error = \(e) list(note = conditionMessage(e))))

stack_dose <- \(part) purrr::imap_dfr(dose_res, \(r, nm)
                                      if (!is.null(r[[part]])) bind_cols(tibble(analysis = nm), r[[part]]) else tibble(analysis = nm, note = r$note))
tab_dose_linear <- stack_dose("linear")
tab_dose_quad   <- stack_dose("quad")
tab_dose_linear %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(width = Inf)
tab_dose_quad   %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(width = Inf)

# 11. Additional sensitivity analysis: data source; co-contaminant diagnostics -----
count(dat_pure, data_source)

# SECONDARY dose-response diagnostics (co-contaminant studies), lnRR scale
d_aq <- dat_dose_cred_l

d_aq %>% group_by(contam) %>%
  summarise(k_es = n(), k_studies = n_distinct(study_id))   # how much of the aqueous set is co-contaminant

# Which studies are the co-contaminant ones, and what do their slopes look like?
d_aq %>% group_by(study_id, reference) %>%
  summarise(k_es = n(), n_doses = n_distinct(log_conc),
            contam = paste(unique(contam), collapse = " / "),
            dose_min = min(conc_mg_L), dose_max = max(conc_mg_L),
            slope = if (n_distinct(log_conc) >= 2) coef(lm(yi ~ log_conc))[2] else NA_real_,
            .groups = "drop") %>%
  arrange(contam, desc(abs(slope)))

# Joint model: does the within-study slope survive adjusting for contaminant type?
d_aq2 <- d_aq %>% add_wb("log_conc") %>% mutate(contam_n = as.numeric(contam == "Co-contaminant"))
m_j <- fit_cred(d_aq2, "simple", mods = ~ log_conc_b + log_conc_w + contam_n)
tab_dose_contam_joint <- inf_tbl(m_j, rob(m_j, d_aq2)) %>% mutate(scale = "lnRR", .before = 1)
print(tab_dose_contam_joint %>% mutate(across(where(is.numeric), ~ round(.x, 3))), width = Inf)

# 12. Outlier sensitivity analysis ------------------------------------------
# [V4] Outlier flags are defined on EACH scale's own simple model, fixed in advance:
#   lnRR: studentised residual |z| > 3 (no variance rule: lnRR variances are tiny, so vi > 5 is a g-specific threshold)
#   g   : |z| > 3, or vi > 5, or either
rs_l <- rstudent(fits_P$lnRR$simple)      # refits per row, so it may take a minute
rs_g <- rstudent(fits_P$g$simple)
dat_pure_l <- dat_pure_l %>% mutate(rstud_z = rs_l$z)
dat_pure   <- dat_pure   %>% mutate(rstud_z = rs_g$z)

ids_resid_l <- dat_pure_l$es_id[abs(dat_pure_l$rstud_z) > 3]
ids_resid_g <- dat_pure$es_id[abs(dat_pure$rstud_z) > 3]
ids_var_g   <- dat_pure$es_id[dat_pure$vi > 5]

excl_l <- list("Primary" = integer(0), "Excl. |studentised z| > 3 (lnRR model)" = ids_resid_l)
excl_g <- list("Primary" = integer(0),
               "Excl. |studentised z| > 3 (g model)" = ids_resid_g,
               "Excl. vi > 5 (g)"                    = ids_var_g,
               "Excl. either rule (g)"               = union(ids_resid_g, ids_var_g))

# What gets removed (check this before reading any results)
tab_outliers_removed <- bind_rows(
  dat_pure_l %>% filter(es_id %in% ids_resid_l) %>% mutate(scale = "lnRR", rule = "|z| > 3 (lnRR model)"),
  dat_pure %>% filter(es_id %in% union(ids_resid_g, ids_var_g)) %>%
    mutate(scale = "Hedges' g",
           rule = case_when(es_id %in% ids_resid_g & es_id %in% ids_var_g ~ "both",
                            es_id %in% ids_resid_g                        ~ "residual (g model)",
                            TRUE                                          ~ "variance (vi > 5)"))) %>%
  select(scale, study_id, reference, species, conc_mg_L, yi, vi, rstud_z, rule) %>%
  arrange(factor(scale, levels = scale_levels), desc(abs(yi))) %>%
  tibble::as_tibble()
print(tab_outliers_removed, n = Inf, width = Inf)

# Whole-study influence on the overall mean (lnRR)
sort(cooks.distance(fits_P$lnRR$simple, cluster = dat_pure_l$study_id), decreasing = TRUE)[1:5]

# Overall effect and dose slope after each exclusion, on the scale the flags were defined on
tab_outlier_overall <- bind_rows(
  purrr::imap_dfr(excl_l, \(ids, nm) red(filter(dat_pure_l, !es_id %in% ids), nm) %>% mutate(rows_removed = length(ids))),
  purrr::imap_dfr(excl_g, \(ids, nm) red(filter(dat_pure,   !es_id %in% ids), nm) %>% mutate(rows_removed = length(ids))))
tab_outlier_overall %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(width = Inf)

tab_outlier_dose <- bind_rows(
  purrr::imap_dfr(excl_l, \(ids, nm) ov_dose(filter(dat_pure_l, !es_id %in% ids)) %>% mutate(analysis = nm, rows_removed = length(ids))),
  purrr::imap_dfr(excl_g, \(ids, nm) ov_dose(filter(dat_pure,   !es_id %in% ids)) %>% mutate(analysis = nm, rows_removed = length(ids))))
tab_outlier_dose %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(width = Inf)

# zero-hatching arms in the primary (lnRR needs a floor for these; see section 13)
sum(dat_pure$treat_mhs <= 0)

# 13. Sensitivity complete --------------------------------------------------
pct_rank <- \(d) d %>% mutate(pc = 100 * (treat_mhs / control_mhs - 1)) %>%
  group_by(study_id, reference) %>% summarise(med_pct = median(pc), .groups = "drop") %>% arrange(med_pct)

# One code path for the sensitivity sets, shared by the primary (here) and secondary (section 16.4). Sets are defined on the
# g data and converted to lnRR with make_lnrr(); the outlier sets differ by scale (each scale's own flags).
make_sens <- function(base, base_r3, big) {
  base_l <- make_lnrr(base)
  z_g <- rstudent(fit_cred(base,   "simple"))$z
  z_l <- rstudent(fit_cred(base_l, "simple"))$z
  id_zg <- base$es_id[abs(z_g) > 3]; id_vg <- base$es_id[base$vi > 5]; id_zl <- base$es_id[abs(z_l) > 3]
  common <- list(
    "Base dataset (CRED)"                                    = base,
    "R3 retained"                                            = base_r3,
    "Excl. PlotDigitiser"                                    = filter(base, data_source != "PlotDigitiser"),
    "Excl. imputed SDs"                                      = filter(base, !sd_imputed),
    "Constant class-mean SD (old rule, all imputed studies)" = swap_sd(base),
    "Excl. 3 largest-effect studies (post hoc)"              = filter(base, !study_id %in% big),
    "Excl. zero-hatching arms"                               = filter(base, !zero_treat))
  list(lnRR = c(purrr::map(common, make_lnrr),
                list("Excl. |z| > 3 (lnRR model)" = make_lnrr(filter(base, !es_id %in% id_zl)))),
       g    = c(common,
                list("Excl. |z| > 3 (g model)"      = filter(base, !es_id %in% id_zg),
                     "Excl. vi > 5 (g)"              = filter(base, !es_id %in% id_vg),
                     "Excl. either rule (g)"         = filter(base, !es_id %in% union(id_zg, id_vg)))))
}
big_pri  <- pct_rank(dat_pure)$study_id[1:3]
tab_top3 <- pct_rank(dat_pure) %>% slice(1:5)
print(tab_top3)        # the three largest reductions (studies change with the data)
sens_pri <- make_sens(dat_pure, sets_g[["R3 retained (pure plastic)"]], big_pri)

sens_run <- \(sens, ds) purrr::imap_dfr(sens, \(sets, sc) purrr::imap_dfr(sets, \(d, nm)
                                                                          ov_dose(d) %>% mutate(dataset = ds, analysis = nm)))
sens_tbl <- bind_rows(
  sens_run(sens_pri, "Primary"),
  ov_dose(dat_sec_l) %>% mutate(dataset = "Secondary", analysis = "Secondary (co-contaminants incl.)"),
  ov_dose(dat_sec)   %>% mutate(dataset = "Secondary", analysis = "Secondary (co-contaminants incl.)")) %>%
  mutate(analysis = factor(analysis, levels = unique(analysis))) %>%
  arrange(factor(scale, levels = scale_levels), dataset, analysis)
sens_tbl %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# [V4] lnRR floor for complete failure: floors of 0.1, 0.5 and 1 (% hatching) and with zero-hatching arms excluded,
# primary and secondary. The zero-hatching arms are listed so they can be checked against the sources.
tab_zero_arms <- dat_all %>% filter(zero_treat) %>%
  select(study_id, reference, contam, cred_r, conc_mg_L, control_mhs, treat_mhs, treat_sd) %>% tibble::as_tibble()
print(tab_zero_arms, n = Inf, width = Inf)
floor_run <- \(d, ds) bind_rows(
  purrr::map_dfr(c(0.1, 0.5, 1), \(fl) ov_dose(make_lnrr(d, fl)) %>% mutate(floor_pct = fl, note = "zero-hatching arms floored")),
  ov_dose(make_lnrr(filter(d, !zero_treat))) %>% mutate(floor_pct = NA_real_, note = "zero-hatching arms excluded")) %>%
  mutate(dataset = ds, .before = 1)
tab_lnrr_floor <- bind_rows(floor_run(dat_pure, "Primary"), floor_run(dat_sec, "Secondary"))
tab_lnrr_floor %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# Random-effects structure x sensitivity set (both scales)
sens_all <- c(
  purrr::imap(sens_pri$lnRR, \(d, nm) list(analysis = nm, data = d)),
  purrr::imap(sens_pri$g,    \(d, nm) list(analysis = nm, data = d)),
  list(list(analysis = "Secondary (co-contaminants incl.)", data = dat_sec_l),
       list(analysis = "Secondary (co-contaminants incl.)", data = dat_sec)))
grid_s <- tidyr::crossing(idx = seq_along(sens_all), struct = c("simple", "species", "phylo", "phylo_sp"))
safe_struct <- function(d, s) tryCatch(ov_dose(d, s), error = function(e) tibble(note = conditionMessage(e)))

sens_struct <- purrr::pmap_dfr(grid_s, \(idx, struct)
                               bind_cols(tibble(analysis = sens_all[[idx]]$analysis, struct = struct),
                                         safe_struct(sens_all[[idx]]$data, struct))) %>%
  arrange(factor(scale, levels = scale_levels), analysis, struct)
sens_struct %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

sens_struct %>% filter(!is.na(est)) %>% group_by(scale, analysis) %>%
  summarise(n_fit = n(), est_range = diff(range(est)), within_range = diff(range(within, na.rm = TRUE)), .groups = "drop")

# 14. Collect and export all results ----------------------------------------
## 14a. Additional checks (moderator p-value grid, aqueous vs non-aqueous, SD/CV screens)
tab_mod_p <- omni %>% filter(struct == "simple") %>%
  mutate(cell = case_when(!is.na(note) ~ "NE", is.na(p_robust) ~ "NE (df<4)", TRUE ~ sprintf("%.3f", p_robust)),
         set  = paste(str_extract(analysis, "^[A-Za-z0-9]+"), if_else(scale == "lnRR", "lnRR", "g"))) %>%
  select(moderator, set, cell) %>% pivot_wider(names_from = set, values_from = cell)
print(tab_mod_p)

dat_nonaq   <- filter(dat_pure,   conc_class != "aqueous")
dat_nonaq_l <- filter(dat_pure_l, conc_class != "aqueous")
tab_medium_ov <- bind_rows(ov(fit_cred(dat_pure_l_dose, "simple"), "aqueous",     dat_pure_l_dose),
                           ov(fit_cred(dat_nonaq_l,     "simple"), "non-aqueous", dat_nonaq_l),
                           ov(fit_cred(dat_pure_dose,   "simple"), "aqueous",     dat_pure_dose),
                           ov(fit_cred(dat_nonaq,       "simple"), "non-aqueous", dat_nonaq))
print(tab_medium_ov)

medium_sets <- list("Primary lnRR" = dat_pure_l, "Primary g" = dat_pure) %>%
  purrr::map(\(d) mutate(d, medium = if_else(conc_class == "aqueous", "aqueous", "non-aqueous")))
tab_medium_test <- purrr::imap_dfr(medium_sets, \(d, nm) bind_cols(tibble(analysis = nm),
                                                                   run_mod("medium", d)$omnibus %>% select(k_es, k_studies, F_robust, df1, df2, p_robust, p_model, reliable, R2_study)))
tab_medium_levels <- purrr::imap_dfr(medium_sets, \(d, nm) bind_cols(tibble(analysis = nm), run_mod("medium", d)$levels))
print(tab_medium_test); print(tab_medium_levels)

# SDs vs % change, by contaminant type (observed SDs only: imputed arms would just return the class CV)
tab_cv_contam <- dat_cred %>% filter(!sd_imputed) %>%
  mutate(cv_ctrl = control_sd / control_mhs, cv_trt = treat_sd / pmax(treat_mhs, 0.5),
         pct_change = 100 * (treat_mhs / control_mhs - 1)) %>%
  group_by(contam) %>%
  summarise(k_es = n(), k_studies = n_distinct(study_id),
            med_cv_ctrl = median(cv_ctrl), med_cv_trt = median(cv_trt),
            med_pct_change = median(pct_change), med_g = median(yi), .groups = "drop")
print(tab_cv_contam)

tab_cv_screen <- dat_cred %>% filter(!sd_imputed) %>% mutate(cv_ctrl = control_sd / control_mhs) %>%
  group_by(study_id, reference) %>% summarise(median_cv = median(cv_ctrl), .groups = "drop") %>%
  arrange(median_cv)

## 14b. Collect (objects that do not exist yet are skipped) and export
nm <- c(flow = "flow", arm_log = "arm_log", class_cv = "class_cv", imputed_arms = "imp_log",
        sd_cap_check = "tab_sd_cap_check", control_check = "tab_control_check", zero_arms = "tab_zero_arms",
        overall = "tab_overall", varcomp = "tab_varcomp", aic = "tab_aic", lrt = "tab_lrt",
        class_test = "tab_class", phylosig = "tab_phylosig", rho = "tab_rho", rho_secondary = "tab_rho_sec",
        loo_overall = "tab_loo_overall", class_descriptive = "tab_class_desc",
        moderators_omnibus = "omni", moderator_levels = "levs", moderators_by_structure = "omni_struct",
        level_support = "tab_level_support", moderator_p_grid = "tab_mod_p",
        taxa_adjusted = "tab_adj", within_study_variation = "tab_within_var",
        within_between = "tab_wb", dose_linear = "tab_dose_linear", dose_quadratic = "tab_dose_quad",
        dose_contam_joint = "tab_dose_contam_joint",
        loo_dose = "tab_loo_dose", outliers_removed = "tab_outliers_removed",
        outlier_overall = "tab_outlier_overall", outlier_dose = "tab_outlier_dose",
        lnrr_floor = "tab_lnrr_floor", cv_screen = "tab_cv_screen", cv_by_contam = "tab_cv_contam",
        medium_overall = "tab_medium_ov", medium_test = "tab_medium_test", medium_levels = "tab_medium_levels",
        sensitivity = "sens_tbl", sensitivity_by_structure = "sens_struct", top3_primary = "tab_top3",
        table1 = "tab1", table1_study_pct = "study_pct", study_pct_summary = "tab_study_pct_summary",
        median_vs_mean = "tab_median_mean", cr_comparison = "cmp",
        table2 = "tab2", table2_contrast = "tab2_contrast", structures_secondary = "tab2s",
        composition = "comp", study_types = "study_types", contrast_relative = "tab_contrast_rel",
        weight_share_secondary = "tab_wshare", flagged_rows_secondary = "tab_flagged_sec",
        sensitivity_primary_vs_secondary = "sens_cmp", loo_secondary = "loo2",
        loo_secondary_summary = "tab_loo2_summary", loo_secondary_top = "tab_loo2_top",
        top3_secondary = "tab_big3_sec", top3_secondary_studies = "tab_big3_sec_desc",
        top3_secondary_composition = "tab_big3_sec_comp", top3_secondary_effect = "tab_big3_sec_red",
        random_draws_secondary = "tab_draws_sec", typical_study_secondary = "tab_typical_sec",
        typical_by_contam = "tab_typical_contam", contam_sensitivity = "tab_contam_sens",
        s38_profile = "tab_s38_profile", s38_vs_others = "tab_s38_vs_others", s38_weight = "tab_s38_weight",
        s38_comparison = "tab_s38_cmp", s38_downsample = "tab_s38_downsample",
        s38_moderators = "tab_s38_moderators", s38_rho = "tab_s38_rho", sd_scenarios = "tab_sd_scenarios",
        imputed_n = "tab_imp_n", imputed_overall = "tab_imp_overall", imputed_dose = "tab_imp_dose",
        imputed_moderators = "tab_imp_mod", imputed_contam = "tab_imp_contam",
        contrast_floor_structure = "tab_contrast_floor", contrast_residual_structure = "tab_contrast_resid",
        secondary_residual_variance = "tab_hetero_overall", polymer_contrasts = "tab_polymer_contrasts",
        polymer_within_study = "tab_polymer_within", polymer_between = "tab_polymer_between",
        polymer_between_n = "tab_polymer_between_n",
        secondary_resid_vs_common = "tab_secondary_resid", contrast_common_variance = "tab2_contrast_common",
        secondary_sensitivity_bytype = "tab_sec_sens", loo_secondary_contrast = "loo2_contrast",
        loo_secondary_contrast_summary = "tab_loo2_contrast_summary",
        contrast_sensitivity_common_variance = "tab_contam_sens_common")

collect_results <- function(print_all = TRUE) {
  res <- purrr::map(nm, \(x) if (exists(x, envir = globalenv())) get(x, envir = globalenv()) else NULL)
  res <- res[!purrr::map_lgl(res, is.null)]
  if (print_all) {
    cat("Not created yet (sections 15-18 run later):", paste(setdiff(names(nm), names(res)), collapse = ", "), "\n")
    purrr::iwalk(res, \(x, n) { cat("\n====", n, "====\n"); print(tibble::as_tibble(x), n = Inf, width = Inf) })
  }
  purrr::iwalk(res, \(x, n) write_csv(tibble::as_tibble(x), file.path(out_dir, paste0(n, ".csv"))))
  saveRDS(res, file.path(out_dir, "CRED_results.rds"))
  invisible(res)
}
results <- collect_results(print_all = TRUE)

# 15. Results section 1: initial models (pure plastic) -----------------------
# lnRR is the primary scale and comes first; Hedges' g is the labelled secondary scale.
d_l <- dat_pure_l; d_g <- dat_pure
fits_l <- fits_P$lnRR; fits_g <- fits_P$g                      # simple, species, phylo, phylo_sp (section 8)
lab_of <- c(simple = "Simple", species = "Species", phylo = "Phylo", phylo_sp = "Phylo + species")

# --- heterogeneity: typical sampling variance, I2, tau (sampling covariances from build_V) ---
het_stats <- function(m, d) {
  V <- build_V(d)
  W <- solve(V); X <- model.matrix(m)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  s2_samp <- (m$k - m$p) / sum(diag(P))
  s2_tot  <- s2_total(m, d)
  tibble(s2_sampling = s2_samp, tau = sqrt(s2_tot), I2_total = 100 * s2_tot / (s2_tot + s2_samp))
}

row_stats <- function(m, d, label) {
  rb  <- rob_safe(m, d$study_id); it <- inf_tbl(m, rb, 1)
  s2  <- \(nm) if (nm %in% m$s.names) m$sigma2[m$s.names == nm] else NA_real_
  pint <- it$est + c(-1, 1) * qt(0.975, m$k - m$p) * sqrt(m$se^2 + s2_total(m, d))   # new effect from a new study (record-averaged residual variance)
  tibble(scale = sc_lab(d), model = label, records = m$k, studies = n_distinct(d$study_id),
         est = it$est, ci_lb = it$ci_lb, ci_ub = it$ci_ub, df = it$df, p_robust = it$p_robust, p_model = it$p_model,
         infer = attr(rb, "infer"),
         resid = resid_of(m),
         s2_study = s2("study_id"), s2_obs = mean(obs_var_vec(m, d)),
         s2_obs_pure = obs_by_type(m)[1], s2_obs_co = obs_by_type(m)[2],
         s2_species = s2("species"), s2_phylo = s2("species_phylo"),
         QE = m$QE, QE_df = m$k - m$p, QE_p = m$QEp,
         pi_lb = pint[1], pi_ub = pint[2], AIC = as.numeric(AIC(m))) %>%
    bind_cols(het_stats(m, d))
}

lrt_vs_simple <- \(fits) c(NA_real_, purrr::map_dbl(fits[-1], \(m) anova(m, fits[[1]])$pval))   # each structure vs Simple (first fit)

tab1 <- bind_rows(purrr::imap_dfr(fits_l, \(m, s) row_stats(m, d_l, lab_of[[s]])),
                  purrr::imap_dfr(fits_g, \(m, s) row_stats(m, d_g, lab_of[[s]])))
tab1$LRT_p <- c(lrt_vs_simple(fits_l), lrt_vs_simple(fits_g))
write_csv(tab1, file.path(out_dir, "table1_initial_models.csv"))

# pooled estimates used by the figures
ov_l <- filter(tab1, scale == "lnRR",      model == "Simple")
ov_g <- filter(tab1, scale == "Hedges' g", model == "Simple")
pct_breaks <- c(-90, -75, -50, -25, 0, 50, 100)
pct_scale  <- scale_x_continuous(breaks = log(1 + pct_breaks / 100), labels = paste0(pct_breaks, "%"))

# --- phylogenetic signal on species means (primary, lnRR): computed in section 8 ---
phy_row <- filter(tab_phylosig, dataset == "Primary", scale == "lnRR")

# --- formatted table ---
dgt <- ifelse(tab1$scale == "lnRR", 3L, 2L)
tab1_fmt <- tab1 %>% transmute(
  scale, Model = model, Records = records, Studies = studies,
  Pct = if_else(scale == "lnRR",
                if_else(is.na(ci_lb), sprintf("%.1f (NE)", pct(est)), sprintf("%.1f (%.1f, %.1f)", pct(est), pct(ci_lb), pct(ci_ub))),
                NA_character_),
  Est = f(est, dgt), CI = fmt_ci(ci_lb, ci_ub, dgt), p_rob = fmt_pr(p_robust), p_mod = fmt_p(p_model),
  s2_study = f(s2_study, dgt), s2_obs = f(s2_obs, dgt), s2_species = f(s2_species, dgt), s2_phylo = f(s2_phylo, dgt),
  tau = f(tau, dgt), I2 = f(I2_total, 1L), QE = sprintf("%.1f (%d)", QE, as.integer(QE_df)),
  PI = if_else(scale == "lnRR", sprintf("%.0f%%, %.0f%%", pct(pi_lb), pct(pi_ub)),
               paste0(f(pi_lb, 2L), ", ", f(pi_ub, 2L))),
  AIC = f(AIC, 1L), LRT = fmt_p(LRT_p))
print(tab1_fmt, width = Inf)

fn <- paste0(
  "Pure plastic exposures, CRED R1/R2 only (", tab1$records[1], " records, ", tab1$studies[1], " studies). ",
  "lnRR (log ratio of treated to control hatching, shown as % change) is the primary effect size; Hedges' g is the secondary scale ",
  "and is not on a common scale across studies because it standardises by the SD between replicates (single embryos in some studies, tanks in others). ",
  "Sampling covariances between records that share a control are computed from sample sizes and variances. ",
  "Robust CIs and p-values are cluster-robust on study (CR2 with Satterthwaite df for Simple; CR1 for the ",
  "models with a species or phylogenetic term, for which CR2 is unavailable) and are reported only where df >= ", df_min,
  " (otherwise NE); the model-based p (residual df) is shown alongside. ",
  if (all(tab1$QE_p < 0.001)) "QE p < 0.001 in every model. " else "See QE p-values in the supplement. ",
  "PI = 95% prediction interval for the effect in a new study (all variance components included). ",
  "LRT p compares each structure with Simple (the only nested comparison common to all). Phylogenetic signal in species mean lnRR: Pagel's lambda = ",
  sprintf("%.3f (p = %.2f)", phy_row$lambda, phy_row$lambda_p), "; Blomberg's K = ", sprintf("%.2f (p = %.2f).", phy_row$K, phy_row$K_p))

tab1_fmt %>% select(-scale) %>%
  kbl(col.names = c("Model", "Records", "Studies", "% change (95% CI)", "Estimate", "95% CI", "p (robust)", "p (model)",
                    "\u03c3\u00b2 study", "\u03c3\u00b2 obs.", "\u03c3\u00b2 species", "\u03c3\u00b2 phylo",
                    "\u03c4", "I\u00b2 (%)", "QE (df)", "95% PI", "AIC", "LRT p"),
      align = c("l", rep("c", 17)), booktabs = TRUE, na = "\u2014",
      caption = paste("Table 1. Overall effect of pure-plastic microplastic exposure on hatching success under four",
                      "random-effects structures. lnRR (primary): log response ratio of treated to control hatching, shown with its",
                      "% change. Hedges' g (secondary): negative values indicate reduced hatching.")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = TRUE, font_size = 12) %>%
  add_header_above(c(" " = 3, "Effect" = 5, "Variance components" = 4, "Heterogeneity" = 4, "Model fit" = 2)) %>%
  row_spec(0, bold = TRUE) %>% row_spec(c(1, 5), bold = TRUE) %>%
  pack_rows("Log response ratio (lnRR; primary scale)", 1, 4) %>% pack_rows("Hedges' g (secondary scale)", 5, 8) %>%
  footnote(general = fn, general_title = "", footnote_as_chunk = TRUE)

# --- Fig. 1: study-level forest, % change in hatching (lnRR, panel a) and g (panel b) ---
study_est <- function(d) {       # inverse-variance mean of a study's records; SE uses the shared-control covariances (build_V)
  V <- build_V(d); w <- (1 / d$vi) / sum(1 / d$vi)
  tibble(est = sum(w * d$yi), se = sqrt(drop(t(w) %*% V %*% w)), k = nrow(d), reference = d$reference[1])
}
lab_fix <- \(x) str_replace(x, "Gonz.{1,3}lez-Doncel", "Gonz\u00e1lez-Doncel")   # repair garbled accent in the source file

se_l <- d_l %>% group_by(study_id) %>% group_modify(~ study_est(.x)) %>% ungroup() %>%
  mutate(label = paste0(lab_fix(reference), " (k=", k, ")"))
se_g <- d_g %>% group_by(study_id) %>% group_modify(~ study_est(.x)) %>% ungroup() %>%
  mutate(label = se_l$label[match(study_id, se_l$study_id)])

# --- class per study ---
cls_tab <- d_g %>% group_by(study_id) %>%
  summarise(n_class = n_distinct(taxa_class),
            class = names(sort(table(taxa_class), decreasing = TRUE))[1], .groups = "drop") %>%
  mutate(label = se_l$label[match(study_id, se_l$study_id)])
filter(cls_tab, n_class > 1)       # should be empty (one class per study)

class_pal <- c(Actinopterygii = "#0072B2", Crustacea = "#E69F00", Mollusca = "#009E73",
               Amphibia = "#CC79A7", Nematoda = "#D55E00", Oligochaeta = "#56B4E9")   # colour-blind safe
stopifnot(all(unique(dat$taxa_class) %in% names(class_pal)))                          # a new class needs a colour

# --- order: class blocks (top to bottom as listed), then lnRR within each block (most negative at the top) ---
class_order <- c("Actinopterygii", "Crustacea", "Mollusca", "Amphibia", "Nematoda", "Oligochaeta")
stopifnot(all(unique(cls_tab$class) %in% class_order))

se_cls <- se_l %>% left_join(select(cls_tab, study_id, class), by = "study_id")

# lev runs bottom to top (Overall = 0), so the last class and the least negative studies come first
lev <- c("Overall",
         se_cls %>% arrange(desc(match(class, class_order)), desc(est)) %>% pull(label))

rec_l <- d_l %>% mutate(label = se_l$label[match(study_id, se_l$study_id)]) %>% select(label, yi)
rec_g <- d_g %>% mutate(label = se_l$label[match(study_id, se_l$study_id)]) %>% select(label, yi)

ypos <- tibble(label = lev, y = seq_along(lev) - 1)     # Overall = 0 (bottom); studies stack upward

# horizontal separators between class blocks (and between the studies and "Overall")
block_breaks <- ypos %>% left_join(distinct(se_cls, label, class), by = "label") %>%
  arrange(y) %>% mutate(new_block = class != lag(class)) %>%
  filter(new_block) %>% pull(y) - 0.5
block_breaks <- c(0.5, block_breaks)

mk_df <- function(se_tab, ov, rec) bind_rows(
  se_tab %>% transmute(label, kind = "study", est, lo = est - 1.96 * se, hi = est + 1.96 * se,
                       w = pmin(1 / se^2, quantile(1 / se^2, 0.9))),
  rec %>% transmute(label, kind = "record", est = yi),
  tibble(label = "Overall", kind = "overall", est = ov$est, lo = ov$ci_lb, hi = ov$ci_ub,
         pi_lo = ov$pi_lb, pi_hi = ov$pi_ub)) %>%
  left_join(select(cls_tab, label, class), by = "label") %>%
  left_join(ypos, by = "label")

make_forest <- function(df, xlab, xlim, x_scale = NULL, lab_fun = \(x) sprintf("%.1f", x)) {
  st <- filter(df, kind == "study"); ov <- filter(df, kind == "overall"); rc <- filter(df, kind == "record")
  gap    <- 0.03 * diff(xlim)
  off_lo <- filter(st, est < xlim[1]); off_hi <- filter(st, est > xlim[2])
  p <- ggplot(df, aes(y = y)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(data = rc, aes(x = est, colour = class), shape = 1, size = 1.3, alpha = 0.8, show.legend = FALSE,
               position = position_jitter(height = 0.2, width = 0, seed = 1)) +
    geom_linerange(data = st, aes(xmin = lo, xmax = hi, colour = class), linewidth = 0.5, show.legend = FALSE) +
    geom_point(data = st, aes(x = est, size = w, colour = class), shape = 15) +
    geom_linerange(data = ov, aes(xmin = pi_lo, xmax = pi_hi), colour = "black", linewidth = 0.4) +
    geom_linerange(data = ov, aes(xmin = lo, xmax = hi), colour = "black", linewidth = 1.6) +
    geom_point(data = ov, aes(x = est), shape = 18, size = 4.5, colour = "black") +
    geom_segment(data = off_lo, aes(x = xlim[1] + 2.5 * gap, xend = xlim[1] + 0.3 * gap, yend = y, colour = class), show.legend = FALSE,
                 arrow = arrow(length = unit(0.08, "in"), type = "closed")) +
    geom_text(data = off_lo, aes(x = xlim[1] + 3 * gap, label = lab_fun(est), colour = class), hjust = 0, size = 2.7, show.legend = FALSE) +
    geom_segment(data = off_hi, aes(x = xlim[2] - 2.5 * gap, xend = xlim[2] - 0.3 * gap, yend = y, colour = class), show.legend = FALSE,
                 arrow = arrow(length = unit(0.08, "in"), type = "closed")) +
    geom_text(data = off_hi, aes(x = xlim[2] - 3 * gap, label = lab_fun(est), colour = class), hjust = 1, size = 2.7, show.legend = FALSE) +
    scale_colour_manual(values = class_pal, name = "Taxonomic class", drop = TRUE) +
    scale_size_continuous(range = c(1, 2.4), guide = "none") +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    scale_y_continuous(breaks = ypos$y, labels = ypos$label, expand = expansion(add = 0.7)) +
    coord_cartesian(xlim = xlim) + labs(x = xlab, y = NULL) + theme_classic(base_size = 10)
  if (!is.null(x_scale)) p <- p + x_scale
  p
}

p_l <- make_forest(mk_df(se_l, ov_l, rec_l), "Change in hatching relative to control (log scale)",
                   xlim = log(1 + c(-80, 50) / 100), x_scale = pct_scale,
                   lab_fun = \(x) sprintf("%.0f%%", 100 * (exp(x) - 1)))
p_g <- make_forest(mk_df(se_g, ov_g, rec_g), "Hedges' g (negative = reduced hatching; secondary scale)", c(-8, 4)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  guides(colour = "none")          # one legend (from panel a)

fig1 <- (p_l + p_g + plot_layout(widths = c(1.35, 1), guides = "collect") +
           plot_annotation(tag_levels = "a")) & theme(legend.position = "bottom")
fig1
ggsave(file.path(out_dir, "fig1_primary_forest.pdf"), fig1, width = 10, height = 7.5)
ggsave(file.path(out_dir, "fig1_primary_forest.png"), fig1, width = 10, height = 7.5, dpi = 150)

# 15b. Fill-ins and checks for the Section 1 paragraph ----------------------

## 1. Bracketed values: QE, I2, tau, prediction interval, and the % version (lnRR)
tab1 %>% filter(model == "Simple") %>%
  transmute(scale, records, studies, est = round(est, 3), QE = round(QE, 1), QE_df, QE_p = signif(QE_p, 2),
            I2_total = round(I2_total, 1), tau = round(tau, 2),
            PI = if_else(scale == "lnRR", sprintf("%.0f%% to %.0f%%", pct(pi_lb), pct(pi_ub)),
                         sprintf("%.2f to %.2f", pi_lb, pi_ub))) %>%
  print(width = Inf)

## 2. Phylogenetic signal (all datasets and scales) and the table footnote
print(tab_phylosig)
cat(fn, "\n")

## 3. CR1 vs CR2: same estimate, three kinds of p, for each structure and scale
cmp <- purrr::imap_dfr(list(lnRR = list(fits_l, d_l), `Hedges' g` = list(fits_g, d_g)), \(x, sc)
                       purrr::imap_dfr(x[[1]], \(m, s) {
                         cr1 <- robust(m, cluster = x[[2]]$study_id)
                         tibble(scale = sc, model = lab_of[[s]], est = as.numeric(coef(m)), p_model = m$pval, p_CR1 = cr1$pval,
                                CR1_lb = cr1$ci.lb, CR1_ub = cr1$ci.ub,
                                p_CR2 = if (s == "simple") rob(m, x[[2]])$pval else NA_real_)
                       }))
cmp %>% mutate(across(where(is.numeric), ~ round(.x, 4))) %>% print(width = Inf)

## 4. Study-level median change in hatching (raw scale; independent of the effect-size scale)
study_pct <- d_g %>%
  mutate(pct_change = 100 * (treat_mhs / control_mhs - 1)) %>%
  group_by(study_id, reference) %>%
  summarise(k_es = n(), med_pct = median(pct_change), .groups = "drop") %>% arrange(med_pct)
print(study_pct, n = Inf)

tab_study_pct_summary <- study_pct %>% summarise(n_studies = n(), median = median(med_pct),
                                                 q25 = quantile(med_pct, 0.25), q75 = quantile(med_pct, 0.75),
                                                 n_below_m10 = sum(med_pct < -10), n_below_m25 = sum(med_pct < -25),
                                                 n_positive = sum(med_pct > 0))
print(tab_study_pct_summary)

## 5. Median vs mean on the same raw scale, and how much a few large reductions move the mean
tab_median_mean <- d_g %>% mutate(pct = 100 * (treat_mhs / control_mhs - 1)) %>%
  summarise(level = "records", n = n(), median = median(pct), mean = mean(pct)) %>%
  bind_rows(study_pct %>% summarise(level = "studies", n = n(), median = median(med_pct), mean = mean(med_pct)),
            study_pct %>% filter(med_pct > -25) %>%
              summarise(level = "studies, excl. those below -25%", n = n(), median = median(med_pct), mean = mean(med_pct)))
print(tab_median_mean)

# 16. Results section 2: co-contaminant rows (secondary) ---------------------
d2_l     <- dat_sec_l;  d2_g   <- dat_sec
dat_co_l <- filter(d2_l, contam == "Co-contaminant");  dat_co <- filter(d2_g, contam == "Co-contaminant")

# composition: records and studies per type, and how many studies contribute both
comp <- d2_g %>% group_by(contam) %>% summarise(records = n(), studies = n_distinct(study_id), .groups = "drop")
study_types <- d2_g %>% group_by(study_id) %>%
  summarise(type = case_when(n_distinct(contam) == 2 ~ "both types",
                             first(contam) == "Pure plastic" ~ "pure plastic only",
                             TRUE ~ "co-contaminant only"), .groups = "drop") %>%
  count(type, name = "studies")

# ---- Table 2 [V5]: lnRR first. Rows: the two contaminant types from the SECONDARY joint model (shared study variance, separate
# observation-level variance per type), then the pooled, within-study and between-study contrasts from that same model, so the by-type
# effects, the contrasts and Fig. 2 all come from one model. The standalone pure-plastic (Table 1) and co-contaminant-only fits are not
# Table 2 rows. For the contrast rows, Records / Studies are the data informing each estimate (within-study: only studies with both types).
t2_levels <- c("Secondary joint model: Pure plastic", "Secondary joint model: Co-contaminant")
tab2 <- purrr::map_dfr(c("lnRR", "g"), \(sc) {
  d <- if (sc == "lnRR") d2_l else d2_g
  bytype(d) %>% mutate(model = paste0("Secondary joint model: ", type))
}) %>%
  mutate(model = factor(model, levels = t2_levels)) %>%
  arrange(factor(scale, levels = scale_levels), model)

# ---- contrasts: pooled, within-study and between-study (co-contaminant minus pure plastic); lnRR first ----
# Default = separate residual variance by type (the secondary model); the common-variance version is saved as a labelled sensitivity.
tab2_contrast        <- bind_rows(contam_contrast(d2_l), contam_contrast(d2_g)) %>% rename(model = contrast)
tab2_contrast_common <- bind_rows(contam_contrast(d2_l, resid = "common"), contam_contrast(d2_g, resid = "common")) %>% rename(model = contrast)
write_csv(tab2, file.path(out_dir, "table2_datasets.csv"))
write_csv(tab2_contrast, file.path(out_dir, "table2_contrasts.csv"))
write_csv(tab2_contrast_common, file.path(out_dir, "table2_contrasts_common_variance.csv"))

# ---- format and render ----
d3 <- ifelse(tab2$scale == "lnRR", 3L, 2L)
t2_main <- tab2 %>% transmute(scale, Model = as.character(model), Records = records, Studies = studies,
                              Pct = if_else(scale == "lnRR",
                                            if_else(is.na(ci_lb), sprintf("%.1f (NE)", pct(est)),
                                                    sprintf("%.1f (%.1f, %.1f)", pct(est), pct(ci_lb), pct(ci_ub))), NA_character_),
                              Est = f(est, d3), CI = fmt_ci(ci_lb, ci_ub, d3), p_rob = fmt_pr(p_robust), p_mod = fmt_p(p_model),
                              s2_study = f(s2_study, d3), s2_obs = f(s2_obs, d3), tau = f(tau, d3),
                              PI = if_else(scale == "lnRR", sprintf("%.0f%%, %.0f%%", pct(pi_lb), pct(pi_ub)),
                                           paste0(f(pi_lb, 2L), ", ", f(pi_ub, 2L))))
dc <- ifelse(tab2_contrast$scale == "lnRR", 3L, 2L)
t2_con <- tab2_contrast %>% transmute(scale, Model = model, Records = records, Studies = studies,
                                      Pct = if_else(scale == "lnRR",
                                                    if_else(is.na(ci_lb), sprintf("%.1f (NE)", pct(est)),
                                                            sprintf("%.1f (%.1f, %.1f)", pct(est), pct(ci_lb), pct(ci_ub))), NA_character_),
                                      Est = f(est, dc), CI = fmt_ci(ci_lb, ci_ub, dc), p_rob = fmt_pr(p_robust), p_mod = fmt_p(p_model),
                                      s2_study = NA_character_, s2_obs = NA_character_, tau = NA_character_, PI = NA_character_)
tab2_fmt <- bind_rows(t2_main, t2_con) %>% arrange(factor(scale, levels = scale_levels))
print(tab2_fmt, width = Inf)

n_l <- sum(tab2_fmt$scale == "lnRR")
tab2_fmt %>% select(-scale) %>%
  kbl(col.names = c("Dataset / contrast", "Records", "Studies", "% change (95% CI)", "Estimate", "95% CI", "p (robust)", "p (model)",
                    "\u03c3\u00b2 study", "\u03c3\u00b2 obs.", "\u03c4", "95% PI"),
      align = c("l", rep("c", 11)), booktabs = TRUE, na = "\u2014",
      caption = paste("Table 2. Effect of microplastic exposure on hatching success for each contaminant type in the secondary joint model,",
                      "with the contrast between types.")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = TRUE, font_size = 12) %>%
  add_header_above(c(" " = 3, "Effect" = 5, "Variance" = 2, "Heterogeneity" = 2)) %>%
  row_spec(0, bold = TRUE) %>%
  pack_rows("Log response ratio (lnRR; primary scale)", 1, n_l) %>%
  pack_rows("Hedges' g (secondary scale)", n_l + 1, nrow(tab2_fmt)) %>%
  footnote(general = paste0("Simple multilevel model; sampling covariances between records that share a control computed from sample sizes and variances. ",
                            "The secondary joint model estimates each contaminant type with a shared study variance and a separate observation-level variance for pure-plastic ",
                            "and co-contaminant records (the sigma-squared obs. column is the type-specific value; tau is the SD of true effects for that type). The model was fitted to all ",
                            nrow(d2_g), " records from ", n_distinct(d2_g$study_id), " studies. ",
                            "Robust (CR2, cluster = study) CIs and p are shown only where df >= ", df_min, " (otherwise NE); model-based p alongside. ",
                            "Contrast rows are co-contaminant minus pure plastic from a model with contaminant type as a moderator and the same separate residual variances ",
                            "(negative = larger reduction with co-contaminants); on the lnRR scale the estimate is a difference in log response ratio and the % column is the ",
                            "relative difference in the hatching ratio. Records / Studies show the data informing each estimate: the pooled and between-study contrasts use all ",
                            "records, whereas the within-study contrast is informed only by the studies that tested both contaminant types. ",
                            "PI = 95% prediction interval for a new record of that type. ",
                            "Hedges' g is not on a common scale across studies (see Table 1)."),
           general_title = "", footnote_as_chunk = TRUE)

# ---- supplementary: the four random-effects structures on the secondary dataset ----
fmt_struct <- function(tab) {
  dgt <- ifelse(tab$scale == "lnRR", 3L, 2L)
  tab %>% transmute(scale, Model = model, Records = records, Studies = studies,
                    Pct = if_else(scale == "lnRR",
                                  if_else(is.na(ci_lb), sprintf("%.1f (NE)", pct(est)),
                                          sprintf("%.1f (%.1f, %.1f)", pct(est), pct(ci_lb), pct(ci_ub))), NA_character_),
                    Est = f(est, dgt), CI = fmt_ci(ci_lb, ci_ub, dgt), p_rob = fmt_pr(p_robust), p_mod = fmt_p(p_model),
                    s2_study = f(s2_study, dgt),
                    s2_obs = if_else(is.na(s2_obs_co), f(s2_obs, dgt), paste0(f(s2_obs_pure, dgt), " / ", f(s2_obs_co, dgt))),
                    s2_species = f(s2_species, dgt), s2_phylo = f(s2_phylo, dgt),
                    tau = f(tau, dgt), I2 = f(I2_total, 1L), QE = sprintf("%.1f (%d)", QE, as.integer(QE_df)),
                    PI = if_else(scale == "lnRR", sprintf("%.0f%%, %.0f%%", pct(pi_lb), pct(pi_ub)),
                                 paste0(f(pi_lb, 2L), ", ", f(pi_ub, 2L))),
                    AIC = f(AIC, 1L), LRT = fmt_p(LRT_p))
}
fits2_l <- fits_S$lnRR; fits2_g <- fits_S$g
tab2s <- bind_rows(purrr::imap_dfr(fits2_l, \(m, s) row_stats(m, d2_l, lab_of[[s]])),
                   purrr::imap_dfr(fits2_g, \(m, s) row_stats(m, d2_g, lab_of[[s]])))
tab2s$LRT_p <- c(lrt_vs_simple(fits2_l), lrt_vs_simple(fits2_g))
write_csv(tab2s, file.path(out_dir, "table2s_structures_secondary.csv"))

# 16.3 Fig. 2 ----------------------------------------------------------------
# Layout: pure plastic (top block), co-contaminant (bottom block); studies ordered by lnRR within each block
st_type <- \(d) d %>% group_by(study_id, contam) %>% group_modify(~ study_est(.x)) %>% ungroup()
set_l <- st_type(d2_l); set_g <- st_type(d2_g)

n_co    <- n_distinct(set_l$study_id[set_l$contam == "Co-contaminant"])
both_id <- intersect(set_l$study_id[set_l$contam == "Pure plastic"],
                     set_l$study_id[set_l$contam == "Co-contaminant"])

lay <- bind_rows(
  set_l %>% filter(contam == "Co-contaminant") %>% arrange(desc(est)) %>% mutate(y = row_number()),            # least negative lowest
  set_l %>% filter(contam == "Pure plastic")   %>% arrange(desc(est)) %>% mutate(y = n_co + 2 + row_number())
) %>%
  mutate(label = paste0(lab_fix(reference), " (k=", k, ")", if_else(study_id %in% both_id, " \u2020", ""))) %>%
  select(study_id, contam, y, label)

y_pool <- c("Co-contaminant" = 0, "Pure plastic" = n_co + 2)
pal2   <- c("Pure plastic" = "#0072B2", "Co-contaminant" = "#D55E00", "All records" = "black")

pool_df <- \(sc) tab2 %>% filter(scale == sc, grepl("^Secondary joint", model)) %>%     # by-type effects from the secondary joint model
  transmute(kind = "pooled",
            contam = if_else(grepl("Pure", model), "Pure plastic", "Co-contaminant"),
            y = unname(y_pool[contam]),
            est, lo = ci_lb, hi = ci_ub, pi_lo = pi_lb, pi_hi = pi_ub)

add_y <- \(d) d %>% mutate(contam = as.character(contam)) %>%
  left_join(select(lay, study_id, contam, y), by = c("study_id", "contam"))

mk_df2 <- \(st, d, sc) bind_rows(
  add_y(st) %>% transmute(kind = "study", contam, y, est, lo = est - 1.96 * se, hi = est + 1.96 * se,
                          w = pmin(1 / se^2, quantile(1 / se^2, 0.9))),
  add_y(d)  %>% transmute(kind = "record", contam, y, est = yi),
  pool_df(sc))

ybr  <- c(y_pool[["Co-contaminant"]], y_pool[["Pure plastic"]], lay$y)
ylab <- c("Co-contaminant (joint model)", "Pure plastic (joint model)", lay$label)

make_forest2 <- function(df, xlab, xlim, x_scale = NULL, lab_fun = \(x) sprintf("%.1f", x)) {
  st <- filter(df, kind == "study"); rc <- filter(df, kind == "record"); po <- filter(df, kind == "pooled")
  gap <- 0.03 * diff(xlim)
  off_lo <- filter(st, est < xlim[1]); off_hi <- filter(st, est > xlim[2])
  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 2, colour = "grey85") +
    geom_point(data = rc, show.legend = FALSE, aes(x = est, y = y, colour = contam), shape = 1, size = 1.1, alpha = 0.7,
               position = position_jitter(height = 0.06, width = 0, seed = 1)) +
    geom_linerange(data = st, show.legend = FALSE, aes(y = y, xmin = lo, xmax = hi, colour = contam), linewidth = 0.45) +
    geom_point(data = st, aes(x = est, y = y, size = w, colour = contam), shape = 15) +
    geom_linerange(data = po, show.legend = FALSE, aes(y = y, xmin = pi_lo, xmax = pi_hi, colour = contam), linewidth = 0.4) +
    geom_linerange(data = po, show.legend = FALSE, aes(y = y, xmin = lo, xmax = hi, colour = contam), linewidth = 1.6) +
    geom_point(data = po, show.legend = FALSE, aes(x = est, y = y, colour = contam), shape = 18, size = 4.5) +
    geom_segment(data = off_lo, show.legend = FALSE, aes(x = xlim[1] + 2.5 * gap, xend = xlim[1] + 0.3 * gap, y = y, yend = y, colour = contam),
                 arrow = arrow(length = unit(0.07, "in"), type = "closed")) +
    geom_text(data = off_lo, show.legend = FALSE, aes(x = xlim[1] + 3 * gap, y = y, label = lab_fun(est), colour = contam), hjust = 0, size = 2.5) +
    geom_segment(data = off_hi, show.legend = FALSE, aes(x = xlim[2] - 2.5 * gap, xend = xlim[2] - 0.3 * gap, y = y, yend = y, colour = contam),
                 arrow = arrow(length = unit(0.07, "in"), type = "closed")) +
    geom_text(data = off_hi, show.legend = FALSE, aes(x = xlim[2] - 3 * gap, y = y, label = lab_fun(est), colour = contam), hjust = 1, size = 2.5) +
    scale_colour_manual(values = pal2, breaks = c("Pure plastic", "Co-contaminant"), name = NULL) +
    scale_size_continuous(range = c(1, 2.4), guide = "none") +
    scale_y_continuous(breaks = ybr, labels = ylab, expand = expansion(add = 0.7)) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    coord_cartesian(xlim = xlim) + labs(x = xlab, y = NULL) + theme_classic(base_size = 9)
  if (!is.null(x_scale)) p <- p + x_scale
  p
}

p_l2 <- make_forest2(mk_df2(set_l, d2_l, "lnRR"), "Change in hatching relative to control (log scale)",
                     xlim = log(1 + c(-95, 100) / 100), x_scale = pct_scale,
                     lab_fun = \(x) sprintf("%.0f%%", 100 * (exp(x) - 1)))
p_g2 <- make_forest2(mk_df2(set_g, d2_g, "Hedges' g"), "Hedges' g (negative = reduced hatching; secondary scale)", c(-28, 4)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  guides(colour = "none")

fig2 <- (p_l2 + p_g2 + plot_layout(widths = c(1.45, 1), guides = "collect") +
           plot_annotation(tag_levels = "a")) & theme(legend.position = "bottom")
fig2
ggsave(file.path(out_dir, "fig2_secondary_forest.pdf"), fig2, width = 10, height = 9)
ggsave(file.path(out_dir, "fig2_secondary_forest.png"), fig2, width = 10, height = 9, dpi = 150)

# 16.4 Sensitivity sets (pct_rank and make_sens are defined in section 13) -----
big_sec  <- pct_rank(dat_sec)$study_id[1:3]
sens_sec <- make_sens(dat_sec, dat_all, big_sec)

run_sens <- \(sens, ds) purrr::imap_dfr(sens, \(sets, sc) purrr::imap_dfr(sets, \(d, nm)
                                                                          red(d, nm) %>% mutate(dataset = ds)))
sens_cmp <- bind_rows(run_sens(sens_pri, "Primary: pure plastic"),
                      run_sens(sens_sec, "Secondary: overall mean (separate residual variances)")) %>%
  mutate(analysis = factor(analysis, levels = unique(analysis)))
write_csv(sens_cmp, file.path(out_dir, "sensitivity_primary_vs_secondary.csv"))

# cell formatter shared by the two sensitivity tables
cell <- \(est, lb, ub, as_pct) case_when(
  is.na(est) ~ NA_character_,
  is.na(lb)  ~ if (as_pct) sprintf("%.1f (NE)", pct(est)) else sprintf("%.2f (NE)", est),
  TRUE       ~ if (as_pct) sprintf("%.1f (%.1f, %.1f)", pct(est), pct(lb), pct(ub)) else sprintf("%.2f (%.2f, %.2f)", est, lb, ub))

# ---- PRIMARY sensitivity table: lnRR (primary) then g; sets that exist for one scale only show a dash on the other ----
# sens table 1 (primary) ------------------------------------------------------
res <- "results_39studies_v5"
rd  <- function(f) read_csv(file.path(res, f), show_col_types = FALSE)

pct <- function(x) 100 * (exp(x) - 1)
fp  <- function(p) ifelse(is.na(p), "NE", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
ci_pct <- function(e, l, u) ifelse(is.na(l), "NE", sprintf("%.1f (%.1f, %.1f)", pct(e), pct(l), pct(u)))
ci_g   <- function(e, l, u) ifelse(is.na(l), "NE", sprintf("%.2f (%.2f, %.2f)", e, l, u))

sens  <- rd("sensitivity.csv") %>% filter(dataset == "Primary")    # one row per analysis and scale
sdscn <- rd("sd_scenarios.csv")                                    # CV multipliers, Study 38 (Prata) variants
sbs   <- rd("sensitivity_by_structure.csv")                        # random-effect structure grid
rho   <- rd("rho.csv") %>% filter(dataset == "Primary")            # fixed-rho covariance assumptions
loo   <- rd("loo_overall.csv")                                     # leave-one-study-out refits

dash <- "\u2013"
cell_pct <- \(x) if (nrow(x) == 0) dash else ci_pct(x$est, x$ci_lb, x$ci_ub)
cell_g   <- \(x) if (nrow(x) == 0) dash else ci_g(x$est, x$ci_lb, x$ci_ub)
cell_n   <- \(x) if (nrow(x) == 0) dash else paste0(x$k_es, " / ", x$k_studies)
cell_p   <- \(x) if (nrow(x) == 0) dash else fp(x$p_robust)

# one row from any table with (scale, analysis, k_es, k_studies, est, ci_lb, ci_ub, p_robust)
row_s <- function(lnRR_name, label, group, g_name = lnRR_name, src = sens) {
  sl <- src %>% filter(scale == "lnRR",      analysis == lnRR_name)
  sg <- src %>% filter(scale == "Hedges' g", analysis == g_name)
  tibble(Group = group, Analysis = label,
         N_l = cell_n(sl), pc = cell_pct(sl), p_l = cell_p(sl),
         N_g = cell_n(sg), g  = cell_g(sg),   p_g = cell_p(sg))
}

# fixed within-study correlation (model-based CIs only, so no robust p)
base_l <- sens %>% filter(scale == "lnRR",      analysis == "Base dataset (CRED)")
base_g <- sens %>% filter(scale == "Hedges' g", analysis == "Base dataset (CRED)")
row_rho <- function(r, label, group) {
  rl <- rho %>% filter(scale == "lnRR",      method == "fixed rho", rho == r)
  rg <- rho %>% filter(scale == "Hedges' g", method == "fixed rho", rho == r)
  tibble(Group = group, Analysis = label,
         N_l = cell_n(base_l), pc = cell_pct(rl), p_l = dash,
         N_g = cell_n(base_g), g  = cell_g(rg),   p_g = dash)
}

# alternative random-effect structure (base dataset)
row_st <- function(struct, label, group) {
  sl <- sbs %>% filter(analysis == "Base dataset (CRED)", struct == !!struct, scale == "lnRR")
  sg <- sbs %>% filter(analysis == "Base dataset (CRED)", struct == !!struct, scale == "Hedges' g")
  tibble(Group = group, Analysis = label,
         N_l = cell_n(sl), pc = cell_pct(sl), p_l = cell_p(sl),
         N_g = cell_n(sg), g  = cell_g(sg),   p_g = cell_p(sg))
}

# leave-one-study-out: range of the pooled estimate and how many refits stay significant
lo <- function(sc) loo %>% filter(scale == sc) %>%
  summarise(lo = min(est), hi = max(est), n_sig = sum(p_robust < 0.05, na.rm = TRUE),
            n = n(), ks = max(k_studies))
lo_l <- lo("lnRR"); lo_g <- lo("Hedges' g")

rows <- bind_rows(
  row_s("Base dataset (CRED)", "Primary analysis (pure plastic)", "Primary"),
  
  row_rho(0.2, "Fixed correlation \u03c1 = 0.2", "Covariance assumption"),
  row_rho(0.5, "Fixed correlation \u03c1 = 0.5", "Covariance assumption"),
  row_rho(0.8, "Fixed correlation \u03c1 = 0.8", "Covariance assumption"),
  
  row_s("R3 retained",         "Retain CRED R3 records",          "Data inclusion"),
  row_s("Excl. PlotDigitiser", "Exclude digitised data",          "Data inclusion"),
  
  row_s("Excl. imputed SDs",   "Exclude studies with imputed SDs", "SD handling"),
  row_s("Constant class-mean SD (old rule, all imputed studies)",
        "Constant class-mean SD instead of CV rule",             "SD handling"),
  row_s("CV x0.5", "Imputed CV \u00d7 0.5",                        "SD handling", src = sdscn),
  row_s("CV x2",   "Imputed CV \u00d7 2",                          "SD handling", src = sdscn),
  row_s("Excl. Study 38", "Exclude Prata et al. (2022)",           "SD handling", src = sdscn),
  row_s("Study 38 SD from class CV (earlier rule)",
        "Prata et al. (2022): class-CV SD instead of binomial SD", "SD handling", src = sdscn),
  
  row_st("species", "Add species random effect (CR1)",           "Model structure"),
  row_st("phylo",   "Add phylogenetic random effect (CR1)",      "Model structure"),
  
  row_s("Excl. |z| > 3 (lnRR model)", "Exclude |studentised residual| > 3", "Outliers and influence",
        g_name = "Excl. |z| > 3 (g model)"),
  tibble(Group = "Outliers and influence", Analysis = "Leave-one-study-out (range)",
         N_l = paste0(dash, " / ", lo_l$ks),
         pc  = sprintf("%.1f to %.1f", pct(lo_l$lo), pct(lo_l$hi)),
         p_l = sprintf("%d/%d significant", lo_l$n_sig, lo_l$n),
         N_g = paste0(dash, " / ", lo_g$ks),
         g   = sprintf("%.2f to %.2f", lo_g$lo, lo_g$hi),
         p_g = sprintf("%d/%d significant", lo_g$n_sig, lo_g$n)),
  row_s("Excl. 3 largest-effect studies (post hoc)",
        "Exclude the three largest-effect studies (post hoc)",   "Outliers and influence")
)

grp <- rle(rows$Group)

kb <- kbl(rows %>% select(-Group), format = "html", booktabs = TRUE, align = c("l", rep("r", 6)),
          caption = "Table S[X]. Sensitivity of the overall effect of pure-plastic exposure on hatching success to analytical decisions.",
          col.names = c("Analysis", "Records / studies", "% change (95% CI)", "p",
                        "Records / studies", "Hedges' g (95% CI)", "p")) %>%
  add_header_above(c(" " = 1, "lnRR (primary)" = 3, "Hedges' g (secondary)" = 3)) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE, position = "left") %>%
  row_spec(0, bold = TRUE) %>% row_spec(1, bold = TRUE)

start <- 1
for (i in seq_along(grp$lengths)) {
  kb <- pack_rows(kb, grp$values[i], start, start + grp$lengths[i] - 1); start <- start + grp$lengths[i]
}

kb %>% footnote(general = paste(
  "Estimates are from the primary (pure-plastic) multilevel model with study and observation random effects unless stated.",
  "lnRR is shown as % change in hatching relative to controls; Hedges' g is the secondary scale.",
  "p values are cluster-robust (CR2, Satterthwaite df); NE = not estimable (df < 4); the species and phylogenetic models use CR1, which is anti-conservative.",
  "Fixed-\u03c1 rows show model-based confidence intervals only and are not given a robust p.",
  "Outlier rules were applied separately on each scale (lnRR: 99 records; g: 98 records).",
  "The primary dataset contains no zero-hatching records, so the zero-hatching floor and exclusion analyses leave it unchanged and are not shown.",
  "Leave-one-study-out gives the range of pooled estimates across the refits, each omitting one study."),
  general_title = "Note: ", footnote_as_chunk = TRUE)

# ---- [V5] SECONDARY sensitivity table: effect of each contaminant type and the within-study contrast, under the secondary model
#      (separate residual variance by type). The last rows repeat the base dataset under the earlier common-variance model. ----
sec_row <- function(d, nm, resid = "auto") {
  bt <- bytype(d, resid = resid); cc <- contam_contrast(d, resid = resid)
  w <- filter(cc, grepl("^Within", contrast)); pu <- filter(bt, type == "Pure plastic"); co <- filter(bt, type == "Co-contaminant")
  tibble(scale = sc_lab(d), analysis = nm, records = nrow(d), studies = n_distinct(d$study_id), resid = bt$resid[1],
         pure_est = pu$est, pure_lb = pu$ci_lb, pure_ub = pu$ci_ub, pure_p = pu$p_robust,
         co_est = co$est, co_lb = co$ci_lb, co_ub = co$ci_ub, co_p = co$p_robust,
         within_est = w$est, within_lb = w$ci_lb, within_ub = w$ci_ub, within_df = w$df,
         within_p = w$p_robust, within_p_model = w$p_model, infer = w$infer)
}
safe_sec_row <- function(d, nm, resid = "auto") {
  r <- purrr::possibly(sec_row, otherwise = NULL)(d, nm, resid)
  if (is.null(r)) tibble(scale = sc_lab(d), analysis = nm, note = "not estimable") else r
}
tab_sec_sens <- bind_rows(
  purrr::imap_dfr(sens_sec$lnRR, \(d, nm) safe_sec_row(d, nm)),
  purrr::imap_dfr(sens_sec$g,    \(d, nm) safe_sec_row(d, nm)),
  safe_sec_row(dat_sec_l, "Base dataset, common residual variance (earlier model)", "common"),
  safe_sec_row(dat_sec,   "Base dataset, common residual variance (earlier model)", "common")) %>%
  mutate(analysis = factor(analysis, levels = unique(analysis)))
if (!"note" %in% names(tab_sec_sens)) tab_sec_sens$note <- NA_character_      # no "not estimable" rows -> no note column
write_csv(tab_sec_sens, file.path(out_dir, "secondary_sensitivity_bytype.csv"))

sec_wide <- tab_sec_sens %>% filter(is.na(note)) %>% mutate(sc = if_else(scale == "lnRR", "l", "g")) %>%
  select(analysis, sc, records, studies, pure_est, pure_lb, pure_ub, co_est, co_lb, co_ub,
         within_est, within_lb, within_ub, within_p, within_p_model) %>%
  pivot_wider(names_from = sc, values_from = c(records, studies, pure_est, pure_lb, pure_ub, co_est, co_lb, co_ub,
                                               within_est, within_lb, within_ub, within_p, within_p_model)) %>%
  arrange(analysis)
for (cc in c("records_l", "studies_l", "pure_est_l", "pure_lb_l", "pure_ub_l", "co_est_l", "co_lb_l", "co_ub_l", "within_est_l",
             "within_lb_l", "within_ub_l", "within_p_l", "within_p_model_l", "records_g", "studies_g", "pure_est_g", "pure_lb_g",
             "pure_ub_g", "co_est_g", "co_lb_g", "co_ub_g", "within_est_g", "within_lb_g", "within_ub_g", "within_p_g", "within_p_model_g"))
  if (!cc %in% names(sec_wide)) sec_wide[[cc]] <- NA_real_
sens_sec_fmt <- sec_wide %>%
  transmute(Analysis = as.character(analysis),
            N_l = if_else(is.na(pure_est_l), NA_character_, paste0(records_l, " / ", studies_l)),
            pu_l = cell(pure_est_l, pure_lb_l, pure_ub_l, TRUE), co_l = cell(co_est_l, co_lb_l, co_ub_l, TRUE),
            wi_l = cell(within_est_l, within_lb_l, within_ub_l, TRUE),
            pw_l = if_else(is.na(within_est_l), NA_character_, paste0(fmt_pr(within_p_l), " / ", fmt_p(within_p_model_l))),
            N_g = if_else(is.na(pure_est_g), NA_character_, paste0(records_g, " / ", studies_g)),
            pu_g = cell(pure_est_g, pure_lb_g, pure_ub_g, FALSE), co_g = cell(co_est_g, co_lb_g, co_ub_g, FALSE),
            wi_g = cell(within_est_g, within_lb_g, within_ub_g, FALSE),
            pw_g = if_else(is.na(within_est_g), NA_character_, paste0(fmt_pr(within_p_g), " / ", fmt_p(within_p_model_g))))
sens_sec_fmt %>%
  kbl(col.names = c("Analysis", "Records / studies", "Pure plastic", "Co-contaminant", "Within-study contrast", "p (robust / model)",
                    "Records / studies", "Pure plastic", "Co-contaminant", "Within-study contrast", "p (robust / model)"),
      align = c("l", rep("c", 10)), booktabs = TRUE, na = "\u2014",
      caption = paste("Sensitivity of the secondary analysis: % change in hatching for each contaminant type and the within-study contrast",
                      "(relative difference in the hatching ratio, co-contaminant minus pure plastic), with a separate residual variance per type.",
                      "Hedges' g (secondary scale) alongside; outlier sets are defined on each scale's own model.")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), font_size = 11) %>%
  add_header_above(c(" " = 1, "lnRR (primary): % change" = 5, "Hedges' g (secondary)" = 5)) %>%
  row_spec(0, bold = TRUE)

# 16.5 Checks: run top to bottom and compare against the paragraphs ----------
cat("\n== A. composition ==\n"); print(comp); print(study_types)
cat("\n== B. Table 2 ==\n"); print(tab2_fmt, width = Inf)
cat("\n== C. contrast on the lnRR scale as a relative difference in hatching ratio (%) ==\n")
tab_contrast_rel <- tab2_contrast %>% filter(scale == "lnRR") %>%
  mutate(rel = pct(est), rel_lb = pct(ci_lb), rel_ub = pct(ci_ub)) %>%
  select(model, est, p_robust, p_model, df, rel, rel_lb, rel_ub)
print(tab_contrast_rel, width = Inf)

cat("\n== D. structures, secondary ==\n"); print(fmt_struct(tab2s), width = Inf)

cat("\n== E. weight share: how much of the pooled weight do co-contaminant records carry? ==\n")
w_share <- \(d, m) d %>% mutate(w = 1 / (vi + extra_var(m) + obs_var_vec(m, d))) %>% group_by(contam) %>%
  summarise(records = n(), med_vi = median(vi), w = sum(w), .groups = "drop") %>%
  mutate(weight_share = round(w / sum(w), 3)) %>% select(-w)
tab_wshare <- bind_rows(w_share(d2_l, fit_cred(d2_l, "simple")) %>% mutate(scale = "lnRR"),
                        w_share(d2_g, fit_cred(d2_g, "simple")) %>% mutate(scale = "Hedges' g"))
print(tab_wshare)

cat("\n== F. sensitivity: primary, and secondary by contaminant type ==\n"); print(sens_fmt, n = Inf, width = Inf)
print(sens_sec_fmt, n = Inf, width = Inf)

cat("\n== G. flagged rows in the secondary (each scale's own flags) ==\n")
z_sec_l <- rstudent(fit_cred(d2_l, "simple"))$z
z_sec_g <- rstudent(fit_cred(d2_g, "simple"))$z
tab_flagged_sec <- bind_rows(
  d2_l %>% mutate(rstud_z = z_sec_l, scale = "lnRR",
                  rule = if_else(abs(z_sec_l) > 3, "residual (lnRR model)", NA_character_)),
  d2_g %>% mutate(rstud_z = z_sec_g, scale = "Hedges' g",
                  rule = case_when(abs(z_sec_g) > 3 & vi > 5 ~ "both", abs(z_sec_g) > 3 ~ "residual (g model)",
                                   vi > 5 ~ "variance (vi > 5)"))) %>%
  filter(!is.na(rule)) %>% select(scale, study_id, reference, contam, conc_mg_L, yi, vi, rstud_z, rule) %>%
  arrange(factor(scale, levels = scale_levels), desc(abs(yi))) %>% tibble::as_tibble()
print(tab_flagged_sec, n = Inf, width = Inf)

cat("\n== H. leave-one-study-out, both scales ==\n")
loo2 <- purrr::map_dfr(list(d2_l, d2_g), \(d) purrr::map_dfr(unique(d$study_id), \(s)
                                                             red(filter(d, study_id != s), "drop") %>% mutate(dropped = s))) %>%
  left_join(distinct(d2_g, study_id, reference), by = c("dropped" = "study_id"))
# [V5] the contrast is the main secondary quantity, so leave-one-study-out is repeated on the within-study contrast
loo2_contrast <- purrr::map_dfr(list(d2_l, d2_g), \(d) purrr::map_dfr(unique(d$study_id), \(s) {
  r <- purrr::possibly(contam_contrast, otherwise = NULL)(filter(d, study_id != s))
  if (is.null(r)) tibble(scale = sc_lab(d), dropped = s, note = "not estimable")
  else filter(r, grepl("^Within", contrast)) %>% mutate(dropped = s) })) %>%
  left_join(distinct(d2_g, study_id, reference), by = c("dropped" = "study_id"))
tab_loo2_contrast_summary <- loo2_contrast %>% filter(!is.na(est)) %>% group_by(scale) %>%
  summarise(min_est = min(est), max_est = max(est), n_sig = sum(p_robust < 0.05, na.rm = TRUE),
            n_reportable = sum(reportable), n = n(), .groups = "drop") %>%
  mutate(min_pct = if_else(scale == "lnRR", pct(min_est), NA_real_), max_pct = if_else(scale == "lnRR", pct(max_est), NA_real_)) %>%
  arrange(factor(scale, levels = scale_levels))
print(tab_loo2_contrast_summary)
loo2_contrast %>% filter(!is.na(est)) %>% group_by(scale) %>% arrange(desc(est), .by_group = TRUE) %>% slice_head(n = 4) %>%
  select(scale, dropped, reference, est, ci_lb, ci_ub, p_robust, p_model) %>% ungroup() %>% print(width = Inf)

tab_loo2_summary <- loo2 %>% group_by(scale) %>%
  summarise(min_est = min(est), max_est = max(est), n_sig = sum(p_robust < 0.05, na.rm = TRUE),
            n_reportable = sum(reportable), n = n(), .groups = "drop") %>%
  mutate(min_pct = if_else(scale == "lnRR", pct(min_est), NA_real_), max_pct = if_else(scale == "lnRR", pct(max_est), NA_real_)) %>%
  arrange(factor(scale, levels = scale_levels))
print(tab_loo2_summary)
tab_loo2_top <- loo2 %>% group_by(scale) %>% arrange(desc(est), .by_group = TRUE) %>% slice_head(n = 4) %>%
  select(scale, dropped, reference, est, ci_lb, ci_ub, p_robust, p_model) %>% ungroup() %>%
  arrange(factor(scale, levels = scale_levels))
print(tab_loo2_top, width = Inf)

cat("\n== I. the three largest-effect studies in the secondary ==\n")
tab_big3_sec <- pct_rank(dat_sec) %>% slice(1:5)
print(tab_big3_sec)
tab_big3_sec_desc <- d2_g %>% filter(study_id %in% big_sec) %>% count(study_id, reference, contam, taxa_class, data_source)
print(tab_big3_sec_desc, width = Inf)
tab_big3_sec_comp <- d2_g %>% mutate(group = if_else(study_id %in% big_sec, "large-effect", "other")) %>% count(group, contam)
print(tab_big3_sec_comp)
tab_big3_sec_red <- bind_rows(sec_row(d2_l, "all"), sec_row(filter(d2_l, !study_id %in% big_sec), "excl. big 3"),
                              sec_row(d2_g, "all"), sec_row(filter(d2_g, !study_id %in% big_sec), "excl. big 3")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(tab_big3_sec_red, width = Inf)

set.seed(1)
draws2 <- replicate(300, sample(unique(d2_g$study_id), 3), simplify = FALSE)
est3  <- \(d, drop) as.numeric(coef(fit_cred(filter(d, !study_id %in% drop), "simple")))
safe3 <- purrr::possibly(est3, otherwise = NA_real_)
rl <- purrr::map_dbl(draws2, \(x) safe3(d2_l, x)); rg <- purrr::map_dbl(draws2, \(x) safe3(d2_g, x))
ol <- est3(d2_l, big_sec); og <- est3(d2_g, big_sec)
tab_draws_sec <- tibble(scale = c("lnRR", "Hedges' g"), observed = c(ol, og),
                        rand_q05 = c(quantile(rl, .05, na.rm = TRUE), quantile(rg, .05, na.rm = TRUE)),
                        rand_median = c(median(rl, na.rm = TRUE), median(rg, na.rm = TRUE)),
                        rand_q95 = c(quantile(rl, .95, na.rm = TRUE), quantile(rg, .95, na.rm = TRUE)),
                        prop_as_close_to_zero = c(mean(rl >= ol, na.rm = TRUE), mean(rg >= og, na.rm = TRUE))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(tab_draws_sec)

cat("\n== J. typical study: change in hatching ==\n")
sp2 <- pct_rank(dat_sec)
tab_typical_sec <- sp2 %>% summarise(n = n(), median = median(med_pct), q25 = quantile(med_pct, .25), q75 = quantile(med_pct, .75),
                                     n_below_m10 = sum(med_pct < -10), n_no_reduction = sum(med_pct >= 0))
print(tab_typical_sec)
tab_typical_contam <- d2_g %>% mutate(pc = 100 * (treat_mhs / control_mhs - 1)) %>% group_by(contam) %>%
  summarise(records = n(), studies = n_distinct(study_id), median = median(pc), mean = mean(pc), .groups = "drop")
print(tab_typical_contam)

cat("\n== K. sample-size covariances vs fixed rho (secondary) ==\n"); print(tab_rho_sec)

# 17. Study 38 influence and SD-imputation scenarios -------------------------
s38 <- 38

## a) What is Study 38, and how does it differ from the other studies?
tab_s38_profile <- dat_cred %>% filter(study_id == s38) %>%
  summarise(reference = first(reference), records = n(), control_groups = n_distinct(control_id),
            classes = paste(unique(taxa_class), collapse = "/"), species = paste(unique(species), collapse = ", "),
            polymers = paste(unique(polymer), collapse = "/"), contam = paste(unique(as.character(contam)), collapse = "/"),
            routes = paste(unique(route), collapse = "/"), lifestages = paste(unique(lifestage), collapse = "/"),
            units = paste(unique(concentration_unit), collapse = "/"),
            conc_mg_L = if (all(is.na(conc_mg_L))) NA_character_ else paste(range(conc_mg_L, na.rm = TRUE), collapse = " to "),
            source = paste(unique(data_source), collapse = "/"),
            med_treat_reps = median(treat_reps), med_ctrl_reps = median(control_reps),
            share_sd_imputed = mean(sd_imputed), g_min = min(yi), g_max = max(yi),
            prop_negative = mean(yi < 0), med_vi = median(vi))
print(tab_s38_profile, width = Inf)

tab_s38_vs_others <- dat_cred %>%
  group_by(group = if_else(study_id == s38, "Study 38", "other studies")) %>%
  summarise(records = n(), studies = n_distinct(study_id), med_vi = median(vi),
            med_treat_reps = median(treat_reps), med_g = median(yi),
            share_sd_imputed = mean(sd_imputed), .groups = "drop")
print(tab_s38_vs_others)

## b) Share of the records vs share of the statistical weight (intercept-only model, exact sampling covariances)
study_weight <- function(m, d) {
  V   <- build_V(d)
  ids <- d$study_id
  M   <- V + extra_var(m) * outer(ids, ids, "==") + diag(obs_var_vec(m, d), nrow = nrow(d))   # marginal covariance
  u   <- unique(ids)
  w   <- vapply(u, \(s) { i <- which(ids == s); sum(solve(M[i, i, drop = FALSE])) }, numeric(1))
  tibble(study_id = u, records = as.integer(table(ids)[as.character(u)]), weight_share = w / sum(w)) %>%
    mutate(record_share = records / sum(records))
}
tab_s38_weight <- bind_rows(
  study_weight(fit_cred(dat_pure_l, "simple"), dat_pure_l) %>% mutate(set = "Primary, lnRR"),
  study_weight(fit_cred(dat_pure, "simple"), dat_pure)     %>% mutate(set = "Primary, g"),
  study_weight(fit_cred(dat_sec_l, "simple"), dat_sec_l)   %>% mutate(set = "Secondary, lnRR"),
  study_weight(fit_cred(dat_sec, "simple"), dat_sec)       %>% mutate(set = "Secondary, g"))
tab_s38_weight %>% filter(study_id == s38) %>% print()
tab_s38_weight %>% group_by(set) %>% arrange(desc(weight_share), .by_group = TRUE) %>% slice_head(n = 4) %>% print()
w38_of <- \(d) { w <- study_weight(fit_cred(d, "simple"), d); x <- w$weight_share[w$study_id == s38]
if (length(x) == 0) NA_real_ else x }

## c) Overall effect, variance components and dose slope with and without Study 38 (lnRR first)
sets38 <- list("Primary" = dat_pure, "Primary excl. 38" = filter(dat_pure, study_id != s38),
               "Secondary" = dat_sec, "Secondary excl. 38" = filter(dat_sec, study_id != s38))
tab_s38_cmp <- purrr::imap_dfr(sets38, \(d, nm) {
  one <- \(dd) { m <- fit_cred(dd, "simple")
  ov_dose(dd) %>% mutate(s2_study = m$sigma2[1], s2_obs = mean(obs_var_vec(m, dd))) }
  bind_rows(one(make_lnrr(d)), one(d)) %>% mutate(analysis = nm) }) %>%
  mutate(pct = if_else(scale == "lnRR", 100 * (exp(est) - 1), NA_real_))
tab_s38_cmp %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(width = Inf)

## d) Downweighting by record count: keep a typical number of Study 38's records, 200 random draws
n38   <- sum(dat_pure$study_id == s38)
k_typ <- round(median(table(dat_pure$study_id[dat_pure$study_id != s38])))
tab_s38_downsample <- if (n38 > k_typ) {
  down <- \(d, seed) { set.seed(seed)
    ids  <- d$es_id[d$study_id == s38]
    keep <- ids[sample.int(length(ids), k_typ)]
    filter(d, study_id != s38 | es_id %in% keep) }
  est1 <- \(d) as.numeric(coef(fit_cred(d, "simple")))
  ds <- purrr::map_dfr(1:200, \(i) tibble(draw = i, lnRR = est1(down(dat_pure_l, i)), g = est1(down(dat_pure, i))))
  tibble(scale = c("lnRR", "Hedges' g"), records_38 = n38, records_kept = k_typ,
         full    = c(est1(dat_pure_l), est1(dat_pure)),
         excl_38 = c(est1(filter(dat_pure_l, study_id != s38)), est1(filter(dat_pure, study_id != s38))),
         down_q05 = c(quantile(ds$lnRR, .05), quantile(ds$g, .05)),
         down_median = c(median(ds$lnRR), median(ds$g)),
         down_q95 = c(quantile(ds$lnRR, .95), quantile(ds$g, .95)))
} else tibble(note = "Study 38 has no more records in the primary set than a typical study")
tab_s38_downsample %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(width = Inf)

## e) Moderator tests with and without Study 38 (robust p where df >= df_min, model-based p alongside)
tab_s38_moderators <- bind_rows(
  run_grid(purrr::map(sets38[1:2], make_lnrr), specs, "lnRR")$omnibus,
  run_grid(sets38[1:2], specs, "Hedges' g")$omnibus) %>%
  filter(moderator != "contam") %>%
  select(scale, analysis, moderator, k_studies, p_robust, p_model, reliable)
tab_s38_moderators %>%
  tidyr::pivot_wider(names_from = analysis, values_from = c(k_studies, p_robust, p_model, reliable)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

## f) Sample-size covariances vs fixed rho, with and without Study 38
tab_s38_rho <- purrr::map_dfr(names(sets38)[1:2], \(set)
                              bind_rows(rho_one(make_lnrr(sets38[[set]]), set), rho_one(sets38[[set]], set)))
tab_s38_rho %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(width = Inf)

## g) SD-imputation scenarios: scale the class CV, or swap in the old rule, or drop the imputed studies
cv_scen <- function(d, mult, incl_bin = FALSE) {      # incl_bin = TRUE also applies the class CV to binomial_ids (the earlier rule)
  d %>% select(-yi, -vi) %>%
    left_join(select(class_cv, taxa_class, cv_ctrl, cv_trt), by = "taxa_class") %>%
    mutate(use = sd_imputed & (incl_bin | !binom),
           control_sd = if_else(use, pmax(pmin(mult * cv_ctrl * control_mhs, sd_cap(control_mhs)), sd_floor), control_sd),
           treat_sd   = if_else(use, pmax(pmin(mult * cv_trt  * treat_mhs,  sd_cap(treat_mhs)),  sd_floor), treat_sd)) %>%
    select(-cv_ctrl, -cv_trt, -use) %>%
    escalc(measure = "SMD", m1i = treat_mhs, sd1i = treat_sd, n1i = treat_reps,
           m2i = control_mhs, sd2i = control_sd, n2i = control_reps, data = .) %>%
    tibble::as_tibble()
}
stopifnot(isTRUE(all.equal(as.numeric(cv_scen(dat_pure, 1)$yi), as.numeric(dat_pure$yi))))   # x1 must reproduce the primary

scen <- c("CV x0.5" = 0.5, "CV x1 (primary)" = 1, "CV x1.5" = 1.5, "CV x2" = 2)
scen_sets <- c(purrr::map(scen, \(m) cv_scen(dat_pure, m)),
               list("Study 38 SD from class CV (earlier rule)"          = cv_scen(dat_pure, 1, incl_bin = TRUE),
                    "Constant class-mean SD (old rule, all imputed studies)" = swap_sd(dat_pure),
                    "Excl. imputed SDs"                                  = filter(dat_pure, !sd_imputed),
                    "Excl. Study 38"                                     = filter(dat_pure, study_id != s38)))

tab_sd_scenarios <- purrr::imap_dfr(scen_sets, \(d, nm) {
  one <- \(dd) ov_dose(dd) %>% mutate(w38 = w38_of(dd))
  bind_rows(one(make_lnrr(d)), one(d)) %>% mutate(analysis = nm) }) %>%
  mutate(pct = if_else(scale == "lnRR", 100 * (exp(est) - 1), NA_real_))
tab_sd_scenarios %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# 18. Primary with vs without the imputed-SD studies --------------------------
# Needs sections 1-14 only. with_imputed = current primary and secondary; excl_imputed = the imputed-SD studies
# (3, 27, 38 and 39) removed from both, i.e. imputed SDs used as a sensitivity set only.
imp_sets <- list(
  with_imputed = list(Primary = dat_pure,                      Secondary = dat_sec),
  excl_imputed = list(Primary = filter(dat_pure, !sd_imputed), Secondary = filter(dat_sec, !sd_imputed)))

tab_imp_n <- purrr::imap_dfr(imp_sets, \(x, nm)
                             purrr::imap_dfr(x, \(d, ds) tibble(version = nm, dataset = ds, records = nrow(d), studies = n_distinct(d$study_id))))
print(tab_imp_n)

# overall effect and variance components, both datasets and scales (lnRR first)
tab_imp_overall <- purrr::imap_dfr(imp_sets, \(x, nm) {
  purrr::imap_dfr(x, \(d, ds) bind_rows(red(make_lnrr(d), ds), red(d, ds))) %>% mutate(version = nm)
}) %>%
  rename(dataset = analysis) %>%
  mutate(pct    = if_else(scale == "lnRR", 100 * (exp(est) - 1),   NA_real_),
         pct_lb = if_else(scale == "lnRR", 100 * (exp(ci_lb) - 1), NA_real_),
         pct_ub = if_else(scale == "lnRR", 100 * (exp(ci_ub) - 1), NA_real_)) %>%
  select(version, dataset, scale, everything())
tab_imp_overall %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# dose-response (within-study slope), primary and secondary, both scales
tab_imp_dose <- purrr::imap_dfr(imp_sets, \(x, nm) {
  purrr::imap_dfr(x, \(d, ds) bind_rows(ov_dose(make_lnrr(d)), ov_dose(d)) %>% mutate(dataset = ds)) %>%
    mutate(version = nm)
}) %>% select(version, dataset, scale, k_es, k_studies, k_dose, within, w_lb, w_ub, p_within, p_within_model, df_within)
tab_imp_dose %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# moderators: robust p (df-gated) and model-based p side by side
tab_imp_mod <- purrr::imap_dfr(imp_sets, \(x, nm) {
  bind_rows(run_grid(purrr::map(x, make_lnrr), specs, "lnRR")$omnibus,
            run_grid(x, specs, "Hedges' g")$omnibus) %>% mutate(version = nm)
}) %>%
  filter(!(moderator == "contam" & analysis == "Primary")) %>%
  select(version, analysis, scale, moderator, k_studies, p_robust, p_model) %>%
  tidyr::pivot_wider(names_from = version, values_from = c(k_studies, p_robust, p_model))
tab_imp_mod %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# contaminant type (secondary): pooled, within- and between-study contrasts
tab_imp_contam <- purrr::imap_dfr(imp_sets, \(x, nm)
                                  bind_rows(contam_contrast(make_lnrr(x$Secondary)), contam_contrast(x$Secondary)) %>% mutate(version = nm)) %>%
  select(version, scale, contrast, records, studies, est, ci_lb, ci_ub, p_robust, p_model, df)
tab_imp_contam %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# contaminant-type contrast under every secondary sensitivity set (is it carried by a few extreme records?), under the secondary model
# (separate residual variance by type) and, for comparison, the earlier common-variance model
cs_row <- \(d, nm, resid = "auto") { r <- purrr::possibly(contam_contrast, otherwise = NULL)(d, resid = resid)
if (is.null(r)) tibble(analysis = nm, note = "not estimable") else mutate(r, analysis = nm) }
tab_contam_sens <- bind_rows(purrr::imap_dfr(sens_sec$lnRR, \(d, nm) cs_row(d, nm)),
                             purrr::imap_dfr(sens_sec$g,    \(d, nm) cs_row(d, nm)))
tab_contam_sens %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)
tab_contam_sens_common <- bind_rows(purrr::imap_dfr(sens_sec$lnRR, \(d, nm) cs_row(d, nm, "common")),
                                    purrr::imap_dfr(sens_sec$g,    \(d, nm) cs_row(d, nm, "common")))

# 18b. Contaminant contrast: robustness to the lnRR floor, the random-effects structure and the residual variance -----
# (a) lnRR floor for complete failure (0.1, 0.5, 1 % hatching) x structure (simple, species) x residual variance (separate by type, common);
#     the lnRR zero-hatching arm is Bellas & Gil 2020
contrast_grid <- tidyr::crossing(structure = c("simple", "species"), residual = c("by_contam", "common"), floor_pct = c(0.1, 0.5, 1))
safe_contrast <- function(d, s, r = "auto") {
  x <- purrr::possibly(contam_contrast, otherwise = NULL)(d, s, r)
  if (is.null(x)) tibble(scale = sc_lab(d), note = "not estimable") else x
}
tab_contrast_floor <- purrr::pmap_dfr(contrast_grid, \(structure, residual, floor_pct)
                                      safe_contrast(make_lnrr(dat_sec, floor_pct), structure, residual) %>%
                                        mutate(structure = structure, residual = residual, floor_pct = floor_pct, .before = 1))
tab_contrast_floor %>% mutate(pct = if_else(scale == "lnRR", 100 * (exp(est) - 1), NA_real_),
                              pct_lb = if_else(scale == "lnRR", 100 * (exp(ci_lb) - 1), NA_real_),
                              pct_ub = if_else(scale == "lnRR", 100 * (exp(ci_ub) - 1), NA_real_)) %>%
  select(structure, residual, floor_pct, contrast, records, studies, est, ci_lb, ci_ub, df, p_robust, p_model, infer, pct, pct_lb, pct_ub) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# (b) structure x residual variance (common vs separate for pure plastic and co-contaminant records), both scales
resid_grid <- tidyr::crossing(structure = c("simple", "species"), resid = c("common", "by_contam"))
tab_contrast_resid <- purrr::pmap_dfr(resid_grid, \(structure, resid)
                                      bind_rows(safe_contrast(dat_sec_l, structure, resid), safe_contrast(dat_sec, structure, resid)) %>%
                                        mutate(structure = structure, residual = resid, .before = 1))
tab_contrast_resid %>% mutate(pct = if_else(scale == "lnRR", 100 * (exp(est) - 1), NA_real_),
                              pct_lb = if_else(scale == "lnRR", 100 * (exp(ci_lb) - 1), NA_real_),
                              pct_ub = if_else(scale == "lnRR", 100 * (exp(ci_ub) - 1), NA_real_)) %>%
  select(scale, structure, residual, contrast, est, ci_lb, ci_ub, df, p_robust, p_model, infer, pct, pct_lb, pct_ub) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# (c) overall secondary effect and variance components under each structure x residual variance
het_overall <- function(d, struct, resid) {
  m <- fit_cred(d, struct, resid = resid); rb <- rob_safe(m, d$study_id); it <- inf_tbl(m, rb, 1)
  s <- m$sigma2
  tibble(scale = sc_lab(d), structure = struct, residual = resid, k_es = m$k, est = it$est, ci_lb = it$ci_lb, ci_ub = it$ci_ub,
         df = it$df, p_robust = it$p_robust, p_model = it$p_model, infer = attr(rb, "infer"),
         s2_study   = s[1],
         s2_species = if (struct == "species") s[if (resid == "by_contam") 2 else 3] else NA_real_,
         s2_obs_pure = if (resid == "by_contam") m$tau2[1] else s[2],
         s2_obs_co   = if (resid == "by_contam") m$tau2[2] else s[2],
         AIC = as.numeric(AIC(m)))
}
tab_hetero_overall <- purrr::pmap_dfr(resid_grid, \(structure, resid)
                                      bind_rows(het_overall(dat_sec_l, structure, resid), het_overall(dat_sec, structure, resid)))
tab_hetero_overall %>% mutate(pct = if_else(scale == "lnRR", 100 * (exp(est) - 1), NA_real_)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# (d) polymer: PE and Other each compared with PS (own df per contrast; the 2-df omnibus is NE), and the in-study comparison
poly_contrasts <- function(d) {
  d <- d %>% filter(!is.na(polymer_test)) %>% mutate(polymer_test = factor(as.character(polymer_test), levels = c("PS", "PE", "Other")))
  m <- fit_cred(d, "simple", mods = ~ polymer_test); rb <- rob_safe(m, d$study_id)
  inf_tbl(m, rb, 2:3) %>%
    mutate(contrast = paste(sub("polymer_test", "", term), "minus PS"), scale = sc_lab(d), k_es = m$k,
           k_studies = n_distinct(d$study_id), infer = attr(rb, "infer"), p_holm = p.adjust(p_robust, "holm")) %>%
    select(-term)
}
tab_polymer_contrasts <- bind_rows(
  purrr::map_dfr(list(Primary = dat_pure_l, Secondary = dat_sec_l, Primary = dat_pure, Secondary = dat_sec), poly_contrasts, .id = "dataset"))
tab_polymer_contrasts %>% mutate(pct = if_else(scale == "lnRR", 100 * (exp(est) - 1), NA_real_)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)

# studies whose arms differ in polymer: mean effect per polymer within each study (descriptive; df cannot support a test)
tab_polymer_within <- bind_rows(Primary = dat_pure_l, Secondary = dat_sec_l, .id = "dataset") %>%
  filter(!is.na(polymer_test)) %>% group_by(dataset, study_id) %>% filter(n_distinct(polymer_test) >= 2) %>%
  group_by(dataset, study_id, reference, polymer_test) %>%
  summarise(k = n(), mean_lnRR = mean(yi), pct = 100 * (exp(mean(yi)) - 1), .groups = "drop")
print(tab_polymer_within, n = Inf)

# [V5] polymer, BETWEEN studies only: keep studies that tested a single polymer, so there is no within-study leverage (Study 38 etc.).
# Of the three tested levels only PS and PE have >= 3 such studies (Other is mostly carried by the mixed-polymer studies).
single_poly <- \(d) d %>% filter(!is.na(polymer_test)) %>% group_by(study_id) %>%
  filter(n_distinct(polymer_test) == 1) %>% ungroup()
tab_polymer_between_n <- purrr::map_dfr(list(Primary = dat_pure_l, Secondary = dat_sec_l), \(d) single_poly(d) %>%
                                          group_by(polymer_test) %>% summarise(k_es = n(), k_studies = n_distinct(study_id), .groups = "drop"), .id = "dataset")
print(tab_polymer_between_n)

poly_between <- function(d) {
  dd <- single_poly(d) %>% filter(polymer_test %in% c("PS", "PE")) %>%
    mutate(polymer_test = factor(as.character(polymer_test), levels = c("PS", "PE")))
  m <- fit_cred(dd, "simple", mods = ~ polymer_test); rb <- rob_safe(m, dd$study_id)
  inf_tbl(m, rb, 2) %>%
    mutate(contrast = "PE minus PS (single-polymer studies)", scale = sc_lab(dd), resid = resid_of(m), k_es = m$k,
           k_studies = n_distinct(dd$study_id), n_PS = n_distinct(dd$study_id[dd$polymer_test == "PS"]),
           n_PE = n_distinct(dd$study_id[dd$polymer_test == "PE"]), infer = attr(rb, "infer")) %>% select(-term)
}
tab_polymer_between <- purrr::map_dfr(list(Primary = dat_pure_l, Secondary = dat_sec_l, Primary = dat_pure, Secondary = dat_sec),
                                      poly_between, .id = "dataset")
tab_polymer_between %>% mutate(pct = if_else(scale == "lnRR", 100 * (exp(est) - 1), NA_real_)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = Inf, width = Inf)


# Final export: re-collect everything, now including sections 15-18b -----------
results <- collect_results(print_all = FALSE)

#-tidying tables
#flow -----------------------------------------
flow <- read_csv("results_39studies_v5/flow.csv", show_col_types = FALSE)

flow_tab <- flow %>%
  mutate(
    # Show changes relative to the previous step (blank for the first row
    # and for the final secondary-only row, which is a subset, not a step)
    `Records removed` = c(NA, -diff(records[1:3]), NA, NA),
    `Studies removed` = c(NA, -diff(studies[1:3]), NA, NA)
  ) %>%
  transmute(
    Step = step,
    Records = records,
    Studies = studies,
    `Records removed`,
    `Studies removed`
  )

kbl(
  flow_tab,
  format    = "html",            # use "latex" for PDF/Word-via-LaTeX output
  caption   = "Table S1. Flow of records and studies through the screening and analysis stages.",
  align     = c("l", "r", "r", "r", "r"),
  col.names = c("Step", "Records", "Studies", "Records removed", "Studies removed"),
  booktabs  = TRUE,
  na        = "–"
) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position   = "left"
  ) %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(4, bold = TRUE) %>%                 # primary analysis dataset
  pack_rows("Screening", 1, 3) %>%
  pack_rows("Analysis datasets", 4, 5) %>%
  footnote(
    general = paste(
      "Primary analysis = pure-plastic records only.",
      "Co-contaminant records were used only in the secondary analysis",
      "(primary + co-contaminant = 134 records, 34 studies).",
      "R3 = CRED reliability criterion 3."
    ),
    general_title = "Note: ",
    footnote_as_chunk = TRUE
  )

#sensitivity table -----------------------------------------------------
res <- "results_39studies_v5"            # folder written by CRED_39studies_v5.R
rd  <- function(f) read_csv(file.path(res, f), show_col_types = FALSE)

pct   <- function(x) 100 * (exp(x) - 1)
fp    <- function(p) ifelse(is.na(p), "–", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
f_pct <- function(e, l, u) sprintf("%.1f (%.1f, %.1f)", pct(e), pct(l), pct(u))
f_g   <- function(e, l, u) sprintf("%.2f (%.2f, %.2f)", e, l, u)

# ---- pull the pieces ---------------------------------------------------------
sens <- rd("sensitivity.csv")    %>% filter(dataset == "Primary")
rho  <- rd("rho.csv")            %>% filter(dataset == "Primary", method == "fixed rho")
rho_l <- filter(rho, scale == "lnRR"); rho_g <- filter(rho, scale == "Hedges' g")
sdsc <- rd("sd_scenarios.csv")
s38  <- rd("s38_comparison.csv") %>% filter(analysis == "Primary excl. 38")
loo  <- rd("loo_overall.csv")    %>% filter(analysis == "drop")

# one-row-per-analysis from a table with scale / est / ci / p_robust / k columns
pair <- function(d, label, group) {
  dl <- d %>% filter(scale == "lnRR"); dg <- d %>% filter(scale == "Hedges' g")
  tibble(Group = group, Analysis = label,
         Records = dl$k_es, Studies = dl$k_studies,
         lnRR = f_pct(dl$est, dl$ci_lb, dl$ci_ub), p_l = fp(dl$p_robust),
         g    = f_g(dg$est, dg$ci_lb, dg$ci_ub),   p_g = fp(dg$p_robust))
}
S <- function(nm, scl = NULL) sens %>% filter(analysis == nm)

rows <- bind_rows(
  pair(S("Base dataset (CRED)"), "Primary analysis", "Primary"),
  
  # covariance assumption (fixed-rho rows report model-based CIs only; no robust p)
  tibble(Group = "Covariance assumption",
         Analysis = paste0("Fixed correlation ρ = ", rho_l$rho),
         Records = 102L, Studies = 32L,
         lnRR = f_pct(rho_l$est, rho_l$ci_lb, rho_l$ci_ub), p_l = "–",
         g    = f_g(rho_g$est, rho_g$ci_lb, rho_g$ci_ub),   p_g = "–"),
  
  pair(S("R3 retained"),            "Retain CRED R3 records",                    "Data inclusion"),
  pair(S("Excl. PlotDigitiser"),    "Exclude digitised data",                    "Data inclusion"),
  
  pair(S("Excl. imputed SDs"),      "Exclude studies with imputed SDs",          "SD handling"),
  pair(S("Constant class-mean SD (old rule, all imputed studies)"),
       "Constant class-mean SD instead of CV rule", "SD handling"),
  pair(sdsc %>% filter(analysis == "CV x0.5"), "Imputed CV × 0.5",               "SD handling"),
  pair(sdsc %>% filter(analysis == "CV x2"),   "Imputed CV × 2",                 "SD handling"),
  pair(s38 %>% mutate(analysis = "x"),         "Exclude Prata et al. (2022)",     "SD handling"),
  pair(sdsc %>% filter(analysis == "Study 38 SD from class CV (earlier rule)"),
       "Prata et al. (2022): class-CV SD instead of binomial SD", "SD handling"),
  
  # outliers: lnRR and g use their own rules (see footnote)
  { ol <- S("Excl. |z| > 3 (lnRR model)") %>% filter(scale == "lnRR")
  og <- S("Excl. |z| > 3 (g model)")    %>% filter(scale == "Hedges' g")
  tibble(Group = "Outliers and influence", Analysis = "Exclude |studentised residual| > 3",
         Records = ol$k_es, Studies = ol$k_studies,
         lnRR = f_pct(ol$est, ol$ci_lb, ol$ci_ub), p_l = fp(ol$p_robust),
         g    = f_g(og$est, og$ci_lb, og$ci_ub),   p_g = fp(og$p_robust)) },
  tibble(Group = "Outliers and influence", Analysis = "Leave-one-study-out (range)",
         Records = NA_integer_, Studies = 31,
         lnRR = {l <- loo %>% filter(scale == "lnRR");
         sprintf("%.1f to %.1f", pct(min(l$est)), pct(max(l$est)))},
         p_l  = {l <- loo %>% filter(scale == "lnRR"); paste0("≤ ", fp(max(l$p_robust, na.rm = TRUE)))},
         g    = {l <- loo %>% filter(scale == "Hedges' g");
         sprintf("%.2f to %.2f", min(l$est), max(l$est))},
         p_g  = {l <- loo %>% filter(scale == "Hedges' g"); paste0("≤ ", fp(max(l$p_robust, na.rm = TRUE)))}),
  
  pair(S("Excl. 3 largest-effect studies (post hoc)"),
       "Exclude the three largest-effect studies (post hoc)", "Outliers and influence")
)

tab <- rows %>% mutate(Records = ifelse(is.na(Records), "–", as.character(Records)))
grp <- rle(tab$Group)

kbl(tab %>% select(-Group), format = "html", booktabs = TRUE, align = c("l", rep("r", 7)),
    caption = "Table 2. Sensitivity of the pure-plastic effect on hatching success to analytical decisions.",
    col.names = c("Analysis", "Records", "Studies", "lnRR, % change (95% CI)", "p", "Hedges' g (95% CI)", "p")) %>%
  add_header_above(c(" " = 3, "Primary scale (lnRR)" = 2, "Secondary scale (g)" = 2)) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE, position = "left") %>%
  row_spec(0, bold = TRUE) %>% row_spec(1, bold = TRUE) %>%
  { x <- .; start <- 1
  for (i in seq_along(grp$lengths)) {
    x <- pack_rows(x, grp$values[i], start, start + grp$lengths[i] - 1); start <- start + grp$lengths[i] }
  x } %>%
  footnote(general = paste(
    "p values are cluster-robust (CR2, Satterthwaite df); all estimates have df ≥ 4.",
    "Fixed-ρ rows show model-based confidence intervals only and are not given a robust p.",
    "Outlier rules were applied separately to each scale (lnRR: 99 records, 31 studies; g: 98 records, 31 studies).",
    "Leave-one-study-out rows show the range of estimates and the largest p across the 32 models."),
    general_title = "Note: ", footnote_as_chunk = TRUE)


# sens table 2 ------------------------------------------------------
res <- "results_39studies_v5"
rd  <- function(f) read_csv(file.path(res, f), show_col_types = FALSE)

pct <- function(x) 100 * (exp(x) - 1)
fp  <- function(p) ifelse(is.na(p), "NE", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
ci_pct <- function(e, l, u) ifelse(is.na(l), "NE", sprintf("%.1f (%.1f, %.1f)", pct(e), pct(l), pct(u)))
ci_g   <- function(e, l, u) ifelse(is.na(l), "NE", sprintf("%.2f (%.2f, %.2f)", e, l, u))

bt  <- rd("secondary_sensitivity_bytype.csv")         # by-type effects + within-study contrast
cf  <- rd("contrast_floor_structure.csv")             # floor and species-structure grid
loo <- rd("loo_secondary_contrast_summary.csv")

# one row from the by-type file; lnRR_name / g_name allow different outlier rules per scale
row_bt <- function(lnRR_name, label, group, g_name = lnRR_name) {
  bl <- bt %>% filter(scale == "lnRR",     analysis == lnRR_name)
  bg <- bt %>% filter(scale == "Hedges' g", analysis == g_name)
  tibble(Group = group, Analysis = label,
         Records = bl$records, Studies = bl$studies,
         Pure = ci_pct(bl$pure_est, bl$pure_lb, bl$pure_ub),
         Co   = ci_pct(bl$co_est,   bl$co_lb,   bl$co_ub),
         Contrast = ci_pct(bl$within_est, bl$within_lb, bl$within_ub), p_l = fp(bl$within_p),
         g  = ci_g(bg$within_est, bg$within_lb, bg$within_ub),         p_g = fp(bg$within_p))
}

# rows that only vary the lnRR floor or the random-effect structure (contrast only)
row_cf <- function(struct, floor, label, group) {
  w <- cf %>% filter(structure == struct, residual == "by_contam", floor_pct == floor,
                     scale == "lnRR", grepl("^Within", contrast))
  tibble(Group = group, Analysis = label, Records = w$records, Studies = w$studies,
         Pure = "–", Co = "–",
         Contrast = ci_pct(w$est, w$ci_lb, w$ci_ub), p_l = fp(w$p_robust),
         g = "–", p_g = "–")
}

lo_l <- loo %>% filter(scale == "lnRR"); lo_g <- loo %>% filter(scale == "Hedges' g")

rows <- bind_rows(
  row_bt("Base dataset (CRED)", "Secondary analysis (primary specification)", "Primary"),
  
  row_bt("R3 retained",                    "Retain CRED R3 records",                     "Data inclusion"),
  row_bt("Excl. PlotDigitiser",            "Exclude digitised data",                     "Data inclusion"),
  
  row_bt("Excl. imputed SDs",              "Exclude studies with imputed SDs",           "SD handling"),
  row_bt("Constant class-mean SD (old rule, all imputed studies)",
         "Constant class-mean SD instead of CV rule",  "SD handling"),
  
  row_cf("simple", 0.1, "lnRR floor of 0.1% (primary: 0.5%)",  "Zero-hatching record"),
  row_cf("simple", 1.0, "lnRR floor of 1%",                    "Zero-hatching record"),
  row_bt("Excl. zero-hatching arms",       "Exclude the zero-hatching record (Bellas & Gil 2020)", "Zero-hatching record"),
  
  row_cf("species", 0.5, "Add species random effect (CR1)",    "Model structure"),
  row_bt("Base dataset, common residual variance (earlier model)",
         "Common residual variance", "Model structure"),
  
  row_bt("Excl. |z| > 3 (lnRR model)",     "Exclude |studentised residual| > 3", "Outliers and influence",
         g_name = "Excl. |z| > 3 (g model)"),
  tibble(Group = "Outliers and influence", Analysis = "Leave-one-study-out (range)",
         Records = NA_integer_, Studies = 33L, Pure = "–", Co = "–",
         Contrast = sprintf("%.1f to %.1f", lo_l$min_pct, lo_l$max_pct),
         p_l = sprintf("%d/%d significant", lo_l$n_sig, lo_l$n),
         g   = sprintf("%.2f to %.2f", lo_g$min_est, lo_g$max_est),
         p_g = sprintf("%d/%d significant", lo_g$n_sig, lo_g$n)),
  row_bt("Excl. 3 largest-effect studies (post hoc)",
         "Exclude the three largest-effect studies (post hoc)", "Outliers and influence")
) %>% mutate(Records = ifelse(is.na(Records), "–", as.character(Records)))

grp <- rle(rows$Group)

kb <- kbl(rows %>% select(-Group), format = "html", booktabs = TRUE, align = c("l", rep("r", 8)),
          caption = "Table S[Y]. Sensitivity of the by-type effects and the within-study contrast (co-contaminant minus pure plastic) to analytical decisions.",
          col.names = c("Analysis", "Records", "Studies", "Pure plastic", "Co-contaminant",
                        "Contrast (95% CI)", "p", "Contrast (95% CI)", "p")) %>%
  add_header_above(c(" " = 3, "By-type effect, lnRR % change (95% CI)" = 2,
                     "Within-study contrast, lnRR %" = 2, "Within-study contrast, g" = 2)) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE, position = "left") %>%
  row_spec(0, bold = TRUE) %>% row_spec(1, bold = TRUE)

start <- 1
for (i in seq_along(grp$lengths)) {
  kb <- pack_rows(kb, grp$values[i], start, start + grp$lengths[i] - 1); start <- start + grp$lengths[i]
}

kb %>% footnote(general = paste(
  "Estimates are from the secondary joint model with separate residual variances by contaminant type unless stated.",
  "The contrast is the within-study difference in the 14 studies testing both types (14–15 depending on the dataset).",
  "p values are cluster-robust (CR2, Satterthwaite df); NE = not estimable (df < 4); the species model uses CR1, which is anti-conservative.",
  "Floor, species and common-variance rows affect lnRR only; g is shown only where it changes.",
  "Outlier rules were applied separately on each scale (lnRR: 129 records; g: 124 records).",
  "Leave-one-study-out gives the range of contrasts across 34 refits (33 studies each)."),
  general_title = "Note: ", footnote_as_chunk = TRUE)


#moderator table ----------------------------------------
res <- "results_39studies_v5"
rd  <- function(f) read_csv(file.path(res, f), show_col_types = FALSE)

fp  <- function(p) ifelse(is.na(p), "NE", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fpc <- function(e, l, u) ifelse(is.na(l), sprintf("%.1f (NE)", e), sprintf("%.1f (%.1f, %.1f)", e, l, u))

lev <- rd("moderator_levels.csv") %>%
  filter(scale == "lnRR", struct == "simple", moderator != "contam",
         analysis %in% c("Primary: pure plastic, CRED", "Secondary: co-contaminants included")) %>%
  mutate(set = ifelse(grepl("^Primary", analysis), "P", "S"))
omn <- rd("moderators_omnibus.csv") %>%
  filter(scale == "lnRR", struct == "simple", moderator != "contam",
         analysis %in% c("Primary: pure plastic, CRED", "Secondary: co-contaminants included")) %>%
  mutate(set = ifelse(grepl("^Primary", analysis), "P", "S"))

mod_lab <- c(taxa_class = "Taxonomic class", polymer_test = "Polymer", env_relevant = "Environmentally relevant",
             lifestage_test = "Life stage exposed", route_test = "Exposure route",
             size_bin = "Particle size class", log_size = "Particle size (continuous)",
             log_conc = "Concentration (continuous)")
mod_ord <- names(mod_lab)

# level rows: records / studies / % change (95% CI), for each analysis
wide <- lev %>%
  mutate(est = fpc(pct_change, pct_lb, pct_ub)) %>%
  select(moderator, level, set, k_es, k_studies, est) %>%
  pivot_wider(names_from = set, values_from = c(k_es, k_studies, est), names_vary = "slowest")

# keep a stable level order as in the CSV
ord <- lev %>% filter(set == "P") %>% distinct(moderator, level) %>% bind_rows(
  lev %>% filter(set == "S") %>% distinct(moderator, level)) %>% distinct() %>%
  mutate(o = row_number())
wide <- wide %>% left_join(ord, by = c("moderator", "level")) %>%
  mutate(m = match(moderator, mod_ord)) %>% arrange(m, o)

# continuous moderators: label the slope row clearly
wide <- wide %>% mutate(level = case_when(
  moderator == "log_size" ~ "% change per tenfold increase in size",
  moderator == "log_conc" ~ "% change per tenfold increase in concentration",
  TRUE ~ level))

# omnibus tests go on a header row for each moderator
# Holm is recomputed here across the moderator tests shown in each analysis (contaminant type is already filtered out above)
om <- omn %>%
  group_by(set) %>% mutate(p_holm = p.adjust(p_robust, "holm")) %>% ungroup() %>%
  select(moderator, set, p_robust, p_holm) %>%
  pivot_wider(names_from = set, values_from = c(p_robust, p_holm), names_vary = "slowest")

hdr <- om %>% mutate(m = match(moderator, mod_ord)) %>% arrange(m) %>%
  transmute(moderator,
            omni_P = ifelse(is.na(p_robust_P), "NE", sprintf("%s (Holm %s)", fp(p_robust_P), fp(p_holm_P))),
            omni_S = ifelse(is.na(p_robust_S), "NE", sprintf("%s (Holm %s)", fp(p_robust_S), fp(p_holm_S))))

tab <- wide %>% left_join(hdr, by = "moderator") %>%
  transmute(moderator, Level = level,
            recP = ifelse(is.na(k_es_P), "–", as.character(k_es_P)), stP = ifelse(is.na(k_studies_P), "–", as.character(k_studies_P)),
            estP = ifelse(is.na(est_P), "–", est_P),
            recS = ifelse(is.na(k_es_S), "–", as.character(k_es_S)), stS = ifelse(is.na(k_studies_S), "–", as.character(k_studies_S)),
            estS = ifelse(is.na(est_S), "–", est_S),
            omni_P, omni_S)

# one block per moderator; the omnibus p is shown in the block heading
kb <- kbl(tab %>% select(Level, recP, stP, estP, recS, stS, estS), format = "html", booktabs = TRUE,
          align = c("l", rep("r", 6)),
          caption = "Table 3. Effects of microplastic exposure on hatching success by moderator level (lnRR, % change relative to control).",
          col.names = c("Level", "Records", "Studies", "% change (95% CI)", "Records", "Studies", "% change (95% CI)")) %>%
  add_header_above(c(" " = 1, "Primary (pure plastic)" = 3, "Secondary (pure plastic + co-contaminant)" = 3)) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"), full_width = FALSE, position = "left") %>%
  row_spec(0, bold = TRUE)

grp <- rle(tab$moderator); start <- 1
for (i in seq_along(grp$lengths)) {
  m <- grp$values[i]; h <- hdr %>% filter(moderator == m)
  kb <- pack_rows(kb, sprintf("%s — omnibus p: primary %s; secondary %s", mod_lab[[m]], h$omni_P, h$omni_S),
                  start, start + grp$lengths[i] - 1, label_row_css = "text-align:left; font-weight:bold;")
  start <- start + grp$lengths[i]
}

kb %>% footnote(general = paste(
  "Estimates are from separate multilevel models per moderator (the moderator replaces the single intercept); levels are shown as group means.",
  "Omnibus p is the cluster-robust (CR2) F test, with Holm adjustment across the moderator tests shown in each analysis; NE = not estimable (df < 4 or fewer than 3 studies per level).",
  "Confidence intervals are shown only where robust df ≥ 4; otherwise the point estimate is given with (NE).",
  "Power is low: several levels are supported by fewer than 5 studies. Continuous moderators show the slope per tenfold increase.",
  "The planned contrast between co-contaminant and pure-plastic records is not a moderator test and is reported in Table 2."),
  general_title = "Note: ", footnote_as_chunk = TRUE)


#citations
library(grateful)
cite_packages(
  pkgs    = "Session",
  omit    = c("base", "maps", "Matrix", "metadat", "numDeriv"),
  out.dir = "."
)

