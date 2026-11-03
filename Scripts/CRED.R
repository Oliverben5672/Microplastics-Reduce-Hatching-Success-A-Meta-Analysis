library(tidyverse)
library(janitor)
library(metafor)
library(readr)
CRED_FINAL <- read_csv("Data/CRED_FINAL.csv")
View(CRED_FINAL)

# 1. Clean names, drop empty Excel rows, tidy text -------------------------
dat <- CRED_FINAL %>%
  clean_names() %>%
  filter(!is.na(study_id)) %>%
  rename(contam    = contaminant_category,
         pathway   = exposure_pathway,
         lifestage = exposure_lifestage) %>%
  select(-any_of("x25")) %>%
  mutate(across(tidyselect::where(is.character),
                ~ .x %>% str_replace_all("\u00a0", " ") %>%
                  str_replace_all("\u00c2", "") %>% str_squish()),
         es_id = row_number())

# 2. Exclusions and SD imputation -------------------------------------------
exclude_ids <- c(2, 28)   # n = 1 per arm
impute_ids  <- c(3, 27)   # flagged IMPUTE in Excel

excluded_log <- dat %>% filter(study_id %in% exclude_ids) %>%
  distinct(study_id, reference, control_reps, treat_reps)
print(excluded_log)

dat <- dat %>% filter(!study_id %in% exclude_ids)
stopifnot(!any(dat$control_reps < 2 | dat$treat_reps < 2))

class_sd <- dat %>%
  select(taxa_class, control_sd, treat_sd) %>%
  pivot_longer(c(control_sd, treat_sd), values_to = "sd") %>%
  filter(!is.na(sd)) %>%
  group_by(taxa_class) %>%
  summarise(mean_sd = mean(sd), n_sd = n(), .groups = "drop")
print(class_sd)

dat <- dat %>%
  left_join(class_sd, by = "taxa_class") %>%
  mutate(sd_imputed = study_id %in% impute_ids & (is.na(control_sd) | is.na(treat_sd)),
         control_sd = if_else(sd_imputed & is.na(control_sd), mean_sd, control_sd),
         treat_sd   = if_else(sd_imputed & is.na(treat_sd),   mean_sd, treat_sd)) %>%
  select(-mean_sd, -n_sd)

stopifnot(!anyNA(dat$control_sd), !anyNA(dat$treat_sd))

# 3. Effect sizes (negative g = reduced hatching) ---------------------------
dat <- escalc(measure = "SMD",
              m1i = treat_mhs,   sd1i = treat_sd,   n1i = treat_reps,
              m2i = control_mhs, sd2i = control_sd, n2i = control_reps,
              data = dat)
stopifnot(!anyNA(dat$yi), !anyNA(dat$vi), all(dat$vi > 0))

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
      polymer %in% c("Mix","CTR","AFH")  ~ polymer,
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

# Testable levels only (>= 3 studies under CRED); the rest stay descriptive
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
library(rotl) 
library(ape)
library(phytools)   # load after cleaning: ape masks dplyr::where

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

dat <- dat %>%
  left_join(lookup, by = "species") %>%
  mutate(control_id = paste(study_id, control_mhs, control_sd, control_reps, sep = "_"))

stopifnot(!anyNA(dat$species_phylo), all(dat$species_phylo %in% rownames(phylo_cor)))
dat %>% count(control_id) %>% count(n, name = "n_control_groups")

# 7. Analysis datasets (always last, after every column exists) -------------
dat_all       <- dat
dat_cred      <- dat %>% filter(cred_r != "R3")
dat_dose      <- dat %>% filter(conc_class == "aqueous")
dat_dose_cred <- dat_dose %>% filter(cred_r != "R3")

stopifnot(all(c("control_id", "species_phylo", "conc_mg_L") %in% names(dat_cred)))

# 8. Model function and the three base CRED models --------------------------
fit_cred <- function(data, struct = c("simple","phylo"), mods = ~ 1, rho = 0.5, method = "REML") {
  struct <- match.arg(struct)
  V <- vcalc(vi, cluster = control_id, obs = es_id, rho = rho, data = data)
  random <- list(~ 1 | study_id, ~ 1 | es_id)
  if (struct == "phylo") random <- c(random, list(~ 1 | species_phylo))
  rma.mv(yi, V, mods = mods, random = random,
         R = if (struct == "phylo") list(species_phylo = phylo_cor),
         data = data, method = method, test = "t")
}

# Class support: >= 3 STUDIES
class_support <- dat_cred %>% group_by(taxa_class) %>%
  summarise(k_es = n(), k_studies = n_distinct(study_id), .groups = "drop")
print(class_support)
dat_tax <- dat_cred %>% filter(taxa_class %in% class_support$taxa_class[class_support$k_studies >= 3])

m_simple <- fit_cred(dat_cred, "simple")
m_phylo  <- fit_cred(dat_cred, "phylo")

m_tax0     <- fit_cred(dat_tax, "simple")                          # same rows, no class
m_tax      <- fit_cred(dat_tax, "simple", mods = ~ 0 + taxa_class) # class-specific means
m_tax_omni <- fit_cred(dat_tax, "simple", mods = ~ taxa_class)     # QM: do classes differ?

AIC(m_simple, m_phylo)
anova(m_simple, m_phylo)
m_tax_omni

# 1. Are the variance components identifiable? (a clear peak = yes, flat = no)
profile(m_phylo)

# 2. Phylogenetic signal, as in the old paper
mean_es <- dat_cred %>% group_by(species_phylo) %>% summarise(m = mean(yi)) %>% deframe()
tr <- keep.tip(tree, names(mean_es)); mean_es <- mean_es[tr$tip.label]
phylosig(tr, mean_es, method = "lambda", test = TRUE)
phylosig(tr, mean_es, method = "K", test = TRUE, nsim = 1000)

# 3. Sensitivity of the overall effect to rho
purrr::map_dfr(c(0.2, 0.5, 0.8), \(r) {
  m <- fit_cred(dat_cred, "simple", rho = r)
  tibble(rho = r, est = as.numeric(coef(m)), ci_lb = m$ci.lb, ci_ub = m$ci.ub,
         s2_study = m$sigma2[1], s2_es = m$sigma2[2])
})

# 4. Robust SEs, clustered on study
robust(m_simple, cluster = dat_cred$study_id)

# 5. Influence: leave one study out
purrr::map_dfr(unique(dat_cred$study_id), \(s) {
  m <- fit_cred(filter(dat_cred, study_id != s), "simple")
  tibble(dropped = s, est = as.numeric(coef(m)), ci_lb = m$ci.lb, ci_ub = m$ci.ub)
}) %>% arrange(est)

# Descriptive summary for the classes that can't be tested
dat_cred %>%
  filter(!taxa_class %in% c("Actinopterygii","Crustacea")) %>%
  group_by(taxa_class) %>%
  summarise(k_es = n(), k_studies = n_distinct(study_id),
            mean_g = mean(yi), min_g = min(yi), max_g = max(yi))


#9 single moderators -------------------------------------------------
dat_dose_cred <- dat_dose_cred %>% mutate(log_conc = log10(conc_mg_L))

fit_mod <- function(var, data, struct = "simple", min_k = 3, continuous = FALSE) {
  d <- data %>% filter(!is.na(.data[[var]]))
  
  if (!continuous) {
    ok <- d %>% group_by(lvl = as.character(.data[[var]])) %>%
      summarise(ks = n_distinct(study_id), .groups = "drop") %>%
      filter(ks >= min_k) %>% pull(lvl)
    d <- d %>% filter(as.character(.data[[var]]) %in% ok)
    d[[var]] <- droplevels(if (is.factor(d[[var]])) d[[var]] else factor(d[[var]]))
  }
  
  m_om <- fit_cred(d, struct, mods = as.formula(paste("~", var)))   # intercept + contrasts: omnibus test
  omni <- tibble(moderator = var, struct = struct, k_es = m_om$k,
                 k_studies = n_distinct(d$study_id),
                 QM = as.numeric(m_om$QM), df1 = m_om$QMdf[1], df2 = m_om$QMdf[2],
                 p = m_om$QMp)
  
  if (continuous) {
    lev <- tibble(moderator = var, struct = struct, level = "slope per unit",
                  k_es = m_om$k, k_studies = n_distinct(d$study_id),
                  g = as.numeric(coef(m_om)[2]), ci_lb = m_om$ci.lb[2], ci_ub = m_om$ci.ub[2])
  } else {
    m_lv <- fit_cred(d, struct, mods = as.formula(paste("~ 0 +", var)))   # level-specific means
    n_tab <- d %>% group_by(level = as.character(.data[[var]])) %>%
      summarise(k_es = n(), k_studies = n_distinct(study_id), .groups = "drop")
    lev <- tibble(moderator = var, struct = struct,
                  level = sub(paste0("^", var), "", rownames(m_lv$beta)),
                  g = as.numeric(m_lv$beta), ci_lb = m_lv$ci.lb, ci_ub = m_lv$ci.ub) %>%
      left_join(n_tab, by = "level")
  }
  list(omnibus = omni, levels = lev)
}

data_list <- list(dat_cred = dat_cred, dat_dose_cred = dat_dose_cred)

specs <- tribble(
  ~var,             ~dset,           ~continuous,
  "taxa_class",     "dat_cred",      FALSE,
  "polymer_test",   "dat_cred",      FALSE,
  "contam",         "dat_cred",      FALSE,
  "env_relevant",   "dat_cred",      FALSE,
  "lifestage_test", "dat_cred",      FALSE,
  "route_test",     "dat_cred",      FALSE,
  "size_bin",       "dat_cred",      FALSE,
  "log_size",       "dat_cred",      TRUE,
  "log_conc",       "dat_dose_cred", TRUE      # aqueous only, mg/L
)

grid <- tidyr::crossing(specs, struct = c("simple", "phylo_sp", "phylo"))
safe_fit <- purrr::possibly(fit_mod, otherwise = NULL)

results <- purrr::pmap(grid, \(var, dset, continuous, struct)
                       safe_fit(var, data_list[[dset]], struct, continuous = continuous))

grid[purrr::map_lgl(results, is.null), ]          # any model that failed to fit

omnibus    <- purrr::map_dfr(results, "omnibus")  %>% group_by(struct) %>%
  mutate(p_holm = p.adjust(p, "holm")) %>% ungroup()
levels_tbl <- purrr::map_dfr(results, "levels")

omnibus %>% arrange(moderator, struct) %>% mutate(across(c(QM, p, p_holm), ~ round(.x, 3)))
levels_tbl %>% filter(struct == "simple") %>% mutate(across(c(g, ci_lb, ci_ub), ~ round(.x, 2)))

library(clubSandwich)
chk <- function(var, data) {
  d <- data %>% filter(!is.na(.data[[var]]))
  m <- fit_cred(d, "simple", mods = as.formula(paste("~", var)))
  anova(robust(m, cluster = d$study_id, clubSandwich = TRUE), btt = 2:length(coef(m)))
}
chk("contam", dat_cred); chk("env_relevant", dat_cred)
chk("lifestage_test", dat_cred); chk("log_conc", dat_dose_cred)

purrr::map_dfr(unique(dat_dose_cred$study_id), \(s) {
  m <- fit_cred(filter(dat_dose_cred, study_id != s), "simple", mods = ~ log_conc)
  tibble(dropped = s, slope = as.numeric(coef(m)[2]), p = m$pval[2])
}) %>% arrange(desc(p))

with(dat_cred, table(contam, env_relevant))
with(dat_cred, table(lifestage_test, contam, useNA = "ifany"))
dat_dose_cred %>% group_by(env_relevant) %>%
  summarise(k_studies = n_distinct(study_id), med = median(conc_mg_L), min = min(conc_mg_L), max = max(conc_mg_L))
dat_cred %>% group_by(contam) %>%
  summarise(k_studies = n_distinct(study_id), prop_env_yes = mean(env_relevant == "Yes"))


dat_tax <- dat_tax %>% mutate(log_conc = log10(conc_mg_L))

# A. Overall estimates for every structure
ov <- function(m, label, d) tibble(model = label, k_es = m$k, k_studies = n_distinct(d$study_id),
                                   g = as.numeric(coef(m)), ci_lb = m$ci.lb, ci_ub = m$ci.ub, p = m$pval)
bind_rows(ov(m_simple, "simple", dat_cred), ov(m_phylo, "phylo", dat_cred),
          ov(m_phylo_sp, "phylo + species", dat_cred), ov(m_tax0, "taxa subset (no class)", dat_tax)) %>%
  mutate(across(c(g, ci_lb, ci_ub), ~ round(.x, 2)), p = signif(p, 2))

# B. Taxa-adjusted moderator models (class + moderator), with robust p
adj <- function(var, data = dat_tax) {
  d <- data %>% filter(!is.na(.data[[var]])) %>% droplevels()
  m <- fit_cred(d, "simple", mods = as.formula(paste("~ taxa_class +", var)))
  idx <- grep(paste0("^", var), rownames(m$beta))
  rb <- robust(m, cluster = d$study_id, clubSandwich = TRUE)
  tibble(moderator = var, model = "taxa + moderator", k_es = m$k, k_studies = n_distinct(d$study_id),
         p_model = anova(m, btt = idx)$QMp, p_robust = anova(rb, btt = idx)$QMp)
}
mods <- c("polymer_test","contam","env_relevant","lifestage_test","route_test","size_bin","log_size","log_conc")
purrr::map_dfr(mods, purrr::possibly(adj, otherwise = NULL)) %>%
  mutate(across(starts_with("p_"), ~ round(.x, 3)))
