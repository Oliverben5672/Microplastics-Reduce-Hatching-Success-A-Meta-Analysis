library(tidyverse)
library(janitor)
library(metafor)
library(readr)
CRED_FINAL <- read_csv("Data/CRED_FINAL.csv")
View(CRED_FINAL)

library(tidyverse)
library(janitor)
library(metafor)

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
library(rotl); library(ape); library(phytools)   # load after cleaning: ape masks dplyr::where

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
fit_cred <- function(data, struct = c("simple","taxonomic","phylo"), mods = ~ 1, rho = 0.5) {
  struct <- match.arg(struct)
  V <- vcalc(vi, cluster = control_id, obs = es_id, rho = rho, data = data)
  random <- list(~ 1 | study_id, ~ 1 | es_id)
  if (struct != "simple") random <- c(random, list(~ 1 | species))
  if (struct == "phylo")  random <- c(random, list(~ 1 | species_phylo))
  rma.mv(yi, V, mods = mods, random = random,
         R = if (struct == "phylo") list(species_phylo = phylo_cor),
         data = data, method = "REML", test = "t")
}

m_simple <- fit_cred(dat_cred, "simple")
m_tax    <- fit_cred(dat_cred, "taxonomic")
m_phylo  <- fit_cred(dat_cred, "phylo")

AIC(m_simple, m_tax, m_phylo)
anova(m_simple, m_tax)
anova(m_tax, m_phylo)
