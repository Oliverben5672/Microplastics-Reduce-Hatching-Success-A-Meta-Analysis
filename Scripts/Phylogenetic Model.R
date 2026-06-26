library(dplyr)
library(metafor)
library(rotl)
library(ape)
library(phytools)

#rename column manually due to persistent encoding issue
MetaAnalysisTest$Species[29] <- "Calanus finmarchicus"
# ── 0. Clean encoding BEFORE anything else ────────────────────────────────
test$Species <- iconv(test$Species, from = "latin1", to = "UTF-8", sub = " ")
test$Species <- gsub("\u00a0", " ", test$Species)  # explicit non-breaking space
test$Species <- gsub("\u00c2", "", test$Species)    # strips the  artifact
test$Species <- trimws(test$Species)

# ── 1. Calculate effect sizes ─────────────────────────────────────────────
escalc_data <- escalc(
  measure = "SMD",
  m1i = Treat_MHS,
  sd1i = Treat_SD,
  n1i = Treat_Reps,
  m2i = Control_MHS,
  sd2i = Control_SD,
  n2i = Control_Reps,
  data = test
)

# Add observation-level ID for overdispersion random effect
escalc_data$obs_id <- seq_len(nrow(escalc_data))

# ── 2. Get phylogenetic tree from Open Tree of Life ───────────────────────
species_list <- unique(escalc_data$Species)
print(species_list)

taxa <- tnrs_match_names(species_list, context_name = "Animals")
print(taxa)

# Retry any NAs without context restriction
failed <- taxa$search_string[is.na(taxa$ott_id)]
failed <- failed[!is.na(failed)]  # exclude NA names before sending to API
if (length(failed) > 0) {
  taxa_retry <- tnrs_match_names(failed)
  taxa$ott_id[is.na(taxa$ott_id) & !is.na(taxa$search_string)] <- taxa_retry$ott_id
}

cat("Matched:", sum(!is.na(taxa$ott_id)), "/", nrow(taxa), "species\n")

tree <- tol_induced_subtree(ott_ids = na.omit(ott_id(taxa)))

# Clean tip labels: strip _ottXXXX suffixes and underscores
tree$tip.label <- gsub("_ott[0-9]+", "", tree$tip.label)
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Strip OTL parenthetical qualifiers e.g. "(species in domain Eukaryota)"
tree$tip.label <- gsub("\\s*\\(.*\\)$", "", tree$tip.label)

# Reconcile OTL accepted synonyms back to dataset names
tree$tip.label[tree$tip.label == "Cornu aspersum"] <- "Cantareus aspersus"

# Diagnostic: should both return character(0) after fixes
cat("\nSpecies in data not in tree:\n")
print(setdiff(unique(escalc_data$Species), tree$tip.label))
cat("\nTips in tree not in data:\n")
print(setdiff(tree$tip.label, unique(escalc_data$Species)))

# ── 3. Prepare tree ───────────────────────────────────────────────────────
if (is.null(tree$edge.length)) {
  tree <- compute.brlen(tree)
}

if (!is.ultrametric(tree)) {
  tree <- force.ultrametric(tree, method = "extend")
}

# ── 4. Match dataset species to tree tips ─────────────────────────────────
escalc_data <- escalc_data %>%
  filter(Species %in% tree$tip.label)

cat("Observations retained:", nrow(escalc_data), "\n")  # should be 29

tree <- keep.tip(tree, unique(escalc_data$Species))

# ── 5. Build phylogenetic correlation matrix ──────────────────────────────
cor_matrix <- vcv(tree, corr = TRUE)

# ── 6. Check phylogenetic signal ──────────────────────────────────────────
mean_es <- tapply(escalc_data$yi, escalc_data$Species, mean)
mean_es <- mean_es[tree$tip.label]

lambda_result <- phylosig(tree, mean_es, method = "lambda", test = TRUE)
K_result      <- phylosig(tree, mean_es, method = "K",      test = TRUE, nsim = 1000)

print(lambda_result)
print(K_result)

# ── 7. Prepare columns for random effects ────────────────────────────────
escalc_data$Species <- as.character(escalc_data$Species)
escalc_data$phylo   <- escalc_data$Species

# ── 8. Fit full phylogenetic meta-analytic model ──────────────────────────
res_phylo <- rma.mv(
  yi, vi,
  random = list(
    ~ 1 | obs_id,    # observation level: overdispersion
    ~ 1 | Paper_ID,  # study level: shared methodology/conditions per paper
    ~ 1 | Species,   # residual species variance not explained by phylogeny
    ~ 1 | phylo      # phylogenetic covariance (structured by cor_matrix)
  ),
  R = list(phylo = cor_matrix),
  data = escalc_data,
  method = "REML"
)

summary(res_phylo)

# ── 9. Compare models ──────────────────────────────────────────────────────
# No-phylogeny version for comparison
res_no_phylo <- rma.mv(
  yi, vi,
  random = list(
    ~ 1 | obs_id,
    ~ 1 | Paper_ID,
    ~ 1 | Species
  ),
  data = escalc_data,
  method = "REML"
)

summary(res_no_phylo)

anova(res_no_phylo, res_phylo)  # likelihood ratio test
AIC(res_no_phylo, res_phylo)    # AIC comparison

# ── 7. Compare with a no-phylogeny model ───────────────────────────────────
res_no_phylo <- rma.mv(
  yi, vi,
  random = list(
    ~ 1 | Species,
    ~ 1 | Paper_ID
  ),
  data = escalc_data,
  method = "REML"
)

# Likelihood ratio test
anova(res_no_phylo, res_phylo)
AIC(res_no_phylo, res_phylo)

#-------------------------------------
phylo_data <- escalc_data %>%
  filter(!is.na(yi),
         !is.na(vi)) %>%
  filter(Species %in% rownames(cor_matrix))  # or however you subsetted

slabs <- res_phylo$data$Paper_ID
species <- res_phylo$data$Species

forest(res_phylo,
       xlim = c(-30, 10),
       slab = slabs,
       ilab = species,
       ilab.xpos = -22,
       mlab = "Overall effect")
