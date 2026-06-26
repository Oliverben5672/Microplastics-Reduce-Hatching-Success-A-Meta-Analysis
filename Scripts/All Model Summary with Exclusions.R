library(dplyr)
library(metafor)
library(rotl)
library(ape)
library(phytools)
library(kableExtra)

# -- fix name encoding issue
MetaAnalysisTest$Species[29] <- "Calanus finmarchicus"

test <- MetaAnalysisTest[!is.na(MetaAnalysisTest$Control_MHS), ]

# ── Build tree and cor_matrix ─────────────────────────────────────────────
species_list <- unique(test$Species)
taxa         <- tnrs_match_names(species_list, context_name = "Animals")

failed <- taxa$search_string[is.na(taxa$ott_id)]
failed <- failed[!is.na(failed)]
if (length(failed) > 0) {
  taxa_retry <- tnrs_match_names(failed)
  taxa$ott_id[is.na(taxa$ott_id) & !is.na(taxa$search_string)] <- taxa_retry$ott_id
}

tree <- tol_induced_subtree(ott_ids = na.omit(ott_id(taxa)))
tree$tip.label <- gsub("_ott[0-9]+",    "",  tree$tip.label)
tree$tip.label <- gsub("_",             " ", tree$tip.label)
tree$tip.label <- gsub("\\s*\\(.*\\)$", "",  tree$tip.label)
tree$tip.label[tree$tip.label == "Cornu aspersum"] <- "Cantareus aspersus"

if (is.null(tree$edge.length)) tree <- compute.brlen(tree)
if (!is.ultrametric(tree))     tree <- force.ultrametric(tree, method = "extend")

cor_matrix <- vcv(tree, corr = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

prepare_phylo_data <- function(escalc_input, cor_matrix, tree) {
  df         <- escalc_input %>% filter(Species %in% tree$tip.label)
  df$obs_id  <- seq_len(nrow(df))
  df$Species <- as.character(df$Species)
  df$phylo   <- df$Species
  return(df)
}

prepare_class_data <- function(escalc_input) {
  class_counts  <- table(escalc_input$Class)
  valid_classes <- names(class_counts[class_counts >= 3])
  df            <- escalc_input %>% filter(Class %in% valid_classes)
  df$Class      <- as.factor(df$Class)
  return(df)
}

fit_simple_model <- function(data) {
  rma(yi, vi, data = data, method = "REML")
}

fit_class_model <- function(data) {
  rma(yi, vi, mods = ~ Class - 1, data = data, method = "REML")
}

fit_phylo_model <- function(data, cor_matrix) {
  rma.mv(
    yi, vi,
    random = list(
      ~ 1 | obs_id,
      ~ 1 | Paper_ID,
      ~ 1 | Species,
      ~ 1 | phylo
    ),
    R      = list(phylo = cor_matrix),
    data   = data,
    method = "REML"
  )
}

# ── Detect outlier flagged by both Cook's distance and studentised residuals ──
# Returns the Paper_ID of the observation flagged by both methods, or NULL
detect_outlier <- function(model, data) {
  cooks    <- cooks.distance(model)
  rstudent <- rstudent(model)
  
  cook_flag  <- which(cooks > mean(cooks, na.rm = TRUE) + 3 * sd(cooks, na.rm = TRUE))
  resid_flag <- which(abs(rstudent$z) > 1.96)
  both       <- intersect(cook_flag, resid_flag)
  
  if (length(both) == 0) {
    cat("  No observation flagged by both methods.\n")
    return(NULL)
  }
  
  flagged_papers <- data$Paper_ID[both]
  cat("  Outlier(s) flagged by both methods: Paper_ID(s)", paste(flagged_papers, collapse = ", "), "\n")
  return(flagged_papers)
}

# ══════════════════════════════════════════════════════════════════════════════
# EXTRACT FUNCTIONS 
# ══════════════════════════════════════════════════════════════════════════════

extract_simple_row <- function(model, label) {
  s      <- summary(model)
  est    <- round(s$b[1], 3)
  lb     <- round(s$ci.lb, 3)
  ub     <- round(s$ci.ub, 3)
  pv     <- s$pval
  pv_fmt <- ifelse(pv < 0.001, "<0.001***",
                   ifelse(pv < 0.01,  paste0(round(pv,3),"**"),
                          ifelse(pv < 0.05,  paste0(round(pv,3),"*"),
                                 as.character(round(pv,3)))))
  data.frame(
    Label    = label,
    k        = model$k,
    Estimate = as.character(est),
    CI_95    = paste0("[", lb, ", ", ub, "]"),
    p_value  = pv_fmt,
    AIC      = round(AIC(model), 2),
    Q_stat   = paste0(round(s$QE, 2), " (", model$k - 1, ")"),
    Het_var  = paste0("I\u00b2=", round(s$I2, 2), "%"),
    stringsAsFactors = FALSE
  )
}

extract_class_rows <- function(model, label) {
  s       <- summary(model)
  classes <- gsub("^Class", "", rownames(s$b))
  est     <- round(s$b[,1], 3)
  lb      <- round(s$ci.lb, 3)
  ub      <- round(s$ci.ub, 3)
  pv      <- s$pval
  pv_fmt  <- ifelse(pv < 0.001, "<0.001***",
                    ifelse(pv < 0.01,  paste0(round(pv,3),"**"),
                           ifelse(pv < 0.05,  paste0(round(pv,3),"*"),
                                  as.character(round(pv,3)))))
  tau2    <- round(s$tau2, 3)
  
  main <- data.frame(
    Label    = label,
    k        = model$k,
    Estimate = "\u2014",
    CI_95    = "Varies by class",
    p_value  = "\u2014",
    AIC      = round(AIC(model), 2),
    Q_stat   = paste0(round(s$QE, 2), " (", model$k - 1, ")"),
    Het_var  = paste0("I\u00b2=", round(s$I2, 2), "%"),
    stringsAsFactors = FALSE
  )
  
  subs <- data.frame(
    Label    = paste0("  \u21b3 ", classes),
    k        = NA,
    Estimate = as.character(est),
    CI_95    = paste0("[", lb, ", ", ub, "]"),
    p_value  = pv_fmt,
    AIC      = NA,
    Q_stat   = NA,
    Het_var  = c(paste0("\u03c4\u00b2=", tau2), rep("", length(classes)-1)),
    stringsAsFactors = FALSE
  )
  
  bind_rows(main, subs)
}

extract_phylo_row <- function(model, label) {
  s      <- summary(model)
  est    <- round(s$b[1], 3)
  lb     <- round(s$ci.lb, 3)
  ub     <- round(s$ci.ub, 3)
  pv     <- s$pval
  pv_fmt <- ifelse(pv < 0.001, "<0.001***",
                   ifelse(pv < 0.01,  paste0(round(pv,3),"**"),
                          ifelse(pv < 0.05,  paste0(round(pv,3),"*"),
                                 as.character(round(pv,3)))))
  sp_var <- round(model$sigma2[3], 3)
  data.frame(
    Label    = label,
    k        = model$k,
    Estimate = as.character(est),
    CI_95    = paste0("[", lb, ", ", ub, "]"),
    p_value  = pv_fmt,
    AIC      = round(AIC(model), 2),
    Q_stat   = paste0(round(s$QE, 2), " (", s$QEdf, ")"),
    Het_var  = paste0("\u03c3\u00b2sp=", sp_var),
    stringsAsFactors = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: ALL STUDIES — fit models and detect outliers
# ══════════════════════════════════════════════════════════════════════════════

escalc_all <- escalc(
  measure = "SMD",
  m1i = Treat_MHS, sd1i = Treat_SD, n1i = Treat_Reps,
  m2i = Control_MHS, sd2i = Control_SD, n2i = Control_Reps,
  data = test
)

res_simple_all <- fit_simple_model(escalc_all)
res_class_all  <- fit_class_model(prepare_class_data(escalc_all))
res_phylo_all  <- fit_phylo_model(prepare_phylo_data(escalc_all, cor_matrix, tree), cor_matrix)

cat("\nMODEL 1a: Simple — All studies\n"); print(summary(res_simple_all))
cat("\nMODEL 1b: Class  — All studies\n"); print(summary(res_class_all))
cat("\nMODEL 1c: Phylo  — All studies\n"); print(summary(res_phylo_all))

# ── Detect outliers per model type on full dataset ────────────────────────
cat("\nOutlier detection — Simple (all studies):\n")
outlier_simple <- detect_outlier(res_simple_all, escalc_all)

cat("\nOutlier detection — Class (all studies):\n")
outlier_class  <- detect_outlier(res_class_all,  prepare_class_data(escalc_all))

cat("\nOutlier detection — Phylo (all studies):\n")
outlier_phylo  <- detect_outlier(res_phylo_all,  prepare_phylo_data(escalc_all, cor_matrix, tree))

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: SYNERGISTIC STUDIES REMOVED (Paper IDs: 3a, 3b, 20)
# ══════════════════════════════════════════════════════════════════════════════

test_nosyn   <- test %>% filter(!Paper_ID %in% c("3a", "3b", "20"))
escalc_nosyn <- escalc(
  measure = "SMD",
  m1i = Treat_MHS, sd1i = Treat_SD, n1i = Treat_Reps,
  m2i = Control_MHS, sd2i = Control_SD, n2i = Control_Reps,
  data = test_nosyn
)

res_simple_nosyn <- fit_simple_model(escalc_nosyn)
res_class_nosyn  <- fit_class_model(prepare_class_data(escalc_nosyn))
res_phylo_nosyn  <- fit_phylo_model(prepare_phylo_data(escalc_nosyn, cor_matrix, tree), cor_matrix)

cat("\nMODEL 2a: Simple — Synergistic removed\n"); print(summary(res_simple_nosyn))
cat("\nMODEL 2b: Class  — Synergistic removed\n"); print(summary(res_class_nosyn))
cat("\nMODEL 2c: Phylo  — Synergistic removed\n"); print(summary(res_phylo_nosyn))

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: OUTLIER REMOVED — model-specific exclusions
# ══════════════════════════════════════════════════════════════════════════════

# Simple: remove outlier detected by simple model
escalc_nooutlier_simple <- escalc_all %>% filter(!Paper_ID %in% outlier_simple)
res_simple_nooutlier    <- fit_simple_model(escalc_nooutlier_simple)

simple_outlier_label <- paste("Paper(s)", paste(outlier_simple, collapse = ", "))
cat("\nMODEL 3a: Simple — Outlier removed (", simple_outlier_label, ")\n")
print(summary(res_simple_nooutlier))

# Class: remove outlier detected by class model
escalc_nooutlier_class <- prepare_class_data(escalc_all) %>% filter(!Paper_ID %in% outlier_class)
res_class_nooutlier    <- fit_class_model(escalc_nooutlier_class)

class_outlier_label <- paste("Paper(s)", paste(outlier_class, collapse = ", "))
cat("\nMODEL 3b: Class — Outlier removed (", class_outlier_label, ")\n")
print(summary(res_class_nooutlier))

# Phylo: remove outlier detected by phylo model
phylo_data_nooutlier <- prepare_phylo_data(escalc_all, cor_matrix, tree) %>%
  filter(!Paper_ID %in% outlier_phylo)
res_phylo_nooutlier  <- fit_phylo_model(phylo_data_nooutlier, cor_matrix)

phylo_outlier_label <- paste("Paper(s)", paste(outlier_phylo, collapse = ", "))
cat("\nMODEL 3c: Phylo — Outlier removed (", phylo_outlier_label, ")\n")
print(summary(res_phylo_nooutlier))

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: SYNERGISTIC AND OUTLIER REMOVED — model-specific exclusions
# ══════════════════════════════════════════════════════════════════════════════

# Simple
escalc_clean_simple <- escalc_nosyn %>% filter(!Paper_ID %in% outlier_simple)
res_simple_clean    <- fit_simple_model(escalc_clean_simple)
cat("\nMODEL 4a: Simple — Synergistic + outlier removed\n"); print(summary(res_simple_clean))

# Class
escalc_clean_class <- prepare_class_data(escalc_nosyn) %>% filter(!Paper_ID %in% outlier_class)
res_class_clean    <- fit_class_model(escalc_clean_class)
cat("\nMODEL 4b: Class — Synergistic + outlier removed\n"); print(summary(res_class_clean))

# Phylo
phylo_data_clean <- prepare_phylo_data(escalc_nosyn, cor_matrix, tree) %>%
  filter(!Paper_ID %in% outlier_phylo)
res_phylo_clean  <- fit_phylo_model(phylo_data_clean, cor_matrix)
cat("\nMODEL 4c: Phylo — Synergistic + outlier removed\n"); print(summary(res_phylo_clean))

# ══════════════════════════════════════════════════════════════════════════════
# BUILD TABLE
# ══════════════════════════════════════════════════════════════════════════════

g1 <- bind_rows(
  extract_simple_row(res_simple_all,  "Simple \u2014 All studies"),
  extract_class_rows(res_class_all,   "Class \u2014 All studies"),
  extract_phylo_row( res_phylo_all,   "Phylo \u2014 All studies")
)

g2 <- bind_rows(
  extract_simple_row(res_simple_nosyn, "Simple \u2014 Synergistic removed (3a, 3b, 20)"),
  extract_class_rows(res_class_nosyn,  "Class \u2014 Synergistic removed (3a, 3b, 20)"),
  extract_phylo_row( res_phylo_nosyn,  "Phylo \u2014 Synergistic removed (3a, 3b, 20)")
)

g3 <- bind_rows(
  extract_simple_row(res_simple_nooutlier,
                     paste0("Simple \u2014 Outlier removed (Paper ", paste(outlier_simple, collapse=", "), ")")),
  extract_class_rows(res_class_nooutlier,
                     paste0("Class \u2014 Outlier removed (Paper ", paste(outlier_class, collapse=", "), ")")),
  extract_phylo_row( res_phylo_nooutlier,
                     paste0("Phylo \u2014 Outlier removed (Paper ", paste(outlier_phylo, collapse=", "), ")"))
)

g4 <- bind_rows(
  extract_simple_row(res_simple_clean,
                     paste0("Simple \u2014 All removed (3a, 3b, 20 + Paper ", paste(outlier_simple, collapse=", "), ")")),
  extract_class_rows(res_class_clean,
                     paste0("Class \u2014 All removed (3a, 3b, 20 + Paper ", paste(outlier_class, collapse=", "), ")")),
  extract_phylo_row( res_phylo_clean,
                     paste0("Phylo \u2014 All removed (3a, 3b, 20 + Paper ", paste(outlier_phylo, collapse=", "), ")"))
)

final_table <- bind_rows(g1, g2, g3, g4)

# ── Row index helpers ─────────────────────────────────────────────────────
is_simple <- grepl("^Simple", final_table$Label)
is_class  <- grepl("^Class",  final_table$Label)
is_phylo  <- grepl("^Phylo",  final_table$Label)
is_sub    <- grepl("^  \u21b3", final_table$Label)

simple_rows <- which(is_simple)
class_rows  <- which(is_class)
phylo_rows  <- which(is_phylo)
sub_rows    <- which(is_sub)

g1_end <- nrow(g1)
g2_end <- g1_end + nrow(g2)
g3_end <- g2_end + nrow(g3)
g4_end <- g3_end + nrow(g4)

# ── Dynamic pack_rows labels ──────────────────────────────────────────────
outlier_group_label <- paste0(
  "Outlier Removed (Simple: Paper ", paste(outlier_simple, collapse=", "),
  "; Class: Paper ", paste(outlier_class, collapse=", "),
  "; Phylo: Paper ", paste(outlier_phylo, collapse=", "), ")"
)

all_removed_label <- paste0(
  "All Removed (3a, 3b, 20 + Simple: Paper ", paste(outlier_simple, collapse=", "),
  "; Class: Paper ", paste(outlier_class, collapse=", "),
  "; Phylo: Paper ", paste(outlier_phylo, collapse=", "), ")"
)

# ── Render table ──────────────────────────────────────────────────────────
final_table %>%
  kbl(
    col.names = c("Model / Class", "k", "Estimate (g)", "95% CI",
                  "p-value", "AIC", "Q statistic (df)", "I\u00b2 / \u03c3\u00b2sp"),
    align     = c("l","c","c","l","c","c","c","c"),
    caption   = paste(
      "Sensitivity analysis across meta-analytic model types and dataset configurations.",
      "Hedges' g as effect size (SMD); negative values indicate reduced hatching success",
      "under microplastic exposure. Simple = random-effects model with no moderators;",
      "Class = taxonomic class as moderator; Phylo = phylogenetically structured multilevel model.",
      "Outliers identified per model type by Cook's distance (>mean+3SD) AND studentised",
      "residuals (|z|>1.96); exclusions differ between model types as leverage differs by structure.",
      "\u03c3\u00b2sp = species-level variance; \u03c3\u00b2phylo \u2248 0 across all Phylo models."
    ),
    booktabs  = TRUE,
    na        = "\u2014"
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"),
    full_width        = TRUE,
    font_size         = 13
  ) %>%
  row_spec(0,           background = "#2E75B6", color = "white", bold = TRUE) %>%
  row_spec(simple_rows, background = "#FEF9E7",                  bold = TRUE) %>%
  row_spec(class_rows,  background = "#EAF2FB",                  bold = TRUE) %>%
  row_spec(phylo_rows,  background = "#EBF7EE",                  bold = TRUE) %>%
  row_spec(sub_rows,    background = "#F7FBFF", color = "#444444", italic = TRUE) %>%
  pack_rows("All Studies",              1,          g1_end, bold = TRUE, color = "#2E75B6") %>%
  pack_rows("Synergistic Removed (3a, 3b, 20)", g1_end+1, g2_end, bold = TRUE, color = "#2E75B6") %>%
  pack_rows(outlier_group_label,        g2_end+1,   g3_end, bold = TRUE, color = "#2E75B6") %>%
  pack_rows(all_removed_label,          g3_end+1,   g4_end, bold = TRUE, color = "#2E75B6") %>%
  footnote(
    general           = "Significance: *** p < 0.001  ** p < 0.01  * p < 0.05",
    general_title     = "",
    footnote_as_chunk = TRUE
  )