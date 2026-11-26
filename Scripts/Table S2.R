# ══════════════════════════════════════════════════════════════════════════════
# Supplementary table: every model structure x every dataset configuration
# (replaces the old single-table script; reads the re-run CSVs, no refitting)
# Run after CRED_39studies_v5.R, or point out_dir at the re-run output folder.
# ══════════════════════════════════════════════════════════════════════════════
library(dplyr)
library(readr)
library(kableExtra)

out_dir <- "results_39studies_v5"

sbs <- read_csv(file.path(out_dir, "sensitivity_by_structure.csv"), show_col_types = FALSE)

# ── formatting helpers ───────────────────────────────────────────────────────
fmt_p <- function(p) ifelse(is.na(p), "NE",
                            ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# robust p reported only where df >= 4 (otherwise NE); CR1 (species/phylo) flagged
robust_p <- function(p, df, infer) {
  out <- ifelse(is.na(p) | is.na(df) | df < 4, "NE", fmt_p(p))
  ifelse(infer == "CR1" & out != "NE", paste0(out, "†"), out)
}

est_ci <- function(est, lb, ub, scale) {
  if (scale == "lnRR") {
    f <- function(x) 100 * (exp(x) - 1)
    sprintf("%.1f%%", f(est)) -> e
    ci <- sprintf("[%.1f%%, %.1f%%]", f(lb), f(ub))
  } else {
    e  <- sprintf("%.2f", est)
    ci <- sprintf("[%.2f, %.2f]", lb, ub)
  }
  list(e = e, ci = ci)
}

struct_lab <- c(simple = "Simple (study + record)",
                species = "+ species",
                phylo = "+ phylogeny",
                phylo_sp = "+ phylogeny + species")
struct_ord <- names(struct_lab)

# analysis groups, in presentation order (anything not listed is appended)
grp_ord <- c("Base dataset (CRED)", "R3 retained", "Secondary (co-contaminants incl.)",
             "Excl. PlotDigitiser", "Excl. imputed SDs",
             "Constant class-mean SD (old rule, all imputed studies)",
             "Excl. zero-hatching arms",
             "Excl. 3 largest-effect studies (post hoc)",
             "Excl. |z| > 3 (lnRR model)", "Excl. |z| > 3 (g model)",
             "Excl. vi > 5 (g)", "Excl. either rule (g)")

make_table <- function(sc) {
  d <- sbs %>% filter(scale == sc)
  grp_all <- c(intersect(grp_ord, unique(d$analysis)),
               setdiff(unique(d$analysis), grp_ord))
  d <- d %>%
    mutate(analysis = factor(analysis, levels = grp_all),
           struct   = factor(struct, levels = struct_ord)) %>%
    arrange(analysis, struct)
  
  ec <- Map(est_ci, d$est, d$ci_lb, d$ci_ub, sc)
  tab <- tibble(
    Model    = paste0(struct_lab[as.character(d$struct)],
                      ifelse(!is.na(d$note), "‡", "")),
    k        = paste0(d$k_es, " (", d$k_studies, ")"),
    Estimate = vapply(ec, `[[`, "", "e"),
    CI_95    = vapply(ec, `[[`, "", "ci"),
    p_robust = robust_p(d$p_robust, d$df, d$infer),
    p_model  = fmt_p(d$p_model),
    df       = ifelse(is.na(d$df), "—", sprintf("%.1f", d$df))
  )
  
  grp_n <- table(d$analysis)[grp_all]
  starts <- cumsum(c(1, head(as.integer(grp_n), -1)))
  ends   <- cumsum(as.integer(grp_n))
  
  k <- tab %>%
    kbl(col.names = c("Model structure", "Records (studies)",
                      if (sc == "lnRR") "Effect (% change)" else "Hedges' g",
                      "95% CI", "Robust p", "Model-based p", "Robust df"),
        align = c("l", "c", "c", "l", "c", "c", "c"),
        booktabs = TRUE,
        caption = paste0(
          "All model structures across dataset configurations (",
          if (sc == "lnRR") "lnRR, shown as % change = 100·(e^lnRR − 1)"
          else "Hedges' g",
          "). Robust p = CR2 cluster-robust (Satterthwaite df) for the simple model; ",
          "NE = not estimable (df < 4). † CR1 robust inference for species/phylogeny models ",
          "(anti-conservative; interpret with caution). ‡ Optimiser did not report ",
          "full convergence.")) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = TRUE, font_size = 13) %>%
    row_spec(0, background = "#2E75B6", color = "white", bold = TRUE)
  for (i in seq_along(grp_all))
    k <- pack_rows(k, grp_all[i], starts[i], ends[i], bold = TRUE, color = "#2E75B6")
  k
}

tab_all_lnrr <- make_table("lnRR")
tab_all_g    <- make_table("Hedges' g")

tab_all_lnrr
tab_all_g

# Optional exports
# save_kable(tab_all_lnrr, "Table_S_all_models_lnRR.html")
# save_kable(tab_all_g,    "Table_S_all_models_g.html")