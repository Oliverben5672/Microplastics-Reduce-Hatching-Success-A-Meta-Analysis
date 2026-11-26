# Publication bias assessment (v5) ---------------------------------------------------------------------------
# Run AFTER CRED_39studies_v5.R in the SAME R session: it uses that script's data (dat_pure_l, dat_pure, dat_sec_l, dat_sec),
# model function (fit_cred), helpers (rob_safe, inf_tbl, sc_lab, pct, study_est, st_type), palettes and out_dir.
#
# Design (replaces the old metafor::regtest / trimfill on a single univariate model):
#   * The TEST is run on the RECORD level in the multilevel model (study + observation random effects, exact sampling
#     covariances, CR2 cluster-robust inference clustered on study, robust CIs/p only where df >= df_min), so non-independence
#     of records from the same study is handled. regtest() on records would ignore it.
#   * Precision is sqrt(1 / n_eff), n_eff = 4 nT nC / (nT + nC), not the sampling SE. For lnRR and g the SE is mechanically related
#     to the effect size, which creates spurious funnel asymmetry (Pustejovsky & Rodgers 2019; Nakagawa et al. 2022).
#   * Two bias terms: small-study effect (Egger-type, se_eff) and time-lag / decline effect (publication year, per decade).
#   * PRIMARY (pure plastic) and SECONDARY (plastic + co-contaminant; contaminant type is kept in the model so the small-study test is
#     not confounded by the very different effects of the two types). lnRR first (primary scale), then g.
#   * Funnel plots: (a) RECORD level, precision sqrt(1 / n_eff) with the regression line from the multilevel test (matches the test; points
#     from one study cluster); (b) STUDY level, the classic funnel (study effect vs its SE with a grey pseudo-95% region). Panel b is for
#     display only: the lnRR / g SE mostly reflects replicate variability and between-study heterogeneity is far larger than sampling error,
#     so the funnel shape is not diagnostic and inference rests on the multilevel regression.
#   * A study-level Egger-type regression (study effect on mean sqrt(1/n_eff)) is reported as a sensitivity check. Trim-and-fill is
#     OFF by default (run_trimfill <- FALSE): it assumes independent studies, one true effect and SE-based precision, so it is unreliable
#     with this much heterogeneity (tau = 0.26 vs sampling SEs of 0.01-0.05) and invalid on a mixture of contaminant types.

stopifnot(exists("fit_cred"), exists("dat_pure_l"), exists("dat_sec_l"), exists("study_est"), exists("cls_tab"), exists("pal2"))

## 0. Precision and publication year ---------------------------------------------------------------------------
pb_prep <- function(d) {
  d %>% mutate(n_eff  = 4 * treat_reps * control_reps / (treat_reps + control_reps),
               se_eff = sqrt(1 / n_eff),
               year   = as.integer(stringr::str_extract(reference, "(19|20)[0-9]{2}"))) %>%
    mutate(year_c = (year - mean(year, na.rm = TRUE)) / 10)         # years centred, expressed per decade
}
pb_sets <- list("Primary lnRR"    = pb_prep(dat_pure_l), "Primary g"    = pb_prep(dat_pure),
                "Secondary lnRR"  = pb_prep(dat_sec_l),  "Secondary g"  = pb_prep(dat_sec))
if (anyNA(pb_sets[[1]]$year)) print(pb_sets[[1]] %>% filter(is.na(year)) %>% distinct(study_id, reference))   # fix any reference without a year
stopifnot(!anyNA(pb_sets[[1]]$year), !anyNA(pb_sets[[3]]$year))

tab_pb_cor <- purrr::imap_dfr(pb_sets, \(d, nm) tibble(dataset = nm, records = nrow(d), year_min = min(d$year), year_max = max(d$year),
                                                       r_se_year = cor(d$se_eff, d$year), r_se_effect = cor(d$se_eff, d$yi),
                                                       r_year_effect = cor(d$year, d$yi)))
print(tab_pb_cor %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

## 1. Multilevel bias tests (record level) ---------------------------------------------------------------------
# A negative se_eff slope on lnRR / g means that less precise (small) studies report LARGER reductions in hatching, which is the
# pattern publication bias towards significant reductions would produce; a positive slope means the opposite. The intercept is the
# effect extrapolated to se_eff = 0 (a very large study), i.e. the effect adjusted for small-study effects. year_c is the change
# in effect per decade of publication (a negative slope = later studies report larger reductions).
pb_terms <- list("Reference (no bias terms)" = character(0),
                 "Small-study effect"        = "se_eff",
                 "Time-lag effect"           = "year_c",
                 "Small-study + time-lag"    = c("se_eff", "year_c"))

pb_one <- function(d, label, terms, secondary) {
  tm <- c(if (secondary) "contam", terms)
  f  <- if (length(tm) == 0) ~ 1 else reformulate(tm)
  m  <- fit_cred(d, "simple", mods = f); rb <- rob_safe(m, d$study_id)
  inf_tbl(m, rb) %>%
    mutate(model = label, scale = sc_lab(d), records = m$k, studies = n_distinct(d$study_id),
           resid = resid_of(m), infer = attr(rb, "infer"))
}
tab_pb_models <- purrr::imap_dfr(pb_sets, \(d, nm)
                                 purrr::imap_dfr(pb_terms, \(tm, lab) pb_one(d, lab, tm, startsWith(nm, "Secondary"))) %>%
                                   mutate(dataset = sub(" .*", "", nm), .before = 1))

# also adjusting for taxonomic class (classes with >= 3 studies), as the earlier class-level bias check did
pb_class <- function(d, nm) {
  cs <- d %>% group_by(taxa_class) %>% summarise(k = n_distinct(study_id), .groups = "drop")
  dd <- d %>% filter(taxa_class %in% cs$taxa_class[cs$k >= 3]) %>% mutate(taxa_class = droplevels(factor(as.character(taxa_class))))
  tm <- c(if (startsWith(nm, "Secondary")) "contam", "se_eff", "taxa_class")
  m  <- fit_cred(dd, "simple", mods = reformulate(tm)); rb <- rob_safe(m, dd$study_id)
  inf_tbl(m, rb) %>%
    mutate(model = "Small-study effect + taxonomic class", scale = sc_lab(dd), records = m$k, studies = n_distinct(dd$study_id),
           resid = resid_of(m), infer = attr(rb, "infer"), dataset = sub(" .*", "", nm), .before = 1)
}
tab_pb_class <- purrr::imap_dfr(pb_sets, \(d, nm) purrr::possibly(pb_class, otherwise = NULL)(d, nm))
tab_pb_models <- bind_rows(tab_pb_models, tab_pb_class)

# key rows: the two bias slopes, and the adjusted (se_eff = 0) effect
tab_pb_key <- tab_pb_models %>%
  filter(term %in% c("se_eff", "year_c") | (term == "intrcpt" & model %in% c("Small-study effect", "Small-study + time-lag"))) %>%
  mutate(pct = if_else(scale == "lnRR", 100 * (exp(est) - 1), NA_real_),
         pct_lb = if_else(scale == "lnRR", 100 * (exp(ci_lb) - 1), NA_real_),
         pct_ub = if_else(scale == "lnRR", 100 * (exp(ci_ub) - 1), NA_real_)) %>%
  select(dataset, scale, model, term, records, studies, est, ci_lb, ci_ub, df, p_robust, p_model, infer, pct, pct_lb, pct_ub) %>%
  arrange(dataset, factor(scale, levels = scale_levels), model)
print(tab_pb_key %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = Inf, width = Inf)

## 2. Study-level sensitivity (one point per study) -------------------------------------------------------------
run_trimfill <- FALSE      # set TRUE to add trim-and-fill (study level, SE-based, primary only; see the header note before using it)

pb_study_est <- function(d, by_contam = FALSE) {
  g   <- if (by_contam) c("study_id", "contam") else "study_id"
  est <- d %>% group_by(across(all_of(g))) %>% group_modify(~ study_est(.x)) %>% ungroup()
  pr  <- d %>% group_by(across(all_of(g))) %>% summarise(se_eff = mean(se_eff), .groups = "drop")      # study-mean sqrt(1 / n_eff)
  est %>% left_join(pr, by = g) %>% left_join(select(cls_tab, study_id, class), by = "study_id")
}

pb_study_stats <- function(d) {
  st <- pb_study_est(d)
  r0 <- metafor::rma(yi = est, sei = se, data = st, method = "REML")                                       # pooled study-level mean
  r1 <- metafor::rma(yi = est, sei = se, mods = ~ se_eff, data = st, method = "REML", test = "knha")      # Egger-type: study effect on mean sqrt(1/n_eff)
  out <- tibble(scale = sc_lab(d), studies = r0$k, pooled = as.numeric(r0$b),
                slope = as.numeric(r1$beta[2]), slope_lb = r1$ci.lb[2], slope_ub = r1$ci.ub[2],
                slope_df = r1$ddf, slope_p = r1$pval[2])
  if (run_trimfill) {
    tf <- tryCatch(metafor::trimfill(r0), error = function(e) NULL)
    out <- out %>% mutate(tf_k0 = if (!is.null(tf)) tf$k0 else NA_integer_, tf_side = if (!is.null(tf)) tf$side else NA_character_,
                          tf_adjusted = if (!is.null(tf)) as.numeric(tf$b) else NA_real_,
                          tf_lb = if (!is.null(tf)) tf$ci.lb else NA_real_, tf_ub = if (!is.null(tf)) tf$ci.ub else NA_real_)
  }
  out
}
tab_pb_study <- tryCatch(
  bind_rows(pb_study_stats(pb_sets[["Primary lnRR"]]), pb_study_stats(pb_sets[["Primary g"]])) %>%
    mutate(pct_pooled = if_else(scale == "lnRR", 100 * (exp(pooled) - 1), NA_real_)),
  error = function(e) { message("[pubbias] study-level statistics failed: ", conditionMessage(e)); tibble() })
print(tab_pb_study %>% mutate(across(where(is.numeric), ~ round(.x, 3))), width = Inf)

## 3. Funnel plots ---------------------------------------------------------------------------------------------
# y axis: precision sqrt(1 / n_eff) (study mean at the study level), reversed so the most precise records are at the top and the axis
# spans only the observed range; a small vertical jitter separates records that share the same replicate numbers.
x_scale <- \(sc) if (sc == "lnRR") pct_scale else NULL
x_lab   <- \(sc) if (sc == "lnRR") "Change in hatching relative to control (log scale)" else "Hedges' g (negative = reduced hatching)"
y_lab   <- \(txt) if (grepl("study", txt)) expression("Precision (study mean): " * sqrt(1 / n[eff])) else expression("Precision: " * sqrt(1 / n[eff]))
y_pad   <- 0.02

pb_record_fig <- function(d, nm) {
  secondary <- startsWith(nm, "Secondary"); sc <- sc_lab(d)
  m <- fit_cred(d, "simple", mods = reformulate(c(if (secondary) "contam", "se_eff"))); b <- coef(m)
  grid <- tibble(se_eff = seq(min(d$se_eff), max(d$se_eff), length.out = 50))        # observed range only
  d <- d %>% mutate(col = if (secondary) as.character(contam) else as.character(taxa_class))
  p <- ggplot(d, aes(x = yi, y = se_eff)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = col), alpha = 0.7, size = 1.8, position = position_jitter(width = 0, height = 0.006, seed = 1))
  if (secondary) {
    line <- bind_rows(grid %>% mutate(grp = "Pure plastic",   est = b[["intrcpt"]] + b[["se_eff"]] * se_eff),
                      grid %>% mutate(grp = "Co-contaminant", est = b[["intrcpt"]] + b[["contamCo-contaminant"]] + b[["se_eff"]] * se_eff))
    p <- p + geom_path(data = line, aes(x = est, y = se_eff, colour = grp), linewidth = 0.9, inherit.aes = FALSE) +
      scale_colour_manual(values = pal2, breaks = c("Pure plastic", "Co-contaminant"), name = NULL)
  } else {
    line <- grid %>% mutate(est = b[["intrcpt"]] + b[["se_eff"]] * se_eff)
    p <- p + geom_path(data = line, aes(x = est, y = se_eff), colour = "black", linewidth = 0.9, inherit.aes = FALSE) +
      scale_colour_manual(values = class_pal, name = "Taxonomic class")
  }
  p <- p + scale_y_reverse(limits = c(max(d$se_eff) + y_pad, min(d$se_eff) - y_pad)) +
    labs(x = x_lab(sc), y = y_lab("Precision: ")) + theme_classic(base_size = 10)
  if (!is.null(x_scale(sc))) p <- p + x_scale(sc)
  p
}

pb_study_fig <- function(d, nm) {          # classic SE-based funnel with a grey pseudo-95% region (study level, one point per study)
  secondary <- startsWith(nm, "Secondary"); sc <- sc_lab(d)
  if (!secondary) {
    st  <- pb_study_est(d); res <- metafor::rma(yi = est, sei = se, data = st, method = "REML")
    mu  <- as.numeric(res$b); smax <- max(st$se) * 1.05
    tri <- tibble(x = c(mu, mu - 1.96 * smax, mu + 1.96 * smax), y = c(0, smax, smax))
    p <- ggplot(st, aes(x = est, y = se)) +
      geom_polygon(data = tri, aes(x = x, y = y), inherit.aes = FALSE, fill = "grey93", colour = "grey70") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_vline(xintercept = mu, colour = "black") +
      geom_point(aes(colour = class), size = 2.2) +
      scale_colour_manual(values = class_pal, name = "Taxonomic class")
    if (run_trimfill) {
      tf <- tryCatch(metafor::trimfill(res), error = function(e) NULL)
      if (!is.null(tf) && tf$k0 > 0) {
        fill <- tibble(est = tf$yi[tf$fill], se = sqrt(tf$vi[tf$fill]))
        p <- p + geom_point(data = fill, aes(x = est, y = se), shape = 2, size = 2.4, inherit.aes = FALSE) +
          geom_vline(xintercept = as.numeric(tf$b), linetype = "dotted")
      }
    }
  } else {
    st  <- pb_study_est(d, by_contam = TRUE)
    mus <- st %>% group_by(contam) %>% group_modify(~ tibble(mu = { dd <- .x; as.numeric(metafor::rma(yi = est, sei = se, data = dd, method = "REML")$b) }))
    p <- ggplot(st, aes(x = est, y = se, colour = as.character(contam))) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_vline(data = mus, aes(xintercept = mu, colour = as.character(contam)), inherit.aes = FALSE) +
      geom_point(size = 2.2) +
      scale_colour_manual(values = pal2, breaks = c("Pure plastic", "Co-contaminant"), name = NULL)
  }
  p <- p + scale_y_reverse() + labs(x = x_lab(sc), y = "Standard error") + theme_classic(base_size = 10)
  if (!is.null(x_scale(sc))) p <- p + x_scale(sc)
  p
}

for (nm in names(pb_sets)) {
  d <- pb_sets[[nm]]
  message("[pubbias] drawing: ", nm)
  pa <- tryCatch(pb_record_fig(d, nm), error = function(e) { message("[pubbias] RECORD panel failed (", nm, "): ", conditionMessage(e)); NULL })
  pb <- tryCatch(pb_study_fig(d, nm),  error = function(e) { message("[pubbias] STUDY panel failed (", nm, "): ", conditionMessage(e)); NULL })
  panels <- Filter(Negate(is.null), list(pa, pb))
  if (length(panels) == 0) next
  fig <- (patchwork::wrap_plots(panels, nrow = 1) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")) &
    theme(legend.position = "bottom")
  tag <- tolower(gsub("[^A-Za-z]+", "_", nm))
  print(fig)
  ggsave(file.path(out_dir, paste0("pubbias_", tag, ".pdf")), fig, width = 5 * length(panels), height = 5)
  ggsave(file.path(out_dir, paste0("pubbias_", tag, ".png")), fig, width = 5 * length(panels), height = 5, dpi = 150)
}

## 4. Export ---------------------------------------------------------------------------------------------------
write_csv(tab_pb_cor,    file.path(out_dir, "pubbias_correlations.csv"))
write_csv(tab_pb_models, file.path(out_dir, "pubbias_models.csv"))
write_csv(tab_pb_key,    file.path(out_dir, "pubbias_key.csv"))
write_csv(tab_pb_study,  file.path(out_dir, "pubbias_study_level.csv"))