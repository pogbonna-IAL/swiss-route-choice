# 05_lc_2class.R -- latent class with covariate-driven class allocation
# Swiss route choice
#
# 06_lc_multiclass.R answers "how many classes and how big are they". This
# script answers the question that follows: WHO is in each class. The
# constant-only allocation is replaced by a linear index in the respondent
# covariates, so class membership probability varies by person:
#
#   V[class 1] = 0                                    (reference)
#   V[class k] = delta_k + g_inc_k  * log_income
#                        + g_car_k  * car_availability
#                        + g_com_k  * commute
#                        + g_shop_k * shopping
#                        + g_bus_k  * business        (leisure = reference)
#
# Estimated for K = 2, 3 and 4, each nesting its constant-only counterpart from
# 06, so every comparison is a clean likelihood-ratio test on (K-1) x 5 df.
#
# Despite the file name this covers K = 2 to 4; the 2-class model is reported
# in most detail because it is the only one whose membership story is simple
# enough to read directly off the coefficients.
#
# Classes are reported in canonical order (sorted by b_tt, most negative
# first), the same rule 06 uses, so "class 1" means the same thing in both
# scripts. Because class 1's whole allocation index is normalised to zero, a
# re-ordering turns every allocation coefficient into a contrast against the
# new reference; lc_canonical_params computes their standard errors as such
# rather than just permuting them.
#
# RUN ORDER: after 06_lc_multiclass.R -- this script reads the comparison
# table 06 writes. The file numbering does not reflect the dependency; see
# run_all.R.
#
# Outputs
#   outputs/models/Swiss_LCcov<K>_model.rds
#   outputs/models/Swiss_LCcov<K>_output.txt
#   outputs/tables/05_lccov_comparison.csv    LR tests against constant-only
#   outputs/tables/05_lccov_parameters.csv    all estimates
#   outputs/tables/05_lccov_allocation.csv    allocation coefficients + odds
#   outputs/tables/05_lccov_profiles.csv      predicted shares by respondent type
#   outputs/tables/05_lccov_runs.csv          per-start log
#   outputs/figures/05_lccov_allocation.png
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))
source(here::here("R", "lc_helpers.R"))

K_SET    <- 2:4
N_STARTS <- 30L

# Completed models are reloaded from outputs/models/ instead of re-estimated.
# Set REFIT=1 in the environment to force a full re-run.
RESUME <- lc_resume_enabled()

apollo_initialise()
data("apollo_swissRouteChoiceData", package = "apollo")

# log_income is centred on its mean so delta_k is the allocation constant for a
# respondent of average income rather than for someone earning 1 CHF.
database <- apollo_swissRouteChoiceData %>%
  mutate(log_income = log(hh_inc_abs) - mean(log(hh_inc_abs)))

# --- Estimation -------------------------------------------------------------
fit_cov <- function(K) {
  cat(sprintf("\n\n=========== LC%d with allocation covariates ===========\n", K))

  lc_install_cov(K)
  par_names <- lc_par_names_cov(K)

  starts <- lapply(seq_len(N_STARTS), function(draw) {
    lc_start_uniform(par_names, SEED + draw)
  })

  # lc_estimate_best re-validates the inputs before the final covariance run.
  # This script previously re-estimated the winning start against whatever
  # apollo_inputs happened to be left in the global environment by the last
  # draw that validated, which is only safe by accident.
  #
  # Wrapped in lc_fit_cached so an interrupted sweep resumes at the next K
  # instead of starting over: three K values at thirty starts each is half an
  # hour, and losing it to a killed process is a poor trade for a cache file.
  fit <- lc_fit_cached(sprintf("Swiss_LCcov%d", K), resume = RESUME, fit_fn = function() {
    lc_estimate_best(
      starts       = starts,
      apollo_fixed = c(),                    # class 1 is fixed at 0 in code
      control = list(
        modelName       = sprintf("Swiss_LCcov%d", K),
        modelDescr      = sprintf("Latent class, %d classes, covariate allocation", K),
        indivID         = "ID",
        outputDirectory = PATH_MODELS,
        panelData       = TRUE,
        seed            = SEED,
        nCores          = 1
      )
    )
  })

  stopifnot(fit$ok)

  # Class shares vary by respondent here, so the degeneracy check uses the
  # share averaged over respondents rather than a single softmax. Averaged
  # over ROWS it would be the same number on this balanced panel but nine
  # times the work, and wrong the moment the panel is unbalanced.
  est <- fit$model$estimate
  resp <- database %>%
    group_by(ID) %>%
    summarise(across(all_of(names(LC_COVARS)), first), .groups = "drop")
  mean_shares <- rowMeans(vapply(seq_len(nrow(resp)), function(i) {
    lc_class_shares_cov(est, K, as.list(resp[i, names(LC_COVARS)]))
  }, numeric(K)))

  fit$verdict <- lc_fit_verdict(fit, K, shares = mean_shares)
  lc_announce_verdict(fit$verdict, sprintf("LCcov%d", K))

  apollo_modelOutput(fit$model)
  # 06 writes a formatted output file per model; this script did not, which
  # left the LCcov models the only ones with no readable record on disk.
  # (The .rds itself is written by lc_fit_cached.)
  writeLines(capture.output(apollo_modelOutput(fit$model)),
             file.path(PATH_MODELS, sprintf("Swiss_LCcov%d_output.txt", K)))

  fit
}

results <- lapply(K_SET, fit_cov)
names(results) <- paste0("K", K_SET)

write_table(
  map_dfr(seq_along(K_SET), function(i) {
    results[[i]]$run_log %>% mutate(classes = K_SET[i], .before = 1)
  }),
  "lccov_runs", prefix = "05")

# --- Comparison against the constant-only models ---------------------------
const_path <- file.path(PATH_TABLES, "06_lc_comparison.csv")
if (!file.exists(const_path)) {
  stop("Run 06_lc_multiclass.R first -- ", const_path, " not found.")
}
constant_only <- read.csv(const_path) %>%
  as_tibble() %>%
  select(classes, LL_const = LL, par_const = parameters, BIC_const = BIC,
         at_best_const = at_best, starts_const = starts)

comparison <- map_dfr(seq_along(K_SET), function(i) {
  K <- K_SET[i]
  fit <- results[[i]]
  m <- fit$model
  n_par <- lc_n_par_cov(K)
  ll <- m$maximum
  tibble(
    classes    = K,
    parameters = n_par,
    LL         = ll,
    AIC        = -2 * ll + 2 * n_par,
    BIC        = -2 * ll + n_par * log(m$nObs),
    at_best    = fit$n_at_best,
    starts     = length(fit$run_log$draw),
    usable     = fit$verdict$usable,
    verdict    = lc_verdict_label(fit$verdict),
    problems   = paste(fit$verdict$flags, collapse = ";")
  )
}) %>%
  left_join(constant_only, by = "classes") %>%
  mutate(
    # The constant-only model is this model with every gamma set to zero, so
    # twice the likelihood gap is chi-squared on the number of gammas.
    lr_stat = 2 * (LL - LL_const),
    lr_df   = parameters - par_const,
    lr_p    = pchisq(lr_stat, df = lr_df, lower.tail = FALSE),
    BIC_improvement = BIC_const - BIC,
    # A likelihood ratio test compares two GLOBAL maxima. If either model's
    # optimum was found by only a handful of starts, the test is being run
    # against a number the search may simply have failed to beat, and the
    # p-value is optimistic. Flag it rather than reporting it bare.
    lr_trustworthy = at_best >= 3L & at_best_const >= 3L
  )

write_table(comparison, "lccov_comparison", prefix = "05")

cat("\n\n=== Covariate allocation vs constant-only ==================\n")
print(as.data.frame(comparison %>%
  transmute(classes, parameters, LL = round(LL, 3), BIC = round(BIC, 2),
            LL_const = round(LL_const, 3), BIC_const = round(BIC_const, 2),
            lr_stat = round(lr_stat, 2), lr_df,
            lr_p = signif(lr_p, 3),
            BIC_better = BIC_improvement > 0,
            at_best = sprintf("%d/%d", at_best, starts),
            at_best_const = sprintf("%d/%d", at_best_const, starts_const),
            lr_ok = lr_trustworthy, verdict)), row.names = FALSE)

if (any(!comparison$usable)) {
  # Be precise about what a singular Hessian does and does not invalidate.
  # The LR STATISTIC uses only log-likelihoods, so it is still computable.
  # What fails is everything else: a singular Hessian means the model is not
  # locally identified, so the "optimum" is a ridge rather than a point, the
  # estimates are arbitrary along it, the degrees of freedom overstate the
  # free parameters, and the BIC penalty is computed on a parameter count
  # the model does not actually have.
  cat(sprintf(paste0(
    "\n!! NOT USABLE: %s.\n",
    "!! Every standard error in those models is NA. The LR statistic is still\n",
    "!! computable from the log-likelihoods, but a singular Hessian means the\n",
    "!! model is not locally identified: the estimates are arbitrary along a\n",
    "!! ridge, the df overstate the free parameters, and BIC is penalising a\n",
    "!! parameter count the model does not have. Do not report coefficients\n",
    "!! from these fits.\n"),
    paste(sprintf("LCcov%d (%s)", comparison$classes[!comparison$usable],
                  comparison$verdict[!comparison$usable]), collapse = ", ")))
}

if (any(!comparison$lr_trustworthy)) {
  cat("\nWARNING: rows with lr_ok = FALSE rest on an optimum that fewer than 3\n")
  cat("starts reproduced, in one model or the other. Treat the LR p-value as an\n")
  cat("upper bound on the evidence, not a result.\n")
}

# --- Parameters -------------------------------------------------------------
parameters <- map_dfr(seq_along(K_SET), function(i) {
  lc_canonical_params(results[[i]]$model, K_SET[i], covariate_alloc = TRUE) %>%
    mutate(classes = K_SET[i], usable = results[[i]]$verdict$usable, .before = 1)
})
write_table(parameters, "lccov_parameters", prefix = "05")

# Allocation coefficients are log-odds relative to class 1, so exponentiating
# gives the multiplicative effect on the odds of membership.
allocation <- parameters %>%
  filter(str_starts(parameter, "g_") | str_starts(parameter, "delta_")) %>%
  mutate(odds_ratio = ifelse(str_starts(parameter, "g_"), exp(estimate), NA_real_),
         significant = !is.na(t_ratio) & abs(t_ratio) > 1.96)
write_table(allocation, "lccov_allocation", prefix = "05")

cat("\n=== Allocation coefficients (log-odds vs class 1) ==========\n")
show_table(allocation %>% select(classes, parameter, estimate, se, t_ratio,
                                 odds_ratio, significant))

# --- Predicted membership profiles -----------------------------------------
# Class shares for archetypal respondents, which is what the allocation model
# is actually for: reading membership off respondent characteristics.
#
# The SD is taken over respondents, not over rows. The panel is balanced at 9
# tasks each so the two coincide here, but that is a property of this dataset
# and not something the code should depend on.
# Shared with 03 via 00_setup.R -- see the note there on why these two
# scripts must describe the same archetypes.
inc_sd   <- respondent_income_sd(database)
profiles <- respondent_profiles(inc_sd)

profile_tbl <- map_dfr(seq_along(K_SET), function(i) {
  K <- K_SET[i]
  # Canonical estimates so the class numbering matches the parameter table.
  est <- lc_canonical_cov(results[[i]]$model$estimate, K)
  map_dfr(names(profiles), function(p) {
    sh <- lc_class_shares_cov(est, K, profiles[[p]])
    tibble(classes = K, profile = p, class = seq_len(K), share = sh)
  })
})
write_table(profile_tbl, "lccov_profiles", prefix = "05")

cat("\n=== Predicted class shares by respondent profile ===========\n")
print(as.data.frame(profile_tbl %>%
        mutate(share = round(share, 3)) %>%
        pivot_wider(names_from = class, values_from = share,
                    names_prefix = "class_")), row.names = FALSE)

# --- Figure -----------------------------------------------------------------
p_alloc <- allocation %>%
  filter(str_starts(parameter, "g_")) %>%
  ggplot(aes(estimate, parameter, colour = significant)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se)) +
  facet_wrap(~ paste0(classes, " classes"), scales = "free_y") +
  labs(title = "Class allocation coefficients (log-odds vs class 1)",
       subtitle = "Intervals crossing zero mean the covariate does not shift membership",
       x = NULL, y = NULL, colour = "p < 0.05")

p_prof <- profile_tbl %>%
  filter(classes == 2) %>%
  ggplot(aes(share, profile, fill = factor(class))) +
  geom_col() +
  labs(title = "Predicted membership, 2-class model",
       x = "share", y = NULL, fill = "class")

fig_path <- write_figure(p_alloc / p_prof, "lccov_allocation", prefix = "05",
                         width = 13, height = 10)

cat(sprintf("\nTables written to: %s\n", PATH_TABLES))
cat(sprintf("Figure written to: %s\n", fig_path))
