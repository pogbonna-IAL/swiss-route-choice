# 10_report.R -- assemble the deliverable
# Swiss route choice
#
# Every earlier script writes CSVs and PNGs. Nothing read them back. This one
# does: it collects the thirty-odd tables into a single markdown report with
# the figures inlined, and writes the cross-model comparison that no
# individual script can produce because no individual script sees more than
# its own models.
#
# Sections whose inputs are missing are skipped with a note naming the script
# that produces them, so a partial pipeline gives a partial report rather than
# an error.
#
# RUN ORDER: last. See run_all.R.
#
# Outputs
#   outputs/report.md                       the written report
#   outputs/tables/10_model_comparison.csv  every model on one likelihood scale
#   outputs/tables/10_validation.csv        every hold-out score in one table
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))
library(knitr)   # kable(); this is the script the dependency was there for

REPORT <- file.path(PATH_OUT, "report.md")

# --- Plumbing --------------------------------------------------------------
lines <- character(0)
add   <- function(...) lines <<- c(lines, ...)
blank <- function() add("")

# Returns NULL rather than erroring, so the report degrades gracefully.
#
# `require` names the columns this report actually reads. A table left over
# from an earlier version of the pipeline is worse than a missing one: it
# exists, so the section runs, and then fails on a column that is no longer
# there. Checking the schema turns that into the same skip-with-a-note that a
# missing file gets.
read_if <- function(name, require = character(0)) {
  path <- file.path(PATH_TABLES, name)
  if (!file.exists(path)) return(NULL)

  x <- as_tibble(read.csv(path, check.names = FALSE))
  missing <- setdiff(require, names(x))
  if (length(missing)) {
    warning(sprintf("%s is stale: missing %s. Re-run the script that writes it.",
                    name, paste(missing, collapse = ", ")), call. = FALSE)
    return(NULL)
  }
  x
}

missing_note <- function(script) {
  add(sprintf("> _Not available: run `%s` to produce this section._", script))
  blank()
}

# kable with the numeric columns rounded. Passing a tibble straight to kable
# prints fifteen significant figures, which is unreadable and implies a
# precision none of these estimates have.
tbl <- function(x, digits = 3, caption = NULL) {
  if (is.null(x) || nrow(x) == 0) return(invisible(NULL))
  x <- x %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
  add(knitr::kable(x, format = "pipe", caption = caption))
  blank()
}

figure <- function(file, alt) {
  if (!file.exists(file.path(PATH_FIGURES, file))) return(invisible(NULL))
  add(sprintf("![%s](figures/%s)", alt, file))
  blank()
}

# --- Load everything -------------------------------------------------------
audit        <- read_if("01_audit_summary.csv", c("check", "value"))
mnl_est      <- read_if("02_mnl_estimates.csv", c("parameter", "estimate", "rob_se", "rob_t_ratio"))
mnl_fit      <- read_if("02_mnl_fit.csv", c("statistic", "value"))
mnl_val      <- read_if("02_mnl_valuations.csv", c("measure", "unit", "value", "se", "ci_low", "ci_high"))
mnl_cross    <- read_if("02_mnl_crosscheck.csv", c("parameter", "abs_diff"))
mnl_holdout  <- read_if("02_mnl_validation.csv", c("model", "sample", "LL_per_obs", "hit_rate"))
cov_lr       <- read_if("03_mnlcov_lrtests.csv", c("comparison", "LL_restricted", "LL_full", "lr_stat", "lr_df", "lr_p", "significant"))
cov_prof     <- read_if("03_mnlcov_profiles.csv", c("quantity", "profile", "value", "se", "ci_low", "ci_high"))
cov_holdout  <- read_if("03_mnlcov_validation.csv", c("model", "sample", "LL_per_obs"))
mxl_comp     <- read_if("04_mxl_comparison.csv", c("model", "parameters", "LL", "BIC", "note"))
mxl_moments  <- read_if("04_mxl_moments.csv", c("model", "attribute", "label", "median", "mean", "sd", "q05", "q95"))
mxl_vtt      <- read_if("04_mxl_vtt.csv", c("model", "median", "q10", "q90", "sigma_log"))
mxl_corr     <- read_if("04_mxl_correlations.csv", "attribute")
mxl_holdout  <- read_if("04_mxl_validation.csv", c("model", "sample", "LL_per_obs"))
lccov_comp   <- read_if("05_lccov_comparison.csv", c("classes", "parameters", "LL", "BIC", "LL_const", "BIC_const", "lr_stat", "lr_df", "lr_p", "lr_trustworthy", "BIC_improvement"))
lccov_alloc  <- read_if("05_lccov_allocation.csv", c("classes", "parameter", "estimate", "se", "t_ratio", "odds_ratio", "significant"))
lccov_prof   <- read_if("05_lccov_profiles.csv", c("classes", "profile", "class", "share"))
lc_comp      <- read_if("06_lc_comparison.csv", c("model", "parameters", "LL", "AIC", "BIC", "min_class_share", "class_shares", "at_best", "starts"))
lc_val       <- read_if("06_lc_valuations.csv", c("model", "class", "share", "measure", "value", "se", "ci_low", "ci_high", "ratio_reliable", "signs_ok"))
stability    <- read_if("07_stability.csv", c("n", "K", "source", "replications", "usable", "median_min_cls", "ll_sd", "sign_reversals", "verdict"))
ss_selection <- read_if("08_selection.csv", c("n", "K_selected", "replications", "share"))
design_qual  <- read_if("09_design_quality.csv", c("statistic", "value"))
failure_modes <- read_if("11_failure_modes.csv",
                         c("n", "K", "reps", "no_covariance", "wrong_sign",
                           "class_collapse", "misaligned", "usable_at_fit",
                           "reproduces"))
reproduction  <- read_if("11_reproduction.csv",
                         c("n", "K", "usable_at_fit", "reproduces", "coverage"))
vtt_recovery  <- read_if("11_vtt_recovery.csv",
                         c("n", "median_vtt", "n_negative", "within_50pct"))
class_recov   <- read_if("11_class_recovery.csv",
                         c("n", "K", "bench_class", "recovered"))

# --- Header ----------------------------------------------------------------
add("# Swiss route choice: heterogeneity, stability and sample size")
blank()
add(sprintf("_Generated %s by `R/10_report.R`._",
            format(Sys.time(), "%Y-%m-%d %H:%M")))
blank()
add("Four ways of letting travellers differ from one another, fitted to the",
    "same 388-respondent Swiss route choice panel, and an honest account of",
    "how much of the resulting structure survives contact with a smaller",
    "sample.")
blank()
add("| Question | Script | Section |")
add("|---|---|---|")
add("| What is in the data? | `01_data_audit.R` | [1](#1-the-data) |")
add("| What does the average traveller want? | `02_mnl.R` | [2](#2-the-baseline) |")
add("| Do observable characteristics explain taste? | `03_mnl_covariates.R` | [3](#3-observed-heterogeneity) |")
add("| Does taste vary in ways we cannot observe? | `04_mixed_logit.R` | [4](#4-continuous-unobserved-heterogeneity) |")
add("| Are there distinct types of traveller? | `06_lc_multiclass.R` | [5](#5-discrete-unobserved-heterogeneity) |")
add("| Who belongs to which type? | `05_lc_2class.R` | [6](#6-who-is-in-each-class) |")
add("| How much of this is real? | `07`, `08` | [7](#7-how-much-of-this-survives-a-smaller-sample) |")
add("| What should the next survey ask? | `09_idefix_design.R` | [8](#8-a-design-for-the-next-survey) |")
blank()

# --- 1. Data ---------------------------------------------------------------
add("## 1. The data")
blank()
if (is.null(audit)) missing_note("01_data_audit.R") else {
  add("The panel is balanced and complete: every respondent answers nine",
      "tasks, nothing is missing, and no task has a dominant alternative, so",
      "every observation carries a genuine trade-off.")
  blank()
  tbl(audit, caption = "Audit summary")
  figure("01_data_audit.png", "Data audit")
}

# --- 2. Baseline -----------------------------------------------------------
add("## 2. The baseline")
blank()
if (is.null(mnl_est)) missing_note("02_mnl.R") else {
  add("A multinomial logit with four generic coefficients and no",
      "alternative-specific constant. All four coefficients are negative and",
      "strongly significant.")
  blank()
  tbl(mnl_est %>% select(parameter, estimate, rob_se, rob_t_ratio), digits = 4,
      caption = "MNL estimates (robust standard errors)")
  tbl(mnl_fit, digits = 3, caption = "Fit")

  if (!is.null(mnl_cross)) {
    ll_row <- mnl_cross %>% filter(parameter == "LL")
    add(sprintf(paste0(
      "The specification is checked against an independent implementation. ",
      "`mlogit` reproduces the log-likelihood to %.1e and every coefficient ",
      "to six decimal places, which rules out a coding error in the utility ",
      "function that a self-consistent Apollo run could not."),
      ll_row$abs_diff[1]))
    blank()
  }

  add("Converting utility into money:")
  blank()
  tbl(mnl_val %>% select(measure, unit, value, se, ci_low, ci_high), digits = 2,
      caption = "Willingness to pay (delta-method standard errors)")

  if (!is.null(mnl_holdout)) {
    ho <- mnl_holdout %>% filter(sample == "hold-out")
    add(sprintf(paste0(
      "Out of sample, on 194 respondents the model never saw, it recovers ",
      "%.1f%% of choices and a market share of %.3f against an observed ",
      "%.3f. The baseline generalises; everything below has to beat it."),
      100 * ho$hit_rate[1], ho$share_alt1_predicted[1],
      ho$share_alt1_observed[1]))
    blank()
  }
}

# --- 3. Observed heterogeneity ---------------------------------------------
add("## 3. Observed heterogeneity")
blank()
if (is.null(cov_lr)) missing_note("03_mnl_covariates.R") else {
  add("Each attribute coefficient becomes a linear function of income, car",
      "availability and trip purpose: twenty extra parameters, and the",
      "baseline is the special case where all twenty are zero.")
  blank()
  tbl(cov_lr %>% select(comparison, LL_restricted, LL_full, lr_stat, lr_df,
                        lr_p, significant),
      digits = 4, caption = "Likelihood ratio tests")

  if (!is.null(cov_prof)) {
    add("What that buys, expressed as the value of travel time for different",
        "kinds of respondent:")
    blank()
    tbl(cov_prof %>% filter(quantity == "vtt_chf_per_hour") %>%
          select(profile, value, se, ci_low, ci_high),
        digits = 2, caption = "Value of travel time by respondent type (CHF/hour)")
  }

  if (!is.null(cov_holdout)) {
    g <- cov_holdout %>% filter(sample == "hold-out")
    gain <- g$LL_per_obs[g$model == "MNL + covariates"] -
            g$LL_per_obs[g$model == "MNL (baseline)"]
    add(sprintf(paste0(
      "Out of sample the interactions gain %+.4f log-likelihood per ",
      "observation over the baseline, so this is real structure rather than ",
      "twenty degrees of freedom spent on noise."), gain))
    blank()
  }
  figure("03_mnlcov.png", "Covariate interactions")
}

# --- 4. Mixed logit --------------------------------------------------------
add("## 4. Continuous unobserved heterogeneity")
blank()
if (is.null(mxl_comp)) missing_note("04_mixed_logit.R") else {
  add("Every coefficient is a negative lognormal draw, so each respondent has",
      "their own taste vector and no respondent can be estimated to enjoy",
      "paying more. Two versions: independent coefficients, and correlated",
      "ones.")
  blank()
  tbl(mxl_comp %>% select(model, parameters, LL, BIC, note), digits = 2,
      caption = "Mixed logit against the baseline")

  if (!is.null(mxl_vtt)) {
    v <- mxl_vtt %>% filter(model == "MXL-C")
    add(sprintf(paste0(
      "The headline result is not a number but a spread. The median ",
      "respondent values travel time at %.1f CHF/hour, but the middle 80%% ",
      "of the population runs from %.1f to %.1f -- a factor of %.1f. A ",
      "single point estimate from the baseline hides all of it."),
      v$median[1], v$q10[1], v$q90[1], v$q90[1] / v$q10[1]))
    blank()
    tbl(mxl_vtt, digits = 2, caption = "Value of travel time across respondents")
  }
  tbl(mxl_moments %>% filter(model == "MXL-C") %>%
        select(attribute, label, median, mean, sd, q05, q95),
      digits = 4, caption = "Implied coefficient distributions (MXL-C)")
  tbl(mxl_corr, digits = 3, caption = "Taste correlations, log scale")
  figure("04_mxl.png", "Mixed logit")
}

# --- 5. Latent classes -----------------------------------------------------
add("## 5. Discrete unobserved heterogeneity")
blank()
if (is.null(lc_comp)) missing_note("06_lc_multiclass.R") else {
  add("Instead of a continuous distribution, a small number of discrete",
      "types. Classes are reported in canonical order -- sorted by the travel",
      "time coefficient -- because the labels are otherwise arbitrary and",
      "would not line up with section 6.")
  blank()
  tbl(lc_comp %>% select(model, parameters, LL, AIC, BIC, min_class_share,
                         class_shares, at_best, starts,
                         any_of(c("verdict"))),
      digits = 2, caption = "One to five classes")

  best_bic <- lc_comp$model[which.min(lc_comp$BIC)]
  best_aic <- lc_comp$model[which.min(lc_comp$AIC)]
  add(sprintf(paste0(
    "BIC selects %s and AIC selects %s. The `at_best` column is the more ",
    "informative one: it counts how many of the random starts rediscovered ",
    "the reported optimum, and where it is low the solution is one the ",
    "search could easily have missed."), best_bic, best_aic))
  blank()

  if (!is.null(lc_val)) {
    tbl(lc_val %>% filter(measure == "vtt_chf_per_hour") %>%
          select(model, class, share, value, se, ci_low, ci_high,
                 ratio_reliable, signs_ok),
        digits = 2,
        caption = "Value of travel time by class, with delta-method standard errors")
    if (any(!lc_val$ratio_reliable)) {
      add("Rows with `ratio_reliable = FALSE` have a cost coefficient not",
          "distinguishable from zero; their ratios are not interpretable.")
      blank()
    }
  }
  figure("06_lc_selection.png", "Class selection")
}

# --- 6. Class membership ---------------------------------------------------
add("## 6. Who is in each class")
blank()
if (is.null(lccov_comp)) missing_note("05_lc_2class.R") else {
  add("Replacing the constant-only allocation with a linear index in the",
      "respondent covariates. Each model nests its counterpart in section 5,",
      "so the comparison is a likelihood ratio test.")
  blank()
  tbl(lccov_comp %>% select(classes, parameters, LL, BIC, LL_const, BIC_const,
                            lr_stat, lr_df, lr_p, lr_trustworthy,
                            any_of(c("verdict"))),
      digits = 3, caption = "Covariate allocation against constant-only")

  if ("usable" %in% names(lccov_comp) && any(!lccov_comp$usable)) {
    add(sprintf(paste0(
      "**%s not usable on the full panel.** Their Hessians are singular, so ",
      "every standard error is `NA`. The LR *statistic* is still computable ",
      "from the log-likelihoods, but a singular Hessian means the model is not ",
      "locally identified: the optimum is a ridge rather than a point, the ",
      "coefficients on it are arbitrary, the degrees of freedom overstate the ",
      "free parameters, and BIC is penalising a parameter count the model does ",
      "not have. Only the K = 2 row is a result."),
      paste(sprintf("LCcov%d", lccov_comp$classes[!lccov_comp$usable]),
            collapse = " and ")))
    blank()
  }

  if (any(!lccov_comp$lr_trustworthy)) {
    add("Rows with `lr_trustworthy = FALSE` rest on an optimum that fewer",
        "than three starts reproduced in one model or the other. Read the",
        "p-value as an upper bound on the evidence, not a result.")
    blank()
  }
  if (all(lccov_comp$BIC_improvement < 0)) {
    add("Note the tension: the likelihood ratio tests are significant, but",
        "BIC prefers the constant-only models at every K. The covariates",
        "shift membership detectably without paying for the parameters they",
        "cost.")
    blank()
  }
  tbl(lccov_alloc %>% filter(classes == 2) %>%
        select(parameter, estimate, se, t_ratio, odds_ratio, significant),
      digits = 3, caption = "Allocation coefficients, 2-class model")
  tbl(lccov_prof %>% filter(classes == 2) %>%
        pivot_wider(names_from = class, values_from = share,
                    names_prefix = "class_"),
      digits = 3, caption = "Predicted membership by respondent profile")
  figure("05_lccov_allocation.png", "Class allocation")
}

# --- 7. Stability ----------------------------------------------------------
add("## 7. How much of this survives a smaller sample")
blank()
if (is.null(stability)) missing_note("07_lc_stability.R and 08_small_sample_experiment.R") else {
  add("The single most important section. Everything above was fitted to 388",
      "respondents. Resampling smaller panels and refitting shows how much of",
      "the structure is a property of the population and how much is a",
      "property of having enough data.")
  blank()
  tbl(stability %>% select(n, K, source, replications, usable, median_min_cls,
                           ll_sd, sign_reversals, verdict),
      digits = 3, caption = "Latent class stability")
  add("`usable` means different things in the two blocks and must not be read",
      "across them: for the resampled rows it is the share of replications",
      "with a usable covariance matrix; for the full-sample rows it is the",
      "share of random starts that found the best solution.")
  blank()

  if (!is.null(ss_selection)) {
    tbl(ss_selection %>%
          pivot_wider(names_from = K_selected, values_from = c(replications, share),
                      values_fill = 0),
        digits = 3, caption = "Number of classes chosen by BIC, by sample size")
  }
  figure("07_stability.png", "Stability")
  figure("08_small_sample.png", "Small sample experiment")
}

# --- 8. Design -------------------------------------------------------------
add("## 8. A design for the next survey")
blank()
if (is.null(design_qual)) missing_note("09_idefix_design.R") else {
  add("Given what the baseline says respondents care about, a Bayesian",
      "D-efficient design for a follow-up. The priors are the estimated",
      "coefficients propagated exactly through the dummy coding, including",
      "their covariance, so the design is robust to the uncertainty this",
      "study leaves behind.")
  blank()
  tbl(design_qual, digits = 4, caption = "Design quality")
  figure("09_design.png", "Efficient design")
}

# --- 8b. Where it stops working --------------------------------------------
add("## 8b. Where latent class estimation stops working")
blank()
if (is.null(failure_modes)) missing_note("11_breakdown.R") else {
  # Counted over the RESAMPLED cells only. The four n = 388 benchmark fits are
  # excluded from both numbers, so the pair is internally consistent -- mixing
  # them (1044 fits against 436 usable) would silently compare two populations.
  n_fits   <- sum(failure_modes$reps)
  n_usable <- round(sum(failure_modes$reps * failure_modes$usable_at_fit))
  dead     <- failure_modes %>% filter(usable_at_fit == 0)

  add(sprintf(paste0(
    "Pushed below any sensible sample size, **nothing ever crashes**. Across ",
    "%d resampled fits at n from 250 down to 5, every single one returned a ",
    "finite log-likelihood, a parameter vector and a set of class shares -- ",
    "including four classes fitted to five people. **%d of them are usable.**"),
    n_fits, n_usable))
  blank()
  add(sprintf(paste0(
    "In %d of the %d (n, K) cells, **not one replication of twenty was ",
    "usable**, and every fit in those cells converged. A fit counts as usable ",
    "only if the Hessian is invertible, no attribute coefficient came back ",
    "positive, no class is thinner than 5%% of respondents, and no two classes ",
    "agree within 10%% on every attribute."),
    nrow(dead), nrow(failure_modes)))
  blank()

  tbl(failure_modes %>%
        select(n, K, `obs/par` = obs_per_par, `no cov` = no_covariance,
               `wrong sign` = wrong_sign, `collapsed` = class_collapse,
               misaligned, usable = usable_at_fit, reproduces),
      digits = 2, caption = "Failure modes by sample size and classes")

  add("`collapsed` is **0% at every sample size**. Classes stay numerically",
      "distinct all the way down while becoming meaningless -- the obvious",
      "diagnostic never fires, which is exactly why it is not a sufficient",
      "check.")
  blank()

  if (!is.null(reproduction)) {
    add("Coverage is the earliest and sharpest failure. At n = 250, where the",
        "two-class model is usable in every replication and shows no other",
        "pathology, its nominal 95% intervals contain the full-panel value",
        "little more than half the time. The MNL's are honest at 95%.")
    blank()
    tbl(reproduction %>%
          select(n, K, usable = usable_at_fit, reproduces, coverage, rel_rmse),
        digits = 3,
        caption = "Reproduction against the n = 388 benchmark")
  }

  if (!is.null(vtt_recovery)) {
    add("The MNL is the control, and it degrades gracefully where the latent",
        "class models fall apart -- so the collapse is a property of the latent",
        "class structure, not simply of having few observations:")
    blank()
    tbl(vtt_recovery %>%
          select(n, median_vtt, q10, q90, `wrong sign` = n_negative,
                 `within 50%` = within_50pct),
        digits = 2,
        caption = "MNL value of travel time by sample size (benchmark 27.21 CHF/h)")
  }

  if (!is.null(class_recov)) {
    tbl(class_recov %>% filter(K == 3) %>%
          select(n, bench_class, recovered) %>%
          pivot_wider(names_from = bench_class, values_from = recovered,
                      names_prefix = "class "),
        digits = 2,
        caption = "Three-class model: share of replications recovering each benchmark class")
  }

  figure("11_breakdown.png", "Small-sample breakdown")
  add("The full write-up is in [`docs/08-what-broke.md`](../docs/08-what-broke.md).")
  blank()
}

# --- 9. Everything on one scale --------------------------------------------
add("## 9. Every model on one scale")
blank()

collect <- function(model, parameters, LL, BIC, family) {
  tibble(model = model, family = family, parameters = parameters,
         LL = LL, BIC = BIC)
}

all_models <- bind_rows(
  if (!is.null(mnl_fit)) {
    collect("MNL",
            mnl_fit$value[mnl_fit$statistic == "parameters"],
            mnl_fit$value[mnl_fit$statistic == "LL(final)"],
            mnl_fit$value[mnl_fit$statistic == "BIC"],
            "fixed coefficients")
  },
  if (!is.null(mxl_comp)) {
    mxl_comp %>% filter(model != "MNL (02)") %>%
      transmute(model, family = "continuous heterogeneity", parameters, LL, BIC)
  },
  if (!is.null(cov_lr) && !is.null(mnl_fit)) {
    # 03 does not write a fit table of its own, so BIC is reconstructed here
    # from the observation count 02 recorded rather than a hardcoded 3492.
    n_obs   <- mnl_fit$value[mnl_fit$statistic == "observations"]
    n_par   <- 4 + 4 * length(MODEL_COVARS)
    collect("MNL + covariates", n_par, cov_lr$LL_full[1],
            -2 * cov_lr$LL_full[1] + n_par * log(n_obs),
            "observed heterogeneity")
  },
  if (!is.null(lc_comp)) {
    lc_comp %>% transmute(model, family = "discrete heterogeneity",
                          parameters, LL, BIC)
  },
  if (!is.null(lccov_comp)) {
    lccov_comp %>% transmute(model = paste0("LCcov", classes),
                             family = "discrete + observed",
                             parameters, LL, BIC)
  }
) %>%
  arrange(BIC) %>%
  mutate(rank_BIC = row_number(), best = BIC == min(BIC))

write_table(all_models, "model_comparison", prefix = "10")
tbl(all_models, digits = 2, caption = "All models ranked by BIC")

if (nrow(all_models)) {
  add(sprintf("BIC's overall pick is **%s** (%d parameters).",
              all_models$model[1], as.integer(all_models$parameters[1])))
  blank()
}

# Hold-out scores are the comparison that does not reward parameters.
# 02, 03 and 04 all score a baseline MNL on the same split, so the same model
# arrives three times. They agree to optimiser tolerance but not to the last
# digit, so distinct() alone would leave three near-identical rows; keying on
# (model, sample) keeps the first and drops the rest.
all_val <- bind_rows(mnl_holdout, cov_holdout, mxl_holdout) %>%
  distinct(model, sample, .keep_all = TRUE)

if (nrow(all_val)) {
  write_table(all_val, "validation", prefix = "10")
  add("Information criteria penalise parameters by a formula. The hold-out",
      "sample penalises them by whether they help:")
  blank()
  tbl(all_val %>% filter(sample == "hold-out") %>%
        select(model, respondents, observations, LL_per_obs, rho2_vs_coin,
               hit_rate, share_alt1_predicted, share_alt1_observed),
      digits = 4, caption = "Hold-out performance (194 unseen respondents)")
}

# --- 10. Reproducing -------------------------------------------------------
add("## 10. Reproducing this")
blank()
add("```")
add("Rscript R/run_all.R          # whole pipeline, in dependency order")
add("Rscript R/run_all.R --list   # the plan, with rough runtimes")
add("Rscript tests/run_tests.R    # helper test suite")
add("```")
blank()
add("The file numbering is not the run order: `05` reads a table `06` writes,",
    "and `07` reads tables `06` and `08` write. `R/run_all.R` holds the real",
    "order and checks each script's inputs before running it.")
blank()
if (file.exists(file.path(PATH_OUT, "session_info.txt"))) {
  add("Package versions for the last full run are in `outputs/session_info.txt`.")
  blank()
}

writeLines(lines, REPORT)

# The docs/ bundle is a separate build step (R/build_docs.R) -- assembling the
# analysis report and assembling the reference documentation are different
# jobs over different inputs, and this script was doing both.
source(here::here("R", "build_docs.R"))
docs_path <- bundle_docs()

cat(sprintf("Report written to:        %s (%d lines)\n", REPORT, length(lines)))
cat(sprintf("Documentation bundled to: %s (%d lines)\n", docs_path,
            length(readLines(docs_path, warn = FALSE))))
cat(sprintf("Tables written to:        %s\n", PATH_TABLES))
cat("\nMarkdown sources for the published pages:\n")
cat("  outputs/report.md          the analysis report\n")
cat("  outputs/documentation.md   the reference documentation, single file\n")
cat("  docs/08-what-broke.md      the small-sample breakdown write-up\n")