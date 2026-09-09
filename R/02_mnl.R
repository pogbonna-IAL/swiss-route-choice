# 02_mnl.R -- baseline multinomial logit
# Swiss route choice
#
# Estimates the baseline MNL on the Swiss route choice data: two unlabelled
# alternatives, four generic coefficients, no alternative-specific constant.
#
# The design justifies each of those restrictions (see 01_data_audit.R):
#   * alternatives are unlabelled and symmetric on every attribute, and the
#     choice split is 49.7/50.3, so an ASC would only pick up position bias
#     that the data says is absent;
#   * coefficients are generic because the same attribute means the same thing
#     in either position.
#
# Three things are checked rather than assumed:
#   1. LL at zero coefficients equals -N log 2.
#   2. An independent implementation (mlogit) reproduces the estimates. Two
#      libraries agreeing on the same likelihood rules out a mis-specified
#      utility function, which an internally consistent Apollo run cannot.
#   3. The model predicts respondents it never saw. In-sample fit statistics
#      cannot distinguish a model that generalises from one that does not.
#
# Outputs
#   outputs/models/Swiss_MNL_model.rds       model object
#   outputs/models/Swiss_MNL_output.txt      formatted output
#   outputs/tables/02_mnl_estimates.csv      tidy estimates
#   outputs/tables/02_mnl_fit.csv            fit statistics
#   outputs/tables/02_mnl_valuations.csv     willingness-to-pay measures
#   outputs/tables/02_mnl_crosscheck.csv     apollo vs mlogit
#   outputs/tables/02_mnl_validation.csv     hold-out prediction
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))
source(here::here("R", "model_helpers.R"))
library(mlogit)   # independent cross-check only; not needed anywhere else

MODEL_NAME <- "Swiss_MNL"
MODEL_RDS  <- file.path(PATH_MODELS, paste0(MODEL_NAME, "_model.rds"))

# Clears Apollo internal state (a stale apollo_inputs or an attached database
# from an earlier model), not the workspace.
apollo_initialise()

# --- Data ------------------------------------------------------------------
data("apollo_swissRouteChoiceData", package = "apollo")
database <- apollo_swissRouteChoiceData

# --- Control ---------------------------------------------------------------
# panelData and seed are set explicitly rather than inherited: Apollo would
# infer panelData = TRUE from the repeated IDs anyway, but the later mixed
# logit and latent class scripts depend on both, so they are stated here.
apollo_control <- list(
  modelName       = MODEL_NAME,
  modelDescr      = "Baseline MNL Swiss route choice",
  indivID         = "ID",
  outputDirectory = PATH_MODELS,
  panelData       = TRUE,
  seed            = SEED,
  nCores          = 1
)

# --- Parameters ------------------------------------------------------------
# Zero starts are safe here: the MNL log-likelihood is globally concave, so
# there is a single maximum and no dependence on the starting point. That
# stops holding from 05_lc_2class.R onward.
apollo_beta <- c(
  b_tt = 0,
  b_tc = 0,
  b_hw = 0,
  b_ch = 0
)

# Nothing to fix: with no ASC and generic coefficients there is no scale or
# location normalisation to impose.
apollo_fixed <- c()

apollo_inputs <- apollo_validateInputs()

# --- Model definition ------------------------------------------------------
apollo_probabilities <- function(apollo_beta, apollo_inputs,
                                 functionality = "estimate") {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()
  V <- list()

  V[["alt1"]] <- b_tt * tt1 + b_tc * tc1 + b_hw * hw1 + b_ch * ch1
  V[["alt2"]] <- b_tt * tt2 + b_tc * tc2 + b_hw * hw2 + b_ch * ch2

  mnl_settings <- list(
    alternatives = c(alt1 = 1, alt2 = 2),
    avail        = 1,          # both routes available in every task
    choiceVar    = choice,
    utilities    = V
  )

  P[["model"]] <- apollo_mnl(mnl_settings, functionality)

  # Multiplies the 9 task probabilities within each respondent. It does not
  # change the MNL point estimates, which factorise over rows either way, but
  # it is the structure the mixed logit in 04 relies on.
  P <- apollo_panelProd(P, apollo_inputs, functionality)

  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# --- Estimation ------------------------------------------------------------
model_mnl <- apollo_estimate(apollo_beta, apollo_fixed,
                             apollo_probabilities, apollo_inputs)

apollo_modelOutput(model_mnl)

# saveRDS rather than apollo_saveOutput: the latter renames any existing file
# to _OLD1, _OLD2, ... on every re-run, so a few re-estimations leave the model
# directory full of stale copies that no script reads and nothing cleans up.
# 06 already avoided it; this script was the last one still generating them.
saveRDS(model_mnl, MODEL_RDS)
writeLines(capture.output(apollo_modelOutput(model_mnl)),
           file.path(PATH_MODELS, paste0(MODEL_NAME, "_output.txt")))

# LL at zero coefficients is -N * log(2) when both alternatives are equally
# likely. A departure here means the utilities were misspecified. LL0 comes
# back named by model component, so drop the name before comparing.
ll_equal_shares <- -nrow(database) * log(2)
ll_zero         <- as.numeric(model_mnl$LL0[["model"]])
stopifnot(isTRUE(all.equal(ll_zero, ll_equal_shares, tolerance = 1e-6)))
cat(sprintf("\nLL(0) check: %.4f == -N*log(2) = %.4f\n",
            ll_zero, ll_equal_shares))

# --- Tidy estimates --------------------------------------------------------
estimates <- tibble(
  parameter    = names(model_mnl$estimate),
  estimate     = as.numeric(model_mnl$estimate),
  se           = as.numeric(sqrt(diag(model_mnl$varcov))),
  rob_se       = as.numeric(sqrt(diag(model_mnl$robvarcov)))
) %>%
  mutate(
    t_ratio     = estimate / se,
    rob_t_ratio = estimate / rob_se,
    # Robust errors exceed classical ones by this factor because the MNL
    # treats a respondent's 9 tasks as independent when they are not.
    se_inflation = rob_se / se
  )

write_table(estimates, "mnl_estimates", prefix = "02")

cat("\n=== Estimates ==============================================\n")
show_table(estimates, digits = 5)

# --- Fit statistics --------------------------------------------------------
n_par <- length(model_mnl$estimate) - length(apollo_fixed)

fit <- tibble(
  statistic = c("observations", "individuals", "parameters",
                "LL(0)", "LL(final)", "rho2", "adj_rho2", "AIC", "BIC"),
  value = c(
    model_mnl$nObs,
    model_mnl$nIndivs,
    n_par,
    ll_zero,
    model_mnl$maximum,
    1 - model_mnl$maximum / ll_zero,
    1 - (model_mnl$maximum - n_par) / ll_zero,
    -2 * model_mnl$maximum + 2 * n_par,
    -2 * model_mnl$maximum + n_par * log(model_mnl$nObs)
  )
)

write_table(fit, "mnl_fit", prefix = "02")

cat("\n=== Fit ====================================================\n")
show_table(fit, digits = 4)

# --- Valuations ------------------------------------------------------------
# Ratios against the cost coefficient convert utility into money. Standard
# errors come from the delta method on the robust covariance matrix, which
# matters because a ratio of two estimates is not normally distributed and its
# uncertainty is far wider than either coefficient suggests on its own.
valuation_expr <- c(
  vtt_chf_per_hour       = "b_tt/b_tc*60",
  headway_chf_per_hour   = "b_hw/b_tc*60",
  interchange_chf        = "b_ch/b_tc",
  interchange_minutes    = "b_ch/b_tt",
  headway_vs_travel_time = "b_hw/b_tt"
)

cat("\n=== Valuations =============================================\n")
valuations <- apollo_deltaMethod(
  model_mnl,
  list(expression = valuation_expr)
)

valuations <- as_tibble(valuations) %>%
  rename(measure = Expression, value = Value, se = `s.e.`,
         t_ratio = `t-ratio (0)`) %>%
  mutate(
    expression = unname(valuation_expr[measure]),
    unit = c("CHF per hour", "CHF per hour", "CHF per interchange",
             "minutes per interchange", "ratio")[match(measure, names(valuation_expr))],
    ci_low  = value - 1.96 * se,
    ci_high = value + 1.96 * se,
    .after = measure
  )

write_table(valuations, "mnl_valuations", prefix = "02")
show_table(valuations, digits = 4)

# --- Independent cross-check against mlogit --------------------------------
# Apollo can only tell us that its own optimiser converged. Re-estimating the
# identical specification in a package written by different people, from a
# differently shaped dataset, is what actually rules out a coding error in the
# utility function: a transposed alternative or a dropped attribute would give
# a self-consistent Apollo run and a mismatch here.
#
# mlogit wants one row per alternative, so the wide design matrix is melted.
# "| 0" suppresses the alternative-specific constants, matching apollo_beta.
mlogit_long <- bind_rows(
  database %>% transmute(obs = row_number(), alt = 1L, chosen = choice == 1L,
                         tt = tt1, tc = tc1, hw = hw1, ch = ch1),
  database %>% transmute(obs = row_number(), alt = 2L, chosen = choice == 2L,
                         tt = tt2, tc = tc2, hw = hw2, ch = ch2)
) %>%
  arrange(obs, alt)

mlogit_idx <- dfidx::dfidx(mlogit_long, idx = list("obs", "alt"),
                           choice = "chosen")
model_mlogit <- mlogit(chosen ~ tt + tc + hw + ch | 0, data = mlogit_idx)

crosscheck <- tibble(
  parameter      = names(model_mnl$estimate),
  apollo         = as.numeric(model_mnl$estimate),
  mlogit         = as.numeric(coef(model_mlogit)[ROUTE_ATTRS]),
  apollo_se      = as.numeric(sqrt(diag(model_mnl$varcov))),
  mlogit_se      = as.numeric(sqrt(diag(vcov(model_mlogit)))[ROUTE_ATTRS])
) %>%
  mutate(abs_diff     = abs(apollo - mlogit),
         rel_diff     = abs_diff / abs(apollo),
         abs_diff_se  = abs(apollo_se - mlogit_se))

ll_diff <- abs(model_mnl$maximum - as.numeric(logLik(model_mlogit)))

write_table(
  crosscheck %>% bind_rows(tibble(parameter = "LL",
                                  apollo = model_mnl$maximum,
                                  mlogit = as.numeric(logLik(model_mlogit)),
                                  abs_diff = ll_diff)),
  "mnl_crosscheck", prefix = "02")

cat("\n=== Cross-check: apollo vs mlogit ==========================\n")
show_table(crosscheck, digits = 6)
cat(sprintf("LL: apollo %.6f | mlogit %.6f | difference %.2e\n",
            model_mnl$maximum, as.numeric(logLik(model_mlogit)), ll_diff))

# Both are maximising the same concave likelihood, so they must agree to
# optimiser tolerance. Anything larger is a specification difference, not
# numerical noise, and should stop the pipeline rather than be written to a
# table nobody reads.
stopifnot(ll_diff < 1e-4, max(crosscheck$rel_diff) < 1e-3)
cat("cross-check passed\n")

# mnl_dv, the hold-out split and the scoring all come from
# model_helpers.R, shared with 03 and 04.

split <- holdout_split(database)
train <- split$train
test  <- split$test

# Re-estimate on the training half only. Nothing from the hold-out touches
# these coefficients.
apollo_control$modelName  <- "Swiss_MNL_train"
apollo_control$modelDescr <- "Baseline MNL, training half"
database      <- train
apollo_inputs <- apollo_validateInputs(silent = TRUE)
model_train   <- apollo_estimate(
  apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
  estimate_settings = list(silent = TRUE, writeIter = FALSE))
database <- apollo_swissRouteChoiceData   # restore

b_train <- model_train$estimate

validation <- bind_rows(
  binary_logit_score(mnl_dv(b_train, train), train, "training", "MNL"),
  binary_logit_score(mnl_dv(b_train, test),  test,  "hold-out", "MNL"),
  # The full-sample model scored on the full sample, as the in-sample
  # reference the hold-out numbers should be read against.
  binary_logit_score(mnl_dv(model_mnl$estimate, database), database,
                     "full (in-sample)", "MNL")
)

write_table(validation, "mnl_validation", prefix = "02")

cat("
=== Hold-out validation (split on respondent) ==============
")
show_table(validation, digits = 4)

cat(sprintf(
  "
hold-out vs training LL per observation: %.4f vs %.4f (gap %.4f)
",
  validation$LL_per_obs[2], validation$LL_per_obs[1],
  validation$LL_per_obs[1] - validation$LL_per_obs[2]))

cat(sprintf("\nModel object: %s\n", MODEL_RDS))
cat(sprintf("Tables written to: %s\n", PATH_TABLES))
