# 03_mnl_covariates.R -- MNL with observed taste heterogeneity
# Swiss route choice
#
# The baseline MNL in 02 gives every respondent the same four coefficients.
# This script relaxes that in the cheapest possible way: each attribute
# coefficient becomes a linear function of the respondent's observed
# characteristics,
#
#   b_a,n = b_a + g_a_inc  * log_income_n     (centred)
#                + g_a_car  * car_availability_n
#                + g_a_com  * commute_n
#                + g_a_shop * shopping_n
#                + g_a_bus  * business_n      (leisure = reference)
#
# for a in {tt, tc, hw, ch}. That is 4 + 20 = 24 parameters against the
# baseline's 4, and the baseline is the special case with all twenty gammas at
# zero, so the comparison is a clean likelihood-ratio test on 20 df.
#
# This is the OBSERVED heterogeneity model. It can only pick up taste
# variation that lines up with a covariate the survey happened to collect.
# 04_mixed_logit.R takes the complementary approach and lets tastes vary in a
# way the covariates do not explain; 05/06 let them vary in discrete classes.
# Reading the three together is the point -- see 10_report.R.
#
# Restricted models are estimated by fixing the relevant gammas at zero with
# apollo_fixed rather than by writing out a second utility function, so every
# model in the comparison is provably the same specification under a
# constraint.
#
# RUN ORDER: after 02_mnl.R (reads the baseline model for the LR test).
#
# Outputs
#   outputs/models/Swiss_MNLcov_model.rds
#   outputs/models/Swiss_MNLcov_output.txt
#   outputs/tables/03_mnlcov_estimates.csv     all 24 estimates
#   outputs/tables/03_mnlcov_lrtests.csv       nested tests, overall and per attribute
#   outputs/tables/03_mnlcov_profiles.csv      implied coefficients and VTT by respondent type
#   outputs/tables/03_mnlcov_validation.csv    hold-out prediction vs the baseline MNL
#   outputs/figures/03_mnlcov.png
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))
source(here::here("R", "model_helpers.R"))

MODEL_NAME <- "Swiss_MNLcov"

apollo_initialise()

# --- Data ------------------------------------------------------------------
data("apollo_swissRouteChoiceData", package = "apollo")

# log_income is centred on its mean so b_a is the coefficient for a respondent
# of average income rather than for someone earning 1 CHF, which would put the
# base coefficient far outside the range of the data and make it
# uninterpretable.
database <- apollo_swissRouteChoiceData %>%
  mutate(log_income = log(hh_inc_abs) - mean(log(hh_inc_abs)))

mnl_path <- file.path(PATH_MODELS, "Swiss_MNL_model.rds")
if (!file.exists(mnl_path)) {
  stop("Run 02_mnl.R first -- ", mnl_path, " not found.")
}
model_base <- readRDS(mnl_path)

# --- Parameters ------------------------------------------------------------
gamma_names <- as.vector(t(outer(paste0("g_", ROUTE_ATTRS),
                                 unname(MODEL_COVARS), paste, sep = "_")))

apollo_beta <- c(
  setNames(as.numeric(model_base$estimate[paste0("b_", ROUTE_ATTRS)]),
           paste0("b_", ROUTE_ATTRS)),   # start at the baseline solution
  setNames(rep(0, length(gamma_names)), gamma_names)
)

apollo_fixed <- c()

apollo_control <- list(
  modelName       = MODEL_NAME,
  modelDescr      = "MNL with covariate-interacted attribute coefficients",
  indivID         = "ID",
  outputDirectory = PATH_MODELS,
  panelData       = TRUE,
  seed            = SEED,
  nCores          = 1
)

apollo_inputs <- apollo_validateInputs()

# --- Model definition ------------------------------------------------------
# Written out in full rather than generated: Apollo parses the source of this
# function, so the parameter names have to appear literally, and at four
# attributes by five covariates the explicit version is still readable.
apollo_probabilities <- function(apollo_beta, apollo_inputs,
                                 functionality = "estimate") {

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  b_tt_n <- b_tt + g_tt_inc * log_income + g_tt_car * car_availability +
                   g_tt_com * commute    + g_tt_shop * shopping +
                   g_tt_bus * business
  b_tc_n <- b_tc + g_tc_inc * log_income + g_tc_car * car_availability +
                   g_tc_com * commute    + g_tc_shop * shopping +
                   g_tc_bus * business
  b_hw_n <- b_hw + g_hw_inc * log_income + g_hw_car * car_availability +
                   g_hw_com * commute    + g_hw_shop * shopping +
                   g_hw_bus * business
  b_ch_n <- b_ch + g_ch_inc * log_income + g_ch_car * car_availability +
                   g_ch_com * commute    + g_ch_shop * shopping +
                   g_ch_bus * business

  P <- list()
  V <- list()

  V[["alt1"]] <- b_tt_n * tt1 + b_tc_n * tc1 + b_hw_n * hw1 + b_ch_n * ch1
  V[["alt2"]] <- b_tt_n * tt2 + b_tc_n * tc2 + b_hw_n * hw2 + b_ch_n * ch2

  mnl_settings <- list(
    alternatives = c(alt1 = 1, alt2 = 2),
    avail        = 1,
    choiceVar    = choice,
    utilities    = V
  )

  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# --- Estimation ------------------------------------------------------------
# One helper for every model in the comparison. `fix` names the gammas held at
# zero; fixing all twenty must reproduce the baseline MNL exactly, which is
# checked below rather than assumed.

# `fix` names the gammas held at zero. apollo_estimate needs its four
# arguments to be visible by name in the calling frame, not just passed as
# values, which is why everything is assigned into the global environment
# before the call rather than handed over as local variables.
estimate_variant <- function(fix, name) {
  beta <- apollo_beta
  beta[fix] <- 0

  assign("apollo_control", modifyList(apollo_control, list(modelName = name)),
         envir = globalenv())
  assign("apollo_beta",  beta, envir = globalenv())
  assign("apollo_fixed", fix,  envir = globalenv())
  assign("apollo_inputs", quietly(apollo_validateInputs(silent = TRUE)),
         envir = globalenv())

  quietly(apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities,
                          apollo_inputs,
                          estimate_settings = list(silent = TRUE,
                                                   writeIter = FALSE)))
}

cat("\n=== Estimating the full interaction model ==================\n")
model_cov <- apollo_estimate(apollo_beta, apollo_fixed,
                             apollo_probabilities, apollo_inputs)
apollo_modelOutput(model_cov)

saveRDS(scrub_local_paths(model_cov),
        file.path(PATH_MODELS, paste0(MODEL_NAME, "_model.rds")))
writeLines(capture.output(apollo_modelOutput(model_cov)),
           file.path(PATH_MODELS, paste0(MODEL_NAME, "_output.txt")))

# --- Nested tests ----------------------------------------------------------
# The fully restricted model. If fixing all twenty gammas at zero does not
# reproduce 02's log-likelihood, the interaction terms are not entering the
# utility the way this script claims they are, and every LR test below would
# be measuring the wrong thing.
cat("\n=== Restricted models ======================================\n")
model_null <- estimate_variant(gamma_names, "Swiss_MNLcov_null")
stopifnot(isTRUE(all.equal(model_null$maximum, model_base$maximum,
                           tolerance = 1e-6)))
cat(sprintf("baseline reproduced by fixing all gammas: LL = %.6f (02 gave %.6f)\n",
            model_null$maximum, model_base$maximum))

# Each attribute's five interactions, dropped as a block: does travel time
# sensitivity vary with who the respondent is, holding the other three
# attributes' interactions in the model?
per_attribute <- map_dfr(ROUTE_ATTRS, function(a) {
  fix <- grep(paste0("^g_", a, "_"), gamma_names, value = TRUE)
  m   <- estimate_variant(fix, paste0("Swiss_MNLcov_no_", a))
  lr  <- 2 * (model_cov$maximum - m$maximum)
  tibble(comparison  = paste0("no ", a, " interactions"),
         description = sprintf("%s coefficient held constant across respondents",
                               ROUTE_ATTR_LABELS[[a]]),
         LL_restricted = m$maximum, LL_full = model_cov$maximum,
         lr_stat = lr, lr_df = length(fix),
         lr_p = pchisq(lr, df = length(fix), lower.tail = FALSE))
})

lr_all <- 2 * (model_cov$maximum - model_null$maximum)
lrtests <- bind_rows(
  tibble(comparison  = "baseline MNL",
         description = "all twenty interactions zero (02_mnl.R)",
         LL_restricted = model_null$maximum, LL_full = model_cov$maximum,
         lr_stat = lr_all, lr_df = length(gamma_names),
         lr_p = pchisq(lr_all, df = length(gamma_names), lower.tail = FALSE)),
  per_attribute
) %>%
  mutate(significant = lr_p < 0.05)

write_table(lrtests, "mnlcov_lrtests", prefix = "03")

cat("\n=== Likelihood ratio tests =================================\n")
show_table(lrtests %>% select(-description), digits = 4)

# --- Estimates -------------------------------------------------------------
se_cov  <- sqrt(diag(model_cov$varcov))
rse_cov <- sqrt(diag(model_cov$robvarcov))

estimates <- tibble(
  parameter = names(model_cov$estimate),
  estimate  = as.numeric(model_cov$estimate),
  se        = as.numeric(se_cov[names(model_cov$estimate)]),
  rob_se    = as.numeric(rse_cov[names(model_cov$estimate)])
) %>%
  mutate(
    attribute   = str_extract(parameter, "(?<=^[bg]_)[a-z]{2}"),
    covariate   = ifelse(str_starts(parameter, "g_"),
                         str_extract(parameter, "(?<=_)[a-z]+$"), "(base)"),
    t_ratio     = estimate / se,
    rob_t_ratio = estimate / rob_se,
    significant = abs(rob_t_ratio) > 1.96
  )

write_table(estimates, "mnlcov_estimates", prefix = "03")

cat("\n=== Estimates ==============================================\n")
show_table(estimates %>% select(parameter, attribute, covariate, estimate,
                                rob_se, rob_t_ratio, significant), digits = 5)

# --- Implied coefficients and valuations by respondent type ----------------
# The gammas are only interesting through what they do to the coefficients,
# and the coefficients are only interpretable as money. Standard errors come
# from the delta method with the profile's covariate values substituted in as
# numeric constants, so a profile VTT carries the uncertainty of the whole
# linear combination rather than of b_tt alone.
# Archetypes and the income SD both come from 00_setup.R, so 03 and 05
# describe the SAME respondents. They had drifted to different profile sets,
# which made the two tables in the report silently non-comparable.
inc_sd   <- respondent_income_sd(database)
profiles <- respondent_profiles(inc_sd)

# "b_tt + g_tt_inc*0.42 + g_tt_car*1 + ..." for a given attribute and profile.
coef_expression <- function(a, x) {
  terms <- c(paste0("b_", a),
             sprintf("%s*%.10g", paste0("g_", a, "_", unname(MODEL_COVARS)),
                     unlist(x[names(MODEL_COVARS)])))
  paste(terms, collapse = " + ")
}

# Profiles are keyed prof1..profN in the expression labels rather than by
# their display names: apollo_deltaMethod round-trips those labels through a
# data frame and names with spaces and brackets do not survive the trip.
profile_key <- tibble(profile_id = paste0("prof", seq_along(profiles)),
                      profile    = names(profiles))

profile_expr <- unlist(lapply(seq_along(profiles), function(i) {
  x  <- profiles[[i]]
  id <- profile_key$profile_id[i]
  c(
    setNames(lapply(ROUTE_ATTRS, function(a) coef_expression(a, x)),
             paste0("b_", ROUTE_ATTRS, "__", id)),
    setNames(list(sprintf("(%s)/(%s)*60", coef_expression("tt", x),
                          coef_expression("tc", x))),
             paste0("vtt_chf_per_hour__", id)),
    setNames(list(sprintf("(%s)/(%s)", coef_expression("ch", x),
                          coef_expression("tc", x))),
             paste0("interchange_chf__", id))
  )
}))

profile_dm <- quietly(apollo_deltaMethod(
  model_cov, list(expression = unlist(profile_expr))))

profile_tbl <- as_tibble(profile_dm) %>%
  rename(quantity = Expression, value = Value, se = `s.e.`,
         t_ratio = `t-ratio (0)`) %>%
  separate(quantity, into = c("quantity", "profile_id"), sep = "__") %>%
  left_join(profile_key, by = "profile_id") %>%
  select(-profile_id) %>%
  mutate(ci_low = value - 1.96 * se, ci_high = value + 1.96 * se)

write_table(profile_tbl, "mnlcov_profiles", prefix = "03")

cat("\n=== Implied coefficients and valuations by profile =========\n")
show_table(profile_tbl %>% filter(str_starts(quantity, "vtt|interchange")),
           digits = 3)

# --- Hold-out validation ---------------------------------------------------
# Twenty extra parameters will always improve in-sample fit. The question this
# answers is whether they improve PREDICTION for respondents the model never
# saw, which is the only thing that distinguishes real observed heterogeneity
# from twenty degrees of freedom spent on noise.
cov_dv <- function(b, db) {
  X <- attr_diff_matrix(db)
  # Row-varying coefficient for each attribute, then the usual inner product.
  coefs <- vapply(ROUTE_ATTRS, function(a) {
    b[[paste0("b_", a)]] +
      rowSums(vapply(names(MODEL_COVARS), function(cv) {
        b[[paste0("g_", a, "_", MODEL_COVARS[[cv]])]] * db[[cv]]
      }, numeric(nrow(db))))
  }, numeric(nrow(db)))
  rowSums(coefs * X)
}

split <- holdout_split(database)   # same split as 02 and 04: same seed, same rule
train <- split$train
test  <- split$test

cat("\n=== Hold-out estimation ====================================\n")
assign("database", train, envir = globalenv())
model_cov_train  <- estimate_variant(character(0), "Swiss_MNLcov_train")
model_base_train <- estimate_variant(gamma_names,  "Swiss_MNL_train_from03")
assign("database", apollo_swissRouteChoiceData %>%
         mutate(log_income = log(hh_inc_abs) - mean(log(hh_inc_abs))),
       envir = globalenv())

validation <- bind_rows(
  binary_logit_score(mnl_dv(model_base_train$estimate, train), train,
                     "training", "MNL"),
  binary_logit_score(mnl_dv(model_base_train$estimate, test), test,
                     "hold-out", "MNL"),
  binary_logit_score(cov_dv(model_cov_train$estimate, train), train,
                     "training", "MNL + covariates"),
  binary_logit_score(cov_dv(model_cov_train$estimate, test), test,
                     "hold-out", "MNL + covariates")
)

write_table(validation, "mnlcov_validation", prefix = "03")

cat("\n=== Hold-out validation ====================================\n")
show_table(validation, digits = 4)

gain <- validation %>%
  filter(sample == "hold-out") %>%
  summarise(gain = LL_per_obs[model == "MNL + covariates"] -
                   LL_per_obs[model == "MNL"]) %>%
  pull(gain)

cat(sprintf(
  "\nhold-out LL per observation, covariates minus baseline: %+.5f\n", gain))
cat(if (gain > 0) {
  "The interactions pay for themselves out of sample.\n"
} else {
  "The interactions do NOT generalise: they buy in-sample fit and lose it\nagain on respondents the model never saw.\n"
})

# --- Figure ----------------------------------------------------------------
p_gamma <- estimates %>%
  filter(str_starts(parameter, "g_")) %>%
  mutate(attribute = factor(attribute, ROUTE_ATTRS,
                            unname(ROUTE_ATTR_LABELS[ROUTE_ATTRS]))) %>%
  ggplot(aes(estimate, covariate, colour = significant)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = estimate - 1.96 * rob_se,
                      xmax = estimate + 1.96 * rob_se)) +
  facet_wrap(~ attribute, scales = "free_x") +
  labs(title = "Covariate effects on each attribute coefficient",
       subtitle = "Intervals crossing zero mean that covariate does not shift the taste",
       x = "gamma (robust 95% CI)", y = NULL, colour = "p < 0.05")

p_vtt <- profile_tbl %>%
  filter(quantity == "vtt_chf_per_hour") %>%
  ggplot(aes(value, fct_reorder(profile, value))) +
  geom_vline(xintercept = profile_tbl$value[
    profile_tbl$quantity == "vtt_chf_per_hour" &
      profile_tbl$profile == "average respondent"],
    linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high)) +
  labs(title = "Value of travel time by respondent type",
       subtitle = "Dashed line is the average respondent; delta-method 95% CI",
       x = "CHF per hour", y = NULL)

p_lr <- lrtests %>%
  ggplot(aes(lr_stat, fct_reorder(comparison, lr_stat), fill = significant)) +
  geom_col() +
  labs(title = "Likelihood ratio statistics",
       subtitle = "Each bar drops one block of interactions from the full model",
       x = "LR statistic", y = NULL, fill = "p < 0.05")

p_val <- validation %>%
  ggplot(aes(LL_per_obs, sample, fill = model)) +
  geom_col(position = "dodge") +
  coord_cartesian(xlim = c(min(validation$LL_per_obs) - 0.01, 0)) +
  labs(title = "Fit per observation, training vs hold-out",
       subtitle = "Closer to zero is better",
       x = "LL per observation", y = NULL, fill = NULL)

fig <- (p_gamma | p_vtt) / (p_lr | p_val)
fig_path <- write_figure(fig, "mnlcov", prefix = "03", width = 15, height = 10)

cat(sprintf("\nModel object: %s\n",
            file.path(PATH_MODELS, paste0(MODEL_NAME, "_model.rds"))))
cat(sprintf("Tables written to: %s\n", PATH_TABLES))
cat(sprintf("Figure written to: %s\n", fig_path))
