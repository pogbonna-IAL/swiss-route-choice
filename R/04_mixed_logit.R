# 04_mixed_logit.R -- random parameters logit
# Swiss route choice
#
# 03 lets tastes vary with covariates the survey collected. This script lets
# them vary in ways it did not: each respondent draws their own coefficient
# vector from a continuous distribution, and only that distribution's
# parameters are estimated. It is the third member of the heterogeneity family
# -- observed (03), continuous unobserved (here), discrete unobserved (05/06)
# -- and the one the panel structure in 02 was set up for.
#
# SPECIFICATION. Every coefficient is a NEGATIVE LOGNORMAL,
#
#   b_a,n = -exp(mu_a + (L z_n)_a),      z_n ~ N(0, I)
#
# rather than a normal. A normal puts mass on both sides of zero, so a share
# of respondents would be estimated to prefer longer, dearer journeys with
# more interchanges -- 06 already treats a wrong-signed coefficient as a
# symptom rather than a finding, and a distribution that guarantees them by
# construction is worse than a class that stumbles into one. The lognormal
# also makes willingness to pay a ratio of two lognormals, which is itself
# lognormal and has closed-form quantiles.
#
# Two models are estimated:
#   MXL-I   L diagonal: four independent taste distributions      (8 par)
#   MXL-C   L lower triangular: tastes correlated across attributes (14 par)
#
# MXL-I nests the MNL only in the limit (all sigmas to zero, which is on the
# boundary of the parameter space), so the MNL comparison is reported on BIC
# and out-of-sample fit rather than as a likelihood ratio test -- the usual
# chi-squared reference distribution does not apply on a boundary. MXL-C
# nests MXL-I in the interior, so that one IS a clean LR test on 6 df.
#
# RUN ORDER: after 02_mnl.R.
#
# Outputs
#   outputs/models/Swiss_MXL_indep_model.rds  / _output.txt
#   outputs/models/Swiss_MXL_corr_model.rds   / _output.txt
#   outputs/tables/04_mxl_comparison.csv      fit against the MNL and each other
#   outputs/tables/04_mxl_estimates.csv       structural parameters
#   outputs/tables/04_mxl_moments.csv         implied coefficient distributions
#   outputs/tables/04_mxl_correlations.csv    taste correlation matrix
#   outputs/tables/04_mxl_vtt.csv             the WTP distribution, not a point
#   outputs/tables/04_mxl_validation.csv      hold-out prediction
#   outputs/figures/04_mxl.png
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))
source(here::here("R", "model_helpers.R"))

# Simulated maximum likelihood is only as accurate as its draws, and the
# estimate is biased for finite R. 200 MLHS draws is enough for a stable
# optimum on four dimensions and keeps a full re-run inside a coffee break;
# published work on this specification would use several thousand. MLHS rather
# than Halton because Halton sequences correlate badly in higher dimensions.
N_DRAWS  <- 200L
DRAW_SET <- c("d_tt", "d_tc", "d_hw", "d_ch")


apollo_initialise()

# --- Data ------------------------------------------------------------------
data("apollo_swissRouteChoiceData", package = "apollo")
database <- apollo_swissRouteChoiceData

mnl_path <- file.path(PATH_MODELS, "Swiss_MNL_model.rds")
if (!file.exists(mnl_path)) {
  stop("Run 02_mnl.R first -- ", mnl_path, " not found.")
}
model_mnl <- readRDS(mnl_path)

# --- Draws -----------------------------------------------------------------
# Inter-respondent draws only: a respondent has ONE taste vector that persists
# across their nine tasks. Intra-respondent draws would say tastes are redrawn
# every task, which is a different and much weaker claim about heterogeneity.
apollo_draws <- list(
  interDrawsType = "mlhs",
  interNDraws    = N_DRAWS,
  interUnifDraws = c(),
  interNormDraws = DRAW_SET,
  intraDrawsType = "mlhs",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)

apollo_control <- list(
  modelName       = "Swiss_MXL_indep",
  modelDescr      = "Mixed logit, independent negative lognormal coefficients",
  indivID         = "ID",
  outputDirectory = PATH_MODELS,
  panelData       = TRUE,
  mixing          = TRUE,
  seed            = SEED,
  nCores          = 1
)

# --- Utility ---------------------------------------------------------------
# Identical for both models; only apollo_randCoeff changes.
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
    avail        = 1,
    choiceVar    = choice,
    utilities    = V
  )

  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  # Product over the respondent's nine tasks FIRST, then average over draws.
  # Doing it the other way round would integrate each task separately and
  # throw away the panel structure, which is the entire point of the model.
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# --- MXL-I : independent coefficients --------------------------------------
apollo_randCoeff <- function(apollo_beta, apollo_inputs) {
  randcoeff <- list()
  randcoeff[["b_tt"]] <- -exp(mu_tt + L_tt_tt * d_tt)
  randcoeff[["b_tc"]] <- -exp(mu_tc + L_tc_tc * d_tc)
  randcoeff[["b_hw"]] <- -exp(mu_hw + L_hw_hw * d_hw)
  randcoeff[["b_ch"]] <- -exp(mu_ch + L_ch_ch * d_ch)
  return(randcoeff)
}

# Start the medians at the MNL point estimates: log|b_a| is the natural
# parameter, so mu_a = log|b_a^MNL| puts the median respondent exactly where
# the fixed-coefficient model put everybody.
mu_start <- setNames(log(abs(as.numeric(model_mnl$estimate[paste0("b_", ROUTE_ATTRS)]))),
                     paste0("mu_", ROUTE_ATTRS))

apollo_beta <- c(
  mu_start,
  L_tt_tt = 0.5, L_tc_tc = 0.5, L_hw_hw = 0.5, L_ch_ch = 0.5
)
apollo_fixed <- c()

cat("\n=== MXL-I: independent lognormal coefficients ==============\n")
apollo_inputs <- apollo_validateInputs()

# Cached: simulated maximum likelihood over 200 draws is the expensive part
# of this script, and it has no reason to be repeated when only the reporting
# or the scoring below has changed. REFIT=1 forces a re-estimation.
model_indep <- cached_model("Swiss_MXL_indep", function() {
  apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
})
apollo_modelOutput(model_indep)

writeLines(capture.output(apollo_modelOutput(model_indep)),
           file.path(PATH_MODELS, "Swiss_MXL_indep_output.txt"))

# --- MXL-C : correlated coefficients ---------------------------------------
# L is the lower-triangular Cholesky factor of the covariance of the LOGGED
# absolute coefficients. Parameterising the factor rather than the covariance
# matrix keeps the estimate positive semi-definite by construction, with no
# constraint for the optimiser to violate.
apollo_randCoeff <- function(apollo_beta, apollo_inputs) {
  randcoeff <- list()
  randcoeff[["b_tt"]] <- -exp(mu_tt + L_tt_tt * d_tt)
  randcoeff[["b_tc"]] <- -exp(mu_tc + L_tc_tt * d_tt + L_tc_tc * d_tc)
  randcoeff[["b_hw"]] <- -exp(mu_hw + L_hw_tt * d_tt + L_hw_tc * d_tc +
                                      L_hw_hw * d_hw)
  randcoeff[["b_ch"]] <- -exp(mu_ch + L_ch_tt * d_tt + L_ch_tc * d_tc +
                                      L_ch_hw * d_hw + L_ch_ch * d_ch)
  return(randcoeff)
}

# Start from the independent solution with zero off-diagonals: the restricted
# model is an interior point of the unrestricted one, so this start is already
# at a good likelihood and the search only has to find the correlations.
apollo_beta <- c(
  model_indep$estimate[paste0("mu_", ROUTE_ATTRS)],
  L_tt_tt = unname(model_indep$estimate[["L_tt_tt"]]),
  L_tc_tt = 0, L_tc_tc = unname(model_indep$estimate[["L_tc_tc"]]),
  L_hw_tt = 0, L_hw_tc = 0, L_hw_hw = unname(model_indep$estimate[["L_hw_hw"]]),
  L_ch_tt = 0, L_ch_tc = 0, L_ch_hw = 0,
  L_ch_ch = unname(model_indep$estimate[["L_ch_ch"]])
)

apollo_control$modelName  <- "Swiss_MXL_corr"
apollo_control$modelDescr <- "Mixed logit, correlated negative lognormal coefficients"

cat("\n=== MXL-C: correlated lognormal coefficients ===============\n")
apollo_inputs <- apollo_validateInputs()
model_corr <- cached_model("Swiss_MXL_corr", function() {
  apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
})
apollo_modelOutput(model_corr)

writeLines(capture.output(apollo_modelOutput(model_corr)),
           file.path(PATH_MODELS, "Swiss_MXL_corr_output.txt"))

# --- Comparison ------------------------------------------------------------
fit_row <- function(m, name, note) {
  n_par <- length(m$estimate)
  tibble(model = name, parameters = n_par, LL = m$maximum,
         AIC = -2 * m$maximum + 2 * n_par,
         BIC = -2 * m$maximum + n_par * log(m$nObs),
         note = note)
}

comparison <- bind_rows(
  fit_row(model_mnl,   "MNL (02)",  "fixed coefficients"),
  fit_row(model_indep, "MXL-I",     "independent lognormals"),
  fit_row(model_corr,  "MXL-C",     "correlated lognormals")
) %>%
  mutate(best_BIC = BIC == min(BIC))

# MXL-C nests MXL-I in the interior of the parameter space (the six
# off-diagonals are zero, not on a boundary), so this test is valid where the
# MNL comparison would not be.
lr_stat <- 2 * (model_corr$maximum - model_indep$maximum)
lr_df   <- length(model_corr$estimate) - length(model_indep$estimate)

comparison <- comparison %>%
  mutate(lr_vs_indep = ifelse(model == "MXL-C", lr_stat, NA_real_),
         lr_df       = ifelse(model == "MXL-C", lr_df, NA_integer_),
         lr_p        = ifelse(model == "MXL-C",
                              pchisq(lr_stat, df = lr_df, lower.tail = FALSE),
                              NA_real_))

write_table(comparison, "mxl_comparison", prefix = "04")

cat("\n=== Model comparison =======================================\n")
show_table(comparison, digits = 4)
cat(paste0(
  "\nThe MNL row is reported for reference only. Testing it against MXL-I by\n",
  "likelihood ratio would put the null on the boundary of the parameter space\n",
  "(all sigmas zero), where the chi-squared reference distribution does not\n",
  "hold. BIC and the hold-out table below are the honest comparisons.\n"))

# --- Estimates -------------------------------------------------------------
estimates <- bind_rows(
  tibble(model = "MXL-I", parameter = names(model_indep$estimate),
         estimate = as.numeric(model_indep$estimate),
         se = as.numeric(sqrt(diag(model_indep$varcov))[names(model_indep$estimate)])),
  tibble(model = "MXL-C", parameter = names(model_corr$estimate),
         estimate = as.numeric(model_corr$estimate),
         se = as.numeric(sqrt(diag(model_corr$varcov))[names(model_corr$estimate)]))
) %>%
  mutate(t_ratio = estimate / se)

write_table(estimates, "mxl_estimates", prefix = "04")

cat("\n=== Estimates ==============================================\n")
show_table(estimates, digits = 4)

# --- Implied distributions -------------------------------------------------
# The structural parameters are on the log scale and nobody can read a taste
# distribution off them. These are the quantities that mean something: where
# the median respondent sits, how much spread there is, and how wide the
# middle 90% of the population is.
chol_matrix <- function(est) {
  L <- matrix(0, 4, 4, dimnames = list(ROUTE_ATTRS, ROUTE_ATTRS))
  for (i in seq_along(ROUTE_ATTRS)) {
    for (j in seq_len(i)) {
      nm <- sprintf("L_%s_%s", ROUTE_ATTRS[i], ROUTE_ATTRS[j])
      if (nm %in% names(est)) L[i, j] <- unname(est[[nm]])
    }
  }
  L
}

moments_of <- function(m, name) {
  est <- m$estimate
  L   <- chol_matrix(est)
  Sig <- L %*% t(L)
  map_dfr(seq_along(ROUTE_ATTRS), function(i) {
    a  <- ROUTE_ATTRS[i]
    mu <- unname(est[[paste0("mu_", a)]])
    s2 <- Sig[i, i]
    s  <- sqrt(s2)
    tibble(
      model      = name,
      attribute  = a,
      label      = unname(ROUTE_ATTR_LABELS[[a]]),
      mu         = mu,
      sigma      = s,
      # Negative lognormal: the coefficient is -exp(mu + s z).
      median     = -exp(mu),
      mean       = -exp(mu + s2 / 2),
      sd         = exp(mu + s2 / 2) * sqrt(exp(s2) - 1),
      q05        = -exp(mu + s * qnorm(0.95)),
      q95        = -exp(mu + s * qnorm(0.05)),
      # How much of the population is within a factor of two of the median.
      within_2x  = pnorm(log(2) / s) - pnorm(-log(2) / s)
    )
  })
}

moments <- bind_rows(moments_of(model_indep, "MXL-I"),
                     moments_of(model_corr,  "MXL-C"))

write_table(moments, "mxl_moments", prefix = "04")

cat("\n=== Implied coefficient distributions ======================\n")
show_table(moments %>% select(-label), digits = 4)

# Taste correlations. This is what MXL-C buys over MXL-I: whether someone who
# hates travel time also hates paying, which is exactly the correlation that
# determines the spread of willingness to pay.
L_corr   <- chol_matrix(model_corr$estimate)
Sig_corr <- L_corr %*% t(L_corr)
D        <- diag(1 / sqrt(diag(Sig_corr)))
Corr     <- D %*% Sig_corr %*% D
dimnames(Corr) <- list(ROUTE_ATTRS, ROUTE_ATTRS)

correlations <- as_tibble(Corr, rownames = "attribute")
write_table(correlations, "mxl_correlations", prefix = "04")

cat("\n=== Taste correlations (log scale, MXL-C) ==================\n")
show_table(correlations, digits = 3)

# --- Willingness to pay as a distribution ----------------------------------
# The MNL in 02 reports one VTT with a confidence interval, which is
# uncertainty about a single population number. This is a different object: a
# distribution ACROSS respondents, whose spread is a finding rather than an
# error bar. log(VTT/60) = mu_tt - mu_tc + (row_tt - row_tc) z is normal, so
# VTT is lognormal and its quantiles are exact.
vtt_dist <- function(m, name) {
  est <- m$estimate
  L   <- chol_matrix(est)
  d   <- L[1, ] - L[2, ]              # row_tt - row_tc
  s   <- sqrt(sum(d^2))
  mu  <- unname(est[["mu_tt"]]) - unname(est[["mu_tc"]]) + log(60)
  tibble(
    model  = name,
    measure = "vtt_chf_per_hour",
    median = exp(mu),
    mean   = exp(mu + s^2 / 2),
    sd     = exp(mu + s^2 / 2) * sqrt(exp(s^2) - 1),
    q10    = exp(mu + s * qnorm(0.10)),
    q25    = exp(mu + s * qnorm(0.25)),
    q75    = exp(mu + s * qnorm(0.75)),
    q90    = exp(mu + s * qnorm(0.90)),
    sigma_log = s
  )
}

vtt <- bind_rows(vtt_dist(model_indep, "MXL-I"), vtt_dist(model_corr, "MXL-C"))
write_table(vtt, "mxl_vtt", prefix = "04")

cat("\n=== Value of travel time: a distribution, not a point ======\n")
show_table(vtt, digits = 3)

mnl_vtt_path <- file.path(PATH_TABLES, "02_mnl_valuations.csv")
if (file.exists(mnl_vtt_path)) {
  mnl_vtt <- read.csv(mnl_vtt_path) %>%
    filter(measure == "vtt_chf_per_hour") %>% pull(value)
  cat(sprintf(
    "\nMNL point estimate: %.2f CHF/h. MXL-C median: %.2f, with the middle 80%%\nof respondents spread from %.2f to %.2f.\n",
    mnl_vtt, vtt$median[2], vtt$q10[2], vtt$q90[2]))
}

# --- Hold-out validation ---------------------------------------------------
# The panel likelihood does not factorise over tasks, so the hold-out score is
# the simulated panel likelihood (panel_logit_score in model_helpers.R): the
# product over a respondent's tasks is taken inside the draw average, not
# outside it. The draws here are independent of the ones used in estimation,
# which is deliberate -- reusing them would flatter the model.
resp_index <- function(db) match(db$ID, sort(unique(db$ID)))

# R here is the number of draws used to SCORE, not to estimate. The log of a
# simulated mean is biased downward and the bias shrinks with R: against this
# model's own in-sample likelihood (-1404.47 from Apollo), 200 iid draws give
# -1421 and 2000 give about -1409. Apollo needs fewer because MLHS draws are
# stratified; plain iid draws here need more to reach the same accuracy.
mxl_dv_draws <- function(est, db, R = 2000L, seed = SEED + 7L) {
  set.seed(seed)
  L    <- chol_matrix(est)
  mu   <- vapply(ROUTE_ATTRS, function(a) unname(est[[paste0("mu_", a)]]),
                 numeric(1))
  X    <- attr_diff_matrix(db)
  ridx <- resp_index(db)
  n    <- max(ridx)

  out <- matrix(0, nrow = nrow(db), ncol = R)
  for (r in seq_len(R)) {
    Z   <- matrix(rnorm(n * 4), nrow = n)
    B   <- -exp(matrix(mu, n, 4, byrow = TRUE) + Z %*% t(L))
    out[, r] <- rowSums(B[ridx, , drop = FALSE] * X)
  }
  out
}

split <- holdout_split(database)   # same split as 02 and 03
train <- split$train
test  <- split$test

cat("\n=== Hold-out estimation (MXL-C on the training half) =======\n")
assign("database", train, envir = globalenv())
apollo_control$modelName  <- "Swiss_MXL_corr_train"
apollo_control$modelDescr <- "Mixed logit, correlated, training half"
assign("apollo_beta", model_corr$estimate, envir = globalenv())
assign("apollo_inputs", quietly(apollo_validateInputs(silent = TRUE)),
       envir = globalenv())
# Cached so the hold-out table can be recomputed -- with more draws, or after
# a change to the scoring code -- without paying for the estimation again.
model_corr_train <- cached_model("Swiss_MXL_corr_train", function() {
  quietly(apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities,
                          apollo_inputs,
                          estimate_settings = list(silent = TRUE,
                                                   writeIter = FALSE)))
})
assign("database", apollo_swissRouteChoiceData, envir = globalenv())
# The baseline has to be fitted on the SAME training half, or the comparison
# is rigged: a full-sample MNL has already seen every hold-out respondent,
# and beating it would prove nothing about the mixed logit. fit_binary_logit
# gives the exact MNL solution without swapping Apollo's global model
# definition mid-script (see model_helpers.R).
b_mnl_train <- fit_binary_logit(train)

validation <- bind_rows(
  binary_logit_score(mnl_dv(b_mnl_train, train), train, "training", "MNL"),
  binary_logit_score(mnl_dv(b_mnl_train, test),  test,  "hold-out", "MNL"),
  panel_logit_score(mxl_dv_draws(model_corr_train$estimate, train), train,
                    "training", "MXL-C"),
  panel_logit_score(mxl_dv_draws(model_corr_train$estimate, test), test,
                    "hold-out", "MXL-C")
)

write_table(validation, "mxl_validation", prefix = "04")

cat("\n=== Hold-out validation ====================================\n")
show_table(validation, digits = 4)

# --- Figure ----------------------------------------------------------------
grid <- tibble(x = seq(0.001, 0.999, length.out = 400))

p_dist <- moments %>%
  filter(model == "MXL-C") %>%
  mutate(label = factor(label, unname(ROUTE_ATTR_LABELS[ROUTE_ATTRS]))) %>%
  crossing(grid) %>%
  mutate(value = -exp(mu + sigma * qnorm(x))) %>%
  ggplot(aes(value, x)) +
  geom_line() +
  geom_vline(xintercept = 0, colour = "grey60") +
  facet_wrap(~ label, scales = "free_x") +
  labs(title = "Where the population sits on each coefficient (MXL-C)",
       subtitle = "Quantile function; every respondent is on the correct side of zero",
       x = "coefficient", y = "population quantile")

p_vtt <- tibble(q = seq(0.005, 0.995, length.out = 400)) %>%
  mutate(MXL_C = exp(log(vtt$median[2]) + vtt$sigma_log[2] * qnorm(q))) %>%
  ggplot(aes(MXL_C, q)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = vtt$median[2], linetype = "dashed", colour = "grey40") +
  scale_x_continuous(trans = "log10") +
  labs(title = "Value of travel time across respondents",
       subtitle = "Dashed line is the median; log scale",
       x = "CHF per hour (log scale)", y = "population quantile")

p_corr <- correlations %>%
  pivot_longer(-attribute, names_to = "with", values_to = "rho") %>%
  ggplot(aes(with, attribute, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.2) +
  scale_fill_gradient2(limits = c(-1, 1), low = "steelblue", high = "firebrick") +
  labs(title = "Taste correlations (MXL-C)",
       subtitle = "On the log-coefficient scale", x = NULL, y = NULL, fill = NULL)

p_fit <- comparison %>%
  ggplot(aes(BIC, fct_reorder(model, -BIC), fill = best_BIC)) +
  geom_col() +
  coord_cartesian(xlim = c(min(comparison$BIC) * 0.97, max(comparison$BIC) * 1.01)) +
  labs(title = "BIC", subtitle = "Lower is better", x = NULL, y = NULL,
       fill = "best")

fig <- (p_dist | p_vtt) / (p_corr | p_fit)
fig_path <- write_figure(fig, "mxl", prefix = "04", width = 15, height = 10)

cat(sprintf("\nModels written to: %s\n", PATH_MODELS))
cat(sprintf("Tables written to: %s\n", PATH_TABLES))
cat(sprintf("Figure written to: %s\n", fig_path))
