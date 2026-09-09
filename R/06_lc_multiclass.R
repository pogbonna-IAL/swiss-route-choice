# 06_lc_multiclass.R -- latent class models, 1 to 5 classes
# Swiss route choice
#
# Estimates LC1 through LC5: latent class models in which every class has its
# own full set of four route attribute coefficients, with a constant-only
# class allocation. LC1 is the plain MNL, so the sequence nests the baseline
# from 02_mnl.R and the whole comparison rests on a common likelihood.
#
# Latent class likelihoods are multi-modal. A single set of starting values
# routinely lands on a local optimum, so each K is estimated from N_STARTS
# perturbed starts and the best log-likelihood is kept. The share of starts
# that reach that best value is recorded as a first read on how fragile the
# solution is; 07_lc_stability.R probes this properly.
#
# Class labels are identified only up to permutation, so every model is put in
# canonical order (classes sorted by b_tt, most negative first) before
# anything is reported. Without that, "class 1" here and "class 1" in
# 05_lc_2class.R are not the same class and the two scripts cannot be read
# together.
#
# RUN ORDER: after 02_mnl.R, before 05_lc_2class.R (which reads the comparison
# table written here). See run_all.R.
#
# Outputs
#   outputs/models/Swiss_LC<K>_model.rds     model object per K
#   outputs/models/Swiss_LC<K>_output.txt    formatted output per K
#   outputs/tables/06_lc_comparison.csv      LL / AIC / BIC / class shares
#   outputs/tables/06_lc_parameters.csv      class-specific estimates
#   outputs/tables/06_lc_valuations.csv      class-specific WTP with delta-method s.e.
#   outputs/tables/06_lc_runs.csv            per-start log for every K
#   outputs/figures/06_lc_selection.png      information criteria vs K
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))
source(here::here("R", "lc_helpers.R"))

K_MAX    <- 5L
N_STARTS <- 50L    # perturbed starting values per K (K >= 2)

# Completed models are reloaded from outputs/models/ instead of re-estimated.
# Set REFIT=1 in the environment to force a full re-run.
RESUME <- lc_resume_enabled()

# --- Data ------------------------------------------------------------------
apollo_initialise()
data("apollo_swissRouteChoiceData", package = "apollo")
database <- apollo_swissRouteChoiceData

# Baseline MNL estimates anchor the starting values: classes are drawn as
# perturbations around them rather than from arbitrary constants.
mnl_path <- file.path(PATH_MODELS, "Swiss_MNL_model.rds")
if (!file.exists(mnl_path)) {
  stop("Run 02_mnl.R first -- ", mnl_path, " not found.")
}
b_mnl <- readRDS(mnl_path)$estimate

# --- Estimation loop -------------------------------------------------------
# The code generators (lc_make_lcPars / lc_make_probabilities), the starting
# value scheme and the multi-start driver all live in lc_helpers.R. This
# script used to carry its own copy of all four, which meant two definitions
# of the function that writes Apollo's likelihood could drift apart silently.
results  <- list()
all_runs <- list()
verdicts <- list()

for (K in seq_len(K_MAX)) {

  cat(sprintf("\n\n================ LC%d : %d class%s ================\n",
              K, K, if (K == 1L) "" else "es"))

  lc_install(K)

  n_draws <- if (K == 1L) 1L else N_STARTS
  starts  <- lapply(seq_len(n_draws), function(draw) {
    lc_start_values(K, draw, b_mnl, seed = SEED + 1000L * K + draw)
  })

  # Cached per K: five models at fifty starts each runs close to an hour, and
  # an interruption used to discard all of it. See lc_fit_cached.
  fit <- lc_fit_cached(paste0("Swiss_LC", K), resume = RESUME, fit_fn = function() {
    lc_estimate_best(
      starts       = starts,
      apollo_fixed = if (K == 1L) c() else c("delta_1"),
      control = list(
        modelName       = paste0("Swiss_LC", K),
        modelDescr      = sprintf("Latent class route choice, %d class(es)", K),
        indivID         = "ID",
        outputDirectory = PATH_MODELS,
        panelData       = TRUE,
        seed            = SEED,
        nCores          = 1
      )
    )
  })

  stopifnot(fit$ok)

  # A warning() here was too quiet: R buffers warnings to the end of the run,
  # where they sit under a "There were N warnings" line a reader scrolls past.
  # The verdict is printed inline, next to the model it condemns.
  verdicts[[K]] <- lc_fit_verdict(
    fit, K, shares = if (K > 1L) lc_class_shares(fit$model$estimate, K) else NULL)
  lc_announce_verdict(verdicts[[K]], sprintf("LC%d", K))

  model <- fit$model
  apollo_modelOutput(model)

  # The .rds is written by lc_fit_cached. Neither uses apollo_saveOutput,
  # which renames any existing file to _OLD1, _OLD2, ... on every re-run and
  # clutters a five-model sweep with copies nothing reads.
  writeLines(capture.output(apollo_modelOutput(model)),
             file.path(PATH_MODELS, sprintf("Swiss_LC%d_output.txt", K)))

  results[[K]]  <- fit
  all_runs[[K]] <- fit$run_log %>% mutate(model = paste0("LC", K), .before = 1)
}

write_table(bind_rows(all_runs), "lc_runs", prefix = "06")

# --- Comparison ------------------------------------------------------------
comparison <- map_dfr(seq_len(K_MAX), function(K) {
  fit   <- results[[K]]
  m     <- fit$model
  n_par <- lc_n_par(K)
  ll    <- m$maximum
  runs  <- fit$run_log

  # Class shares reported in canonical order so the string lines up with the
  # class numbering used in the parameter and valuation tables.
  ord  <- lc_class_order(m$estimate, K)
  pi_k <- lc_class_shares(m$estimate, K)[ord]

  tibble(
    model           = paste0("LC", K),
    classes         = K,
    parameters      = n_par,
    LL              = ll,
    AIC             = -2 * ll + 2 * n_par,
    BIC             = -2 * ll + n_par * log(m$nObs),
    min_class_share = min(pi_k),
    max_class_share = max(pi_k),
    class_shares    = paste(sprintf("%.3f", pi_k), collapse = " "),
    starts          = nrow(runs),
    converged       = sum(runs$converged),
    # How often the search rediscovered the best solution. A low share means
    # the likelihood surface is rough and the optimum is easy to miss.
    at_best         = fit$n_at_best,
    hessian_ok      = fit$hessian_ok,
    # Carried into the table so a reader of the CSV sees what a reader of the
    # console saw. A comparison table that ranks models on BIC without saying
    # which of them are estimable is a trap.
    usable          = verdicts[[K]]$usable,
    verdict         = lc_verdict_label(verdicts[[K]]),
    problems        = paste(verdicts[[K]]$flags, collapse = ";")
  )
})

# lc_n_par is the arithmetic definition (4K free coefficients plus K-1 free
# deltas). Cross-check it against what Apollo actually estimated so a change
# to the specification cannot silently desynchronise the two.
for (K in seq_len(K_MAX)) {
  apollo_free <- length(results[[K]]$model$estimate) - (if (K == 1L) 0L else 1L)
  stopifnot(identical(as.integer(comparison$parameters[K]), as.integer(apollo_free)))
}

comparison <- comparison %>%
  mutate(
    delta_LL   = LL - lag(LL),
    lr_stat    = 2 * delta_LL,
    delta_par  = parameters - lag(parameters),
    delta_BIC  = BIC - lag(BIC),
    best_BIC   = BIC == min(BIC),
    best_AIC   = AIC == min(AIC)
  )

write_table(comparison, "lc_comparison", prefix = "06")

cat("\n\n=== Model comparison =======================================\n")
show_table(comparison %>%
             select(model, classes, parameters, LL, AIC, BIC,
                    min_class_share, at_best, starts, verdict))

# The selection verdict must not be read off BIC alone. If the model BIC
# prefers is not estimable, saying so IS the result.
i_best <- which.min(comparison$BIC)
if (!comparison$usable[i_best]) {
  cat(sprintf(paste0(
    "\n!! BIC selects %s, and that model is NOT USABLE (%s).\n",
    "!! Reporting a BIC ranking without this line would name a model that\n",
    "!! cannot be estimated as the best in the sweep.\n"),
    comparison$model[i_best], comparison$verdict[i_best]))
} else {
  cat(sprintf("\nBIC selects %s, which is usable.\n", comparison$model[i_best]))
}
if (any(!comparison$usable)) {
  cat(sprintf("unusable models in this sweep: %s\n",
              paste(comparison$model[!comparison$usable], collapse = ", ")))
}

# --- Class-specific parameters ---------------------------------------------
# Canonical order, with the standard errors permuted alongside the estimates
# and the deltas re-referenced as contrasts against the new class 1 (see
# lc_canonical_params).
parameters <- map_dfr(seq_len(K_MAX), function(K) {
  lc_canonical_params(results[[K]]$model, K, covariate_alloc = FALSE) %>%
    mutate(model = paste0("LC", K), classes = K, .before = 1)
})

write_table(parameters, "lc_parameters", prefix = "06")

# --- Valuations -------------------------------------------------------------
# Willingness-to-pay per class WITH delta-method standard errors on the robust
# covariance matrix. A ratio of two estimates is not normally distributed and
# its uncertainty is far wider than either coefficient suggests alone, so
# reporting the point ratio on its own -- as this script used to -- says
# nothing about whether two classes actually value time differently.
valuations <- map_dfr(seq_len(K_MAX), function(K) {
  fit <- results[[K]]
  if (!fit$hessian_ok) return(tibble())
  m   <- fit$model
  ord <- lc_class_order(m$estimate, K)
  est <- lc_canonical(m$estimate, K)

  shares <- lc_class_shares(m$estimate, K)[ord]
  coefs  <- map_dfr(seq_len(K), function(k) {
    b <- vapply(ROUTE_ATTRS, function(a) unname(est[[paste0("b_", a, "_", k)]]),
                numeric(1))
    tibble(class = k, share = shares[k],
           b_tt = b[["tt"]], b_tc = b[["tc"]],
           b_hw = b[["hw"]], b_ch = b[["ch"]],
           signs_ok = all(b < 0))
  })

  lc_valuations(m, K, ord = ord) %>%
    left_join(coefs, by = "class") %>%
    mutate(model = paste0("LC", K), classes = K, .before = 1)
})

write_table(valuations, "lc_valuations", prefix = "06")

cat("\n=== Class-specific valuations (delta-method s.e.) ==========\n")
show_table(valuations %>%
             select(model, class, share, measure, value, se, t_ratio,
                    ci_low, ci_high, ratio_reliable, signs_ok))

if (any(!valuations$ratio_reliable)) {
  cat("\nNOTE: rows with ratio_reliable = FALSE have a cost coefficient that is\n")
  cat("not distinguishable from zero. Their WTP ratios are not interpretable.\n")
}

# --- Selection figure ------------------------------------------------------
ic_long <- comparison %>%
  select(classes, AIC, BIC) %>%
  pivot_longer(-classes, names_to = "criterion", values_to = "value")

p_ic <- ggplot(ic_long, aes(classes, value, colour = criterion)) +
  geom_line() +
  geom_point(size = 2) +
  geom_point(data = ic_long %>% group_by(criterion) %>% slice_min(value, n = 1),
             size = 4, shape = 21, stroke = 1.2, fill = NA) +
  scale_x_continuous(breaks = seq_len(K_MAX)) +
  labs(title = "Information criteria by number of classes",
       subtitle = "Circled point is the minimum for each criterion",
       x = "classes", y = NULL, colour = NULL)

p_ll <- ggplot(comparison, aes(classes, LL)) +
  geom_line() + geom_point(size = 2) +
  scale_x_continuous(breaks = seq_len(K_MAX)) +
  labs(title = "Log-likelihood by number of classes",
       subtitle = "Always improves with K; the criteria above penalise that",
       x = "classes", y = "LL")

p_share <- ggplot(comparison, aes(classes, min_class_share)) +
  geom_col() +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq_len(K_MAX)) +
  labs(title = "Smallest class share",
       subtitle = "Below the dashed 5% line a class is too thin to interpret",
       x = "classes", y = "share")

p_stab <- ggplot(comparison, aes(classes, at_best / starts)) +
  geom_col() +
  scale_x_continuous(breaks = seq_len(K_MAX)) +
  ylim(0, 1) +
  labs(title = "Share of starts reaching the best solution",
       subtitle = "Lower means a rougher likelihood surface",
       x = "classes", y = "share of starts")

fig <- (p_ic | p_ll) / (p_share | p_stab)
fig_path <- write_figure(fig, "lc_selection", prefix = "06", width = 14, height = 9)

cat(sprintf("\nModels written to:  %s\n", PATH_MODELS))
cat(sprintf("Tables written to:  %s\n", PATH_TABLES))
cat(sprintf("Figure written to:  %s\n", fig_path))
