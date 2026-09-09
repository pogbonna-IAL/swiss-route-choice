# 07_lc_stability.R -- when is a latent class solution trustworthy?
# Swiss route choice
#
# Builds the stability table: for each sample size and number of classes, how
# often estimation succeeds, how small the smallest class gets, how much the
# log-likelihood moves across resamples, and how often a coefficient comes back
# with the wrong sign.
#
# Two different kinds of instability are measured, and they must not be
# confused:
#
#   * SAMPLING instability (n < FULL_N): the same model refitted to 20
#     different draws of respondents. Variation here is what a researcher with
#     one dataset of that size would face without ever knowing it.
#   * STARTING VALUE instability (n = FULL_N): the full panel has no sampling
#     variation at all, since it is the whole population of respondents. What
#     varies there is only the optimiser's starting point, taken from the
#     50-start sweep in 06_lc_multiclass.R.
#
# Reading a single "convergence rate" across both would be meaningless, so the
# full-sample rows are labelled separately and carry their own verdict scale.
#
# RUN ORDER: after 06_lc_multiclass.R and 08_small_sample_experiment.R. The
# file numbering does not reflect the dependency; see run_all.R.
#
# Inputs   outputs/tables/08_fits.csv, 08_parameters.csv,
#          06_lc_comparison.csv, 06_lc_parameters.csv
# Outputs  outputs/tables/07_stability.csv
#          outputs/figures/07_stability.png
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

# The smallest sample this table describes. Below it every answer is the
# same answer -- 24 of 57 rows read "meaningless" once the grid was extended
# to n = 5 -- which tells a reader nothing and buries the range where the
# verdict scale actually discriminates. 11_breakdown.R owns the tail and asks
# a different question of it: not "is this trustworthy" but "how does it fail".
MIN_N <- 30L

fits        <- read_required("08_fits.csv",          "08_small_sample_experiment.R")
params      <- read_required("08_parameters.csv",    "08_small_sample_experiment.R")
comparison  <- read_required("06_lc_comparison.csv", "06_lc_multiclass.R")
full_params <- read_required("06_lc_parameters.csv", "06_lc_multiclass.R")

# --- Sign reversals --------------------------------------------------------
# Every route attribute coefficient must be negative: more travel time, cost,
# headway or interchanges can only reduce utility. A positive estimate in any
# class is not a finding, it is a symptom that the class is being fitted to
# noise. Allocation constants (delta_*) are unrestricted and excluded.
is_attribute_par <- function(x) str_starts(x, "b_")

sign_check <- params %>%
  filter(is_attribute_par(parameter)) %>%
  group_by(n, seed, K) %>%
  summarise(wrong_sign_pars = sum(estimate > 0),
            any_wrong_sign  = as.integer(any(estimate > 0)),
            .groups = "drop")

# The same check on the full-sample models. This used to be hardcoded to zero,
# which happened to be the right answer but asserted a result instead of
# measuring one -- and would have gone on asserting it after a respecification.
full_sign_check <- full_params %>%
  filter(is_attribute_par(parameter)) %>%
  group_by(K = classes) %>%
  summarise(wrong_sign_pars = sum(estimate > 0),
            any_wrong_sign  = as.integer(any(estimate > 0)),
            .groups = "drop")

# --- Sampling instability (n < FULL_N) -------------------------------------
sampled_K <- sort(unique(fits$K[fits$n != FULL_N & fits$n >= MIN_N]))

sampling <- fits %>%
  filter(n != FULL_N, n >= MIN_N) %>%
  left_join(sign_check, by = c("n", "seed", "K")) %>%
  group_by(n, K) %>%
  summarise(
    replications   = n(),
    # "At least one start returned a finite likelihood" is 1.000 in every cell
    # of this experiment, so it is kept only as a check that no cell failed
    # outright -- it is not a convergence rate and must not be read as one.
    any_fit        = mean(ok),
    # This is the honest optimiser-reliability measure: of the starts actually
    # attempted in a replication, how many converged.
    mean_starts_ok = mean(n_starts_ok),
    # A covariance matrix is necessary but not sufficient. `usable` is the
    # full fit-time verdict recorded by 08 (lc_fit_verdict): converged, an
    # invertible Hessian, no wrong-signed coefficient, no class too thin to
    # read, and no two classes collapsed onto each other. Both are kept
    # because the gap between them is itself informative -- a cell where the
    # Hessian exists but the fit is still unusable failed for a different
    # reason than one where it does not.
    has_covariance = mean(hessian_ok),
    usable         = mean(usable),
    median_min_cls = median(min_class_share, na.rm = TRUE),
    min_min_cls    = min(min_class_share, na.rm = TRUE),
    ll_sd          = sd(LL, na.rm = TRUE),
    sign_reversals = sum(any_wrong_sign, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "resampling", resampled = TRUE)

# --- Starting value instability (n = FULL_N) -------------------------------
# 06 sweeps K = 1..5 while the resampling experiment covers K = 1..4, so the
# K = 5 row has no sampling counterpart. It is kept -- the starting-value
# evidence is real -- but flagged, so a reader does not compare a row that has
# resampling evidence against one that does not.
full_row <- comparison %>%
  transmute(
    n              = FULL_N,
    K              = classes,
    replications   = starts,
    any_fit        = converged / starts,
    has_covariance = NA_real_,
    # Deliberately NA: for the resampling rows this is the mean number of
    # converged starts within one replication, and the full sample has no
    # replications. Its counterpart here is `usable` -- the share of starts
    # that reached the best solution -- which is already reported.
    mean_starts_ok = NA_real_,
    # At the full sample the meaningful reliability figure is how often a
    # random start reaches the best solution found.
    usable         = at_best / starts,
    median_min_cls = min_class_share,
    min_min_cls    = min_class_share,
    ll_sd          = NA_real_,
    source         = "starting values",
    resampled      = classes %in% sampled_K
  ) %>%
  left_join(full_sign_check %>% select(K, sign_reversals = any_wrong_sign),
            by = "K")

# --- Verdict ---------------------------------------------------------------
# Thresholds are deliberately blunt. A solution is only called stable when it
# converges essentially always, keeps every class large enough to interpret,
# and never produces a wrong-signed coefficient.
stability <- bind_rows(full_row, sampling) %>%
  mutate(
    verdict = case_when(
      # The full-sample rows measure a different quantity -- how often a random
      # start finds the optimum -- so the sampling thresholds do not apply to
      # them. What they imply is how many starts the model needs, not whether
      # its solution is trustworthy.
      source == "starting values" & usable >= 0.30 ~ "20 starts enough",
      source == "starting values" & usable >= 0.10 ~ "50+ starts",
      source == "starting values"                  ~ "100+ starts",
      usable >= 0.95 & median_min_cls >= 0.10 & sign_reversals == 0 ~ "Yes",
      usable >= 0.80 & median_min_cls >= 0.05 & sign_reversals <= 1 ~ "Mostly",
      usable >= 0.50 & sign_reversals <= 3                          ~ "Weak",
      usable >= 0.25                                                ~ "No",
      TRUE                                                          ~ "meaningless"
    )
  ) %>%
  arrange(desc(n), K)

write_table(stability, "stability", prefix = "07")

cat("\n=== Latent class stability ==================================\n")
print(as.data.frame(
  stability %>%
    transmute(
      N = n, Classes = K,
      `Starts ok`   = ifelse(is.na(mean_starts_ok), "--",
                             sprintf("%.1f", mean_starts_ok)),
      `Has cov`     = ifelse(is.na(has_covariance), "--",
                             sprintf("%.0f%%", 100 * has_covariance)),
      `Usable`      = sprintf("%.0f%%", 100 * usable),
      `Median min class` = sprintf("%.1f%%", 100 * median_min_cls),
      `LL SD`       = ifelse(is.na(ll_sd), "--", sprintf("%.1f", ll_sd)),
      `Sign rev.`   = sign_reversals,
      `Resampled`   = ifelse(resampled, "", "no resampling"),
      `Stable?`     = verdict
    )), row.names = FALSE)

cat(sprintf(paste0(
  "\nn = %d is the whole panel and measures something different. There is\n",
  "no sampling variation to measure -- it is every respondent there is -- so\n",
  "'Usable' is the share of random starts that reached the best solution,\n",
  "not a convergence rate; 'Sign rev.' is a 0/1 flag on the one\n",
  "full-sample model rather than a count over replications; and 'Starts ok'\n",
  "is blank because there are no replications to average over. Rows marked\n",
  "'no resampling' have no small-sample counterpart because 06 sweeps\n",
  "five classes and 08 only four.\n"), FULL_N))

# --- Figure ----------------------------------------------------------------
plot_dat <- stability %>% filter(source == "resampling")

# The log axis itself is log_size_axis() in 00_setup.R; only the set of
# breaks differs between scripts, so this is a binding, not a second copy.
size_axis <- function(p) log_size_axis(p, plot_dat$n)

p_usable <- size_axis(
  ggplot(plot_dat, aes(n, usable, colour = factor(K))) +
    geom_line() + geom_point(size = 1.8) +
    geom_hline(yintercept = c(0.8, 0.95), linetype = "dashed", colour = "grey60") +
    ylim(0, 1)) +
  labs(title = "Share of replications with a usable covariance matrix",
       x = "respondents (log scale)", y = NULL, colour = "classes")

p_cls <- size_axis(
  ggplot(plot_dat, aes(n, median_min_cls, colour = factor(K))) +
    geom_line() + geom_point(size = 1.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey60")) +
  labs(title = "Median smallest class share",
       subtitle = "Below 5% a class cannot be interpreted",
       x = "respondents (log scale)", y = NULL, colour = "classes")

p_ll <- size_axis(
  ggplot(plot_dat, aes(n, ll_sd, colour = factor(K))) +
    geom_line() + geom_point(size = 1.8)) +
  labs(title = "Log-likelihood SD across resamples",
       x = "respondents (log scale)", y = NULL, colour = "classes")

p_sign <- size_axis(
  ggplot(plot_dat, aes(n, sign_reversals, colour = factor(K))) +
    geom_line() + geom_point(size = 1.8)) +
  labs(title = "Replications with a wrong-signed coefficient",
       subtitle = sprintf("out of %d", max(plot_dat$replications)),
       x = "respondents (log scale)", y = NULL, colour = "classes")

fig <- (p_usable | p_cls) / (p_ll | p_sign)
fig_path <- write_figure(fig, "stability", prefix = "07", width = 14, height = 9)

cat(sprintf("\nTable written to:  %s\n",
            file.path(PATH_TABLES, "07_stability.csv")))
cat(sprintf("Figure written to: %s\n", fig_path))
