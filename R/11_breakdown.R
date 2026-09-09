# 11_breakdown.R -- where latent class estimation stops working
# Swiss route choice
#
# 07_lc_stability.R asks whether a solution is trustworthy at a plausible
# study size. This script asks a different and blunter question: pushed below
# any sensible sample size, HOW does latent class estimation fail, and at what
# n does each failure mode appear?
#
# The distinction matters because the failures are not one failure. A model
# can converge, report standard errors, and still have produced nothing: two
# classes collapsed onto each other, or a class holding one respondent, or a
# coefficient with the wrong sign. Each of those is a different pathology with
# a different diagnostic, and they appear at different sample sizes.
#
# Reproduction is judged against the n = 388 benchmark, which is not a "true"
# parameter vector but is the whole population of respondents this study has.
# The question is therefore precisely: how small a subsample still recovers
# what the full panel says?
#
# RUN ORDER: after 08_small_sample_experiment.R.
#
# Inputs   outputs/tables/08_fits.csv, 08_parameters.csv
# Outputs  outputs/tables/11_failure_modes.csv    rate of each pathology by (n, K)
#          outputs/tables/11_reproduction.csv     what fraction recovers the benchmark
#          outputs/tables/11_separation.csv       class collapse diagnostics
#          outputs/tables/11_vtt_recovery.csv     MNL willingness-to-pay recovery
#          outputs/figures/11_breakdown.png
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))
source(here::here("R", "lc_helpers.R"))

# Thresholds. Deliberately blunt, and stated once here rather than buried in
# the case_when below.
# The class-share and class-collapse thresholds are NOT redefined here. They
# live in lc_helpers.R (LC_MIN_CLASS_SHARE, LC_COLLAPSE_REL) and are applied
# at fit time; this script reads the resulting flags. A second copy of a
# threshold is a second thing that can disagree with the first.
EXPLODE_AT <- 10   # a coefficient this many times the benchmark is not an estimate

# All permutations of 1..K, without a dependency for K <= 4.
combinat_perms <- function(K) {
  if (K == 1L) return(list(1L))
  out <- list()
  for (i in seq_len(K)) {
    rest <- combinat_perms(K - 1L)
    for (r in rest) out[[length(out) + 1L]] <- c(i, setdiff(seq_len(K), i)[r])
  }
  out
}

fits   <- read_required("08_fits.csv",       "08_small_sample_experiment.R")
params <- read_required("08_parameters.csv", "08_small_sample_experiment.R")

# The fit-time verdict (lc_fit_verdict) is recorded by 08. This script does
# NOT recompute those flags: two definitions of "usable" that can disagree is
# exactly the failure this project keeps running into. What it adds are the
# checks only a benchmark can support -- whether the classes still line up
# with the full-panel classes, and whether a coefficient has left the plausible
# range entirely.
if (!all(c("usable", "verdict", "problems") %in% names(fits))) {
  stop("08_fits.csv predates the usability verdict. Re-run ",
       "08_small_sample_experiment.R once -- it migrates the file in place.",
       call. = FALSE)
}

has_flag <- function(problems, flag) {
  vapply(strsplit(ifelse(is.na(problems), "", problems), ";"),
         function(f) as.integer(flag %in% f), integer(1))
}

is_attr_par <- function(x) str_starts(x, "b_")

# --- Benchmark -------------------------------------------------------------
# The full panel, one fit per K, in canonical class order (08 canonicalises
# before recording, so class k here means the same thing as class k there).
benchmark <- params %>%
  filter(n == FULL_N, is_attr_par(parameter)) %>%
  select(K, parameter, truth = estimate)

# Per-attribute scale for normalising distances: the benchmark magnitude of
# that coefficient, averaged over classes. Coefficients differ by two orders
# of magnitude across attributes (b_ch is ~1.2, b_hw is ~0.04), so an
# unnormalised Euclidean distance between class vectors would be entirely
# driven by interchanges.
attr_scale <- benchmark %>%
  mutate(attribute = str_extract(parameter, "(?<=^b_)[a-z]{2}")) %>%
  group_by(K, attribute) %>%
  summarise(scale = mean(abs(truth)), .groups = "drop")

# --- Per-replication diagnostics -------------------------------------------
# One row per (n, seed, K), each column a distinct way the fit can be useless.
attr_long <- params %>%
  filter(is_attr_par(parameter)) %>%
  mutate(attribute = str_extract(parameter, "(?<=^b_)[a-z]{2}"),
         class     = as.integer(str_extract(parameter, "(?<=_)[0-9]+$")))

# Class separation: for every pair of classes within a fit, the largest
# per-attribute gap measured in benchmark units. If the closest pair is under
# LC_COLLAPSE_REL the model has fitted one class twice and the extra parameters
# bought nothing.
separation <- attr_long %>%
  inner_join(attr_scale, by = c("K", "attribute")) %>%
  mutate(scaled = estimate / scale) %>%
  select(n, seed, K, attribute, class, scaled) %>%
  group_by(n, seed, K) %>%
  group_modify(function(d, key) {
    if (key$K < 2L || anyNA(d$scaled)) {
      return(tibble(min_pair_gap = NA_real_))
    }
    w <- d %>% pivot_wider(names_from = attribute, values_from = scaled) %>%
      arrange(class)
    m <- as.matrix(w %>% select(-class))
    if (nrow(m) < 2L) return(tibble(min_pair_gap = NA_real_))
    gaps <- combn(nrow(m), 2, function(ij) max(abs(m[ij[1], ] - m[ij[2], ])))
    tibble(min_pair_gap = min(gaps))
  }) %>%
  ungroup()

# Does canonical ordering still ALIGN classes with the benchmark?
#
# Sorting by b_tt fixes the order within a fit, but it only aligns two fits if
# the ordering statistic is estimated well enough to put the same population
# class in the same position. At small n it is not. Here the fit's classes are
# matched to the benchmark's by the permutation minimising total distance; if
# the winner is not the identity, canonical ordering has silently compared
# different classes -- which makes every bias, RMSE and coverage number for
# that replication a comparison between mismatched classes.
bench_vec <- benchmark %>%
  mutate(attribute = str_extract(parameter, "(?<=^b_)[a-z]{2}"),
         class     = as.integer(str_extract(parameter, "(?<=_)[0-9]+$"))) %>%
  inner_join(attr_scale, by = c("K", "attribute")) %>%
  mutate(scaled = truth / scale) %>%
  select(K, class, attribute, scaled)

alignment <- attr_long %>%
  inner_join(attr_scale, by = c("K", "attribute")) %>%
  mutate(scaled = estimate / scale) %>%
  select(n, seed, K, class, attribute, scaled) %>%
  group_by(n, seed, K) %>%
  group_modify(function(d, key) {
    K <- key$K
    if (K < 2L || anyNA(d$scaled)) return(tibble(aligned = NA_integer_))
    fit_m <- d %>% pivot_wider(names_from = attribute, values_from = scaled) %>%
      arrange(class) %>% select(-class) %>% as.matrix()
    b <- bench_vec %>% filter(K == !!K) %>%
      pivot_wider(names_from = attribute, values_from = scaled) %>%
      arrange(class) %>% select(-class, -K) %>% as.matrix()
    # Both matrices must be K x 4 in the same attribute order before they can
    # be differenced; the shape check is the guard, not an assumption.
    if (!identical(dim(fit_m), dim(b)) ||
        !identical(colnames(fit_m), colnames(b))) {
      return(tibble(aligned = NA_integer_))
    }
    perms <- combinat_perms(K)
    cost <- vapply(perms, function(p) sum((fit_m - b[p, , drop = FALSE])^2),
                   numeric(1))
    tibble(aligned = as.integer(identical(perms[[which.min(cost)]], seq_len(K))))
  }) %>%
  ungroup()

explosion <- attr_long %>%
  inner_join(attr_scale, by = c("K", "attribute")) %>%
  group_by(n, seed, K) %>%
  summarise(max_rel_coef = max(abs(estimate) / scale, na.rm = TRUE),
            .groups = "drop")

diag_tbl <- fits %>%
  filter(n != FULL_N) %>%
  # separation is written out as a deliverable below; f_collapse itself comes
  # from the fit-time flag, so min_pair_gap is not joined in here.
  left_join(alignment,  by = c("n", "seed", "K")) %>%
  left_join(explosion,  by = c("n", "seed", "K")) %>%
  mutate(
    # Recorded at fit time by lc_fit_verdict().
    f_no_fit      = has_flag(problems, "no_fit"),
    f_no_cov      = has_flag(problems, "no_covariance"),
    f_wrong_sign  = has_flag(problems, "wrong_sign"),
    f_thin_class  = has_flag(problems, "degenerate_class"),
    f_collapse    = has_flag(problems, "classes_collapsed"),
    # Added here, because both need the n = 388 benchmark to be defined.
    f_explode     = as.integer(ok & !is.na(max_rel_coef) & max_rel_coef > EXPLODE_AT),
    f_misaligned  = as.integer(ok & K > 1L & !is.na(aligned) & aligned == 0L),
    # The composite: usable at fit time AND still comparable to the benchmark.
    # A fit can be perfectly well-behaved on its own terms and still have found
    # a different set of classes, which makes every bias and coverage number
    # computed from it a comparison between different things.
    reproduces = as.integer(
      usable &
        (is.na(max_rel_coef) | max_rel_coef <= EXPLODE_AT) &
        (K == 1L | is.na(aligned) | aligned == 1L)
    )
  )

# --- Failure modes by cell -------------------------------------------------
failure_modes <- diag_tbl %>%
  group_by(n, K) %>%
  summarise(
    reps            = n(),
    obs_per_par     = round(mean(n * 9 / n_par), 1),
    no_fit          = mean(f_no_fit),
    no_covariance   = mean(f_no_cov),
    wrong_sign      = mean(f_wrong_sign),
    thin_class      = mean(f_thin_class),
    class_collapse  = mean(f_collapse),
    exploded        = mean(f_explode),
    misaligned      = mean(f_misaligned),
    usable_at_fit   = mean(usable),
    reproduces      = mean(reproduces),
    .groups = "drop"
  ) %>%
  arrange(desc(n), K)

write_table(failure_modes, "failure_modes", prefix = "11")

cat("\n=== Failure modes by sample size and classes ================\n")
cat("(share of 20 replications exhibiting each pathology)\n\n")
print(as.data.frame(failure_modes %>%
  transmute(N = n, K, `obs/par` = obs_per_par,
            `no fit`      = sprintf("%.0f%%", 100 * no_fit),
            `no cov`      = sprintf("%.0f%%", 100 * no_covariance),
            `wrong sign`  = sprintf("%.0f%%", 100 * wrong_sign),
            `thin class`  = sprintf("%.0f%%", 100 * thin_class),
            `collapsed`   = sprintf("%.0f%%", 100 * class_collapse),
            `exploded`    = sprintf("%.0f%%", 100 * exploded),
            `misaligned`  = sprintf("%.0f%%", 100 * misaligned),
            `usable`      = sprintf("%.0f%%", 100 * usable_at_fit),
            `REPRODUCES`  = sprintf("%.0f%%", 100 * reproduces))),
  row.names = FALSE)

# --- Where each mode first appears -----------------------------------------
# The headline: reading down the sample sizes, at what n does each pathology
# first affect more than a fifth of replications?
onset <- failure_modes %>%
  select(n, K, no_covariance, wrong_sign, thin_class, class_collapse, exploded,
         misaligned) %>%
  pivot_longer(-c(n, K), names_to = "mode", values_to = "rate") %>%
  filter(rate > 0.20) %>%
  group_by(K, mode) %>%
  summarise(first_n = max(n), .groups = "drop") %>%
  pivot_wider(names_from = mode, values_from = first_n)

cat("\n=== Largest N at which each failure exceeds 20% of replications ===\n")
print(as.data.frame(onset), row.names = FALSE)

# --- Quantitative reproduction ---------------------------------------------
# Coverage is the sharpest single number: does the small-sample 95% interval
# contain the full-panel value? At the nominal level it should be 0.95, and a
# collapse below that means the reported uncertainty is a fiction.
accuracy <- params %>%
  filter(n != FULL_N, is_attr_par(parameter)) %>%
  inner_join(benchmark, by = c("K", "parameter")) %>%
  mutate(attribute = str_extract(parameter, "(?<=^b_)[a-z]{2}")) %>%
  inner_join(attr_scale, by = c("K", "attribute")) %>%
  mutate(covered = !is.na(se) & se > 0 &
           truth >= estimate - 1.96 * se & truth <= estimate + 1.96 * se) %>%
  group_by(n, K) %>%
  summarise(coverage = mean(covered, na.rm = TRUE),
            # Normalised by the ATTRIBUTE scale, not by each parameter's own
            # benchmark value. Dividing by the individual truth explodes
            # whenever a class coefficient sits near zero, which produced
            # relative RMSEs in the hundreds that described the denominator
            # rather than the estimate.
            rel_rmse = sqrt(mean(((estimate - truth) / scale)^2, na.rm = TRUE)),
            .groups = "drop")

reproduction <- failure_modes %>%
  select(n, K, reps, obs_per_par, usable_at_fit, reproduces) %>%
  left_join(accuracy, by = c("n", "K")) %>%
  arrange(desc(n), K)

write_table(reproduction, "reproduction", prefix = "11")

cat("\n=== Reproduction against the n = 388 benchmark ==============\n")
print(as.data.frame(reproduction %>%
  transmute(N = n, K, `obs/par` = obs_per_par,
            `usable` = sprintf("%.0f%%", 100 * usable_at_fit),
            `reproduces` = sprintf("%.0f%%", 100 * reproduces),
            `95% coverage` = sprintf("%.2f", coverage),
            `rel RMSE` = sprintf("%.2f", rel_rmse))), row.names = FALSE)

write_table(separation %>% filter(!is.na(min_pair_gap)), "separation", prefix = "11")

# --- Which individual classes survive? --------------------------------------
# "The classes did not reproduce" is too coarse. Some do. For each benchmark
# class, this measures the distance to the CLOSEST class the subsample found,
# in benchmark units, irrespective of labelling. A class that is genuinely
# recoverable will have some fitted class sitting near it in almost every
# replication; one that is an artefact of the full panel will not.
class_recovery <- attr_long %>%
  filter(n != FULL_N) %>%
  inner_join(attr_scale, by = c("K", "attribute")) %>%
  mutate(scaled = estimate / scale) %>%
  select(n, seed, K, class, attribute, scaled) %>%
  group_by(n, seed, K) %>%
  group_modify(function(d, key) {
    Kx <- key$K
    if (Kx < 2L || anyNA(d$scaled)) return(tibble(bench_class = integer(0),
                                                  nearest = numeric(0)))
    fit_m <- d %>% pivot_wider(names_from = attribute, values_from = scaled) %>%
      arrange(class) %>% select(-class) %>% as.matrix()
    b <- bench_vec %>% filter(K == Kx) %>%
      pivot_wider(names_from = attribute, values_from = scaled) %>%
      arrange(class) %>% select(-class, -K) %>% as.matrix()
    if (!identical(dim(fit_m), dim(b))) return(tibble(bench_class = integer(0),
                                                      nearest = numeric(0)))
    # Distance RELATIVE to the benchmark class's own magnitude. An absolute
    # threshold would make the near-zero class a much bigger target -- any
    # badly estimated class shrunk toward the origin lands beside it -- and
    # would manufacture exactly the asymmetry this table is meant to test.
    tibble(bench_class = seq_len(Kx),
           nearest = vapply(seq_len(Kx), function(j) {
             bj <- matrix(b[j, ], nrow(fit_m), ncol(b), byrow = TRUE)
             min(sqrt(rowSums((fit_m - bj)^2))) / sqrt(sum(b[j, ]^2))
           }, numeric(1)))
  }) %>%
  ungroup() %>%
  group_by(n, K, bench_class) %>%
  summarise(median_dist = median(nearest, na.rm = TRUE),
            # Within 50% of the benchmark class's own size.
            recovered   = mean(nearest < 0.5, na.rm = TRUE),
            .groups = "drop")

write_table(class_recovery, "class_recovery", prefix = "11")

cat("
=== Which benchmark classes are recoverable? ================
")
cat("(share of replications with SOME fitted class within 0.5 benchmark units)

")
print(as.data.frame(class_recovery %>%
  filter(K == 3) %>%
  mutate(recovered = sprintf("%3.0f%%", 100 * recovered)) %>%
  pivot_wider(names_from = bench_class, values_from = c(recovered, median_dist),
              names_prefix = "class") %>%
  select(n, starts_with("recovered"))), row.names = FALSE)

# --- What DOES survive: the MNL ---------------------------------------------
# The MNL is the control. If it degrades gracefully where the latent class
# models fall apart, the failure is a property of the latent class structure
# and not simply of having few observations.
vtt <- params %>%
  filter(K == 1, parameter %in% c("b_tt_1", "b_tc_1")) %>%
  select(n, seed, parameter, estimate) %>%
  pivot_wider(names_from = parameter, values_from = estimate) %>%
  mutate(vtt = 60 * b_tt_1 / b_tc_1)

vtt_bench <- vtt %>% filter(n == FULL_N) %>% pull(vtt)

vtt_recovery <- vtt %>%
  filter(n != FULL_N) %>%
  group_by(n) %>%
  summarise(
    reps       = n(),
    median_vtt = median(vtt, na.rm = TRUE),
    q10        = quantile(vtt, 0.10, na.rm = TRUE),
    q90        = quantile(vtt, 0.90, na.rm = TRUE),
    # A ratio estimator with a noisy denominator produces wild values once the
    # cost coefficient approaches zero; the sign flip is the visible symptom.
    n_negative = sum(vtt < 0, na.rm = TRUE),
    within_50pct = mean(abs(vtt / vtt_bench - 1) < 0.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

write_table(vtt_recovery, "vtt_recovery", prefix = "11")

cat(sprintf("\n=== MNL value of travel time (benchmark %.2f CHF/h) =========\n",
            vtt_bench))
print(as.data.frame(vtt_recovery %>%
  transmute(N = n, median = round(median_vtt, 1),
            `10th` = round(q10, 1), `90th` = round(q90, 1),
            `wrong sign` = n_negative,
            `within 50%` = sprintf("%.0f%%", 100 * within_50pct))),
  row.names = FALSE)

# --- Figure ----------------------------------------------------------------
sizes <- sort(unique(diag_tbl$n))
# The log axis itself is log_size_axis() in 00_setup.R; only the set of
# breaks differs between scripts, so this is a binding, not a second copy.
size_axis <- function(p) log_size_axis(p, sizes)

p_repro <- size_axis(
  ggplot(failure_modes, aes(n, reproduces, colour = factor(K))) +
    geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed", colour = "grey70") +
    geom_line() + geom_point(size = 1.7) + ylim(0, 1)) +
  labs(title = "Replications reproducing the full-panel result",
       subtitle = "Right signs, interpretable distinct classes, usable inference",
       x = "respondents (log scale)", y = NULL, colour = "classes")

mode_long <- failure_modes %>%
  select(n, K, `no covariance` = no_covariance, `wrong sign` = wrong_sign,
         `thin class` = thin_class, `class collapse` = class_collapse,
         `misaligned` = misaligned) %>%
  pivot_longer(-c(n, K), names_to = "mode", values_to = "rate")

p_modes <- size_axis(
  ggplot(mode_long %>% filter(K > 1), aes(n, rate, colour = mode)) +
    geom_line() + geom_point(size = 1.4) + ylim(0, 1) +
    facet_wrap(~ paste0(K, " classes"), nrow = 1)) +
  labs(title = "How it fails, not just whether",
       subtitle = "Each pathology appears at a different sample size",
       x = "respondents (log scale)", y = "share of replications", colour = NULL)

p_cov <- size_axis(
  ggplot(accuracy, aes(n, coverage, colour = factor(K))) +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
    geom_line() + geom_point(size = 1.7) + ylim(0, 1)) +
  labs(title = "Coverage of the nominal 95% interval",
       subtitle = "Below the dashed line the reported uncertainty understates the error",
       x = "respondents (log scale)", y = NULL, colour = "classes")

p_vtt <- size_axis(
  ggplot(vtt %>% filter(n != FULL_N), aes(n, vtt, group = n)) +
    geom_hline(yintercept = vtt_bench, linetype = "dashed", colour = "grey30") +
    geom_boxplot(outlier.size = 0.6, width = 0.06) +
    coord_cartesian(ylim = c(-50, 150))) +
  labs(title = "MNL value of travel time by sample size",
       subtitle = sprintf("Dashed line is the full-panel %.1f CHF/h; axis clipped",
                          vtt_bench),
       x = "respondents (log scale)", y = "CHF per hour")

fig <- (p_repro | p_cov) / p_modes / p_vtt + plot_layout(heights = c(1, 1, 1))
fig_path <- write_figure(fig, "breakdown", prefix = "11", width = 15, height = 13)

cat(sprintf("\nTables written to: %s\n", PATH_TABLES))
cat(sprintf("Figure written to: %s\n", fig_path))
