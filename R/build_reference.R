# build_reference.R -- freeze the headline results as a reproduction target
# Swiss route choice
#
# A pipeline that runs to completion is not the same thing as a pipeline that
# reproduces. Without a reference, every script here could execute cleanly and
# silently produce different numbers -- a changed seed, a package update that
# moves an optimiser, a specification edit nobody noticed.
#
# This writes tests/reference/expected_results.csv: a small set of quantities
# that a correct run must recover, with the tolerance each is checked at.
# tests/testthat/test-reproduction.R does the checking.
#
# Run this ONLY when you have deliberately changed what the pipeline should
# produce, and commit the diff so the change is reviewable:
#
#   Rscript R/build_reference.R
#
# Tolerances are per-quantity because the quantities are not alike. Log-
# likelihoods come from an optimiser and agree to about 1e-4; counts are
# exact; a simulated hold-out score depends on draws and is looser.
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

PATH_REF <- here::here("tests", "reference")
if (!dir.exists(PATH_REF)) dir.create(PATH_REF, recursive = TRUE)

get1 <- function(tbl, filter_expr, column) {
  d <- read_required(tbl, "the pipeline")
  v <- d %>% filter(!!rlang::parse_expr(filter_expr)) %>% pull(!!sym(column))
  if (length(v) != 1L) {
    stop(sprintf("%s [%s]$%s matched %d rows, expected 1",
                 tbl, filter_expr, column, length(v)), call. = FALSE)
  }
  as.numeric(v)
}

ref <- function(quantity, tbl, filter_expr, column, tol, why) {
  tibble(quantity = quantity, table = tbl, filter = filter_expr,
         column = column, value = get1(tbl, filter_expr, column),
         tol = tol, why = why)
}

reference <- bind_rows(
  # --- The baseline. If these move, everything downstream has moved. --------
  ref("mnl_ll",        "02_mnl_fit.csv",       'statistic == "LL(final)"', "value",
      1e-4, "globally concave likelihood; must be reproducible exactly"),
  ref("mnl_bic",       "02_mnl_fit.csv",       'statistic == "BIC"', "value",
      1e-3, "derived from LL and a fixed parameter count"),
  ref("mnl_b_tt",      "02_mnl_estimates.csv", 'parameter == "b_tt"', "estimate",
      1e-6, "single maximum, no dependence on starting values"),
  ref("mnl_b_tc",      "02_mnl_estimates.csv", 'parameter == "b_tc"', "estimate",
      1e-6, "as above"),
  ref("mnl_vtt",       "02_mnl_valuations.csv", 'measure == "vtt_chf_per_hour"', "value",
      1e-3, "the headline willingness-to-pay figure"),
  ref("mnl_crosscheck","02_mnl_crosscheck.csv", 'parameter == "LL"', "abs_diff",
      1e-6, "apollo vs mlogit; near zero or the specification differs"),

  # --- Latent class sweep. Multi-start, so tolerance is optimiser-level. ----
  ref("lc2_ll", "06_lc_comparison.csv", 'model == "LC2"', "LL", 1e-3,
      "best of 50 seeded starts; reproducible given the seed"),
  ref("lc3_ll", "06_lc_comparison.csv", 'model == "LC3"', "LL", 1e-3, "as above"),
  ref("lc4_ll", "06_lc_comparison.csv", 'model == "LC4"', "LL", 1e-3, "as above"),
  ref("lc5_ll", "06_lc_comparison.csv", 'model == "LC5"', "LL", 1e-3, "as above"),
  ref("lc4_bic", "06_lc_comparison.csv", 'model == "LC4"', "BIC", 1e-2,
      "LC4 is BIC's pick among the constant-only models"),

  # --- Mixed logit. Simulated ML, so looser. -------------------------------
  ref("mxl_indep_ll", "04_mxl_comparison.csv", 'model == "MXL-I"', "LL", 1e-2,
      "simulated maximum likelihood over 200 MLHS draws"),
  ref("mxl_corr_ll",  "04_mxl_comparison.csv", 'model == "MXL-C"', "LL", 1e-2, "as above"),
  ref("mxl_corr_bic", "04_mxl_comparison.csv", 'model == "MXL-C"', "BIC", 1e-1,
      "the best BIC in the whole project"),

  # --- Covariate models ----------------------------------------------------
  ref("mnlcov_lr",  "03_mnlcov_lrtests.csv", 'comparison == "baseline MNL"', "lr_stat",
      1e-2, "20 df test of observed heterogeneity"),
  ref("lccov2_ll",  "05_lccov_comparison.csv", 'classes == 2', "LL", 1e-3,
      "the only covariate-allocation model with an invertible Hessian"),

  # --- Experimental design -------------------------------------------------
  ref("design_d_error", "09_design_quality.csv",
      'statistic == "Bayesian D-error (chosen design)"', "value", 1e-3,
      "modified Fedorov is serial and seeded, so this is deterministic"),

  # --- The small-sample experiment. Counts, so exact. -----------------------
  ref("small_n_usable_250_k2", "11_reproduction.csv", 'n == 250 & K == 2',
      "usable_at_fit", 1e-9, "share of usable fits; a count over 20 replications"),
  ref("small_n_coverage_250_k2", "11_reproduction.csv", 'n == 250 & K == 2',
      "coverage", 1e-9, "the 54% coverage finding"),
  ref("small_n_usable_30_k4", "11_reproduction.csv", 'n == 30 & K == 4',
      "usable_at_fit", 1e-9, "zero usable fits at 30 respondents, four classes")
)

write.csv(reference, file.path(PATH_REF, "expected_results.csv"), row.names = FALSE)

cat("\n=== Reference results ======================================\n")
show_table(reference %>% select(quantity, value, tol), digits = 6)
cat(sprintf("\nWritten to: %s\n", file.path(PATH_REF, "expected_results.csv")))
cat("Commit the diff so any change to expected output is reviewable.\n")
