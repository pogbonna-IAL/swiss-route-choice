# run_tests.R -- entry point for the test suite
# Swiss route choice
#
#   Rscript tests/run_tests.R
#
# The tests cover the pure helpers only: canonical class ordering, the
# contrast standard errors that ordering implies, the parameter-count
# arithmetic every comparison table depends on, and the hold-out scoring.
# They deliberately estimate nothing -- a suite that takes an hour is a suite
# nobody runs, and the estimation itself is checked inside the scripts by the
# LL(0) identity in 02, the mlogit cross-check in 02, the parameter-count
# assertion in 06 and the restricted-model assertion in 03.
# ---------------------------------------------------------------------------

library(testthat)
res <- test_dir(here::here("tests", "testthat"), stop_on_failure = TRUE)
