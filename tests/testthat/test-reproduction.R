# Does a re-run actually reproduce?
#
# Every other test in this suite checks a pure helper. This one checks the
# PIPELINE: it reads the tables the scripts wrote and compares the headline
# numbers against tests/reference/expected_results.csv.
#
# Without it, every script here could execute cleanly and silently produce
# different numbers -- a changed seed, a package update that moves an
# optimiser, a specification edit nobody noticed. "It ran" and "it reproduced"
# are different claims and only one of them is worth anything.
#
# The tests skip when the outputs are absent, so the suite still passes on a
# fresh clone before the pipeline has been run. Regenerate the reference with
# R/build_reference.R, deliberately, and commit the diff.

ref_path <- here::here("tests", "reference", "expected_results.csv")

test_that("the reference file exists and is well formed", {
  skip_if_not(file.exists(ref_path), "no reference file")
  ref <- read.csv(ref_path)
  expect_true(all(c("quantity", "table", "filter", "column", "value", "tol", "why")
                  %in% names(ref)))
  expect_gt(nrow(ref), 15)
  expect_false(any(duplicated(ref$quantity)))
  # Every reference must carry a stated reason for its tolerance, so a future
  # reader can tell a loose tolerance from a lazy one.
  expect_true(all(nzchar(ref$why)))
})

test_that("the pipeline reproduces its reference results", {
  skip_if_not(file.exists(ref_path), "no reference file")
  ref <- read.csv(ref_path)

  present <- file.exists(file.path(PATH_TABLES, unique(ref$table)))
  skip_if_not(all(present),
              paste("pipeline outputs missing; run Rscript R/run_all.R first --",
                    paste(unique(ref$table)[!present], collapse = ", ")))

  actual <- vapply(seq_len(nrow(ref)), function(i) {
    d <- as_tibble(read.csv(file.path(PATH_TABLES, ref$table[i])))
    v <- d %>%
      dplyr::filter(!!rlang::parse_expr(ref$filter[i])) %>%
      dplyr::pull(!!dplyr::sym(ref$column[i]))
    if (length(v) != 1L) return(NA_real_)
    as.numeric(v)
  }, numeric(1))

  drift <- tibble(quantity = ref$quantity, expected = ref$value,
                  actual = actual, tol = ref$tol) %>%
    mutate(diff = abs(actual - expected),
           ok = !is.na(actual) & diff <= tol)

  # Report every failure at once. Finding out about them one re-run at a time
  # is how a small drift becomes an afternoon.
  if (any(!drift$ok)) {
    bad <- drift %>% filter(!ok)
    fail(paste0(
      sprintf("%d of %d reference values did not reproduce:\n", nrow(bad), nrow(drift)),
      paste(sprintf("  %-26s expected %14.6f  got %14.6f  (diff %.3g > tol %.3g)",
                    bad$quantity, bad$expected, bad$actual, bad$diff, bad$tol),
            collapse = "\n")))
  }
  expect_true(all(drift$ok))
})

test_that("the fit ledger is intact", {
  # These are counts, not estimates: any change means the experiment itself
  # changed, not that an optimiser wandered.
  path <- file.path(PATH_TABLES, "08_fits.csv")
  skip_if_not(file.exists(path), "run 08_small_sample_experiment.R first")
  f <- read.csv(path)

  expect_equal(nrow(f), 1044)
  expect_true(all(f$ok))                       # nothing ever fails to converge
  expect_equal(sum(f$usable), 440)             # but fewer than half are usable
  expect_equal(sum(f$n == FULL_N), 4)          # one benchmark fit per K
  expect_setequal(unique(f$K), 1:4)
  expect_setequal(unique(f$n),
                  c(FULL_N, 250, 150, 100, 75, 50, 30, 25, 20, 15, 12, 10, 8, 5))
})
