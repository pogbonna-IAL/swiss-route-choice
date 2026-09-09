# The usability verdict.
#
# These exist because the pipeline used to print "converged" over fits whose
# Hessian was singular in every replication. Convergence is not usability, and
# the distinction has to be enforced somewhere testable.

mk_fit <- function(est, hessian_ok = TRUE, ok = TRUE) {
  list(ok = ok, hessian_ok = hessian_ok,
       model = if (ok) list(estimate = unlist(est)) else NULL)
}

good_2class <- list(
  b_tt_1 = -0.05, b_tt_2 = -0.20, b_tc_1 = -0.30, b_tc_2 = -0.60,
  b_hw_1 = -0.01, b_hw_2 = -0.03, b_ch_1 = -0.40, b_ch_2 = -1.10
)

test_that("a healthy fit is usable", {
  v <- lc_fit_verdict(mk_fit(good_2class), 2L, shares = c(0.45, 0.55))
  expect_true(v$usable)
  expect_equal(v$verdict, "ok")
  expect_length(v$flags, 0)
})

test_that("convergence without a covariance matrix is NOT usable", {
  # The whole point. The optimiser succeeded; every standard error is
  # undefined; the old code called this "converged" and moved on.
  v <- lc_fit_verdict(mk_fit(good_2class, hessian_ok = FALSE), 2L,
                      shares = c(0.45, 0.55))
  expect_false(v$usable)
  expect_equal(v$verdict, "NO COVARIANCE")
  expect_true("no_covariance" %in% v$flags)
  expect_match(v$detail, "UNDEFINED")
})

test_that("a wrong-signed attribute coefficient is caught", {
  bad <- good_2class; bad$b_hw_2 <- 0.02
  v <- lc_fit_verdict(mk_fit(bad), 2L, shares = c(0.45, 0.55))
  expect_false(v$usable)
  expect_true("wrong_sign" %in% v$flags)
})

test_that("a class too thin to interpret is caught", {
  v <- lc_fit_verdict(mk_fit(good_2class), 2L, shares = c(0.02, 0.98))
  expect_false(v$usable)
  expect_true("degenerate_class" %in% v$flags)
  # And is not raised when the class is merely small but readable.
  expect_true(lc_fit_verdict(mk_fit(good_2class), 2L,
                             shares = c(0.12, 0.88))$usable)
})

test_that("two classes collapsed onto one are caught", {
  coll <- list(b_tt_1 = -0.100, b_tt_2 = -0.101, b_tc_1 = -0.300, b_tc_2 = -0.302,
               b_hw_1 = -0.020, b_hw_2 = -0.0201, b_ch_1 = -0.800, b_ch_2 = -0.805)
  v <- lc_fit_verdict(mk_fit(coll), 2L, shares = c(0.45, 0.55))
  expect_false(v$usable)
  expect_true("classes_collapsed" %in% v$flags)
})

test_that("class separation is scale-free", {
  # b_ch is ~30x b_hw in magnitude. A raw distance would be dominated by
  # interchanges and would call two well-separated classes identical whenever
  # their interchange coefficients happened to agree.
  m <- lc_min_class_gap(unlist(good_2class), 2L)
  expect_true(m > LC_COLLAPSE_REL)

  scaled_up <- lapply(good_2class, function(x) x * 1000)
  expect_equal(lc_min_class_gap(unlist(scaled_up), 2L), m)
})

test_that("a fit that never converged reports NO FIT", {
  v <- lc_fit_verdict(mk_fit(NULL, ok = FALSE), 2L)
  expect_false(v$usable)
  expect_equal(v$verdict, "NO FIT")
})

test_that("the most severe problem leads, and every problem is listed", {
  bad <- good_2class; bad$b_hw_2 <- 0.02
  v <- lc_fit_verdict(mk_fit(bad, hessian_ok = FALSE), 2L, shares = c(0.01, 0.99))
  expect_equal(v$verdict, "NO COVARIANCE")
  expect_setequal(v$flags, c("no_covariance", "wrong_sign", "degenerate_class"))
})

test_that("a one-class model skips the class checks", {
  v <- lc_fit_verdict(mk_fit(list(b_tt_1 = -0.06, b_tc_1 = -0.13,
                                  b_hw_1 = -0.04, b_ch_1 = -1.15)), 1L)
  expect_true(v$usable)
  expect_true(is.na(lc_min_class_gap(c(b_tt_1 = -0.06), 1L)))
})

test_that("the announcement is greppable and unmissable", {
  v <- lc_fit_verdict(mk_fit(good_2class, hessian_ok = FALSE), 2L,
                      shares = c(0.45, 0.55))
  out <- capture.output(lc_announce_verdict(v, "LC2"))
  expect_true(any(grepl("UNUSABLE", out, fixed = TRUE)))
  expect_true(any(grepl("!!!!!!", out, fixed = TRUE)))
  expect_true(any(grepl("not a result", out, fixed = TRUE)))

  quiet <- capture.output(lc_announce_verdict(
    lc_fit_verdict(mk_fit(good_2class), 2L, shares = c(0.45, 0.55)), "LC2"))
  expect_false(any(grepl("UNUSABLE", quiet, fixed = TRUE)))
})
