# Parameter counts, names, and the standard errors that canonical ordering
# implies.

test_that("parameter counts match the names actually generated", {
  for (K in 1:5) {
    expect_equal(lc_n_par(K), length(lc_par_names(K)) - (K > 1L),
                 info = paste("K =", K))     # delta_1 is fixed, not free
    expect_equal(lc_n_par_cov(K), length(lc_par_names_cov(K)),
                 info = paste("K =", K))     # class 1 has no parameters at all
  }
})

test_that("the covariate model costs five extra parameters per added class", {
  for (K in 2:5) {
    expect_equal(lc_n_par_cov(K) - lc_n_par(K), (K - 1L) * length(LC_COVARS))
  }
})

test_that("allocation parameters become contrasts, not permutations", {
  # Reordering changes the reference class, so a delta is no longer V_kk but
  # V_kk + V_jj - 2 V_kj. Reporting the permuted s.e. would be wrong whenever
  # the two deltas covary, which they always do.
  V <- matrix(c(0.04, 0.01,
                0.01, 0.09), 2, 2,
              dimnames = list(c("delta_2", "delta_3"), c("delta_2", "delta_3")))
  expect_equal(lc_contrast_se(V, "delta_2", "delta_3"),
               sqrt(0.04 + 0.09 - 2 * 0.01))
  # A structurally-zero parameter contributes nothing.
  expect_equal(lc_contrast_se(V, "delta_2", "delta_1"), sqrt(0.04))
  expect_equal(lc_contrast_se(V, "delta_2", NA), sqrt(0.04))
  # And a contrast with itself is exactly zero, never a small negative under
  # the square root.
  expect_equal(lc_contrast_se(V, "delta_2", "delta_2"), 0)
})

test_that("lc_canonical_params permutes beta s.e. and contrasts the deltas", {
  raw  <- fake_lc3()
  nms  <- names(raw)
  free <- setdiff(nms, "delta_1")
  set.seed(42)
  A <- matrix(rnorm(length(free)^2), length(free))
  V <- crossprod(A) / length(free)
  dimnames(V) <- list(free, free)

  model <- list(estimate = unlist(raw), varcov = V)
  out   <- lc_canonical_params(model, 3L, covariate_alloc = FALSE)
  ord   <- lc_class_order(raw, 3L)

  # Within-class coefficients: pure relabelling, so the s.e. travels with it.
  row <- out[out$parameter == "b_tt_1", ]
  expect_equal(row$se, sqrt(V[paste0("b_tt_", ord[1]), paste0("b_tt_", ord[1])]))

  # Allocation constants: contrast against the new reference class.
  row <- out[out$parameter == "delta_2", ]
  expect_equal(row$se, lc_contrast_se(V, paste0("delta_", ord[2]),
                                      paste0("delta_", ord[1])))
  expect_equal(row$estimate,
               lc_par_value(raw, paste0("delta_", ord[2])) -
                 lc_par_value(raw, paste0("delta_", ord[1])))
})

test_that("uniform starts are reproducible and pin delta_1 at zero", {
  nm <- lc_par_names(3L)
  a <- lc_start_uniform(nm, seed = 7L)
  b <- lc_start_uniform(nm, seed = 7L)
  expect_equal(a, b)
  expect_equal(unname(a[["delta_1"]]), 0)
  expect_false(isTRUE(all.equal(a, lc_start_uniform(nm, seed = 8L))))
})
