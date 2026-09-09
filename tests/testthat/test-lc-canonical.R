# Canonical class ordering.
#
# Latent class labels are identified only up to permutation. 05, 06 and 08 all
# reorder classes before reporting, and the whole point is that reordering
# must not change the model -- only the labels. These tests pin that down.

test_that("classes come back sorted by travel time coefficient", {
  est <- lc_canonical(fake_lc3(), 3L)
  b_tt <- unlist(est[c("b_tt_1", "b_tt_2", "b_tt_3")])
  expect_equal(unname(b_tt), c(-0.20, -0.10, -0.05))
  expect_false(is.unsorted(unname(b_tt)))    # ascending = most negative first
})

test_that("every attribute moves with its class, not independently", {
  raw <- fake_lc3()
  est <- lc_canonical(raw, 3L)
  ord <- lc_class_order(raw, 3L)
  for (a in LC_ATTRS) {
    for (k in 1:3) {
      expect_equal(est[[paste0("b_", a, "_", k)]],
                   raw[[paste0("b_", a, "_", ord[k])]],
                   info = sprintf("attribute %s, class %d", a, k))
    }
  }
})

test_that("constant-only reordering leaves the class shares a permutation", {
  raw <- fake_lc3()
  ord <- lc_class_order(raw, 3L)
  expect_equal(lc_class_shares(lc_canonical(raw, 3L), 3L),
               lc_class_shares(raw, 3L)[ord])
})

test_that("delta_1 is renormalised to zero after reordering", {
  est <- lc_canonical(fake_lc3(), 3L)
  expect_equal(est[["delta_1"]], 0)
})

test_that("covariate reordering leaves membership probabilities unchanged", {
  # This is the substantive claim: re-referencing the allocation index must
  # preserve the likelihood for EVERY respondent, not just the average one.
  raw <- fake_lc3(covariate_alloc = TRUE)
  can <- lc_canonical_cov(raw, 3L)
  ord <- lc_class_order(raw, 3L)

  set.seed(1)
  for (i in 1:25) {
    x <- list(log_income = rnorm(1), car_availability = rbinom(1, 1, 0.5),
              commute = rbinom(1, 1, 0.3), shopping = rbinom(1, 1, 0.1),
              business = rbinom(1, 1, 0.1))
    expect_equal(lc_class_shares_cov(can, 3L, x),
                 lc_class_shares_cov(raw, 3L, x)[ord])
  }
})

test_that("a single class is returned untouched", {
  est <- list(b_tt_1 = -0.06, b_tc_1 = -0.13, b_hw_1 = -0.04, b_ch_1 = -1.15)
  expect_equal(lc_canonical(est, 1L), est)
  expect_equal(lc_canonical_cov(est, 1L), est)
  expect_equal(lc_class_shares(est, 1L), 1)
})

test_that("helpers accept Apollo's named numeric vector, not just a list", {
  # model$estimate is an atomic vector, and est[["absent"]] on an atomic
  # vector is an ERROR rather than NULL. The covariate allocation model has no
  # delta_1 and no g_*_1 at all -- class 1's whole index is normalised to zero
  # -- so looking up an absent name is the normal case, not a mistake.
  vec <- unlist(fake_lc3(covariate_alloc = TRUE))
  expect_type(vec, "double")

  expect_equal(lc_par_value(vec, "delta_1"), 0)
  expect_equal(lc_par_value(vec, "g_inc_1"), 0)
  expect_equal(lc_par_value(vec, "delta_2"), 0.5)
  expect_equal(lc_par_value(list(a = 3), "absent"), 0)

  ord <- lc_class_order(vec, 3L)
  expect_equal(lc_class_shares_cov(lc_canonical_cov(vec, 3L), 3L,
                                   list(log_income = 0.4, car_availability = 1,
                                        commute = 0, shopping = 0, business = 1)),
               lc_class_shares_cov(vec, 3L,
                                   list(log_income = 0.4, car_availability = 1,
                                        commute = 0, shopping = 0, business = 1))[ord])
})

test_that("lc_canonical_params works on a covariate model with no class-1 parameters", {
  raw  <- fake_lc3(covariate_alloc = TRUE)
  vec  <- unlist(raw)
  free <- names(vec)
  set.seed(11)
  A <- matrix(rnorm(length(free)^2), length(free))
  V <- crossprod(A) / length(free)
  dimnames(V) <- list(free, free)

  out <- lc_canonical_params(list(estimate = vec, varcov = V), 3L,
                             covariate_alloc = TRUE)

  # 12 betas + 2 classes x (1 delta + 5 gammas)
  expect_equal(nrow(out), 12 + 2 * (1 + length(LC_COVARS)))
  expect_false(any(is.na(out$estimate)))
  expect_true(all(out$se >= 0))
  # No class-1 allocation parameter is ever reported: it does not exist.
  expect_false(any(out$parameter %in% c("delta_1", "g_inc_1")))
})
