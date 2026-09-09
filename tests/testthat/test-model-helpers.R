# Hold-out splitting and scoring.

test_that("the hold-out split is on respondent, never on task", {
  db <- swiss_data()
  s  <- holdout_split(db)

  # No respondent may appear on both sides: nine tasks split across training
  # and hold-out would let the model see the person it is being tested on.
  expect_length(intersect(unique(s$train$ID), unique(s$test$ID)), 0)
  expect_equal(nrow(s$train) + nrow(s$test), nrow(db))
  expect_equal(n_distinct(s$train$ID) + n_distinct(s$test$ID), n_distinct(db$ID))

  # Whole respondents only.
  tasks <- db %>% count(ID) %>% pull(n) %>% unique()
  expect_equal(s$train %>% count(ID) %>% pull(n) %>% unique(), tasks)
})

test_that("the split is reproducible and seed-dependent", {
  db <- swiss_data()
  expect_equal(holdout_split(db)$train_ids, holdout_split(db)$train_ids)
  expect_false(identical(sort(holdout_split(db, seed = 1L)$train_ids),
                         sort(holdout_split(db, seed = 2L)$train_ids)))
})

test_that("attribute differences are alternative 1 minus alternative 2", {
  db <- swiss_data()
  X  <- attr_diff_matrix(db)
  expect_equal(colnames(X), ROUTE_ATTRS)
  expect_equal(X[, "tt"], db$tt1 - db$tt2)
  expect_equal(X[, "ch"], db$ch1 - db$ch2)
})

test_that("scoring a zero utility difference gives exactly the coin flip", {
  db <- swiss_data()
  s  <- binary_logit_score(rep(0, nrow(db)), db, "null")
  expect_equal(s$LL, -nrow(db) * log(2))
  expect_equal(s$rho2_vs_coin, 0)
  expect_equal(s$share_alt1_predicted, 0.5)
})

test_that("hit rate and recovered share are computed from the same probabilities", {
  db <- data.frame(ID = c(1, 1, 2, 2), choice = c(1L, 2L, 1L, 1L))
  dv <- c(10, 10, -10, 10)          # strongly predicts 1, 1, 2, 1
  s  <- binary_logit_score(dv, db, "toy")
  expect_equal(s$hit_rate, 0.5)     # right on rows 1 and 4 only
  expect_equal(s$observations, 4)
  expect_equal(s$respondents, 2)
  expect_true(s$LL < 0)
})

test_that("the panel score collapses to the closed form when draws are identical", {
  # With no taste variation the simulated panel likelihood must reproduce the
  # analytic binary logit exactly. If it does not, the log-sum-exp averaging
  # is wrong.
  db <- swiss_data() %>% filter(ID %in% head(unique(ID), 20))
  dv <- as.vector(attr_diff_matrix(db) %*% c(-0.06, -0.13, -0.04, -1.15))

  flat <- panel_logit_score(matrix(rep(dv, 5), ncol = 5), db, "flat")
  exact <- binary_logit_score(dv, db, "exact")

  expect_equal(flat$LL, exact$LL)
  expect_equal(flat$hit_rate, exact$hit_rate)
  expect_equal(flat$share_alt1_predicted, exact$share_alt1_predicted)
})

test_that("scoring survives utility differences that overflow the naive formula", {
  # exp(745) is Inf in double precision, so 1/(1+exp(-dv)) returns EXACTLY
  # zero below about dv = -745 and log(0) is -Inf. 04 draws lognormal
  # coefficients whose tails reach well past that, so this is the real
  # operating range, not a contrived one.
  expect_identical(1 / (1 + exp(-(-800))), 0)          # the naive form fails
  expect_true(is.finite(plogis(-800, log.p = TRUE)))   # the used form does not

  db <- swiss_data() %>% filter(ID %in% head(unique(ID), 5))

  for (dv_val in c(-100, -800, -1e5)) {
    dv <- rep(dv_val, nrow(db))
    b  <- binary_logit_score(dv, db, "extreme")
    p  <- panel_logit_score(matrix(rep(dv, 3), ncol = 3), db, "extreme")
    expect_true(is.finite(b$LL), info = paste("binary, dv =", dv_val))
    expect_true(is.finite(p$LL), info = paste("panel, dv =", dv_val))
    expect_true(b$LL < 0)
    expect_true(p$LL < 0)
  }
})

test_that("a respondent with a near-zero sequence probability stays finite", {
  # Nine tasks at p = 1e-300 multiply to 1e-2700, which is zero as a double.
  # Summing logs within respondent before averaging over draws is the only
  # reason the panel score returns a number at all.
  db <- swiss_data() %>% filter(ID %in% head(unique(ID), 5))
  dv <- ifelse(db$choice == 1L, -700, 700)
  s  <- panel_logit_score(matrix(rep(dv, 3), ncol = 3), db, "underflow")
  expect_true(is.finite(s$LL))
  expect_equal(s$LL, sum(rep(-700, nrow(db))), tolerance = 1e-6)
})

test_that("the direct binary logit reproduces the estimated MNL", {
  # 04 uses fit_binary_logit rather than re-plumbing Apollo to get a
  # training-half baseline. If it did not agree with the model 02 estimates,
  # the hold-out comparison in 04 would be against a different model.
  path <- file.path(PATH_MODELS, "Swiss_MNL_model.rds")
  skip_if_not(file.exists(path), "run 02_mnl.R first")

  apollo_b <- readRDS(path)$estimate[paste0("b_", ROUTE_ATTRS)]
  direct   <- fit_binary_logit(swiss_data())

  expect_equal(unname(direct), unname(as.numeric(apollo_b)), tolerance = 1e-5)
})

test_that("the panel score actually uses every draw, not just the first", {
  # Regression test. ifelse() returns a result shaped like its TEST argument,
  # so ifelse(db$choice == 1L, p, 1 - p) on an observations-by-draws matrix
  # silently collapses to a vector holding only draw 1 -- no error, no
  # warning, and a log-likelihood that does not move when draws are added.
  # The identical-draws test above cannot catch it, because there draw 1 is
  # the right answer.
  db <- data.frame(ID = rep(1:2, each = 3), choice = c(1L,1L,1L, 2L,2L,2L))

  # Draw A: p(chosen) = 0.5 for everyone. Draw B: p(chosen) = plogis(4).
  dv <- cbind(rep(0, 6), ifelse(db$choice == 1L, 4, -4))
  s  <- panel_logit_score(dv, db, "two draws")

  # Per respondent: log( (0.5^3 + plogis(4)^3) / 2 ), same for both.
  expected <- 2 * log((0.5^3 + plogis(4)^3) / 2)
  expect_equal(s$LL, expected)

  # And the collapsed-to-draw-1 answer must be a DIFFERENT number, or the
  # test would pass against the bug.
  expect_false(isTRUE(all.equal(s$LL, 2 * log(0.5^3))))
})

test_that("adding draws changes the panel likelihood", {
  # The direct symptom of the collapse bug: the score was byte-identical at
  # 200, 1000 and 4000 draws.
  db <- swiss_data() %>% filter(ID %in% head(unique(ID), 30))
  X  <- attr_diff_matrix(db)
  n  <- n_distinct(db$ID)
  ri <- match(db$ID, sort(unique(db$ID)))

  draw <- function(R, seed) {
    set.seed(seed)
    vapply(seq_len(R), function(r) {
      B <- -exp(matrix(c(-1.3, -0.4, -2.4, 1.2), n, 4, byrow = TRUE) +
                  matrix(rnorm(n * 4), n))
      rowSums(B[ri, , drop = FALSE] * X)
    }, numeric(nrow(db)))
  }

  ll_small <- panel_logit_score(draw(5L,  1L), db, "R=5")$LL
  ll_big   <- panel_logit_score(draw(200L, 1L), db, "R=200")$LL
  expect_false(isTRUE(all.equal(ll_small, ll_big)))
  expect_true(is.finite(ll_small) && is.finite(ll_big))
})

test_that("chosen_sign flips only the alternative-2 rows", {
  db <- data.frame(ID = 1:4, choice = c(1L, 2L, 2L, 1L))
  expect_equal(chosen_sign(db), c(1, -1, -1, 1))
})
