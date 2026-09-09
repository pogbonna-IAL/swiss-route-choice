# Code generation. These functions write Apollo's likelihood out as source and
# eval it into the global environment, so a bug here is a silently wrong model
# rather than an error.

test_that("a one-class covariate allocation is refused, not generated", {
  # 2:K runs backwards at K = 1 and would emit code for a class_2 that does
  # not exist. There is nothing to allocate with one class, so the only
  # correct behaviour is to refuse.
  expect_error(lc_make_lcPars_cov(1L))
  expect_error(lc_make_lcPars(1L))
})

test_that("installing a one-class model clears any stale allocation function", {
  lc_install(3L)
  expect_true(exists("apollo_lcPars", envir = globalenv()))
  lc_install(1L)
  expect_false(exists("apollo_lcPars", envir = globalenv()))

  lc_install_cov(3L)
  expect_true(exists("apollo_lcPars", envir = globalenv()))
  lc_install_cov(1L)
  expect_false(exists("apollo_lcPars", envir = globalenv()))
})

test_that("generated allocation code names every class and every covariate", {
  lc_make_lcPars_cov(3L)
  src <- paste(deparse(get("apollo_lcPars", envir = globalenv())), collapse = " ")
  # fixed = TRUE throughout: these are literal code fragments, not patterns.
  expect_true(grepl('V[["class_1"]] <- 0', src, fixed = TRUE))
  for (k in 2:3) {
    expect_match(src, sprintf("delta_%d", k))
    for (cv in unname(LC_COVARS)) expect_match(src, sprintf("g_%s_%d", cv, k))
  }
  rm("apollo_lcPars", envir = globalenv())
})

test_that("the generated one-class model is the plain MNL", {
  lc_install(1L)
  src <- paste(deparse(get("apollo_probabilities", envir = globalenv())),
               collapse = " ")
  expect_match(src, "apollo_mnl")
  expect_false(grepl("apollo_lc(", src, fixed = TRUE))  # no degenerate LC call
  expect_match(src, "b_tt_1")
})
