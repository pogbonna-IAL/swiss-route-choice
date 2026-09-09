# The shared constants in 00_setup.R are read by scripts that never see the
# data themselves (07 sizes its table against FULL_N). If one drifts from the
# data, several tables become quietly wrong rather than failing.

test_that("FULL_N matches the data", {
  expect_equal(n_distinct(swiss_data()$ID), FULL_N)
})

test_that("every route attribute exists for both alternatives", {
  db <- swiss_data()
  for (a in ROUTE_ATTRS) {
    expect_true(paste0(a, 1) %in% names(db))
    expect_true(paste0(a, 2) %in% names(db))
  }
  expect_setequal(names(ROUTE_ATTR_LABELS), ROUTE_ATTRS)
})

test_that("every model covariate is a column, or is derived by the scripts", {
  db <- names(swiss_data())
  # log_income is constructed by 03 and 05 from hh_inc_abs; the rest are raw.
  expect_true(all(setdiff(names(MODEL_COVARS), "log_income") %in% db))
  expect_true("hh_inc_abs" %in% db)
  expect_true(all(ROUTE_COVARS %in% db))
})

test_that("leisure is excluded from MODEL_COVARS as the reference purpose", {
  # The four purpose dummies sum to one for every respondent (01 verifies
  # this), so including all four alongside a constant is exact collinearity.
  expect_false("leisure" %in% names(MODEL_COVARS))
  expect_true(all(c("commute", "shopping", "business") %in% names(MODEL_COVARS)))
})

test_that("write_table and write_figure enforce the prefix convention", {
  expect_equal(basename(write_table(tibble(a = 1), "unit_test_tmp", "99")),
               "99_unit_test_tmp.csv")
  unlink(file.path(PATH_TABLES, "99_unit_test_tmp.csv"))
})
