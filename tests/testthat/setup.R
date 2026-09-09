# Sourced by testthat before any test file.
suppressPackageStartupMessages({
  source(here::here("R", "00_setup.R"))
  source(here::here("R", "lc_helpers.R"))
  source(here::here("R", "model_helpers.R"))
})

# A synthetic 3-class estimate with the classes deliberately OUT of canonical
# order (b_tt is -0.05, -0.20, -0.10), so any test that passes on it would
# fail if the ordering logic were a no-op.
fake_lc3 <- function(covariate_alloc = FALSE) {
  est <- list(
    b_tt_1 = -0.05, b_tt_2 = -0.20, b_tt_3 = -0.10,
    b_tc_1 = -0.30, b_tc_2 = -0.60, b_tc_3 = -0.45,
    b_hw_1 = -0.01, b_hw_2 = -0.02, b_hw_3 = -0.03,
    b_ch_1 = -0.40, b_ch_2 = -0.80, b_ch_3 = -0.60
  )
  if (covariate_alloc) {
    c(est, list(delta_2 = 0.5, delta_3 = -0.2,
                g_inc_2 = 0.3, g_inc_3 = -0.1, g_car_2 = 0.2, g_car_3 = 0.4,
                g_com_2 = -0.5, g_com_3 = 0.1, g_shop_2 = 0.05, g_shop_3 = 0.2,
                g_bus_2 = 0.7, g_bus_3 = -0.3))
  } else {
    c(est, list(delta_1 = 0, delta_2 = 0.5, delta_3 = -0.2))
  }
}

swiss_data <- function() {
  data("apollo_swissRouteChoiceData", package = "apollo",
       envir = environment())
  apollo_swissRouteChoiceData
}
