# 00_setup.R -- setup
# Swiss route choice
# Paths, packages, shared constants and IO helpers. Sourced by every other
# script.
#
# Package versions are pinned by renv (see renv.lock). If a library is missing,
# run renv::restore() rather than install.packages().
#
# Only packages that EVERY script needs are attached here. Script-specific
# dependencies are attached where they are used:
#   mlogit  -> 02_mnl.R          (independent cross-check on the MNL)
#   idefix  -> 09_idefix_design.R (efficient design generation)
#   knitr   -> 10_report.R        (kable tables for the written report)
# idefix in particular attaches shiny, which has no business being on the
# search path of a script that only estimates a model.
# ---------------------------------------------------------------------------

# --- Is the pinned library actually in use? ---------------------------------
# If renv fails to activate -- a failed bootstrap, a path too long for
# Windows, .Rprofile not sourced because R was started with --vanilla -- then
# every library() call below silently resolves against whatever the user
# happens to have installed. The scripts still run. The tests still pass. The
# numbers are no longer produced by the pinned environment and nothing says
# so, which is the worst possible failure for a replication package.
local({
  active <- !is.na(Sys.getenv("RENV_PROJECT", unset = NA))
  if (!active) {
    warning(paste0(
      "renv is NOT active: packages are resolving against your personal
      library, not the versions pinned in renv.lock. Results may differ.
      Start R from the project root so .Rprofile runs, and check that
      renv::status() reports a synchronised library."),
      call. = FALSE, immediate. = TRUE)
  }
})

# --- Environment check ------------------------------------------------------
# renv.lock pins R 4.4.1. renv restores the PACKAGES but cannot change the R
# version, and a mismatch surfaces much later as a confusing package error --
# usually inside apollo, an hour into an estimation. Say it here instead.
LOCK_R_VERSION <- "4.4.1"
local({
  have <- paste(R.version$major, R.version$minor, sep = ".")
  if (!identical(have, LOCK_R_VERSION)) {
    packageStartupMessage(sprintf(paste0(
      "NOTE: this project was built and pinned on R %s; you are on R %s.\n",
      "      Results should still reproduce, but if renv::restore() cannot find\n",
      "      binaries for your version you may need to build from source."),
      LOCK_R_VERSION, have))
  }
})

library(here)
library(tidyverse)   # dplyr / tidyr / ggplot2 / readr
library(apollo)      # MNL, mixed logit and latent class estimation
library(patchwork)   # figure composition

# Pin the tidyverse verb. Nothing on the default search path masks select()
# today, but 09 attaches idefix (and with it shiny) and any future package
# could reintroduce a mask, which fails silently and confusingly.
select <- dplyr::select

# --- Paths -----------------------------------------------------------------
# There is no data/ directory: the estimation data ships with apollo and is
# loaded with data("apollo_swissRouteChoiceData"). Nothing is read from disk.
PATH_OUT     <- here("outputs")
PATH_MODELS  <- here("outputs", "models")
PATH_TABLES  <- here("outputs", "tables")
PATH_FIGURES <- here("outputs", "figures")

for (p in c(PATH_OUT, PATH_MODELS, PATH_TABLES, PATH_FIGURES)) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

# --- Shared constants ------------------------------------------------------
# These were previously redeclared in five scripts, which made it possible for
# a seed or a sample size to drift between the script that wrote a table and
# the script that read it.

SEED   <- 20260903L   # every RNG stream in the project derives from this
FULL_N <- 388L        # respondents in apollo_swissRouteChoiceData

# The four route attributes, in the order they appear in every utility
# function, table and figure.
ROUTE_ATTRS <- c("tt", "tc", "hw", "ch")

ROUTE_ATTR_LABELS <- c(tt = "travel time (min)",
                       tc = "travel cost (CHF)",
                       hw = "headway (min)",
                       ch = "interchanges (n)")

# Respondent-level covariates carried by the data.
ROUTE_COVARS <- c("hh_inc_abs", "car_availability",
                  "commute", "shopping", "business", "leisure")

# The covariates that enter a model, mapped to the suffix their coefficients
# carry. Names are the column in `database`, values are the parameter suffix:
# g_tt_inc in 03, g_inc_2 in 05. leisure is the reference trip purpose --
# the four purpose dummies sum to one for every respondent (verified in
# 01_data_audit.R), so including all four alongside a constant would be exact
# collinearity. log_income is centred on its mean by the scripts that use it,
# so a zero covariate vector means "respondent of average income".
MODEL_COVARS <- c(log_income       = "inc",
                  car_availability = "car",
                  commute          = "com",
                  shopping         = "shop",
                  business         = "bus")

# Unlabelled alternatives.
ALTS <- c(1L, 2L)

# --- Respondent archetypes -------------------------------------------------
# 03 reports willingness to pay for a set of archetypal respondents and 05
# reports predicted class membership for the same idea. They had drifted:
# 03 defined a "shopping trip" archetype, 05 defined "high income + car +
# business" instead, and five of seven matched. The report puts those two
# tables three sections apart, so a reader compares rows that are not the same
# rows. One definition, used by both.
#
# inc_sd is passed in because it is a property of the data, not of the
# archetypes -- both scripts compute it over RESPONDENTS, not rows.
respondent_profiles <- function(inc_sd) {
  base <- setNames(as.list(rep(0, length(MODEL_COVARS))), names(MODEL_COVARS))
  with_profile <- function(...) modifyList(base, list(...))
  list(
    `average respondent`           = base,
    `low income (-1 SD)`           = with_profile(log_income = -inc_sd),
    `high income (+1 SD)`          = with_profile(log_income =  inc_sd),
    `has car`                      = with_profile(car_availability = 1),
    `commuter`                     = with_profile(commute = 1),
    `shopping trip`                = with_profile(shopping = 1),
    `business traveller`           = with_profile(business = 1),
    `high income + car + business` = with_profile(log_income = inc_sd,
                                                  car_availability = 1,
                                                  business = 1)
  )
}

# Standard deviation of centred log income across RESPONDENTS. The panel is
# balanced at nine tasks each so a row-level SD coincides here, but that is a
# property of this dataset rather than something the code should rely on.
respondent_income_sd <- function(db) {
  db %>%
    group_by(ID) %>%
    summarise(log_income = first(log_income), .groups = "drop") %>%
    pull(log_income) %>%
    sd()
}

# --- Shared plumbing -------------------------------------------------------
# Apollo is extremely chatty and several scripts call it in a loop. Errors are
# still raised; only routine output is swallowed. lc_helpers.R aliases this as
# lc_quietly so latent-class code reads consistently.
quietly <- function(expr) {
  suppressWarnings(suppressMessages({
    invisible(capture.output(res <- expr))
    res
  }))
}

# Reads a table another script produced, failing with the name of the script
# that produces it rather than a bare file-not-found.
read_required <- function(name, produced_by) {
  path <- file.path(PATH_TABLES, name)
  if (!file.exists(path)) {
    stop("Run ", produced_by, " first -- ", path, " not found.", call. = FALSE)
  }
  as_tibble(read.csv(path))
}

# Sample sizes span 5 to 388, so every figure that plots against n uses a log
# axis with the actual sizes as breaks.
log_size_axis <- function(p, sizes) {
  p + scale_x_continuous(trans = "log10", breaks = sort(unique(sizes)))
}

# --- IO helpers ------------------------------------------------------------
# Every deliverable table goes through write_table() so the naming convention
# (<script prefix>_<name>.csv) is enforced in one place rather than open-coded
# with file.path() in eight scripts.
write_table <- function(x, name, prefix) {
  path <- file.path(PATH_TABLES, sprintf("%s_%s.csv", prefix, name))
  write.csv(x, path, row.names = FALSE)
  invisible(path)
}

write_figure <- function(plot, name, prefix, width, height, dpi = 150) {
  path <- file.path(PATH_FIGURES, sprintf("%s_%s.png", prefix, name))
  ggsave(path, plot, width = width, height = height, dpi = dpi)
  invisible(path)
}

# Rounds every numeric column for console display without touching the object
# that gets written to disk at full precision.
show_table <- function(x, digits = 3) {
  print(as.data.frame(x %>% mutate(across(where(is.numeric), ~ round(.x, digits)))),
        row.names = FALSE)
  invisible(x)
}

# --- Resuming expensive fits -----------------------------------------------
# 04, 05 and 06 each spend tens of minutes estimating, and an interruption
# part-way through used to discard all of it. They now cache each fitted model
# and reload it on the next run. Set REFIT=1 in the environment to ignore
# every cache and re-estimate from scratch.
resume_enabled <- function() !identical(Sys.getenv("REFIT"), "1")

# Caches one fitted model. `key` is the basename under outputs/models/.
cached_model <- function(key, fit_fn, resume = resume_enabled()) {
  path <- file.path(PATH_MODELS, paste0(key, "_model.rds"))
  if (resume && file.exists(path)) {
    m <- readRDS(path)
    cat(sprintf("  cached model reused: %s (LL = %.4f) -- delete the .rds or set REFIT=1 to refit
",
                key, m$maximum))
    return(m)
  }
  m <- fit_fn()
  saveRDS(m, path)
  m
}

# --- Provenance ------------------------------------------------------------
# A replication that cannot say which package versions produced a number is
# not reproducible. run_all.R calls this once at the end of a full pass.
save_session_info <- function(path = file.path(PATH_OUT, "session_info.txt")) {
  writeLines(c(
    sprintf("Swiss route choice -- session recorded %s",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    strrep("-", 75),
    capture.output(sessionInfo()),
    "",
    strrep("-", 75),
    "renv status:",
    capture.output(try(renv::status(), silent = TRUE))
  ), path)
  invisible(path)
}

# --- Global options --------------------------------------------------------
options(stringsAsFactors = FALSE)
set.seed(SEED)

theme_set(theme_minimal(base_size = 11))
