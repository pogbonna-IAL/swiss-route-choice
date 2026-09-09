# run_all.R -- run the pipeline in dependency order
# Swiss route choice
#
# THE FILE NUMBERING IS NOT THE RUN ORDER. Two scripts read tables that a
# higher-numbered script writes:
#
#   05_lc_2class.R          reads 06_lc_comparison.csv   -> must run after 06
#   07_lc_stability.R       reads 08_fits.csv, 08_parameters.csv,
#                                 06_lc_comparison.csv, 06_lc_parameters.csv
#                                                        -> must run after 06 and 08
#
# Running the scripts in filename order fails at 05. The numbering is kept
# because it is baked into every output filename (05_lccov_*.csv and so on)
# and into 485 latent class fits already on disk; renaming them would orphan
# all of it to fix a problem that one ordered list solves.
#
# Each script runs in its OWN R process. Apollo keeps its model definition in
# the global environment -- apollo_probabilities, apollo_randCoeff,
# apollo_lcPars -- and a leftover apollo_randCoeff from 04 would silently turn
# 05's latent class model into a mixed logit. Process isolation is the only
# reliable guard.
#
# Expensive fits are cached: 04, 05 and 06 reload completed models from
# outputs/models/ and 08 skips completed cells, so re-running after an
# interruption resumes rather than starting over. Set REFIT=1 to force a full
# re-estimation.
#
# Usage
#   Rscript R/run_all.R              # everything, in order
#   Rscript R/run_all.R 06 05 07     # just those, still in dependency order
#   Rscript R/run_all.R --list       # show the plan and exit
#   Rscript R/run_all.R --from 06    # 06 and everything after it
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

PATH_LOGS <- file.path(PATH_OUT, "logs")
if (!dir.exists(PATH_LOGS)) dir.create(PATH_LOGS, recursive = TRUE)

# `needs` is what the script reads that another script wrote. It is the reason
# the order is what it is, and it is checked before each script runs so a
# partial pipeline fails with a useful message instead of a missing-file error
# forty lines into someone else's code.
PIPELINE <- tibble::tribble(
  ~step, ~script,                        ~needs,                          ~minutes, ~note,
  "01",  "01_data_audit.R",              character(0),                    0.5,  "panel structure, missingness, dominance, non-traders",
  "02",  "02_mnl.R",                     character(0),                    0.5,  "baseline MNL + mlogit cross-check + hold-out",
  "03",  "03_mnl_covariates.R",          "Swiss_MNL_model.rds",           1,    "observed heterogeneity: covariate interactions",
  "04",  "04_mixed_logit.R",             "Swiss_MNL_model.rds",           20,   "continuous unobserved heterogeneity",
  "06",  "06_lc_multiclass.R",           "Swiss_MNL_model.rds",           45,   "LC1-LC5, 50 starts each",
  "05",  "05_lc_2class.R",               "06_lc_comparison.csv",          30,   "covariate class allocation, K = 2..4",
  "08",  "08_small_sample_experiment.R", "Swiss_MNL_model.rds",           240,  "resampling experiment (resumes from CSV)",
  "07",  "07_lc_stability.R",            c("08_fits.csv", "08_parameters.csv",
                                            "06_lc_comparison.csv",
                                            "06_lc_parameters.csv"),        0.5,  "stability verdicts, n >= 30",
  "09",  "09_idefix_design.R",           "Swiss_MNL_model.rds",           5,    "Bayesian D-efficient follow-up design",
  "11",  "11_breakdown.R",               c("08_fits.csv",
                                            "08_parameters.csv"),           0.5,  "where latent class estimation stops working",
  "10",  "10_report.R",                  "06_lc_comparison.csv",          0.5,  "assembles the report and bundles docs/ into markdown"
)

# --- Arguments -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

plan <- PIPELINE
if (length(args) > 0 && !identical(args, "--list")) {
  if (args[1] == "--from") {
    if (length(args) < 2 || !args[2] %in% PIPELINE$step) {
      stop("--from needs a step: ", paste(PIPELINE$step, collapse = " "))
    }
    plan <- PIPELINE[which(PIPELINE$step == args[2]):nrow(PIPELINE), ]
  } else {
    unknown <- setdiff(args, PIPELINE$step)
    if (length(unknown)) {
      stop("unknown step(s): ", paste(unknown, collapse = ", "),
           "\nknown steps: ", paste(PIPELINE$step, collapse = " "))
    }
    # Filter, never reorder: the point of this file is the order.
    plan <- PIPELINE[PIPELINE$step %in% args, ]
  }
}

cat("\n=== Pipeline ===============================================\n")
print(as.data.frame(plan %>%
  transmute(step, script, `~min` = minutes, note)), row.names = FALSE)

if (identical(args, "--list")) quit(save = "no", status = 0)

# --- Runner ----------------------------------------------------------------
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows")
  "Rscript.exe" else "Rscript")

find_input <- function(f) {
  if (grepl("\\.rds$", f)) file.path(PATH_MODELS, f) else file.path(PATH_TABLES, f)
}

run_step <- function(step, script, needs, note) {
  needs   <- unlist(needs)   # tribble gives a list column once any cell is a vector
  missing <- needs[!file.exists(vapply(needs, find_input, character(1)))]
  if (length(missing)) {
    stop(sprintf(
      "step %s (%s) needs %s, which does not exist yet.\nRun the earlier steps first: Rscript R/run_all.R",
      step, script, paste(missing, collapse = ", ")), call. = FALSE)
  }

  log_path <- file.path(PATH_LOGS, sub("\\.R$", ".log", script))
  cat(sprintf("\n[%s] %s -- %s\n", step, script, note))
  cat(sprintf("     log: %s\n", log_path))

  t0 <- Sys.time()
  status <- system2(rscript, args = shQuote(here::here("R", script)),
                    stdout = log_path, stderr = log_path)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  if (status != 0) {
    cat(sprintf("     FAILED after %.1f min (exit %d)\n", elapsed, status))
    cat("     last 25 lines of the log:\n")
    tail_lines <- tail(readLines(log_path, warn = FALSE), 25)
    cat(paste0("     | ", tail_lines, collapse = "\n"), "\n")
    stop(sprintf("step %s failed; see %s", step, log_path), call. = FALSE)
  }

  cat(sprintf("     done in %.1f min\n", elapsed))
  tibble(step = step, script = script, minutes = elapsed, status = "ok")
}

# A pipeline run is expensive; an out-of-sync library is cheap to detect. Check
# before spending the time rather than after.
local({
  st <- try(renv::status(), silent = TRUE)
  if (!inherits(st, "try-error") && isFALSE(st$synchronized)) {
    cat("
!! renv reports the library is OUT OF SYNC with renv.lock.
")
    cat("!! Run renv::restore() before trusting these results.

")
  }
})

t_start <- Sys.time()
results <- purrr::pmap_dfr(
  plan %>% select(step, script, needs, note),
  function(step, script, needs, note) run_step(step, script, needs, note)
)

# --- Provenance ------------------------------------------------------------
info_path <- save_session_info()

# --- Did it reproduce? -----------------------------------------------------
# Running to completion is not the same as reproducing. Only meaningful after
# a full pass, so it is skipped when a subset was requested.
if (nrow(plan) == nrow(PIPELINE)) {
  cat("\n=== Reproduction check =====================================\n")
  res <- try(testthat::test_file(
    here::here("tests", "testthat", "test-reproduction.R"),
    reporter = "summary"), silent = TRUE)
  if (inherits(res, "try-error")) {
    cat("could not run the reproduction check: ",
        conditionMessage(attr(res, "condition")), "\n", sep = "")
  }
}

cat("\n=== Complete ===============================================\n")
print(as.data.frame(results %>% mutate(minutes = round(minutes, 2))),
      row.names = FALSE)
cat(sprintf("\ntotal: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat(sprintf("session info: %s\n", info_path))
cat(sprintf("logs:         %s\n", PATH_LOGS))
