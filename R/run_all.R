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
# and into 1,044 latent class fits already on disk; renaming them would orphan
# all of it to fix a problem that one ordered list solves.
#
# Each script runs in its OWN R process. Apollo keeps its model definition in
# the global environment -- apollo_probabilities, apollo_randCoeff,
# apollo_lcPars -- and a leftover apollo_randCoeff from 04 would silently turn
# 05's latent class model into a mixed logit. Process isolation is the only
# reliable guard.
#
# Expensive fits are cached: 04, 05, 06 and 09 reload completed models and 08
# skips completed cells, so re-running after an interruption resumes rather
# than starting over. Set REFIT=1 to force a full re-estimation.
#
# What happens AROUND each step -- run records, marker scanning, notifications
# -- lives in orchestrate.R.
#
# Usage
#   Rscript R/run_all.R              # everything, in order
#   Rscript R/run_all.R 06 05 07     # just those, still in dependency order
#   Rscript R/run_all.R --list       # show the plan and exit
#   Rscript R/run_all.R --from 06    # 06 and everything after it
#   Rscript R/run_all.R --resume     # restart at the first step that did not finish
#
#   NOTIFY_DESKTOP=1        toast when the run ends
#   NOTIFY_WEBHOOK=<url>    POST the summary to Slack/Teams
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))
source(here::here("R", "orchestrate.R"))

PATH_LOGS <- file.path(PATH_OUT, "logs")
if (!dir.exists(PATH_LOGS)) dir.create(PATH_LOGS, recursive = TRUE)

# `needs`  what the step reads that another step wrote -- checked before it
#          runs, so a partial pipeline fails with a useful message instead of
#          a missing-file error inside someone else's code.
# `cache`  the artefacts that let the step skip its estimation. Pre-flight
#          checks these to say what the run is actually about to cost.
PIPELINE <- tibble::tribble(
  ~step, ~script,                        ~needs,                          ~cache, ~minutes, ~note,
  "01",  "01_data_audit.R",              character(0),                    character(0), 0.5,  "panel structure, missingness, dominance, non-traders",
  "02",  "02_mnl.R",                     character(0),                    character(0), 0.5,  "baseline MNL + mlogit cross-check + hold-out",
  "03",  "03_mnl_covariates.R",          "Swiss_MNL_model.rds",           character(0), 1,    "observed heterogeneity: covariate interactions",
  "04",  "04_mixed_logit.R",             "Swiss_MNL_model.rds",
         c("Swiss_MXL_indep_model.rds", "Swiss_MXL_corr_model.rds",
           "Swiss_MXL_corr_train_model.rds"),                             20,   "continuous unobserved heterogeneity",
  "06",  "06_lc_multiclass.R",           "Swiss_MNL_model.rds",
         c("Swiss_LC1_model.rds", "Swiss_LC2_model.rds", "Swiss_LC3_model.rds",
           "Swiss_LC4_model.rds", "Swiss_LC5_model.rds"),                 45,   "LC1-LC5, 50 starts each",
  "05",  "05_lc_2class.R",               "06_lc_comparison.csv",
         c("Swiss_LCcov2_model.rds", "Swiss_LCcov3_model.rds",
           "Swiss_LCcov4_model.rds"),                                     30,   "covariate class allocation, K = 2..4",
  "08",  "08_small_sample_experiment.R", "Swiss_MNL_model.rds",
         c("08_fits.csv", "08_parameters.csv"),                           240,  "resampling experiment (resumes from CSV)",
  "07",  "07_lc_stability.R",            c("08_fits.csv", "08_parameters.csv",
                                            "06_lc_comparison.csv",
                                            "06_lc_parameters.csv"),      character(0), 0.5,  "stability verdicts, n >= 30",
  "09",  "09_idefix_design.R",           "Swiss_MNL_model.rds",           "idefix_design_D.rds", 5, "Bayesian D-efficient follow-up design",
  "11",  "11_breakdown.R",               c("08_fits.csv",
                                            "08_parameters.csv"),         character(0), 0.5,  "where latent class estimation stops working",
  "10",  "10_report.R",                  "06_lc_comparison.csv",          character(0), 0.5,  "assembles the report and bundles docs/ into markdown"
)

# --- Arguments -------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
plan <- PIPELINE

if (length(args) && identical(args[1], "--resume")) {
  last <- orch_latest_run()
  if (is.na(last) || !file.exists(file.path(last, "run.json"))) {
    stop("nothing to resume: no previous run record in ", PATH_RUNS, call. = FALSE)
  }
  prev <- jsonlite::fromJSON(file.path(last, "run.json"))
  done <- prev$steps$step[prev$steps$status == "ok"]
  plan <- PIPELINE[!PIPELINE$step %in% done, ]
  cat(sprintf("resuming from %s: %d of %d steps already completed\n",
              basename(last), length(done), nrow(PIPELINE)))
  if (!nrow(plan)) { cat("nothing left to run.\n"); quit(save = "no", status = 0) }
} else if (length(args) && !identical(args, "--list")) {
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

if (identical(args, "--list")) {
  cat("\n=== Pipeline ===============================================\n")
  print(as.data.frame(plan %>%
    transmute(step, script, `~min` = minutes, note)), row.names = FALSE)
  quit(save = "no", status = 0)
}

# --- Pre-flight ------------------------------------------------------------
local({
  st <- try(renv::status(), silent = TRUE)
  if (!inherits(st, "try-error") && isFALSE(st$synchronized)) {
    cat("\n!! renv reports the library is OUT OF SYNC with renv.lock.\n")
    cat("!! Run renv::restore() before trusting these results.\n")
  }
})
orch_preflight(plan)

run_dir <- orch_new_run()
cat(sprintf("\nrun record: %s\n", run_dir))

# --- Runner ----------------------------------------------------------------
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows")
  "Rscript.exe" else "Rscript")

find_input <- function(f) {
  if (grepl("\\.rds$", f)) file.path(PATH_MODELS, f) else file.path(PATH_TABLES, f)
}

all_markers <- list()

run_step <- function(step, script, needs, note) {
  needs   <- unlist(needs)
  missing <- needs[!file.exists(vapply(needs, find_input, character(1)))]
  if (length(missing)) {
    stop(sprintf(
      "step %s (%s) needs %s, which does not exist yet.\nRun the earlier steps first: Rscript R/run_all.R",
      step, script, paste(missing, collapse = ", ")), call. = FALSE)
  }

  log_path <- file.path(run_dir, sub("\\.R$", ".log", script))
  cat(sprintf("\n[%s] %s -- %s\n", step, script, note))

  started <- Sys.time()
  status  <- system2(rscript, args = shQuote(here::here("R", script)),
                     stdout = log_path, stderr = log_path)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  orch_stamp_log(log_path, step, script, started, elapsed, status)
  # Also keep the conventional per-script location, so the paths quoted in
  # every error message and doc still resolve.
  file.copy(log_path, file.path(PATH_LOGS, basename(log_path)), overwrite = TRUE)

  markers <- orch_scan_log(log_path)
  if (nrow(markers)) {
    all_markers[[step]] <<- markers %>% mutate(step = step, .before = 1)
  }

  if (status != 0) {
    cat(sprintf("     FAILED after %.1f min (exit %d)\n", elapsed, status))
    # The error, not the last 25 lines -- which are usually Apollo's banner.
    lines <- readLines(log_path, warn = FALSE)
    err <- grep("Error|halted|cannot|unable", lines, value = TRUE)
    show <- if (length(err)) tail(err, 6) else tail(lines, 12)
    cat(paste0("     | ", show, collapse = "\n"), "\n")
    cat(sprintf("     full log: %s\n", log_path))
    orch_notify(sprintf("Pipeline FAILED at step %s", step),
                sprintf("%s exited %d after %.1f min", script, status, elapsed),
                level = "fail")
    stop(sprintf("step %s failed; see %s", step, log_path), call. = FALSE)
  }

  badge <- orch_badge(markers)
  cat(sprintf("     done in %.2f min   %s\n", elapsed, badge))
  if (nrow(markers)) {
    crit <- markers[markers$severity %in% c("critical", "fail"), , drop = FALSE]
    for (i in seq_len(nrow(crit))) {
      cat(sprintf("     -> %s (x%d), log line %d\n",
                  crit$label[i], crit$n[i], crit$first_line[i]))
    }
  }

  tibble(step = step, script = script, minutes = elapsed,
         status = "ok", findings = badge)
}

t_start <- Sys.time()
results <- purrr::pmap_dfr(
  plan %>% select(step, script, needs, note),
  function(step, script, needs, note) run_step(step, script, needs, note)
)
elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

markers <- if (length(all_markers)) bind_rows(all_markers) else
  tibble(step = character(), severity = character(), label = character(),
         pattern = character(), n = integer(), first_line = integer(),
         example = character())

# --- Did it reproduce? -----------------------------------------------------
# Running to completion is not the same as reproducing. Unlike a marker, drift
# here DOES fail the run: it means the pipeline no longer produces the numbers
# this repository claims it produces.
repro <- list(ok = NA, detail = "not checked (partial run)")
if (nrow(plan) == nrow(PIPELINE)) {
  cat("\n=== Reproduction check =====================================\n")
  res <- try(testthat::test_file(
    here::here("tests", "testthat", "test-reproduction.R"),
    reporter = "silent"), silent = TRUE)
  if (inherits(res, "try-error")) {
    repro <- list(ok = FALSE,
                  detail = conditionMessage(attr(res, "condition")))
  } else {
    df <- as.data.frame(res)
    nfail <- sum(df$failed) + sum(df$error)
    repro <- list(ok = nfail == 0,
                  detail = if (nfail == 0) "all reference values within tolerance"
                           else sprintf("%d reference check(s) failed", nfail))
  }
  cat(if (isTRUE(repro$ok)) "REPRODUCED: all reference values within tolerance.\n"
      else sprintf("!! REPRODUCTION FAILED -- %s\n", repro$detail))
}

# --- Record, summarise, notify ---------------------------------------------
info_path <- save_session_info()
orch_write_record(run_dir, results, markers, repro, elapsed)
summary_path <- orch_write_summary(run_dir, results, markers, repro, elapsed)
orch_prune_runs()

cat("\n=== Complete ===============================================\n")
print(as.data.frame(results %>% mutate(minutes = round(minutes, 2))),
      row.names = FALSE)

n_crit <- sum(markers$n[markers$severity %in% c("critical", "fail")])
n_warn <- sum(markers$n[markers$severity == "warn"])
cat(sprintf("\ntotal: %.1f min | %d critical finding(s), %d warning(s)\n",
            elapsed, n_crit, n_warn))
if (n_crit) {
  cat("Critical findings do not fail the run -- an inestimable model is a\n")
  cat("result. They are listed in the summary; read it before citing anything.\n")
}
cat(sprintf("summary:      %s\n", summary_path))
cat(sprintf("run record:   %s\n", file.path(run_dir, "run.json")))
cat(sprintf("session info: %s\n", info_path))

orch_notify(
  sprintf("Pipeline %s in %.1f min",
          if (isTRUE(repro$ok) || is.na(repro$ok)) "complete" else "COMPLETE, REPRODUCTION FAILED",
          elapsed),
  sprintf("%d steps, %d critical, %d warnings. %s",
          nrow(results), n_crit, n_warn, repro$detail),
  level = if (isFALSE(repro$ok)) "fail" else if (n_crit) "critical" else "info")

# Reproduction drift is the one thing that fails the run.
if (isFALSE(repro$ok)) quit(save = "no", status = 1)
