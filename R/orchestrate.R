# orchestrate.R -- run records, marker scanning and notifications
# Swiss route choice
#
# run_all.R decides WHAT runs and in what order. This decides what happens
# around each step: what is recorded, what is surfaced, and who is told.
#
# The problem it exists to solve: the orchestrator reported process success and
# called it done. Step 05 printed two UNUSABLE banners -- LCcov3 and LCcov4
# have singular Hessians on the full panel, every standard error NA -- into its
# own log, and the summary line said "ok". The scripts already say the right
# things. Nothing read them back.
#
# That is the same "converged is not usable" confusion the fit-level verdict
# fixed inside the scripts, recurring one level up. A step that produces an
# inestimable model has still SUCCEEDED as a step; the finding has to travel
# without failing the run.
#
# Requires 00_setup.R.
# ---------------------------------------------------------------------------

if (!exists("PATH_OUT")) stop("source 00_setup.R before orchestrate.R")

PATH_RUNS <- file.path(PATH_OUT, "runs")
ORCH_KEEP_RUNS <- 10L

# --- What to look for in a step's log ---------------------------------------
# Matched literally, not as regex: these are the exact strings the scripts
# print, and a pattern that quietly stops matching is worse than no pattern.
#
# "critical" does NOT fail the run. An inestimable model is a result, and a
# pipeline that refuses to finish over one would be useless. It fails only on
# "fail" (the script errored) and on reproduction drift.
# `emitter` is the file that must still contain the pattern. A literal match
# is fast and predictable, but it fails SILENTLY: reword a cat() and the
# marker simply stops firing, the step reports ok, and the finding is lost --
# which is the exact failure this whole layer exists to prevent.
# tests/testthat/test-orchestrate.R asserts every pattern is still present in
# its emitter, so a reworded message breaks a test instead of a report.
#
# Two markers were already dead when the guard was written: "!! WARNING" never
# matched (05 prints "WARNING: rows with lr_ok"), and "no usable covariance
# matrix" was removed from 05 and 06 when the verdict system replaced those
# warning() calls. Both had been silently matching nothing.
#
# NA emitter means R itself prints it, not this project.
ORCH_MARKERS <- tibble::tribble(
  ~pattern,                          ~severity,   ~label,                          ~emitter,
  "UNUSABLE",                        "critical",  "unusable model fit",            "R/lc_helpers.R",
  "NOT USABLE:",                     "critical",  "model not estimable",           "R/05_lc_2class.R",
  "NOT ONE usable replication",      "critical",  "cell with no usable fit",       "R/08_small_sample_experiment.R",
  "renv is NOT active",              "critical",  "packages not from renv.lock",   "R/00_setup.R",
  "did not reproduce",               "critical",  "reproduction drift",            "tests/testthat/test-reproduction.R",
  "WARNING: rows with lr_ok",        "warn",      "LR test on a weak optimum",     "R/05_lc_2class.R",
  "NOTE: rows with ratio_reliable",  "warn",      "uninterpretable WTP ratio",     "R/06_lc_multiclass.R",
  "is stale:",                       "warn",      "stale input table",             "R/10_report.R",
  "Execution halted",                "fail",      "script error",                  NA_character_,
  "Error in ",                       "fail",      "script error",                  NA_character_
)

ORCH_SEVERITY_RANK <- c(info = 0L, warn = 1L, critical = 2L, fail = 3L)

# Counts, not one row per hit: 08 emits hundreds of UNUSABLE banners and a
# summary that lists them all is a summary nobody reads.
orch_scan_log <- function(path, markers = ORCH_MARKERS) {
  empty <- tibble(severity = character(), label = character(),
                  pattern = character(), n = integer(),
                  first_line = integer(), example = character())
  if (!file.exists(path)) return(empty)
  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) return(empty)

  out <- lapply(seq_len(nrow(markers)), function(i) {
    hits <- grep(markers$pattern[i], lines, fixed = TRUE)
    if (!length(hits)) return(NULL)
    tibble(severity = markers$severity[i], label = markers$label[i],
           pattern = markers$pattern[i], n = length(hits),
           first_line = hits[1], example = trimws(lines[hits[1]]))
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(empty)
  bind_rows(out)
}

orch_worst <- function(markers) {
  if (!nrow(markers)) return("info")
  names(ORCH_SEVERITY_RANK)[max(ORCH_SEVERITY_RANK[markers$severity])+ 1L]
}

# A compact badge for the step line, so a marker cannot be scrolled past.
orch_badge <- function(markers) {
  if (!nrow(markers)) return("ok")
  crit <- sum(markers$n[markers$severity %in% c("critical", "fail")])
  warn <- sum(markers$n[markers$severity == "warn"])
  paste(c(if (crit) sprintf("!! %d critical", crit),
          if (warn) sprintf("~ %d warning%s", warn, if (warn > 1) "s" else "")),
        collapse = "  ")
}

# --- Run directories --------------------------------------------------------
orch_new_run <- function() {
  d <- file.path(PATH_RUNS, format(Sys.time(), "%Y%m%d-%H%M%S"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Keeps the most recent ORCH_KEEP_RUNS. These are diagnostic history, not
# results -- every number they describe is reproducible from the pipeline.
orch_prune_runs <- function(keep = ORCH_KEEP_RUNS) {
  if (!dir.exists(PATH_RUNS)) return(invisible(0L))
  ds <- sort(list.dirs(PATH_RUNS, recursive = FALSE), decreasing = TRUE)
  if (length(ds) <= keep) return(invisible(0L))
  unlink(ds[-seq_len(keep)], recursive = TRUE)
  invisible(length(ds) - keep)
}

orch_latest_run <- function() {
  if (!dir.exists(PATH_RUNS)) return(NA_character_)
  ds <- sort(list.dirs(PATH_RUNS, recursive = FALSE), decreasing = TRUE)
  if (!length(ds)) NA_character_ else ds[1]
}

# Prepended after the step finishes: system2() owns the stream while it runs,
# so a header written first would be overwritten.
orch_stamp_log <- function(path, step, script, started, elapsed, status) {
  if (!file.exists(path)) return(invisible(NULL))
  body <- readLines(path, warn = FALSE)
  writeLines(c(
    sprintf("# step %s -- %s", step, script),
    sprintf("# started  %s", format(started, "%Y-%m-%d %H:%M:%S")),
    sprintf("# finished %s (%.2f min)",
            format(started + elapsed * 60, "%Y-%m-%d %H:%M:%S"), elapsed),
    sprintf("# exit     %d", status),
    strrep("-", 75), ""), path)
  cat(body, file = path, sep = "\n", append = TRUE)
  invisible(path)
}

# --- Pre-flight -------------------------------------------------------------
# Says what is about to happen before it happens. A six-hour run should not
# start silently, and neither should one that is about to re-estimate
# everything because a cache was deleted.
# A cache artefact may be a fitted model under outputs/models/ or a resumable
# ledger under outputs/tables/ -- 08 resumes from its CSVs, not from an .rds.
# Resolving everything against PATH_MODELS made 08 look uncached and put a
# 240-minute estimate on a step that takes twelve seconds, which triggered the
# long-run warning on a four-minute run. An estimate that cries wolf gets
# ignored, and then it is not an estimate.
orch_cache_path <- function(f) {
  if (grepl("[.]rds$", f)) file.path(PATH_MODELS, f) else file.path(PATH_TABLES, f)
}

orch_preflight <- function(plan) {
  cached <- vapply(seq_len(nrow(plan)), function(i) {
    art <- unlist(plan$cache[[i]])
    if (!length(art)) return(NA)
    all(file.exists(vapply(art, orch_cache_path, character(1))))
  }, logical(1))

  # Estimate, best source first:
  #   1. what the step actually took last time, from the previous run record
  #   2. a flat 5% of the nominal cost when the cache is warm
  #   3. the nominal cost
  # (1) matters because the heuristic is crude -- it put 12 minutes on a step
  # that takes 10 seconds -- and measuring the same machine doing the same
  # work beats any guess. An estimate that cries wolf gets ignored, and then
  # it is not an estimate.
  refit <- !identical(Sys.getenv("REFIT"), "")
  est <- ifelse(!is.na(cached) & cached, plan$minutes * 0.05, plan$minutes)
  measured <- rep(FALSE, nrow(plan))
  if (!refit) {
    last <- orch_latest_run()
    if (!is.na(last) && file.exists(file.path(last, "run.json"))) {
      prev <- try(jsonlite::fromJSON(file.path(last, "run.json")), silent = TRUE)
      if (!inherits(prev, "try-error") && !is.null(prev$steps$step)) {
        m <- match(plan$step, prev$steps$step)
        measured <- !is.na(m)
        est[measured] <- prev$steps$minutes[m[measured]]
      }
    }
  } else {
    est <- plan$minutes
  }
  n_refit <- sum(!is.na(cached) & !cached)

  cat("\n=== Pre-flight =============================================\n")
  print(as.data.frame(tibble(
    step = plan$step, script = plan$script,
    cache = ifelse(is.na(cached), "-", ifelse(cached, "hit", "MISS")),
    `est min` = round(est, 1),
    source = ifelse(measured, "measured", "nominal"))), row.names = FALSE)

  cat(sprintf("\n%d of %d steps will re-estimate; estimated total %.0f min%s\n",
              n_refit, nrow(plan), sum(est),
              if (refit) "  (REFIT=1: every cache ignored)"
              else "  (measured, where a previous run recorded it)"))

  if (sum(est) > 60) {
    cat("This is a long run. Nothing is lost if it is interrupted -- every\n")
    cat("expensive fit is cached, and 'Rscript R/run_all.R --resume' picks up\n")
    cat("from the first step that did not finish.\n")
  }
  invisible(tibble(step = plan$step, cached = cached, est_minutes = est))
}

# --- Notifications ----------------------------------------------------------
# Terminal always. Desktop and webhook are opt-in through the environment, so
# the default behaviour of a clone is unchanged.
#
#   NOTIFY_DESKTOP=1            a toast when the run ends
#   NOTIFY_WEBHOOK=<url>        POST the summary (Slack / Teams shape)
orch_desktop <- function(title, message) {
  if (!identical(Sys.getenv("NOTIFY_DESKTOP"), "1")) return(invisible(FALSE))
  ok <- try(silent = TRUE, {
    if (.Platform$OS.type == "windows") {
      ps <- sprintf(paste0(
        "[reflection.assembly]::loadwithpartialname('System.Windows.Forms')|Out-Null;",
        "$n=New-Object System.Windows.Forms.NotifyIcon;",
        "$n.Icon=[System.Drawing.SystemIcons]::Information;",
        "$n.BalloonTipTitle='%s';$n.BalloonTipText='%s';",
        "$n.Visible=$true;$n.ShowBalloonTip(10000);Start-Sleep -s 6;$n.Dispose()"),
        gsub("'", "", title), gsub("'", "", message))
      system2("powershell", c("-NoProfile", "-Command", shQuote(ps)),
              stdout = NULL, stderr = NULL, wait = FALSE)
    } else if (Sys.info()[["sysname"]] == "Darwin") {
      system2("osascript", c("-e", shQuote(sprintf(
        'display notification "%s" with title "%s"', message, title))))
    } else {
      system2("notify-send", c(shQuote(title), shQuote(message)))
    }
  })
  invisible(!inherits(ok, "try-error"))
}

orch_webhook <- function(title, message, level = "info") {
  url <- Sys.getenv("NOTIFY_WEBHOOK")
  if (!nzchar(url)) return(invisible(FALSE))
  if (!nzchar(Sys.which("curl"))) {
    cat("NOTIFY_WEBHOOK is set but curl was not found; skipping.\n")
    return(invisible(FALSE))
  }
  payload <- jsonlite::toJSON(list(
    text = paste0("*", title, "*\n", message),
    level = level), auto_unbox = TRUE)
  tmp <- tempfile(fileext = ".json")
  writeLines(payload, tmp)
  ok <- system2("curl", c("-s", "-S", "-X", "POST", "-H",
                          shQuote("Content-Type: application/json"),
                          "-d", shQuote(paste0("@", tmp)), shQuote(url)),
                stdout = NULL, stderr = NULL)
  unlink(tmp)
  invisible(ok == 0)
}

orch_notify <- function(title, message, level = "info") {
  orch_desktop(title, message)
  orch_webhook(title, message, level)
  invisible(NULL)
}

# --- The run record ---------------------------------------------------------
orch_write_record <- function(run_dir, results, markers, repro, elapsed) {
  rec <- list(
    started      = format(Sys.time() - elapsed * 60, "%Y-%m-%d %H:%M:%S"),
    finished     = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    total_minutes = round(elapsed, 2),
    refit        = !identical(Sys.getenv("REFIT"), ""),
    r_version    = paste(R.version$major, R.version$minor, sep = "."),
    steps        = results,
    markers      = markers,
    reproduction = repro
  )
  p <- file.path(run_dir, "run.json")
  writeLines(jsonlite::toJSON(rec, auto_unbox = TRUE, pretty = TRUE, null = "null"), p)
  invisible(p)
}

orch_write_summary <- function(run_dir, results, markers, repro, elapsed) {
  l <- c(
    sprintf("# Pipeline run %s", basename(run_dir)), "",
    sprintf("_%s, %.1f minutes, R %s%s_",
            format(Sys.time(), "%Y-%m-%d %H:%M"), elapsed,
            paste(R.version$major, R.version$minor, sep = "."),
            if (!identical(Sys.getenv("REFIT"), "")) ", REFIT=1" else ""),
    "", "## Steps", "",
    "| Step | Script | Minutes | Status | Findings |",
    "|---|---|---:|---|---|",
    sprintf("| %s | `%s` | %.2f | %s | %s |",
            results$step, results$script, results$minutes,
            results$status, results$findings),
    "")

  crit <- markers[markers$severity %in% c("critical", "fail"), , drop = FALSE]
  if (nrow(crit)) {
    l <- c(l, "## Findings that need reading", "",
           "These do not fail the run. An inestimable model is a result, not a",
           "broken pipeline -- but it must not be silently passed over.", "",
           sprintf("- **%s** x%d in `%s` (line %d)  \n  `%s`",
                   crit$label, crit$n, crit$step, crit$first_line, crit$example),
           "")
  }
  warn <- markers[markers$severity == "warn", , drop = FALSE]
  if (nrow(warn)) {
    l <- c(l, "## Warnings", "",
           sprintf("- %s x%d in `%s`", warn$label, warn$n, warn$step), "")
  }
  if (!is.null(repro) && !is.na(repro$ok)) {
    l <- c(l, "## Reproduction", "",
           if (isTRUE(repro$ok)) {
             "All reference values reproduced within tolerance."
           } else {
             paste0("**Reproduction FAILED.** ", repro$detail)
           }, "")
  }
  l <- c(l, "---", "",
         sprintf("Logs: `%s`  ", run_dir),
         "Regenerate everything with `Rscript R/run_all.R`.")
  p <- file.path(run_dir, "summary.md")
  writeLines(l, p)
  invisible(p)
}
