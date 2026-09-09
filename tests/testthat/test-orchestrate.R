# The orchestration layer.
#
# Its whole job is that a finding printed by a script cannot be scrolled past.
# It does that by matching literal strings in each step's log -- fast and
# predictable, and silent when it breaks. Reword a cat() and the marker simply
# stops firing: the step reports ok, the summary says nothing, and the finding
# is lost exactly as it was before the layer existed.
#
# So the patterns are guarded here. Two of them were already dead when this
# file was written, which is the argument for it.

test_that("every marker is still emitted by the file that should emit it", {
  live <- ORCH_MARKERS[!is.na(ORCH_MARKERS$emitter), ]
  expect_gt(nrow(live), 5)

  missing <- character(0)
  for (i in seq_len(nrow(live))) {
    f <- here::here(live$emitter[i])
    if (!file.exists(f)) {
      missing <- c(missing, sprintf("%s (emitter file not found: %s)",
                                    live$pattern[i], live$emitter[i]))
      next
    }
    src <- paste(readLines(f, warn = FALSE), collapse = "\n")
    if (!grepl(live$pattern[i], src, fixed = TRUE)) {
      missing <- c(missing, sprintf("'%s' no longer appears in %s",
                                    live$pattern[i], live$emitter[i]))
    }
  }
  if (length(missing)) {
    fail(paste0(
      "marker(s) no longer emitted -- the orchestrator would report these\n",
      "steps as clean while the finding went unreported:\n  ",
      paste(missing, collapse = "\n  ")))
  }
  expect_length(missing, 0)
})

test_that("no marker is defined without a severity the runner understands", {
  expect_true(all(ORCH_MARKERS$severity %in% names(ORCH_SEVERITY_RANK)))
  expect_false(any(duplicated(ORCH_MARKERS$pattern)))
  expect_true(all(nzchar(ORCH_MARKERS$label)))
})

test_that("the scanner counts hits rather than listing every one", {
  # 08 emits hundreds of UNUSABLE banners; a summary that lists them all is a
  # summary nobody reads.
  tmp <- tempfile(fileext = ".log")
  writeLines(c("routine output",
               "!! UNUSABLE -- LCcov3 -- NO COVARIANCE",
               "more output",
               "!! UNUSABLE -- LCcov4 -- NO COVARIANCE",
               "!! NOT USABLE: LCcov3 (NO COVARIANCE)."), tmp)

  m <- orch_scan_log(tmp)
  # Two banners. "NOT USABLE" does not contain "UNUSABLE" as a substring, so
  # the summary line counts under its own pattern and not this one -- the two
  # markers are genuinely distinct and must stay that way.
  expect_equal(sum(m$n[m$pattern == "UNUSABLE"]), 2)
  expect_equal(m$n[m$pattern == "NOT USABLE:"], 1)
  expect_equal(m$first_line[m$pattern == "NOT USABLE:"], 5)
  expect_true(all(m$severity == "critical"))
  unlink(tmp)
})

test_that("a clean log produces no markers and an ok badge", {
  tmp <- tempfile(fileext = ".log")
  writeLines(c("=== Estimates ===", "b_tt -0.0598", "done"), tmp)
  m <- orch_scan_log(tmp)
  expect_equal(nrow(m), 0)
  expect_equal(orch_badge(m), "ok")
  unlink(tmp)
})

test_that("the badge distinguishes critical from warning", {
  crit <- tibble(severity = "critical", n = 2L)
  warn <- tibble(severity = "warn", n = 1L)
  expect_match(orch_badge(crit), "2 critical", fixed = TRUE)
  expect_match(orch_badge(warn), "1 warning", fixed = TRUE)
  expect_match(orch_badge(bind_rows(crit, warn)), "critical")
  expect_match(orch_badge(bind_rows(crit, warn)), "warning")
})

test_that("a missing log is not an error", {
  # A step can fail before writing anything; the scanner must not then fail too
  # and mask the real failure.
  expect_equal(nrow(orch_scan_log(tempfile())), 0)
})

test_that("cache artefacts resolve to the right directory", {
  # 08 resumes from CSVs under outputs/tables, not from an .rds under
  # outputs/models. Resolving everything against models made 08 look uncached
  # and put a 240-minute estimate on a 12-second step.
  expect_equal(basename(dirname(orch_cache_path("Swiss_MNL_model.rds"))),
               basename(PATH_MODELS))
  expect_equal(basename(dirname(orch_cache_path("08_fits.csv"))),
               basename(PATH_TABLES))
})

test_that("severity policy is what the pipeline documents", {
  # Critical findings must NOT be in the set that fails a run: an inestimable
  # model is a result, and a pipeline that refused to finish over one would be
  # unusable. Only a script error does, plus reproduction drift, which the
  # runner handles separately.
  expect_true("critical" %in% ORCH_MARKERS$severity)
  expect_true("fail" %in% ORCH_MARKERS$severity)
  expect_lt(ORCH_SEVERITY_RANK[["critical"]], ORCH_SEVERITY_RANK[["fail"]])
})
