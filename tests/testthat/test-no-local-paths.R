# The cached model objects are committed to a public repository.
#
# Apollo records the absolute output directory of the machine that ran the
# estimation in three separate places inside every fitted model, so without
# scrubbing, this repo would publish a local directory tree and every cached
# object would point at a path that exists on exactly one computer. That is
# both untidy and a direct contradiction of the portability the README claims.

test_that("no committed model cache carries this machine's project path", {
  files <- list.files(PATH_MODELS, "[.]rds$", full.names = TRUE)
  skip_if(length(files) == 0, "no cached models present")

  dirty <- Filter(function(f) has_local_path(readRDS(f)), files)
  expect_equal(basename(dirty), character(0))
})

test_that("the output directory fields are relative, not absolute", {
  path <- file.path(PATH_MODELS, "Swiss_MNL_model.rds")
  skip_if_not(file.exists(path), "run 02_mnl.R first")
  m <- readRDS(path)

  for (field in LOCAL_PATH_FIELDS) {
    v <- tryCatch(m[[field]], error = function(e) NULL)
    if (is.null(v)) next
    expect_false(grepl("^[A-Za-z]:", v),
                 info = paste("absolute path at", paste(field, collapse = "$")))
  }
})

test_that("the checker tests for the real root, not a loose pattern", {
  # Regression: "[A-Za-z]:[/\]" matches binary noise in a serialised object
  # and reported nine already-clean files as dirty.
  expect_false(has_local_path(list(junk = "v:/?Bc4NFuVu}xtK@_uD*iz&Z?3")))
  expect_true(has_local_path(list(p = file.path(here::here(), "outputs"))))
})
