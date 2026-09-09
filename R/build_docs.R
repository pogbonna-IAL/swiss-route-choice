# build_docs.R -- bundle docs/ into one markdown file
# Swiss route choice
#
# Every page published as an Artifact has a markdown source in this repo, and
# the markdown is the canonical copy -- the page is a rendering of it, not the
# other way round. The reference documentation lives as separate files under
# docs/ because that is how it is edited; this bundles them into one portable
# document so the published reference has a single-file counterpart rather
# than only a set of fragments.
#
# Sourced by 10_report.R, which calls bundle_docs(). Safe to source on its own:
#   Rscript -e 'source("R/00_setup.R"); source("R/build_docs.R"); bundle_docs()'
# ---------------------------------------------------------------------------

# Every page published as an Artifact has a markdown source in this repo, and
# the markdown is the canonical copy -- the page is a rendering of it. The
# reference documentation lives as eight files under docs/ because that is how
# it is edited; this bundles them into one portable document so the published
# reference has a single-file markdown counterpart rather than only a set of
# fragments.
#
# Cross-links between the docs are rewritten to plain anchors, which resolve
# inside the bundle. The character classes below avoid backslash escapes on
# purpose: [.] is a literal dot, []] a literal bracket, [(] a literal paren.
bundle_docs <- function(out_path = file.path(PATH_OUT, "documentation.md")) {
  doc_dir <- here::here("docs")
  files   <- sort(list.files(doc_dir, "^[0-9][0-9]-.*[.]md$", full.names = TRUE))
  if (!length(files)) return(invisible(NULL))

  # Map "01-models.md" -> the anchor of its own first heading, so a bare file
  # link becomes a working in-document link.
  anchor_of <- function(f) {
    h1 <- grep("^# ", readLines(f, warn = FALSE), value = TRUE)[1]
    h1 <- sub("^# +", "", h1)
    a  <- gsub("[^a-z0-9 -]", "", tolower(h1))
    paste0("#", gsub(" +", "-", trimws(a)))
  }
  anchors <- setNames(vapply(files, anchor_of, character(1)), basename(files))

  fix_links <- function(x) {
    # "](02-estimation.md#23-foo)" -> "](#23-foo)"
    x <- gsub("[]][(][0-9][0-9]-[a-z-]+[.]md#", "](#", x)
    # "](01-models.md)" -> "](#1-models)"
    for (nm in names(anchors)) {
      x <- gsub(paste0("[]][(]", sub("[.]md$", "", nm), "[.]md[)]"),
                paste0("](", anchors[[nm]], ")"), x)
    }
    x
  }

  header <- c(
    "# Swiss route choice -- reference documentation",
    "",
    sprintf("_Bundled %s from `docs/`. The individual files are the canonical source;",
            format(Sys.time(), "%Y-%m-%d %H:%M")),
    "this is a single-file copy for reading or export._",
    ""
  )
  index <- fix_links(readLines(file.path(doc_dir, "README.md"), warn = FALSE))

  body <- unlist(lapply(files, function(f) {
    c("", strrep("-", 75), "", fix_links(readLines(f, warn = FALSE)))
  }))

  writeLines(c(header, index, body), out_path)
  invisible(out_path)
}

