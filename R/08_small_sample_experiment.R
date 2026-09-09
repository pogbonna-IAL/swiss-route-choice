# 08_small_sample_experiment.R -- how estimation degrades with sample size
# Swiss route choice
#
# Resamples respondents from the full 388-person panel at a range of sample
# sizes and re-estimates the MNL and the 2-, 3- and 4-class latent class
# models in every draw. The questions it answers:
#
#   1. How fast does coefficient and WTP precision decay as n falls?
#   2. At what n does latent class estimation stop converging reliably?
#   3. Does BIC under-select the number of classes in small samples?
#
# Design notes
#   * Sampling is on ID, never on rows, so each respondent keeps all 9 tasks.
#   * n = 388 is NOT a sampled cell: drawing 388 from 388 returns the whole
#     panel for every seed, so all seeds would be identical. It is estimated
#     once as the benchmark that bias and coverage are measured against.
#   * Latent class labels are identified only up to permutation, so every
#     fitted model is put in canonical class order (see lc_canonical) before
#     anything is averaged across replications.
#   * Results are appended to CSV after each cell and completed cells are
#     skipped on restart, so a multi-hour run survives interruption.
#
# RUN ORDER: after 02_mnl.R, before 07_lc_stability.R (which reads the two
# per-fit tables written here). See run_all.R.
#
# Outputs
#   outputs/tables/08_fits.csv          one row per (n, seed, K)
#   outputs/tables/08_parameters.csv    long-format estimates and s.e.
#   outputs/tables/08_summary.csv       bias / RMSE / coverage by n and K
#   outputs/tables/08_selection.csv     BIC-selected K by n
#   outputs/figures/08_small_sample.png
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))
source(here::here("R", "lc_helpers.R"))

# FULL_N and SEED come from 00_setup.R; they used to be redeclared here and in
# 07, which made it possible for the script that wrote a table and the script
# that read it to disagree about the size of the full panel.
# Two regimes. The upper block is the plausible-study range; the lower block
# is deliberately past the point of usefulness, to locate where latent class
# estimation stops working rather than merely getting noisy.
#
# For a given seed these samples are NESTED -- sample(ids, 5) is the first
# five elements of sample(ids, 250) under the same seed -- so shrinking n
# removes respondents from a fixed ordering rather than drawing an unrelated
# group. "Where does it break" is therefore a question about the same people,
# progressively fewer, which is much cleaner than comparing disjoint samples.
SAMPLE_SIZES <- c(250L, 150L, 100L, 75L, 50L, 30L,
                  25L, 20L, 15L, 12L, 10L, 8L, 5L)
SEEDS        <- 1:20
K_SET        <- 1:4        # MNL, LC2, LC3, LC4
N_STARTS     <- 8L         # perturbed starts per latent class fit

PATH_FITS   <- file.path(PATH_TABLES, "08_fits.csv")
PATH_PARAMS <- file.path(PATH_TABLES, "08_parameters.csv")

# --- One-time schema migration ---------------------------------------------
# Earlier runs recorded `ok` and `hessian_ok` but no usability verdict, and
# 08_fits.csv represents many hours of estimation that must not be thrown
# away to add three derived columns. This backfills them from the estimates
# already on disk and is a no-op once done. Nothing is invented: the verdict
# is recomputed from the same records a fresh run would use.
migrate_fits_schema <- function(fits_path, params_path) {
  if (!file.exists(fits_path)) return(invisible(FALSE))
  fits <- as_tibble(read.csv(fits_path))
  if (all(c("usable", "verdict", "problems") %in% names(fits))) {
    return(invisible(FALSE))
  }
  cat("migrating 08_fits.csv to the usability schema ... ")
  params <- as_tibble(read.csv(params_path))

  verdicts <- fits %>%
    select(n, seed, K, ok, hessian_ok, min_class_share) %>%
    pmap_dfr(function(n, seed, K, ok, hessian_ok, min_class_share) {
      est <- params %>%
        filter(n == !!n, seed == !!seed, K == !!K) %>%
        select(parameter, estimate) %>%
        deframe()
      fit <- list(ok = ok, hessian_ok = hessian_ok,
                  model = if (ok) list(estimate = est) else NULL)
      v <- lc_fit_verdict(
        fit, K,
        shares = if (K > 1L && !is.na(min_class_share)) {
          c(min_class_share, 1 - min_class_share)   # only the minimum matters
        } else NULL)
      tibble(usable = v$usable, verdict = v$verdict,
             problems = paste(v$flags, collapse = ";"))
    })

  out <- bind_cols(fits, verdicts)
  write.csv(out, fits_path, row.names = FALSE)
  cat(sprintf("done (%d rows, %d now marked unusable)\n",
              nrow(out), sum(!out$usable)))
  invisible(TRUE)
}

migrate_fits_schema(PATH_FITS, PATH_PARAMS)

apollo_initialise()
data("apollo_swissRouteChoiceData", package = "apollo")
full_database <- apollo_swissRouteChoiceData

mnl_path <- file.path(PATH_MODELS, "Swiss_MNL_model.rds")
if (!file.exists(mnl_path)) stop("Run 02_mnl.R first -- ", mnl_path, " not found.")
b_ref <- readRDS(mnl_path)$estimate

# --- Estimation -------------------------------------------------------------
# Apollo reads its inputs from the global environment, so each fit installs its
# objects there rather than passing them down as arguments.
#
# The multi-start loop itself now lives in lc_helpers.R (lc_estimate_best) and
# is shared with 05 and 06. Only the parts specific to this experiment stay
# here: swapping in the resampled database and keeping the throwaway fits out
# of outputs/models. The starting-value seeds are unchanged, so fits already
# recorded in 08_fits.csv remain exactly reproducible.
fit_lc <- function(db, K, seed_base, n_starts = N_STARTS) {
  lc_install(K)
  assign("database", db, envir = globalenv())

  starts <- lapply(seq_len(if (K == 1L) 1L else n_starts), function(draw) {
    lc_start_values(K, draw, b_ref, seed_base + draw)
  })

  lc_estimate_best(
    starts       = starts,
    apollo_fixed = if (K == 1L) c() else c("delta_1"),
    control = list(
      modelName       = sprintf("tmp_LC%d", K),
      modelDescr      = "small sample experiment",
      indivID         = "ID",
      outputDirectory = tempdir(),  # throwaway fits, keep outputs/models clean
      panelData       = TRUE,
      seed            = SEED,
      nCores          = 1
    ),
    verbose = FALSE
  )
}

# Turns one fitted model into the two rows-of-record it contributes.
summarise_fit <- function(res, n, seed, K, elapsed) {
  if (!res$ok) {
    return(list(
      fit = tibble(n = n, seed = seed, K = K, ok = FALSE, n_starts_ok = 0L,
                   LL = NA_real_, n_par = lc_n_par(K), AIC = NA_real_,
                   BIC = NA_real_, min_class_share = NA_real_,
                   hessian_ok = FALSE, secs = elapsed,
                   usable = FALSE, verdict = "NO FIT", problems = "no_fit"),
      par = tibble()
    ))
  }

  m     <- res$model
  est   <- lc_canonical(m$estimate, K)   # fix label switching before recording
  n_par <- lc_n_par(K)
  ll    <- m$maximum
  pi_k  <- lc_class_shares(est, K)

  # NOT named `verdict`: tibble() evaluates its arguments in order and a
  # column defined earlier shadows a like-named object for every argument
  # after it, so `verdict = verdict$verdict` would leave the next line
  # reading `$flags` off a character vector.
  fit_verdict <- lc_fit_verdict(res, K, shares = if (K > 1L) pi_k else NULL)

  se <- rep(NA_real_, length(est))
  names(se) <- names(est)
  if (res$hessian_ok) {
    v <- sqrt(diag(m$varcov))
    # Reorder the standard errors the same way the estimates were reordered.
    ord_names <- names(est)
    common <- intersect(ord_names, names(v))
    se[common] <- v[common]
  }

  list(
    fit = tibble(
      n = n, seed = seed, K = K, ok = TRUE, n_starts_ok = res$n_ok,
      LL = ll, n_par = n_par,
      AIC = -2 * ll + 2 * n_par,
      BIC = -2 * ll + n_par * log(m$nObs),
      min_class_share = min(pi_k),
      hessian_ok = res$hessian_ok, secs = elapsed,
      # The optimiser returning a finite likelihood is NOT the same thing as
      # a fit anyone can use. `ok` records the former; `usable` records the
      # latter, and they diverge sharply below about 30 respondents.
      usable = fit_verdict$usable, verdict = fit_verdict$verdict,
      problems = paste(fit_verdict$flags, collapse = ";")
    ),
    par = tibble(
      n = n, seed = seed, K = K,
      parameter = names(est),
      estimate  = as.numeric(unlist(est)),
      se        = as.numeric(se[names(est)])
    )
  )
}

append_csv <- function(x, path) {
  if (nrow(x) == 0) return(invisible(NULL))
  # Appending a data frame whose columns differ from the existing header
  # writes rows that silently misalign with it -- every later read of the
  # file is then wrong in a way nothing detects. Refuse instead.
  if (file.exists(path)) {
    existing <- names(read.csv(path, nrows = 1))
    if (!identical(existing, names(x))) {
      stop(sprintf(paste0(
        "schema mismatch for %s.\n",
        "  on disk: %s\n",
        "  new:     %s\n",
        "Delete the file to re-run from scratch, or migrate it first."),
        basename(path), paste(existing, collapse = ", "),
        paste(names(x), collapse = ", ")), call. = FALSE)
    }
  }
  write.table(x, path, sep = ",", row.names = FALSE,
              col.names = !file.exists(path), append = file.exists(path),
              qmethod = "double")
}

# --- Cell list --------------------------------------------------------------
# The benchmark carries seed 0 and is estimated exactly once.
cells <- bind_rows(
  tibble(n = FULL_N, seed = 0L),
  expand_grid(n = SAMPLE_SIZES, seed = SEEDS)
)

# Resume is tracked per (n, seed, K), not per cell. Keying on the cell alone
# would mean that widening K_SET -- adding LC4 to an experiment that already
# ran MNL/LC2/LC3 -- marks every existing cell incomplete and re-estimates
# hundreds of models that are already correct and on disk.
done <- if (file.exists(PATH_FITS)) {
  read.csv(PATH_FITS) %>% as_tibble() %>% distinct(n, seed, K)
} else {
  tibble(n = integer(), seed = integer(), K = integer())
}

todo <- cells %>%
  expand_grid(K = K_SET) %>%
  anti_join(done, by = c("n", "seed", "K")) %>%
  arrange(desc(n), seed, K)

todo_cells <- todo %>% distinct(n, seed)

cat(sprintf("fits total: %d | already done: %d | to run: %d (across %d cells)\n",
            nrow(cells) * length(K_SET), nrow(done), nrow(todo),
            nrow(todo_cells)))

# --- Main loop --------------------------------------------------------------
ids <- unique(full_database$ID)
t_start <- Sys.time()

for (i in seq_len(nrow(todo_cells))) {
  n    <- todo_cells$n[i]
  seed <- todo_cells$seed[i]

  if (n == FULL_N) {
    db <- full_database
  } else {
    # Drawn from the same seed as the original run, so a cell re-entered to add
    # a model gets exactly the same respondents as its existing fits.
    set.seed(SEED + seed)
    db <- full_database %>% filter(ID %in% sample(ids, n))
  }

  ks <- todo %>% filter(n == !!n, seed == !!seed) %>% pull(K)

  cat(sprintf("\n[%3d/%3d] n = %3d  seed = %2d  (%d rows)  K: %s\n",
              i, nrow(todo_cells), n, seed, nrow(db),
              paste(ks, collapse = ",")))

  for (K in ks) {
    t0  <- Sys.time()
    res <- fit_lc(db, K, seed_base = SEED + 1000L * K + seed)
    el  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    out <- summarise_fit(res, n, seed, K, el)
    append_csv(out$fit, PATH_FITS)
    append_csv(out$par, PATH_PARAMS)

    cat(sprintf("   K=%d  %-18s LL = %11s  starts ok %d  %5.1fs\n", K,
                out$fit$verdict,
                if (res$ok) sprintf("%.3f", out$fit$LL) else "-",
                res$n_ok, el))

    if (!out$fit$usable) {
      lc_announce_verdict(list(usable = FALSE, verdict = out$fit$verdict,
                               detail = lc_flag_detail(
                                 strsplit(out$fit$problems, ";")[[1]])),
                          sprintf("n=%d seed=%d K=%d", n, seed, K),
                          indent = "     ")
    }
  }
}

cat(sprintf("\nelapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

# --- Analysis ---------------------------------------------------------------
fits   <- read.csv(PATH_FITS)   %>% as_tibble()
params <- read.csv(PATH_PARAMS) %>% as_tibble()

# --- Usability, before anything else ---------------------------------------
# The convergence table below reports that essentially every cell converged,
# which is true and almost meaningless. This block runs first, and loudly,
# because "the optimiser returned a number" and "you may report this number"
# are different claims and the experiment exists to show where they part.
usability <- fits %>%
  filter(n != FULL_N) %>%
  group_by(n, K) %>%
  summarise(reps = n(), usable = mean(usable), .groups = "drop") %>%
  arrange(desc(n), K)

cat("\n=== USABLE FITS (share of replications) ====================\n")
print(as.data.frame(usability %>%
  mutate(usable = sprintf("%3.0f%%", 100 * usable)) %>%
  pivot_wider(names_from = K, values_from = usable, names_prefix = "K=")),
  row.names = FALSE)

worst <- usability %>% filter(usable == 0)
if (nrow(worst)) {
  cat(sprintf(paste0(
    "\n!! NOT ONE usable replication in %d of the %d cells:\n",
    "!! %s\n",
    "!! Every fit in those cells converged. None of them can be reported.\n"),
    nrow(worst), nrow(usability),
    paste(sprintf("n=%d K=%d", worst$n, worst$K), collapse = ", ")))
}

mode_counts <- fits %>%
  filter(n != FULL_N, !usable) %>%
  count(verdict, sort = TRUE)
if (nrow(mode_counts)) {
  cat("\nHow the unusable fits failed:\n")
  print(as.data.frame(mode_counts), row.names = FALSE)
}

benchmark <- params %>%
  filter(n == FULL_N) %>%
  select(K, parameter, truth = estimate)

# Convergence and stability by sample size.
convergence <- fits %>%
  filter(n != FULL_N) %>%
  group_by(n, K) %>%
  summarise(replications = n(),
            converged    = mean(ok),
            hessian_ok   = mean(hessian_ok),
            starts_ok    = mean(n_starts_ok),
            median_secs  = median(secs),
            .groups = "drop")

cat("\n=== Convergence by sample size =============================\n")
print(as.data.frame(convergence %>% mutate(across(where(is.numeric), ~ round(.x, 3)))),
      row.names = FALSE)

# Bias, RMSE and coverage against the full-sample benchmark.
accuracy <- params %>%
  filter(n != FULL_N) %>%
  inner_join(benchmark, by = c("K", "parameter")) %>%
  mutate(err = estimate - truth,
         covered = !is.na(se) & se > 0 &
           truth >= estimate - 1.96 * se & truth <= estimate + 1.96 * se) %>%
  group_by(n, K, parameter) %>%
  summarise(truth    = first(truth),
            mean_est = mean(estimate, na.rm = TRUE),
            bias     = mean(err, na.rm = TRUE),
            rmse     = sqrt(mean(err^2, na.rm = TRUE)),
            rel_rmse = sqrt(mean(err^2, na.rm = TRUE)) / abs(first(truth)),
            coverage = mean(covered, na.rm = TRUE),
            .groups  = "drop")

write_table(accuracy, "summary", prefix = "08")

cat("\n=== MNL accuracy by sample size ============================\n")
print(as.data.frame(accuracy %>% filter(K == 1) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 4)))),
      row.names = FALSE)

# Does BIC pick fewer classes when the sample is small?
selection <- fits %>%
  filter(n != FULL_N, ok) %>%
  group_by(n, seed) %>%
  filter(n_distinct(K) == length(K_SET)) %>%   # only fully estimated cells
  slice_min(BIC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  count(n, K_selected = K, name = "replications") %>%
  group_by(n) %>%
  mutate(share = replications / sum(replications)) %>%
  ungroup()

write_table(selection, "selection", prefix = "08")

cat("\n=== BIC-selected number of classes =========================\n")
print(as.data.frame(selection %>% mutate(share = round(share, 3))), row.names = FALSE)

# Value of travel time recovered per replication, MNL only.
vtt <- params %>%
  filter(K == 1, parameter %in% c("b_tt_1", "b_tc_1")) %>%
  select(n, seed, parameter, estimate) %>%
  pivot_wider(names_from = parameter, values_from = estimate) %>%
  mutate(vtt = 60 * b_tt_1 / b_tc_1)

vtt_benchmark <- vtt %>% filter(n == FULL_N) %>% pull(vtt)

cat(sprintf("\nfull-sample VTT benchmark: %.3f CHF/hour\n", vtt_benchmark))

# --- Figures ----------------------------------------------------------------
# The log axis itself is log_size_axis() in 00_setup.R; only the set of
# breaks differs between scripts, so this is a binding, not a second copy.
size_axis <- function(p) log_size_axis(p, c(SAMPLE_SIZES, FULL_N))

p_vtt <- size_axis(
  ggplot(vtt %>% filter(n != FULL_N), aes(n, vtt, group = n)) +
    geom_boxplot(outlier.size = 0.7) +
    geom_hline(yintercept = vtt_benchmark, linetype = "dashed", colour = "grey30")
) +
  labs(title = "Value of travel time by sample size (MNL)",
       subtitle = "Dashed line is the full-sample estimate",
       x = "respondents (log scale)", y = "CHF per hour")

p_rmse <- size_axis(
  ggplot(accuracy %>% filter(K == 1), aes(n, rel_rmse, colour = parameter)) +
    geom_line() + geom_point(size = 1.6)
) +
  labs(title = "Relative RMSE of MNL coefficients",
       subtitle = "RMSE as a share of the full-sample value",
       x = "respondents (log scale)", y = "RMSE / |benchmark|", colour = NULL)

p_conv <- size_axis(
  ggplot(convergence, aes(n, converged, colour = factor(K))) +
    geom_line() + geom_point(size = 1.6) + ylim(0, 1)
) +
  labs(title = "Convergence rate by sample size",
       x = "respondents (log scale)", y = "share converged", colour = "classes")

p_sel <- ggplot(selection, aes(factor(n), share, fill = factor(K_selected))) +
  geom_col() +
  labs(title = "Number of classes chosen by BIC",
       subtitle = "Small samples should favour fewer classes",
       x = "respondents", y = "share of replications", fill = "K selected")

fig <- (p_vtt | p_rmse) / (p_conv | p_sel)
ggsave(file.path(PATH_FIGURES, "08_small_sample.png"), fig,
       width = 14, height = 9, dpi = 150)

cat(sprintf("\nTables written to: %s\n", PATH_TABLES))
cat(sprintf("Figure written to: %s\n",
            file.path(PATH_FIGURES, "08_small_sample.png")))
