# lc_helpers.R -- shared latent class machinery
# Swiss route choice
#
# Apollo parses the source of apollo_lcPars and apollo_probabilities, so a
# K-class model cannot be assembled with get()/assign() at run time. These
# helpers write the class structure out as literal code for a given K and
# evaluate it into the global environment, which is where Apollo looks.
#
# Sourced by 05_lc_2class.R, 06_lc_multiclass.R and 08_small_sample_experiment.R.
# 06 used to carry its own copy of the four code-generation functions below;
# there is now one copy, here.
#
# Requires 00_setup.R (ROUTE_ATTRS, SEED).
# ---------------------------------------------------------------------------

if (!exists("ROUTE_ATTRS")) stop("source 00_setup.R before lc_helpers.R")

LC_ATTRS <- ROUTE_ATTRS

# Free parameters: four coefficients per class, plus K-1 allocation constants
# (delta_1 is normalised to zero). LC1 has no allocation model at all.
lc_n_par <- function(K) 4L * K + max(K - 1L, 0L)

lc_par_names <- function(K) {
  betas <- as.vector(t(outer(paste0("b_", LC_ATTRS), seq_len(K), paste, sep = "_")))
  if (K == 1L) return(betas)
  c(betas, paste0("delta_", seq_len(K)))
}

# --- Silencing -------------------------------------------------------------
# One definition, in 00_setup.R, aliased here so latent-class code reads
# consistently. 03 and 04 previously carried their own identical copies.
lc_quietly <- quietly

# --- Code generation: constant-only allocation -----------------------------
lc_make_lcPars <- function(K) {
  stopifnot(K >= 2L)

  lists <- vapply(LC_ATTRS, function(a) {
    sprintf('  lcpars[["b_%s"]] <- list(%s)', a,
            paste0("b_", a, "_", seq_len(K), collapse = ", "))
  }, character(1))

  alloc   <- sprintf('  V[["class_%d"]] <- delta_%d', seq_len(K), seq_len(K))
  classes <- paste0("class_", seq_len(K), " = ", seq_len(K), collapse = ", ")

  code <- c(
    "apollo_lcPars <- function(apollo_beta, apollo_inputs) {",
    "  lcpars <- list()",
    lists,
    "  V <- list()",
    alloc,
    sprintf("  classAlloc_settings <- list(classes = c(%s), utilities = V)", classes),
    '  lcpars[["pi_values"]] <- apollo_classAlloc(classAlloc_settings)',
    "  return(lcpars)",
    "}"
  )
  eval(parse(text = paste(code, collapse = "\n")), envir = globalenv())
}

lc_make_probabilities <- function(K) {
  if (K == 1L) {
    # A one-class latent class model is just the MNL; estimating it as such
    # keeps LC1 comparable without asking apollo_lc to handle a degenerate case.
    code <- c(
      'apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {',
      "  apollo_attach(apollo_beta, apollo_inputs)",
      "  on.exit(apollo_detach(apollo_beta, apollo_inputs))",
      "  P <- list()",
      "  V <- list()",
      '  V[["alt1"]] <- b_tt_1 * tt1 + b_tc_1 * tc1 + b_hw_1 * hw1 + b_ch_1 * ch1',
      '  V[["alt2"]] <- b_tt_1 * tt2 + b_tc_1 * tc2 + b_hw_1 * hw2 + b_ch_1 * ch2',
      "  mnl_settings <- list(alternatives = c(alt1 = 1, alt2 = 2), avail = 1,",
      "                       choiceVar = choice, utilities = V)",
      '  P[["model"]] <- apollo_mnl(mnl_settings, functionality)',
      "  P <- apollo_panelProd(P, apollo_inputs, functionality)",
      "  P <- apollo_prepareProb(P, apollo_inputs, functionality)",
      "  return(P)",
      "}"
    )
  } else {
    code <- c(
      'apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {',
      "  apollo_attach(apollo_beta, apollo_inputs)",
      "  on.exit(apollo_detach(apollo_beta, apollo_inputs))",
      "  P <- list()",
      "  mnl_settings <- list(alternatives = c(alt1 = 1, alt2 = 2), avail = 1,",
      "                       choiceVar = choice)",
      sprintf("  for (s in 1:%d) {", K),
      "    V <- list()",
      '    V[["alt1"]] <- b_tt[[s]] * tt1 + b_tc[[s]] * tc1 + b_hw[[s]] * hw1 + b_ch[[s]] * ch1',
      '    V[["alt2"]] <- b_tt[[s]] * tt2 + b_tc[[s]] * tc2 + b_hw[[s]] * hw2 + b_ch[[s]] * ch2',
      "    mnl_settings$utilities     <- V",
      '    mnl_settings$componentName <- paste0("Class_", s)',
      '    P[[paste0("Class_", s)]] <- apollo_mnl(mnl_settings, functionality)',
      '    P[[paste0("Class_", s)]] <- apollo_panelProd(P[[paste0("Class_", s)]],',
      "                                                apollo_inputs, functionality)",
      "  }",
      "  lc_settings <- list(inClassProb = P, classProb = pi_values)",
      '  P[["model"]] <- apollo_lc(lc_settings, apollo_inputs, functionality)',
      "  P <- apollo_prepareProb(P, apollo_inputs, functionality)",
      "  return(P)",
      "}"
    )
  }
  eval(parse(text = paste(code, collapse = "\n")), envir = globalenv())
}

# Installs the pair of user functions for a given K, clearing apollo_lcPars for
# the one-class case so a stale definition cannot leak between models.
lc_install <- function(K) {
  lc_make_probabilities(K)
  if (K > 1L) {
    lc_make_lcPars(K)
  } else if (exists("apollo_lcPars", envir = globalenv())) {
    rm("apollo_lcPars", envir = globalenv())
  }
}

# --- Starting values -------------------------------------------------------
# Class k starts at the reference (full-sample MNL) vector scaled by a
# class-level factor and jittered per coefficient. A positive scale preserves
# the sign of every coefficient, so no start begins somewhere absurd.
lc_start_values <- function(K, draw, b_ref, seed) {
  set.seed(seed)

  beta <- numeric(0)
  # The first draw is deterministic and spreads classes over a fixed range, so
  # at least one start is reproducible independent of the RNG stream.
  scales <- if (draw == 1L) seq(0.5, 1.8, length.out = K) else runif(K, 0.4, 2.0)

  for (a in LC_ATTRS) {
    jitter <- if (draw == 1L) rep(1, K) else exp(rnorm(K, 0, 0.25))
    vals <- unname(b_ref[[paste0("b_", a)]]) * scales * jitter
    names(vals) <- paste0("b_", a, "_", seq_len(K))
    beta <- c(beta, vals)
  }

  # LC1 has no allocation model: a lone delta_1 would be a parameter that does
  # not enter the likelihood, which Apollo rejects outright.
  if (K == 1L) return(beta)

  deltas <- if (draw == 1L) rep(0, K) else c(0, rnorm(K - 1L, 0, 0.5))
  names(deltas) <- paste0("delta_", seq_len(K))
  deltas["delta_1"] <- 0

  c(beta, deltas)
}

# Naive dispersed starts. Tested against perturbing the MNL estimates on LC2:
# uniform starts reached the best optimum in 8 of 20 tries versus 1 of 20 for
# the perturbation scheme. Anchoring on the MNL biases the search toward the
# MNL-like local optimum, which is exactly the one to avoid.
lc_start_uniform <- function(par_names, seed, lo = -0.1, hi = 0.1) {
  set.seed(seed)
  starts <- runif(length(par_names), lo, hi)
  names(starts) <- par_names
  if ("delta_1" %in% par_names) starts["delta_1"] <- 0
  starts
}

# --- Class shares ----------------------------------------------------------
# Class shares under a constant-only allocation: a softmax of the deltas,
# identical for every respondent.
lc_class_shares <- function(est, K) {
  if (K == 1L) return(1)
  d <- vapply(seq_len(K), function(k) unname(est[[paste0("delta_", k)]]), numeric(1))
  exp(d) / sum(exp(d))
}

# --- Canonical class order -------------------------------------------------
# Latent class labels are identified only up to permutation: two runs can find
# the same optimum with classes 1 and 2 swapped. Averaging parameters across
# replications -- or comparing "class 1" between the constant-only model in 06
# and the covariate model in 05 -- without fixing an order would mix distinct
# classes together. Classes are always sorted by their travel time
# coefficient, most negative first.
#
# lc_class_order returns the permutation itself so that anything derived from
# the raw model (standard errors, delta-method valuations) can be relabelled
# the same way the estimates were.
lc_class_order <- function(est, K) {
  if (K == 1L) return(1L)
  b_tt <- vapply(seq_len(K), function(k) unname(est[[paste0("b_tt_", k)]]), numeric(1))
  order(b_tt)
}

# Reorders estimates into canonical order and re-normalises the deltas so the
# new first class is again the reference.
lc_canonical <- function(est, K) {
  if (K == 1L) return(est)

  ord <- lc_class_order(est, K)
  out <- est

  for (a in LC_ATTRS) {
    for (k in seq_len(K)) {
      out[[paste0("b_", a, "_", k)]] <- unname(est[[paste0("b_", a, "_", ord[k])]])
    }
  }

  d  <- vapply(seq_len(K), function(k) unname(est[[paste0("delta_", k)]]), numeric(1))
  d2 <- d[ord] - d[ord][1]
  for (k in seq_len(K)) out[[paste0("delta_", k)]] <- d2[k]

  out
}

# ---------------------------------------------------------------------------
# Class allocation covariates
#
# The constant-only allocation above says how big each class is. These make
# membership depend on the respondent: class k gets a linear index in the
# socio-demographics instead of a bare constant.
#
# Class 1 is the reference and its whole allocation utility is normalised to
# zero, not just its constant -- only differences between classes are
# identified. leisure is the reference trip purpose because the four purpose
# dummies sum to one for every respondent (verified in 01_data_audit.R), so
# including all four alongside a constant would be exact collinearity.
# ---------------------------------------------------------------------------

# Defined in 00_setup.R so 03 (attribute interactions) and 05 (class
# allocation) cannot disagree about which covariates the project uses.
LC_COVARS <- MODEL_COVARS

# Free parameters: 4 per class, plus (K-1) x (1 constant + one per covariate).
lc_n_par_cov <- function(K, covars = LC_COVARS) {
  4L * K + max(K - 1L, 0L) * (1L + length(covars))
}

lc_par_names_cov <- function(K, covars = LC_COVARS) {
  betas <- as.vector(t(outer(paste0("b_", LC_ATTRS), seq_len(K), paste, sep = "_")))
  if (K == 1L) return(betas)
  alloc <- unlist(lapply(2:K, function(k) {
    c(paste0("delta_", k), paste0("g_", covars, "_", k))
  }))
  c(betas, alloc)
}

lc_make_lcPars_cov <- function(K, covars = LC_COVARS) {
  # With one class there is nothing to allocate and 2:K would run backwards,
  # generating code for a class_2 that does not exist. Callers get the
  # constant-only (no-op) path instead.
  stopifnot(K >= 2L)

  lists <- vapply(LC_ATTRS, function(a) {
    sprintf('  lcpars[["b_%s"]] <- list(%s)', a,
            paste0("b_", a, "_", seq_len(K), collapse = ", "))
  }, character(1))

  # Reference class: allocation utility fixed at zero.
  alloc <- '  V[["class_1"]] <- 0'
  for (k in 2:K) {
    terms <- c(sprintf("delta_%d", k),
               sprintf("g_%s_%d * %s", covars, k, names(covars)))
    alloc <- c(alloc, sprintf('  V[["class_%d"]] <- %s', k,
                              paste(terms, collapse = " + ")))
  }

  classes <- paste0("class_", seq_len(K), " = ", seq_len(K), collapse = ", ")

  code <- c(
    "apollo_lcPars <- function(apollo_beta, apollo_inputs) {",
    "  lcpars <- list()",
    lists,
    "  V <- list()",
    alloc,
    sprintf("  classAlloc_settings <- list(classes = c(%s), utilities = V)", classes),
    '  lcpars[["pi_values"]] <- apollo_classAlloc(classAlloc_settings)',
    "  return(lcpars)",
    "}"
  )
  eval(parse(text = paste(code, collapse = "\n")), envir = globalenv())
}

lc_install_cov <- function(K, covars = LC_COVARS) {
  lc_make_probabilities(K)          # the within-class MNL is unchanged
  if (K > 1L) {
    lc_make_lcPars_cov(K, covars)
  } else if (exists("apollo_lcPars", envir = globalenv())) {
    rm("apollo_lcPars", envir = globalenv())
  }
}

# Class shares for one respondent profile x under the covariate allocation.
# x is a named list or vector giving a value for every covariate.
lc_class_shares_cov <- function(est, K, x, covars = LC_COVARS) {
  if (K == 1L) return(1)
  v <- c(0, vapply(2:K, function(k) {
    unname(est[[paste0("delta_", k)]]) +
      sum(vapply(names(covars), function(cv) {
        unname(est[[paste0("g_", covars[[cv]], "_", k)]]) * x[[cv]]
      }, numeric(1)))
  }, numeric(1)))
  exp(v) / sum(exp(v))
}

# Canonical order for a covariate allocation model.
#
# Class 1's allocation utility is identically zero here -- there is no delta_1
# and no g_*_1 -- so re-normalising after a permutation means subtracting the
# new reference class's WHOLE allocation index, coefficient by coefficient:
#
#   V_k' = (delta_k - delta_j) + sum_c (g_ck - g_cj) x_c        j = ord[1]
#
# which leaves every class-share difference, and therefore the likelihood,
# exactly unchanged.
lc_canonical_cov <- function(est, K, covars = LC_COVARS) {
  if (K == 1L) return(est)

  ord <- lc_class_order(est, K)
  out <- est

  for (a in LC_ATTRS) {
    for (k in seq_len(K)) {
      out[[paste0("b_", a, "_", k)]] <- unname(est[[paste0("b_", a, "_", ord[k])]])
    }
  }

  # Pull the allocation coefficients into full K-length vectors with the
  # structural zeros for class 1 made explicit, permute, then re-reference.
  pull <- function(prefix) {
    v <- vapply(seq_len(K), function(k) {
      if (k == 1L) 0 else lc_par_value(est, paste0(prefix, k))
    }, numeric(1))
    v[ord] - v[ord][1]
  }

  d <- pull("delta_")
  g <- lapply(covars, function(cv) pull(paste0("g_", cv, "_")))
  names(g) <- covars

  for (k in 2:K) {
    out[[paste0("delta_", k)]] <- d[k]
    for (cv in covars) out[[paste0("g_", cv, "_", k)]] <- g[[cv]][k]
  }

  out
}

# --- Valuations ------------------------------------------------------------
# Willingness-to-pay per class, with delta-method standard errors on the
# robust covariance matrix. 06 previously reported bare coefficient ratios
# with no uncertainty at all, which is indefensible for a quantity that is a
# ratio of two estimates: its sampling distribution is not normal and its
# spread is far wider than either coefficient suggests on its own.
#
# A class whose cost coefficient sits near zero produces an explosive ratio.
# Those are flagged (ratio_reliable = FALSE) using the t-ratio on b_tc rather
# than silently reported as if they meant something.
LC_VALUATIONS <- c(
  vtt_chf_per_hour     = "b_tt_%1$d/b_tc_%1$d*60",
  headway_chf_per_hour = "b_hw_%1$d/b_tc_%1$d*60",
  interchange_chf      = "b_ch_%1$d/b_tc_%1$d",
  interchange_minutes  = "b_ch_%1$d/b_tt_%1$d"
)

LC_VALUATION_UNITS <- c(
  vtt_chf_per_hour     = "CHF per hour",
  headway_chf_per_hour = "CHF per hour",
  interchange_chf      = "CHF per interchange",
  interchange_minutes  = "minutes per interchange"
)

# model  : a fitted apollo model with a covariance matrix
# K      : number of classes
# ord    : permutation from lc_class_order() applied to the model's estimates,
#          so raw class ord[k] is reported as canonical class k
lc_valuations <- function(model, K, ord = seq_len(K), t_min = 1.96) {
  est <- model$estimate

  expr <- unlist(lapply(seq_len(K), function(k) {
    setNames(sprintf(LC_VALUATIONS, k),
             paste0(names(LC_VALUATIONS), "__class", k))
  }))

  dm <- tryCatch(
    lc_quietly(apollo_deltaMethod(model, list(expression = expr))),
    error = function(e) e
  )
  if (inherits(dm, "error")) {
    warning("delta method failed: ", conditionMessage(dm))
    return(tibble())
  }

  dm <- as_tibble(dm)
  names(dm)[names(dm) == "Expression"]  <- "measure"
  names(dm)[names(dm) == "Value"]       <- "value"
  names(dm)[names(dm) == "s.e."]        <- "se"
  names(dm)[names(dm) == "t-ratio (0)"] <- "t_ratio"

  # Reliability is a property of the denominator, so it is computed once per
  # raw class from the cost coefficient rather than per valuation.
  se_all <- sqrt(diag(model$varcov))
  tc_t <- vapply(seq_len(K), function(k) {
    nm <- paste0("b_tc_", k)
    unname(est[[nm]]) / unname(se_all[[nm]])
  }, numeric(1))

  # Map raw class index -> canonical class index.
  canon_of <- integer(K)
  canon_of[ord] <- seq_len(K)

  dm %>%
    separate(measure, into = c("measure", "raw_class"), sep = "__class",
             convert = TRUE) %>%
    mutate(
      class           = canon_of[raw_class],
      unit            = unname(LC_VALUATION_UNITS[measure]),
      ci_low          = value - 1.96 * se,
      ci_high         = value + 1.96 * se,
      ratio_reliable  = abs(tc_t[raw_class]) >= t_min,
      .after = measure
    ) %>%
    arrange(class, match(measure, names(LC_VALUATIONS)))
}

# --- Multi-start estimation ------------------------------------------------
# 05, 06 and 08 each ran their own copy of this loop, and they had drifted:
# 05 re-estimated the winning start against a STALE apollo_inputs left over
# from whichever draw last validated, while 06 re-validated first. There is
# now one driver and it always re-validates.
#
# starts : list of named numeric vectors, one per draw
# returns: model, run_log (one row per draw), n_at_best, hessian_ok
lc_estimate_best <- function(starts, apollo_fixed, control,
                             label = "", verbose = TRUE, final_hessian = TRUE,
                             tol_at_best = 1e-3) {

  assign("apollo_control", control, envir = globalenv())
  assign("apollo_fixed",   apollo_fixed, envir = globalenv())

  best <- NULL
  run_log <- tibble()

  for (draw in seq_along(starts)) {
    assign("apollo_beta", starts[[draw]], envir = globalenv())

    inputs <- tryCatch(lc_quietly(apollo_validateInputs(silent = TRUE)),
                       error = function(e) e)
    if (inherits(inputs, "error")) {
      run_log <- bind_rows(run_log, tibble(
        draw = draw, converged = FALSE, LL = NA_real_,
        error = paste("input:", conditionMessage(inputs))))
      if (verbose) cat(sprintf("  start %2d/%2d : INPUT ERROR - %s\n",
                               draw, length(starts), conditionMessage(inputs)))
      next
    }
    assign("apollo_inputs", inputs, envir = globalenv())

    # The search pass skips the covariance matrix; only the winning start pays
    # for it. That makes each draw roughly twice as fast. The error message is
    # kept rather than discarded: a start that fails for a specification
    # reason must not look like an ordinary non-convergence.
    m <- tryCatch(
      lc_quietly(apollo_estimate(
        apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
        estimate_settings = list(silent = TRUE, writeIter = FALSE,
                                 hessianRoutine = "none"))),
      error = function(e) e
    )

    failed <- inherits(m, "error")
    ok     <- !failed && is.finite(m$maximum)
    ll     <- if (ok) m$maximum else NA_real_

    run_log <- bind_rows(run_log, tibble(
      draw = draw, converged = ok, LL = ll,
      error = if (failed) conditionMessage(m) else NA_character_))

    if (verbose) {
      cat(sprintf("  start %2d/%2d : %s\n", draw, length(starts),
                  if (ok) sprintf("LL = %12.4f", ll)
                  else if (failed) paste("ERROR -", conditionMessage(m))
                  else "did not converge"))
    }

    if (ok && (is.null(best) || ll > best$maximum)) best <- m
  }

  if (is.null(best)) {
    return(list(model = NULL, run_log = run_log, n_at_best = 0L,
                n_ok = 0L, hessian_ok = FALSE, ok = FALSE))
  }

  n_ok      <- sum(run_log$converged)
  n_at_best <- sum(run_log$converged &
                     abs(run_log$LL - max(run_log$LL, na.rm = TRUE)) < tol_at_best)

  if (!final_hessian) {
    return(list(model = best, run_log = run_log, n_at_best = n_at_best,
                n_ok = n_ok, hessian_ok = FALSE, ok = TRUE))
  }

  # Re-estimate from the winning point with the full covariance matrix, always
  # against freshly validated inputs.
  assign("apollo_beta", best$estimate, envir = globalenv())
  inputs <- tryCatch(lc_quietly(apollo_validateInputs(silent = TRUE)),
                     error = function(e) e)
  if (!inherits(inputs, "error")) assign("apollo_inputs", inputs, envir = globalenv())

  final <- tryCatch(
    lc_quietly(apollo_estimate(apollo_beta, apollo_fixed,
                               apollo_probabilities, apollo_inputs,
                               estimate_settings = list(silent = TRUE,
                                                        writeIter = FALSE))),
    error = function(e) e
  )

  # In small samples the covariance matrix can fail even when the point
  # estimates are fine, so a failure is recorded rather than allowed to lose
  # the fit.
  hessian_ok <- !inherits(final, "error") &&
    !is.null(final$varcov) && all(is.finite(diag(final$varcov)))

  list(model = if (hessian_ok) final else best, run_log = run_log,
       n_at_best = n_at_best, n_ok = n_ok, hessian_ok = hessian_ok, ok = TRUE)
}

# --- Canonical-order parameter tables --------------------------------------
# Reordering classes is only a relabelling for the within-class coefficients,
# so their standard errors are simply permuted alongside them. The allocation
# parameters are different: the reference class changes, so every one of them
# becomes a CONTRAST against the new reference,
#
#   delta_k' = delta_k - delta_j        g_ck' = g_ck - g_cj      j = ord[1]
#
# and its variance is V_kk + V_jj - 2 V_kj, not V_kk. Reporting the permuted
# standard error here would understate or overstate the uncertainty by however
# much the two allocation parameters covary.
#
# A parameter that is structurally zero -- delta_1 in the constant-only model
# (fixed by apollo_fixed) and the whole class-1 index in the covariate model
# (absent from the parameter vector) -- contributes zero variance and zero
# covariance, which is what the is.na branches below encode.
lc_contrast_se <- function(V, a, b) {
  va  <- if (is.na(a) || !a %in% rownames(V)) 0 else V[a, a]
  vb  <- if (is.na(b) || !b %in% rownames(V)) 0 else V[b, b]
  cab <- if (is.na(a) || is.na(b) ||
             !a %in% rownames(V) || !b %in% rownames(V)) 0 else V[a, b]
  sqrt(max(va + vb - 2 * cab, 0))
}

# Apollo returns model$estimate as a NAMED NUMERIC VECTOR, and x[["absent"]]
# on an atomic vector is an error, not NULL -- unlike a list, which is what
# the synthetic estimates in the tests are. Both have to work here: the
# covariate allocation model has no delta_1 and no g_*_1 at all, because
# class 1's entire allocation index is normalised to zero, so looking one up
# is the normal case rather than a mistake.
lc_par_value <- function(est, name) {
  if (!name %in% names(est)) return(0)
  unname(est[[name]])
}

# Returns a tidy tibble of estimates and standard errors with classes in
# canonical order (most negative b_tt first). Used by both 05 and 06 so the
# two scripts cannot disagree about which class is "class 1".
lc_canonical_params <- function(model, K, covariate_alloc = FALSE,
                                covars = LC_COVARS) {
  est <- model$estimate
  V   <- model$varcov
  ord <- lc_class_order(est, K)
  ref <- ord[1]

  se_of <- function(name) {
    if (!name %in% rownames(V)) return(NA_real_)
    sqrt(V[name, name])
  }

  out <- tibble()

  for (a in LC_ATTRS) {
    for (k in seq_len(K)) {
      raw <- paste0("b_", a, "_", ord[k])
      out <- bind_rows(out, tibble(
        parameter = paste0("b_", a, "_", k),
        estimate  = lc_par_value(est, raw),
        se        = se_of(raw),
        raw_parameter = raw
      ))
    }
  }

  if (K > 1L) {
    alloc_prefixes <- if (covariate_alloc) {
      c("delta_", paste0("g_", unname(covars), "_"))
    } else {
      "delta_"
    }
    for (pre in alloc_prefixes) {
      for (k in 2:K) {
        raw_k <- paste0(pre, ord[k])
        raw_j <- paste0(pre, ref)
        out <- bind_rows(out, tibble(
          parameter = paste0(pre, k),
          estimate  = lc_par_value(est, raw_k) - lc_par_value(est, raw_j),
          se        = lc_contrast_se(V, raw_k, raw_j),
          raw_parameter = if (ref == ord[k]) raw_k else paste0(raw_k, " - ", raw_j)
        ))
      }
    }
  }

  out %>% mutate(t_ratio = estimate / se)
}

# --- Resumable fits --------------------------------------------------------
# 06 and 05 each spend tens of minutes on a multi-start sweep, and an
# interruption part-way through -- a killed process, a full disk, a laptop
# closing -- threw all of it away. 08 has always cached its results and
# resumed; these two now do the same, one model at a time.
#
# The cached model and its per-start log are written together, so a partly
# written cache cannot be mistaken for a complete one, and everything the
# comparison tables report (n_at_best, hessian_ok) is recomputed from them
# rather than stored separately and allowed to drift.
lc_fit_result <- function(model, run_log, tol = 1e-3) {
  list(
    model      = model,
    run_log    = run_log,
    ok         = TRUE,
    n_ok       = sum(run_log$converged),
    n_at_best  = sum(run_log$converged &
                       abs(run_log$LL - max(run_log$LL, na.rm = TRUE)) < tol),
    hessian_ok = !is.null(model$varcov) && all(is.finite(diag(model$varcov)))
  )
}

# key    : basename for the cache, e.g. "Swiss_LC3"
# resume : FALSE forces a refit even when a cache exists
# fit_fn : zero-argument function returning an lc_estimate_best() result
lc_fit_cached <- function(key, fit_fn, resume = TRUE) {
  model_path <- file.path(PATH_MODELS, paste0(key, "_model.rds"))
  runs_path  <- file.path(PATH_MODELS, paste0(key, "_runs.rds"))

  if (resume && file.exists(model_path) && file.exists(runs_path)) {
    model   <- readRDS(model_path)
    run_log <- readRDS(runs_path)
    cat(sprintf("  cached fit reused: %s (LL = %.4f, %d starts on record)\n",
                key, model$maximum, nrow(run_log)))
    cat("  delete the _model.rds and _runs.rds to force a refit\n")
    return(lc_fit_result(model, run_log))
  }

  fit <- fit_fn()
  if (isTRUE(fit$ok)) {
    saveRDS(fit$model,   model_path)
    saveRDS(fit$run_log, runs_path)
  }
  fit
}

# The switch itself lives in 00_setup.R (resume_enabled) so 04, which caches a
# single model rather than a multi-start sweep, honours the same REFIT=1.
lc_resume_enabled <- resume_enabled

# --- Is this fit usable? ---------------------------------------------------
# A latent class model can converge, report a log-likelihood, hand back class
# shares, and still be worthless. The optimiser has no opinion about whether
# two classes collapsed onto each other, whether a class holds one respondent,
# whether a coefficient came back with a sign that is economically impossible,
# or whether the covariance matrix it just failed to invert means every
# standard error in the output is undefined.
#
# Left to itself the pipeline printed "converged" for all of these. At n = 20
# that word appeared over fits whose Hessian was singular in EVERY replication.
# This is the single definition of "usable", and everything -- 05, 06, 08 and
# the tables they write -- routes through it so no script can quietly disagree.

LC_MIN_CLASS_SHARE <- 0.05   # thinner than this and a class cannot be read
LC_COLLAPSE_REL    <- 0.10   # closer than this and two classes are one class

# The K x 4 matrix of class coefficient vectors.
lc_class_coefs <- function(est, K) {
  m <- vapply(LC_ATTRS, function(a) {
    vapply(seq_len(K), function(k) lc_par_value(est, paste0("b_", a, "_", k)),
           numeric(1))
  }, numeric(K))
  matrix(m, nrow = K, dimnames = list(NULL, LC_ATTRS))
}

# Smallest relative gap between any pair of classes, scale-free: for each pair
# the largest per-attribute difference divided by the larger magnitude. Needs
# no benchmark, so it can be evaluated at fit time.
lc_min_class_gap <- function(est, K) {
  if (K < 2L) return(NA_real_)
  m <- lc_class_coefs(est, K)
  if (anyNA(m)) return(NA_real_)
  gaps <- utils::combn(K, 2, function(ij) {
    a <- m[ij[1], ]; b <- m[ij[2], ]
    max(abs(a - b) / pmax(abs(a), abs(b), 1e-12))
  })
  min(gaps)
}

# shares: class shares if the caller can compute them (constant-only
# allocation can; covariate allocation needs the data, so 05 passes the
# respondent-averaged shares). NULL skips that one check and says so.
lc_fit_verdict <- function(fit, K, shares = NULL) {
  flags <- character(0)

  if (!isTRUE(fit$ok) || is.null(fit$model)) {
    return(list(usable = FALSE, verdict = "NO FIT", flags = "no_fit",
                detail = "no starting value produced a finite log-likelihood"))
  }

  est <- fit$model$estimate

  if (!isTRUE(fit$hessian_ok)) {
    flags <- c(flags, "no_covariance")
  }
  coefs <- lc_class_coefs(est, K)
  if (any(coefs > 0, na.rm = TRUE)) flags <- c(flags, "wrong_sign")
  if (!is.null(shares) && min(shares, na.rm = TRUE) < LC_MIN_CLASS_SHARE) {
    flags <- c(flags, "degenerate_class")
  }
  gap <- lc_min_class_gap(est, K)
  if (!is.na(gap) && gap < LC_COLLAPSE_REL) flags <- c(flags, "classes_collapsed")

  if (!length(flags)) {
    return(list(usable = TRUE, verdict = "ok", flags = character(0),
                detail = "converged, invertible Hessian, classes distinct and correctly signed"))
  }

  list(usable = FALSE,
       verdict = lc_primary_verdict(flags),
       flags = flags,
       detail = lc_flag_detail(flags))
}

# Ordered by severity: what makes a fit least usable comes first.
LC_FLAG_ORDER <- c("no_fit", "no_covariance", "wrong_sign",
                   "degenerate_class", "classes_collapsed")

lc_flag_detail <- function(flags) {
  detail <- c(
    no_fit            = "no starting value produced a finite log-likelihood",
    no_covariance     = "the Hessian is singular -- every standard error, t-ratio and confidence interval in this fit is UNDEFINED",
    wrong_sign        = "a route attribute has a positive coefficient, which is economically impossible; the class is fitting noise",
    degenerate_class  = sprintf("a class holds less than %.0f%% of respondents and cannot be interpreted", 100 * LC_MIN_CLASS_SHARE),
    classes_collapsed = sprintf("two classes agree to within %.0f%% on every attribute -- the model fitted one class twice", 100 * LC_COLLAPSE_REL)
  )
  flags <- LC_FLAG_ORDER[LC_FLAG_ORDER %in% flags]
  if (!length(flags)) return("")
  paste(unname(detail[flags]), collapse = "; ")
}

lc_primary_verdict <- function(flags) {
  primary <- LC_FLAG_ORDER[LC_FLAG_ORDER %in% flags][1]
  if (is.na(primary)) return("ok")
  toupper(gsub("_", " ", primary))
}

# Prints the verdict so an unusable fit cannot be skimmed past. The marker is
# greppable on purpose: `grep UNUSABLE` over a run log finds every one.
lc_announce_verdict <- function(v, label, indent = "  ") {
  if (isTRUE(v$usable)) {
    cat(sprintf("%s%s: ok\n", indent, label))
    return(invisible(v))
  }
  bar <- strrep("!", 68)
  cat(sprintf("\n%s%s\n", indent, bar))
  cat(sprintf("%s!! UNUSABLE -- %s -- %s\n", indent, label, v$verdict))
  for (line in strwrap(v$detail, width = 62)) {
    cat(sprintf("%s!!   %s\n", indent, line))
  }
  cat(sprintf("%s!! Any number derived from this fit is not a result.\n", indent))
  cat(sprintf("%s%s\n\n", indent, bar))
  invisible(v)
}

# One-line summary for a table column.
lc_verdict_label <- function(v) if (isTRUE(v$usable)) "ok" else v$verdict
