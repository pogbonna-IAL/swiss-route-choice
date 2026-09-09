# model_helpers.R -- shared scoring and hold-out machinery
# Swiss route choice
#
# In-sample fit statistics cannot distinguish a model that generalises from
# one that has fitted noise, and every model in this project is compared to
# every other on likelihood. These helpers give 02 (MNL), 03 (MNL with
# covariates) and 04 (mixed logit) one hold-out split and one scoring
# function, so the numbers in their validation tables are comparable by
# construction rather than by coincidence.
#
# Sourced by 02_mnl.R, 03_mnl_covariates.R and 04_mixed_logit.R.
# Requires 00_setup.R (SEED, ROUTE_ATTRS).
# ---------------------------------------------------------------------------

if (!exists("ROUTE_ATTRS")) stop("source 00_setup.R before model_helpers.R")

# --- Hold-out split --------------------------------------------------------
# Splitting on RESPONDENT, not on task, is what makes this a real test.
# Putting some of a person's nine tasks in training and the rest in the
# hold-out would let the model exploit that person's own revealed taste, and
# the hold-out fit would flatter every model with taste heterogeneity in it --
# which is exactly the class of model 03 and 04 are trying to justify.
holdout_split <- function(db, frac = 0.5, seed = SEED) {
  set.seed(seed)
  ids       <- unique(db$ID)
  train_ids <- sample(ids, floor(frac * length(ids)))
  list(
    train     = db %>% filter(ID %in% train_ids),
    test      = db %>% filter(!ID %in% train_ids),
    train_ids = train_ids
  )
}

# --- Attribute differences -------------------------------------------------
# Both alternatives are always available and the coefficients are generic, so
# only the attribute DIFFERENCE enters the binary logit. Returning it as a
# matrix lets a caller with row-varying coefficients (03) multiply
# element-wise instead of writing the utility out four times.
attr_diff_matrix <- function(db, attrs = ROUTE_ATTRS) {
  m <- vapply(attrs, function(a) db[[paste0(a, 1)]] - db[[paste0(a, 2)]],
              numeric(nrow(db)))
  colnames(m) <- attrs
  m
}

# --- Utility difference for the fixed-coefficient MNL -----------------------
# Defined here rather than three times over: 02, 03 and 04 each carried an
# identical copy, and it belongs beside attr_diff_matrix, which it calls.
mnl_dv <- function(b, db) {
  as.vector(attr_diff_matrix(db) %*% b[paste0("b_", ROUTE_ATTRS)])
}

# --- Scoring ---------------------------------------------------------------
# +1 for a respondent who chose alternative 1, -1 for alternative 2. With two
# alternatives the chosen option's utility advantage is just the signed
# difference, so the chosen probability is plogis(sign * dv) whether dv is a
# vector (02, 03) or an observations-by-draws matrix (04).
chosen_sign <- function(db) ifelse(db$choice == 1L, 1, -1)

# dv is the utility difference V(alt 1) - V(alt 2), one value per row of db.
# For a binary logit that is all the information there is, so the predicted
# probability has a closed form and there is no need to round-trip a second
# database through Apollo's global state to score a hold-out sample.
binary_logit_score <- function(dv, db, label, model = NA_character_) {
  stopifnot(length(dv) == nrow(db))

  # plogis, not 1/(1+exp(-dv)). The naive form overflows: exp(745) is Inf in
  # double precision, so any dv below about -745 gives a probability of
  # EXACTLY zero and a log-likelihood of -Inf. That is not a hypothetical --
  # 04 draws lognormal coefficients whose tail produces utility differences
  # well past that, and the hold-out score came back NA because of it.
  # The chosen alternative's log-probability, obtained by flipping the sign of
  # the utility difference for respondents who picked alternative 2. Written
  # as a multiplication rather than ifelse() on purpose -- see the note in
  # panel_logit_score, where the same expression on a matrix silently drops
  # every column but the first.
  p1      <- plogis(dv)
  log_pch <- plogis(chosen_sign(db) * dv, log.p = TRUE)
  ll <- sum(log_pch)

  tibble(
    model         = model,
    sample        = label,
    respondents   = n_distinct(db$ID),
    observations  = nrow(db),
    LL            = ll,
    LL_per_obs    = ll / nrow(db),
    # Against a coin flip, which is the right null for two unlabelled
    # alternatives with a near 50/50 split.
    rho2_vs_coin  = 1 - ll / (-nrow(db) * log(2)),
    hit_rate      = mean((p1 > 0.5) == (db$choice == 1L)),
    # Recovered market share. A model can hit its overall share while getting
    # every individual task wrong, so this is reported alongside the hit rate,
    # never instead of it.
    share_alt1_observed  = mean(db$choice == 1L),
    share_alt1_predicted = mean(p1)
  )
}

# Panel simulated score for a model with random coefficients (04).
#
# The unconditional probability of a respondent's whole sequence is
#
#   P_n = (1/R) sum_r prod_t P_nt(beta_r)
#
# and the likelihood is the sum of its log over respondents. It does NOT
# factorise over tasks, so unlike the closed-form case above the per-task
# probability is not defined; hit rate and share are therefore computed from
# the draw-averaged per-task probability, which is the quantity a forecaster
# would actually use.
#
# dv_draws: matrix, nrow(db) x R, of utility differences under each draw
panel_logit_score <- function(dv_draws, db, label, model = NA_character_) {
  stopifnot(nrow(dv_draws) == nrow(db))

  # NEVER use ifelse() here. It returns a result shaped like its TEST
  # argument, and db$choice == 1L is a plain vector, so
  # ifelse(db$choice == 1L, p, 1 - p) on an observations-by-draws matrix
  # silently returns a VECTOR holding only the first draw -- no error, no
  # warning, and a log-likelihood that does not change when you add draws.
  # Multiplying by the sign keeps the matrix shape.
  #
  # The log-probability is also computed directly rather than as log(p),
  # because p underflows to exactly zero in the tails of the coefficient
  # distribution and log(0) poisons the respondent's whole sequence.
  p1      <- plogis(dv_draws)
  log_pch <- plogis(chosen_sign(db) * dv_draws, log.p = TRUE)
  stopifnot(identical(dim(log_pch), dim(dv_draws)))

  # Sum log-probabilities within respondent before averaging over draws, and
  # do it on the log scale: a respondent contributes nine probabilities and
  # their product underflows to zero in double precision often enough to
  # matter.
  by_resp <- rowsum(log_pch, group = db$ID, reorder = TRUE)
  mx      <- apply(by_resp, 1, max)
  ll_n    <- mx + log(rowMeans(exp(by_resp - mx)))

  p1_avg  <- rowMeans(p1)

  tibble(
    model         = model,
    sample        = label,
    respondents   = nrow(by_resp),
    observations  = nrow(db),
    LL            = sum(ll_n),
    LL_per_obs    = sum(ll_n) / nrow(db),
    rho2_vs_coin  = 1 - sum(ll_n) / (-nrow(db) * log(2)),
    hit_rate      = mean((p1_avg > 0.5) == (db$choice == 1L)),
    share_alt1_observed  = mean(db$choice == 1L),
    share_alt1_predicted = mean(p1_avg)
  )
}

# --- Exact binary logit ----------------------------------------------------
# With two always-available alternatives and generic coefficients the MNL is a
# plain binary logit in the attribute differences, and its log-likelihood is
# globally concave. Fitting it directly takes milliseconds and, crucially,
# avoids juggling Apollo's global state between two different model families
# inside one script -- 04 estimates a mixed logit, and swapping
# apollo_randCoeff in and out to also fit an MNL is precisely the state leak
# that run_all.R uses separate processes to prevent.
#
# 02 checks this against Apollo (and against mlogit) on the full sample.
fit_binary_logit <- function(db, attrs = ROUTE_ATTRS) {
  X <- attr_diff_matrix(db, attrs)
  y <- as.numeric(db$choice == 1L)

  nll <- function(b) {
    dv <- as.vector(X %*% b)
    # plogis(dv, log.p = TRUE) is log(1/(1+exp(-dv))) computed without
    # overflowing for large |dv|.
    -sum(y * plogis(dv, log.p = TRUE) + (1 - y) * plogis(-dv, log.p = TRUE))
  }
  gr <- function(b) {
    dv <- as.vector(X %*% b)
    -as.vector(crossprod(X, y - plogis(dv)))
  }

  opt <- optim(rep(0, length(attrs)), nll, gr, method = "BFGS",
               control = list(reltol = 1e-12, maxit = 500))
  if (opt$convergence != 0) {
    warning("fit_binary_logit did not converge (code ", opt$convergence, ")")
  }
  setNames(opt$par, paste0("b_", attrs))
}
