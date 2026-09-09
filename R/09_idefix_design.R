# 09_idefix_design.R -- Bayesian D-efficient design for a follow-up survey
# Swiss route choice
#
# The estimation data was collected with a design this project did not choose.
# This script answers the design question the analysis raises: given what the
# MNL in 02 says respondents actually care about, what choice sets should a
# follow-up survey use?
#
# PROVENANCE NOTE. outputs/models/idefix_design_D.rds existed in this project
# before this script did -- it was produced in an ad-hoc session and no code
# generated it, which makes it exactly the kind of artefact a replication must
# not contain. This script regenerates it from a fixed seed and documented
# priors. The recipe matches the orphaned object's structure (12 choice sets,
# 2 alternatives, 4 attributes at 3 levels each, dummy coded, no ASC,
# Bayesian D-optimality, 12 random starts), but the numbers will differ
# because the original seed is unrecoverable. A reproducible design that
# differs is worth more than an irreproducible one that does not.
#
# Priors are the whole point of an efficient design: it maximises information
# ABOUT AN ASSUMED PARAMETER VECTOR, so a design built on the wrong priors can
# be worse than random. They come from 02's estimates, propagated exactly
# through the dummy coding (see prior_mean / prior_cov below) rather than
# guessed.
#
# RUN ORDER: after 02_mnl.R (reads the baseline model for the priors).
#
# Outputs
#   outputs/models/idefix_design_D.rds       full idefix design_list object
#   outputs/tables/09_design_levels.csv      attribute levels used
#   outputs/tables/09_design_priors.csv      dummy-coded priors and their s.d.
#   outputs/tables/09_design.csv             the design in attribute levels
#   outputs/tables/09_design_quality.csv     D-error, balance, overlap vs random
#   outputs/figures/09_design.png
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))
library(idefix)   # attaches shiny as a side effect; needed only here

N_SETS    <- 12L   # choice sets shown to each respondent
N_ALTS    <- 2L    # unlabelled, matching the estimation data
N_LEVELS  <- 3L    # levels per attribute
N_STARTS  <- 12L   # random starts for the modified Fedorov search
N_DRAWS   <- 100L  # prior draws for the Bayesian D-error

set.seed(SEED)

# --- Data and priors -------------------------------------------------------
data("apollo_swissRouteChoiceData", package = "apollo")
database <- apollo_swissRouteChoiceData

mnl_path <- file.path(PATH_MODELS, "Swiss_MNL_model.rds")
if (!file.exists(mnl_path)) {
  stop("Run 02_mnl.R first -- ", mnl_path, " not found.")
}
model_mnl <- readRDS(mnl_path)

b_mnl <- model_mnl$estimate[paste0("b_", ROUTE_ATTRS)]
V_mnl <- model_mnl$robvarcov[paste0("b_", ROUTE_ATTRS), paste0("b_", ROUTE_ATTRS)]

# --- Attribute levels ------------------------------------------------------
# Three levels per attribute spanning the range the estimation data actually
# covers. Going wider would buy statistical efficiency at the cost of asking
# respondents about journeys that do not exist; the quantiles keep every level
# inside observed experience.
level_of <- function(a, probs = c(0.10, 0.50, 0.90)) {
  v <- c(database[[paste0(a, 1)]], database[[paste0(a, 2)]])
  if (a == "ch") {
    # Interchanges are a small count; quantiles would collapse onto repeated
    # values, so the observed integers are used directly.
    return(sort(unique(v))[seq_len(N_LEVELS)])
  }
  round(unname(quantile(v, probs)), if (a == "tc") 1 else 0)
}

levels_list <- setNames(lapply(ROUTE_ATTRS, level_of), ROUTE_ATTRS)

levels_tbl <- map_dfr(ROUTE_ATTRS, function(a) {
  tibble(attribute = a, label = ROUTE_ATTR_LABELS[[a]],
         level_index = seq_len(N_LEVELS), value = levels_list[[a]])
})

write_table(levels_tbl, "design_levels", prefix = "09")

cat("=== Attribute levels =======================================\n")
show_table(levels_tbl)

# --- Dummy-coded priors ----------------------------------------------------
# idefix dummy codes each attribute with the FIRST level as reference, so the
# design parameters are the utility differences between level j and level 1:
#
#   beta_{a,j} = b_a * (level_j - level_1)
#
# That is a linear map A of the four MNL coefficients, so the prior mean is
# A b and the prior covariance is A V A', both exact. Drawing the Bayesian
# priors from that covariance -- rather than from independent guesses -- means
# the design is robust to the uncertainty the estimation data actually left,
# including the correlation between the cost and time coefficients that drives
# every willingness-to-pay ratio.
A <- matrix(0, nrow = length(ROUTE_ATTRS) * (N_LEVELS - 1L),
            ncol = length(ROUTE_ATTRS),
            dimnames = list(NULL, paste0("b_", ROUTE_ATTRS)))

prior_names <- character(0)
row <- 1L
for (i in seq_along(ROUTE_ATTRS)) {
  a  <- ROUTE_ATTRS[i]
  lv <- levels_list[[a]]
  for (j in 2:N_LEVELS) {
    A[row, i] <- lv[j] - lv[1]
    prior_names <- c(prior_names, sprintf("%s_lvl%d", a, j))
    row <- row + 1L
  }
}
rownames(A) <- prior_names

prior_mean <- as.vector(A %*% as.numeric(b_mnl))
prior_cov  <- A %*% V_mnl %*% t(A)
names(prior_mean) <- prior_names

# Draws from the prior. chol() needs a positive definite matrix; A has full
# row rank only if no two level gaps coincide, which is not guaranteed, so a
# tiny ridge keeps the factorisation defined without materially moving the
# draws.
ridge     <- diag(1e-10, nrow(prior_cov))
L_chol    <- t(chol(prior_cov + ridge))
par_draws <- t(prior_mean + L_chol %*% matrix(rnorm(length(prior_mean) * N_DRAWS),
                                              nrow = length(prior_mean)))
colnames(par_draws) <- prior_names

priors_tbl <- tibble(
  parameter  = prior_names,
  attribute  = rep(ROUTE_ATTRS, each = N_LEVELS - 1L),
  # vapply returns levels-in-rows by attributes-in-columns; as.vector reads
  # column-major, which is exactly the (tt2, tt3, tc2, tc3, ...) order that
  # prior_names and A use. Transposing here silently interleaved them.
  level_gap  = as.vector(vapply(ROUTE_ATTRS,
                                function(a) levels_list[[a]][-1] - levels_list[[a]][1],
                                numeric(N_LEVELS - 1L))),
  prior_mean = prior_mean,
  prior_sd   = sqrt(diag(prior_cov))
)

write_table(priors_tbl, "design_priors", prefix = "09")

cat("\n=== Dummy-coded priors (from the 02 MNL) ===================\n")
show_table(priors_tbl, digits = 4)

# --- Search ----------------------------------------------------------------
# Modified Fedorov: swap profiles in and out of the design, keeping whatever
# lowers the Bayesian D-error, restarted from N_STARTS random designs because
# the search is a local one and a single start lands in a local optimum for
# the same reason the latent class likelihood does.
cat(sprintf("\n=== Modified Fedorov search (%d starts) ====================\n",
            N_STARTS))

# Cached like the estimation scripts: the search is deterministic given the
# seed and the priors, so repeating it on every pipeline run spends five
# minutes rediscovering the same design. REFIT=1 forces a fresh search.
# This object is not a fitted model, so the cache is handled here rather
# than through cached_model().
design_path <- file.path(PATH_MODELS, "idefix_design_D.rds")

t0 <- Sys.time()
if (resume_enabled() && file.exists(design_path)) {
  design <- readRDS(design_path)
  cat(sprintf("  cached design reused (D-error %.4f) -- delete the .rds or\n",
              unname(design$BestDesign$DB.error)))
  cat("  set REFIT=1 to re-run the search\n")
} else {
design <- Modfed(
  cand.set  = Profiles(lvls = rep(N_LEVELS, length(ROUTE_ATTRS)),
                       coding = rep("D", length(ROUTE_ATTRS))),
  n.sets    = N_SETS,
  n.alts    = N_ALTS,
  alt.cte   = rep(0, N_ALTS),   # unlabelled alternatives, no ASC
  par.draws = par_draws,
  n.start   = N_STARTS,
  # Modfed parallelises across starts by default, which draws its random
  # start designs from worker RNG streams that set.seed() above does not
  # reach. Serial is slower and reproducible; this script exists because the
  # previous design was not reproducible.
  parallel  = FALSE
)
  saveRDS(design, design_path)
}
cat(sprintf("search took %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))


best <- design$BestDesign

# --- The design in attribute levels ----------------------------------------
# The dummy matrix is what the algorithm optimises; nobody can read a survey
# off it. This decodes it back to the levels a respondent would actually see.
decode_row <- function(r) {
  out <- numeric(length(ROUTE_ATTRS))
  col <- 1L
  for (i in seq_along(ROUTE_ATTRS)) {
    d <- r[col:(col + N_LEVELS - 2L)]
    # All-zero dummies mean the reference level.
    lvl <- if (sum(d) == 0) 1L else which(d == 1)[1] + 1L
    out[i] <- levels_list[[ROUTE_ATTRS[i]]][lvl]
    col <- col + N_LEVELS - 1L
  }
  out
}

decoded <- t(apply(best$design, 1, decode_row))
colnames(decoded) <- ROUTE_ATTRS

design_tbl <- as_tibble(decoded) %>%
  mutate(row = rownames(best$design), .before = 1) %>%
  separate(row, into = c("set", "alt"), sep = "\\.") %>%
  mutate(set = as.integer(str_remove(set, "set")),
         alt = as.integer(str_remove(alt, "alt"))) %>%
  arrange(set, alt)

write_table(design_tbl, "design", prefix = "09")

cat("\n=== Design (attribute levels as shown to respondents) ======\n")
show_table(design_tbl)

# --- Quality ---------------------------------------------------------------
# A D-error means nothing on its own; it only means something against an
# alternative. The comparison is against randomly drawn designs of the same
# shape, which is what a survey built without any optimisation would give.
d_error <- function(des) {
  # idefix reports the Bayesian D-error of a design directly.
  DBerr(par.draws = par_draws, des = des, n.alts = N_ALTS)
}

cand <- Profiles(lvls = rep(N_LEVELS, length(ROUTE_ATTRS)),
                 coding = rep("D", length(ROUTE_ATTRS)))

random_errors <- vapply(seq_len(200L), function(i) {
  rows <- sample(nrow(cand), N_SETS * N_ALTS, replace = TRUE)
  des  <- cand[rows, , drop = FALSE]
  rownames(des) <- rownames(best$design)
  out <- tryCatch(d_error(des), error = function(e) NA_real_)
  if (is.null(out) || length(out) != 1L) NA_real_ else out
}, numeric(1))

overlap_rate <- mean(design_tbl %>%
  group_by(set) %>%
  summarise(across(all_of(ROUTE_ATTRS), ~ .x[1] == .x[2]), .groups = "drop") %>%
  select(all_of(ROUTE_ATTRS)) %>%
  as.matrix())

balance <- map_dfr(seq_along(ROUTE_ATTRS), function(i) {
  a <- ROUTE_ATTRS[i]
  counts <- table(factor(design_tbl[[a]], levels = levels_list[[a]]))
  tibble(attribute = a,
         level = as.numeric(names(counts)),
         n = as.integer(counts),
         share = as.integer(counts) / sum(counts))
})

write_table(balance, "design_balance", prefix = "09")

quality <- tibble(
  statistic = c("Bayesian D-error (chosen design)",
                "Bayesian D-error (random, median of 200)",
                "Bayesian D-error (random, best of 200)",
                "improvement vs median random (%)",
                "A-error", "orthogonality",
                "attribute-level overlap rate",
                "choice sets", "alternatives per set", "prior draws"),
  value = c(
    unname(best$DB.error),
    median(random_errors, na.rm = TRUE),
    min(random_errors, na.rm = TRUE),
    100 * (1 - unname(best$DB.error) / median(random_errors, na.rm = TRUE)),
    unname(best$AB.error),
    unname(best$Orthogonality),
    overlap_rate,
    N_SETS, N_ALTS, N_DRAWS
  )
)

write_table(quality, "design_quality", prefix = "09")

cat("\n=== Design quality =========================================\n")
show_table(quality, digits = 4)

# A design that cannot beat the median random draw is not worth the compute.
stopifnot(unname(best$DB.error) <= median(random_errors, na.rm = TRUE))
cat("efficient design beats the median random design\n")

# --- Figure ----------------------------------------------------------------
p_err <- tibble(d = random_errors[!is.na(random_errors)]) %>%
  ggplot(aes(d)) +
  geom_histogram(bins = 30, fill = "grey70") +
  geom_vline(xintercept = unname(best$DB.error), colour = "firebrick",
             linewidth = 1) +
  labs(title = "Bayesian D-error: chosen design vs 200 random designs",
       subtitle = "Red line is the modified Fedorov solution; lower is better",
       x = "D-error", y = "random designs")

p_bal <- balance %>%
  mutate(attribute = factor(attribute, ROUTE_ATTRS,
                            unname(ROUTE_ATTR_LABELS[ROUTE_ATTRS]))) %>%
  ggplot(aes(factor(level), share)) +
  geom_col(fill = "grey40") +
  geom_hline(yintercept = 1 / N_LEVELS, linetype = "dashed", colour = "grey20") +
  facet_wrap(~ attribute, scales = "free_x") +
  labs(title = "Level balance",
       subtitle = "Dashed line is perfect balance; efficiency does not require it",
       x = NULL, y = "share of profiles")

p_des <- design_tbl %>%
  pivot_longer(all_of(ROUTE_ATTRS), names_to = "attribute", values_to = "value") %>%
  group_by(attribute) %>%
  mutate(scaled = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup() %>%
  mutate(attribute = factor(attribute, ROUTE_ATTRS,
                            unname(ROUTE_ATTR_LABELS[ROUTE_ATTRS]))) %>%
  ggplot(aes(factor(alt), factor(set), fill = scaled)) +
  geom_tile(colour = "white") +
  facet_wrap(~ attribute, nrow = 1) +
  scale_fill_viridis_c(option = "mako", direction = -1) +
  labs(title = "The design", subtitle = "Level within each attribute, rescaled to 0-1",
       x = "alternative", y = "choice set", fill = NULL)

p_prior <- priors_tbl %>%
  ggplot(aes(prior_mean, parameter)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_pointrange(aes(xmin = prior_mean - 1.96 * prior_sd,
                      xmax = prior_mean + 1.96 * prior_sd)) +
  labs(title = "Priors used by the search",
       subtitle = "MNL estimates propagated through the dummy coding, 95% interval",
       x = "utility difference vs the reference level", y = NULL)

fig <- (p_err | p_bal) / (p_des | p_prior)
fig_path <- write_figure(fig, "design", prefix = "09", width = 15, height = 10)

cat(sprintf("\nDesign object: %s\n",
            file.path(PATH_MODELS, "idefix_design_D.rds")))
cat(sprintf("Tables written to: %s\n", PATH_TABLES))
cat(sprintf("Figure written to: %s\n", fig_path))
