# 01_data_audit.R -- data audit
# Swiss route choice
#
# Establishes what the estimation data actually contains before any model is
# fitted: panel structure, missingness, attribute levels, choice shares,
# dominated tasks and non-trading respondents. Every check writes a CSV to
# outputs/tables/ so the audit is reproducible rather than console-only.
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

# --- Load ------------------------------------------------------------------
data("apollo_swissRouteChoiceData", package = "apollo")
database <- apollo_swissRouteChoiceData

# The attribute names, labels and covariate list are defined once in
# 00_setup.R; every downstream script reads the same objects, so a
# respecification cannot leave the audit describing a different model from the
# one that gets estimated.
ATTRS  <- ROUTE_ATTR_LABELS
COVARS <- ROUTE_COVARS

write_audit <- function(x, name) write_table(x, name, prefix = "01")

cat("\n=== 1. Dimensions ==========================================\n")
cat(sprintf("rows: %d   cols: %d\n", nrow(database), ncol(database)))
str(database)

# --- 2. Panel structure ----------------------------------------------------
# Mixed logit and latent class both need the likelihood taken over respondents,
# not rows, so the tasks-per-respondent count has to be established first.
cat("\n=== 2. Panel structure =====================================\n")

tasks_per_id <- database %>%
  count(ID, name = "n_tasks")

n_resp   <- nrow(tasks_per_id)
balanced <- length(unique(tasks_per_id$n_tasks)) == 1L

cat(sprintf("respondents (unique ID): %d\n", n_resp))
# Every downstream script sizes its resampling grid against FULL_N. If the
# data ever changes underneath that constant, the audit is where it should
# surface, not three scripts later in a silently mis-scaled experiment.
stopifnot(identical(as.integer(n_resp), as.integer(FULL_N)))
cat(sprintf("tasks per respondent:    %s\n",
            paste(sort(unique(tasks_per_id$n_tasks)), collapse = ", ")))
cat(sprintf("balanced panel:          %s\n", balanced))
cat(sprintf("rows == sum of tasks:    %s\n",
            identical(nrow(database), sum(tasks_per_id$n_tasks))))

write_audit(tasks_per_id, "tasks_per_respondent")
write_audit(count(tasks_per_id, n_tasks, name = "n_respondents"),
            "panel_balance")

# --- 3. Missingness and duplicates -----------------------------------------
cat("\n=== 3. Missingness / duplicates ============================\n")

missingness <- tibble(
  variable  = names(database),
  n_missing = colSums(is.na(database)),
  n_unique  = sapply(database, function(x) length(unique(x)))
)
print(as.data.frame(missingness), row.names = FALSE)

n_dup_rows <- sum(duplicated(database))
cat(sprintf("\nfully duplicated rows: %d\n", n_dup_rows))

# A constant column would silently break identification of its coefficient.
constant_cols <- missingness$variable[missingness$n_unique == 1L]
cat(sprintf("constant columns: %s\n",
            if (length(constant_cols)) paste(constant_cols, collapse = ", ") else "none"))

write_audit(missingness, "missingness")

# --- 4. Choice distribution ------------------------------------------------
cat("\n=== 4. Choice shares =======================================\n")

choice_shares <- database %>%
  count(choice, name = "n") %>%
  mutate(share = n / sum(n))
print(as.data.frame(choice_shares), row.names = FALSE)

# Unlabelled alternatives: shares near 50/50 mean no strong left/right bias.
cat(sprintf("\nchoice values observed: %s (expected 1, 2)\n",
            paste(sort(unique(database$choice)), collapse = ", ")))

write_audit(choice_shares, "choice_shares")

# --- 5. Attribute levels and descriptives ----------------------------------
# Two kinds of attribute live in this design and they need different treatment:
# hw and ch are fixed design levels shared by every respondent, while tt and tc
# were pivoted off each respondent's own reference trip and so take hundreds of
# distinct values. Enumerating levels is only meaningful for the former.
cat("\n=== 5. Attribute descriptives ==============================\n")

MAX_ENUM <- 12L  # above this an attribute is treated as continuous

attr_summary <- map_dfr(names(ATTRS), function(a) {
  map_dfr(ALTS, function(alt) {
    v <- database[[paste0(a, alt)]]
    n_lv <- length(unique(v))
    tibble(
      attribute   = a,
      label       = ATTRS[[a]],
      alternative = alt,
      type        = if (n_lv <= MAX_ENUM) "design level" else "continuous",
      min         = min(v),
      median      = median(v),
      mean        = round(mean(v), 2),
      max         = max(v),
      sd          = round(sd(v), 2),
      n_levels    = n_lv,
      levels      = if (n_lv <= MAX_ENUM) paste(sort(unique(v)), collapse = " ") else ""
    )
  })
})
print(as.data.frame(attr_summary %>% select(-label)), row.names = FALSE)

write_audit(attr_summary, "attribute_summary")

# Discrete attributes should draw on the same level set for both alternatives
# in an unlabelled design; the continuous ones will not match and need not.
cat("\nlevel sets identical across alternatives (discrete attributes only):\n")
for (a in names(ATTRS)) {
  n_lv <- length(unique(database[[paste0(a, 1)]]))
  if (n_lv > MAX_ENUM) {
    cat(sprintf("  %-3s -- continuous, not applicable\n", a))
    next
  }
  same <- identical(sort(unique(database[[paste0(a, 1)]])),
                    sort(unique(database[[paste0(a, 2)]])))
  cat(sprintf("  %-3s %s\n", a, same))
}

# Which attributes are continuous drives the plotting below.
CONT_ATTRS <- attr_summary %>%
  filter(type == "continuous") %>%
  pull(attribute) %>%
  unique()
DISC_ATTRS <- setdiff(names(ATTRS), CONT_ATTRS)
cat(sprintf("\ncontinuous: %s | discrete: %s\n",
            paste(CONT_ATTRS, collapse = ", "), paste(DISC_ATTRS, collapse = ", ")))

# --- 6. Attribute correlations ---------------------------------------------
# Strong correlation between attributes inflates standard errors and is the
# first thing to check if a coefficient comes out wrong-signed.
cat("\n=== 6. Attribute correlations ==============================\n")

attr_cols <- as.vector(t(outer(names(ATTRS), ALTS, paste0)))
attr_cor  <- cor(database[, attr_cols])
print(round(attr_cor, 3))

write_audit(
  as.data.frame(round(attr_cor, 4)) %>% rownames_to_column("variable"),
  "attribute_correlations"
)

# --- 7. Dominated tasks ----------------------------------------------------
# A task where one alternative is weakly better on all four attributes carries
# no trade-off information; a respondent choosing the dominated option is a
# candidate inattentive response.
cat("\n=== 7. Dominance ===========================================\n")

dom <- database %>%
  transmute(
    ID, choice,
    a1_weakly_better = tt1 <= tt2 & tc1 <= tc2 & hw1 <= hw2 & ch1 <= ch2,
    a2_weakly_better = tt2 <= tt1 & tc2 <= tc1 & hw2 <= hw1 & ch2 <= ch1,
    a1_strictly_any  = tt1 <  tt2 | tc1 <  tc2 | hw1 <  hw2 | ch1 <  ch2,
    a2_strictly_any  = tt2 <  tt1 | tc2 <  tc1 | hw2 <  hw1 | ch2 <  ch1
  ) %>%
  mutate(
    dominant = case_when(
      a1_weakly_better & a1_strictly_any ~ 1L,
      a2_weakly_better & a2_strictly_any ~ 2L,
      TRUE ~ NA_integer_
    ),
    dominated_choice = !is.na(dominant) & choice != dominant
  )

n_dom <- sum(!is.na(dom$dominant))
cat(sprintf("tasks with a dominant alternative: %d (%.1f%% of %d)\n",
            n_dom, 100 * n_dom / nrow(dom), nrow(dom)))
cat(sprintf("dominated alternative chosen:      %d (%.1f%% of dominated tasks)\n",
            sum(dom$dominated_choice),
            if (n_dom > 0) 100 * sum(dom$dominated_choice) / n_dom else 0))

dom_by_id <- dom %>%
  group_by(ID) %>%
  summarise(n_dominant_tasks = sum(!is.na(dominant)),
            n_failed         = sum(dominated_choice), .groups = "drop") %>%
  arrange(desc(n_failed))
write_audit(dom_by_id, "dominance_by_respondent")

# --- 8. Non-traders --------------------------------------------------------
# Respondents who never switch alternative contribute no within-person
# variation; they drive degenerate classes in latent class models.
cat("\n=== 8. Non-trading respondents =============================\n")

traders <- database %>%
  group_by(ID) %>%
  summarise(n_tasks    = n(),
            n_alt1     = sum(choice == 1),
            share_alt1 = mean(choice == 1),
            .groups    = "drop") %>%
  mutate(non_trader = share_alt1 %in% c(0, 1))

cat(sprintf("non-traders: %d of %d respondents (%.1f%%)\n",
            sum(traders$non_trader), nrow(traders),
            100 * mean(traders$non_trader)))
cat(sprintf("  always alt 1: %d | always alt 2: %d\n",
            sum(traders$share_alt1 == 1), sum(traders$share_alt1 == 0)))

write_audit(traders, "choice_by_respondent")

# --- 9. Respondent covariates ----------------------------------------------
cat("\n=== 9. Covariates ==========================================\n")

# Covariates are respondent-level, so summarise one row per ID, not per task.
covar_by_id <- database %>%
  group_by(ID) %>%
  summarise(across(all_of(COVARS), first), .groups = "drop")

# Verify they really are constant within respondent before trusting that.
covar_varies <- database %>%
  group_by(ID) %>%
  summarise(across(all_of(COVARS), ~ length(unique(.x)) > 1L), .groups = "drop") %>%
  select(-ID) %>%
  summarise(across(everything(), sum))
cat("respondents whose covariate changes across tasks (should be 0):\n")
print(as.data.frame(covar_varies), row.names = FALSE)

cat("\nincome (hh_inc_abs) across respondents:\n")
print(summary(covar_by_id$hh_inc_abs))

purpose <- c("commute", "shopping", "business", "leisure")
purpose_tab <- covar_by_id %>%
  summarise(across(all_of(c("car_availability", purpose)), sum)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_respondents") %>%
  mutate(share = n_respondents / nrow(covar_by_id))
cat("\ncovariate counts (respondent level):\n")
print(as.data.frame(purpose_tab), row.names = FALSE)

# The purpose dummies should partition respondents exactly once.
purpose_sum <- rowSums(covar_by_id[, purpose])
cat(sprintf("\npurpose dummies summing to 1: %d of %d respondents\n",
            sum(purpose_sum == 1), nrow(covar_by_id)))
if (any(purpose_sum != 1)) {
  cat("  observed row sums: ",
      paste(sort(unique(purpose_sum)), collapse = ", "), "\n", sep = "")
}

write_audit(covar_by_id, "covariates_by_respondent")
write_audit(purpose_tab, "covariate_shares")

# --- 10. Figures -----------------------------------------------------------
cat("\n=== 10. Figures ============================================\n")

long_attrs <- database %>%
  mutate(.row = row_number()) %>%
  select(.row, all_of(attr_cols)) %>%
  pivot_longer(-.row, names_to = "var", values_to = "value") %>%
  mutate(attribute   = str_sub(var, 1, 2),
         alternative = paste("alt", str_sub(var, 3, 3)))

lab <- function(a) unname(ATTRS[a])

# Continuous attributes need a histogram; the discrete ones read better as
# counts per level, so they are drawn separately rather than force-fitted.
p_cont <- long_attrs %>%
  filter(attribute %in% CONT_ATTRS) %>%
  mutate(attribute = factor(attribute, CONT_ATTRS, lab(CONT_ATTRS))) %>%
  ggplot(aes(value, fill = alternative)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.55) +
  facet_wrap(~ attribute, scales = "free") +
  labs(title = "Continuous attributes (pivoted off each reference trip)",
       x = NULL, y = "tasks", fill = NULL)

p_disc <- long_attrs %>%
  filter(attribute %in% DISC_ATTRS) %>%
  mutate(attribute = factor(attribute, DISC_ATTRS, lab(DISC_ATTRS))) %>%
  ggplot(aes(factor(value), fill = alternative)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ attribute, scales = "free_x") +
  labs(title = "Design levels", x = NULL, y = "tasks", fill = NULL)

# Choice against the attribute difference shows whether the expected trade-off
# is visible in the raw data before any model is estimated. Continuous
# differences are binned into quantiles so each point rests on similar support.
# Quantile breaks collapse when a difference is heavily tied (tc differences
# cluster on a few small integers), so dedupe the breaks and fall back to the
# raw values whenever the attribute has few enough of them to stand alone.
bin_difference <- function(x, n_bins = 9L) {
  if (length(unique(x)) <= n_bins) return(factor(x))
  breaks <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1L),
                            type = 1L, names = FALSE))
  if (length(breaks) < 3L) return(factor(x))
  cut(x, breaks = breaks, include.lowest = TRUE)
}

summarise_difference <- function(a) {
  database %>%
    transmute(chose_alt1 = as.integer(choice == 1),
              difference = .data[[paste0(a, 1)]] - .data[[paste0(a, 2)]]) %>%
    mutate(bin = bin_difference(difference)) %>%
    group_by(bin) %>%
    summarise(share_alt1 = mean(chose_alt1), n = n(), .groups = "drop") %>%
    mutate(attribute = a, .before = 1)
}

# Each attribute is plotted separately rather than faceted: the bin labels of
# hw and ch overlap numerically, and a shared discrete scale would union their
# level sets and scramble the ordering within each panel.
diff_plot <- function(a) {
  ggplot(summarise_difference(a), aes(bin, share_alt1)) +
    geom_col() +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    ylim(0, 1) +
    labs(
      # patchwork drops annotations on a nested assembly, so the explanation
      # rides on the first panel and the y label appears once, on the left.
      title    = if (a == names(ATTRS)[1]) {
        sprintf("Share choosing alt 1 by attribute difference (alt 1 - alt 2): %s",
                ATTRS[[a]])
      } else ATTRS[[a]],
      subtitle = if (a == names(ATTRS)[1]) {
        "Downward slope is the expected trade-off; dashed line is indifference"
      } else NULL,
      x = NULL,
      y = if (a == names(ATTRS)[1]) "share choosing alt 1" else NULL
    )
}

diffs  <- map_dfr(names(ATTRS), summarise_difference)
p_diff <- wrap_plots(map(names(ATTRS), diff_plot), nrow = 1)

write_audit(diffs %>% mutate(bin = as.character(bin)), "choice_by_difference")

p_alt1 <- ggplot(traders, aes(share_alt1)) +
  geom_histogram(bins = 11) +
  labs(title = "Share of tasks where alt 1 was chosen",
       subtitle = "Mass at 0 and 1 marks non-trading respondents",
       x = "share alt 1", y = "respondents")

# Panel balance is reported in the summary table; when it holds, a
# tasks-per-respondent bar is a single column and tells the reader nothing, so
# the slot goes to the income distribution the covariate models will use.
p_second <- if (balanced) {
  ggplot(covar_by_id, aes(hh_inc_abs)) +
    geom_histogram(bins = 8) +
    scale_x_continuous(labels = scales::label_comma()) +
    labs(title = "Household income (respondent level)",
         subtitle = sprintf("Balanced panel: %d respondents x %d tasks",
                            n_resp, unique(tasks_per_id$n_tasks)),
         x = "CHF", y = "respondents")
} else {
  ggplot(tasks_per_id, aes(n_tasks)) +
    geom_bar() +
    labs(title = "Tasks per respondent",
         subtitle = "Unbalanced panel", x = "tasks", y = "respondents")
}

fig <- ((p_cont / p_disc) | (p_second / p_alt1)) / p_diff +
  plot_layout(heights = c(2, 1))
fig_path <- write_figure(fig, "data_audit", prefix = "01",
                         width = 17, height = 11)
cat(sprintf("figure written: %s\n", fig_path))

# --- 11. Audit summary -----------------------------------------------------
audit_summary <- tibble(
  check = c("rows", "columns", "respondents", "tasks per respondent",
            "balanced panel", "missing values", "duplicated rows",
            "constant columns", "share choosing alt 1",
            "tasks with dominant alternative", "dominated alternative chosen",
            "non-trading respondents"),
  value = c(nrow(database), ncol(database), n_resp,
            paste(sort(unique(tasks_per_id$n_tasks)), collapse = "/"),
            balanced, sum(is.na(database)), n_dup_rows,
            length(constant_cols),
            sprintf("%.3f", mean(database$choice == 1)),
            n_dom, sum(dom$dominated_choice), sum(traders$non_trader))
)
write_audit(audit_summary, "audit_summary")

cat("\n=== Audit summary ==========================================\n")
print(as.data.frame(audit_summary), row.names = FALSE)
cat(sprintf("\nTables written to: %s\n", PATH_TABLES))
