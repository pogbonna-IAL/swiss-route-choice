# 3. Post-estimation

Turning a converged model into numbers a reader can use. Three problems:
latent class labels are arbitrary, willingness to pay is a ratio, and both
have standard errors that are easy to get wrong.

---

## 3.1 Label switching

A $K$-class model's likelihood is **invariant to permuting the class labels**.
Swap every parameter of class 1 with class 2 and the log-likelihood is
identical. The labels carry no information; only the partition does.

This is harmless until you compare or average across fits:

- `08` averages parameters over 20 resampled replications per cell. If class 1
  is the time-sensitive class in one replication and the price-sensitive class
  in the next, the average is a blend of two distinct classes and every bias
  and RMSE figure is meaningless.
- `05` and `06` are meant to be read together — `06` says how big each class
  is, `05` says who is in it. If "class 1" denotes different classes in the two
  scripts, the joint reading is nonsense.

### The fix: canonical ordering

Sort classes by the travel time coefficient, most negative first:

```r
lc_class_order <- function(est, K) {
  b_tt <- vapply(1:K, function(k) est[[paste0("b_tt_", k)]], numeric(1))
  order(b_tt)
}
```

Any rule that is a deterministic function of the parameters works; $\beta_{tt}$
is chosen because it is significant in every class of every model fitted here,
so the ordering never has to break a near-tie.

`lc_class_order()` returns the **permutation itself**, not just the reordered
estimates. That matters because anything derived from the raw model — standard
errors, delta-method valuations — has to be relabelled the same way, and
recomputing the ordering separately in each place invites divergence.

Applied in three places, all through the same function:

| Function | Used by |
|---|---|
| `lc_canonical()` | `08`, constant-only |
| `lc_canonical_cov()` | `05`, covariate allocation |
| `lc_canonical_params()` | `05` and `06`, for the reported tables |

---

## 3.2 Re-referencing the allocation, and why the standard errors change

Reordering the within-class coefficients is pure relabelling: $\beta_{a,k}$
becomes $\beta_{a,\text{ord}[k]}$ and its standard error travels with it
unchanged.

The **allocation** parameters are different, because the reference class moves.

### Constant-only allocation

Only differences between the $\delta$'s are identified, and $\delta_1 \equiv 0$
by normalisation. After permuting by $\text{ord}$, class $j = \text{ord}[1]$
becomes the new reference, so every constant must be re-referenced:

$$\delta'_k = \delta_{\text{ord}[k]} - \delta_{j}$$

### Covariate allocation

Here class 1's **entire index** is normalised to zero — there is no $\delta_1$
and no $\gamma_{c1}$ in the parameter vector at all. Re-referencing means
subtracting the new reference class's whole allocation index:

$$W'_{nk} = W_{nk} - W_{nj} = (\delta_k - \delta_j) + \sum_c (\gamma_{ck} - \gamma_{cj})\, z_{nc}$$

Because the index is linear in $z$, this works coefficient by coefficient:

$$\delta'_k = \delta_k - \delta_j, \qquad \gamma'_{ck} = \gamma_{ck} - \gamma_{cj}$$

Every class-share *difference* is unchanged, so the likelihood is untouched.
`tests/testthat/test-lc-canonical.R` verifies this on 25 random covariate
vectors: the canonicalised model's membership probabilities must equal the
raw model's, permuted.

### The standard errors

**This is the part that is easy to get wrong.** After re-referencing, each
allocation parameter is no longer an estimate — it is a **contrast** between
two estimates. Its variance is therefore

$$\operatorname{Var}(\theta_k - \theta_j) = V_{kk} + V_{jj} - 2 V_{kj}$$

not $V_{kk}$. Simply permuting the reported standard errors alongside the
estimates would misstate the uncertainty by however much the two parameters
covary — and allocation parameters in a latent class model always covary,
because they compete for the same respondents.

```r
lc_contrast_se <- function(V, a, b) {
  va  <- if (is.na(a) || !a %in% rownames(V)) 0 else V[a, a]
  vb  <- if (is.na(b) || !b %in% rownames(V)) 0 else V[b, b]
  cab <- if (is.na(a) || is.na(b) ||
             !a %in% rownames(V) || !b %in% rownames(V)) 0 else V[a, b]
  sqrt(max(va + vb - 2 * cab, 0))
}
```

Two details:

- **Structurally zero parameters contribute nothing.** `delta_1` is fixed by
  `apollo_fixed` (and absent from the covariance matrix); the whole class-1
  index is absent in the covariate model. Both cases fall through to the
  `%in% rownames(V)` guard and contribute zero variance and zero covariance,
  which is correct — they are constants, not estimates.
- **`max(..., 0)` before the square root.** When $j = \text{ord}[k]$ the
  contrast is a parameter with itself and the expression is algebraically zero;
  floating point can make it a very small negative number, and `sqrt` of that
  is `NaN`.

`lc_canonical_params()` assembles the whole tidy table — permuted betas with
permuted standard errors, re-referenced allocation parameters with contrast
standard errors — and both `05` and `06` call it, so the two scripts cannot
disagree about what class 1 means or how uncertain its parameters are.

---

## 3.3 The delta method

Coefficients are in utils, which is not a unit anyone can act on. Ratios
against the cost coefficient convert them to money:

$$\text{VTT} = 60 \cdot \frac{\beta_{tt}}{\beta_{tc}} \quad \text{CHF per hour}$$

A ratio of two estimates is **not normally distributed**, and its uncertainty
is far wider than either coefficient suggests alone. The delta method gives a
first-order approximation: for $g(\theta)$,

$$\operatorname{Var}\big(g(\hat\theta)\big) \approx \nabla g(\hat\theta)^{\!\top} \, \mathbf{V} \, \nabla g(\hat\theta)$$

`apollo_deltaMethod` does this symbolically from an expression string.

### Where it is applied

| Script | Quantities |
|---|---|
| `02` | Five population valuations, on the **robust** covariance matrix |
| `03` | Implied coefficients and VTT for seven respondent profiles |
| `06` | Four valuations **per class** for every $K$ |

`06` previously reported bare point ratios with no uncertainty at all — which
cannot support the only interesting claim a multi-class model makes, namely
that two classes value time *differently*. `lc_valuations()` now returns
value, standard error, t-ratio and a 95% interval per class.

### Profile valuations in `03`

The profile covariate values are data, not parameters, so they are substituted
into the expression as numeric literals:

```r
coef_expression <- function(a, x) {
  terms <- c(paste0("b_", a),
             sprintf("%s*%.10g", paste0("g_", a, "_", MODEL_COVARS), unlist(x)))
  paste(terms, collapse = " + ")
}
# -> "b_tt + g_tt_inc*0.4213 + g_tt_car*1 + g_tt_com*0 + ..."
```

The resulting VTT carries the uncertainty of the whole linear combination, not
just of $\beta_{tt}$. Profiles are keyed `prof1..profN` in the expression
labels rather than by display name, because `apollo_deltaMethod` round-trips
those labels through a data frame and names containing spaces and brackets do
not survive.

### When the ratio is meaningless

If a class's cost coefficient is not distinguishable from zero, the ratio
explodes and its delta-method interval is uninformative. `lc_valuations()`
flags this from the t-ratio on $\beta_{tc}$ rather than reporting the number as
if it meant something:

```r
ratio_reliable = abs(tc_t[raw_class]) >= t_min   # default 1.96
```

Reliability is a property of the **denominator**, so it is computed once per
class rather than per valuation.

---

## 3.4 Class shares

Under constant-only allocation the shares are a softmax of the deltas,
identical for every respondent:

$$\pi_k = \frac{e^{\delta_k}}{\sum_l e^{\delta_l}}$$

Under covariate allocation they vary by person, and `lc_class_shares_cov()`
evaluates them for a given covariate vector:

$$\pi_{nk} = \frac{\exp(W_{nk})}{\sum_l \exp(W_{nl})}, \qquad W_{n1} = 0$$

`05` uses this to build **archetype profiles** — the average respondent, ±1 SD
of log income, a car owner, a commuter, a business traveller — because that is
what an allocation model is actually for: reading membership off
characteristics rather than off a coefficient table.

The income SD is taken over **respondents**, not over rows:

```r
inc_sd <- database %>%
  group_by(ID) %>% summarise(log_income = first(log_income)) %>%
  pull(log_income) %>% sd()
```

The panel is balanced at nine tasks each, so the two coincide here — but that
is a property of this dataset, not something the code should depend on.

---

## 3.5 Fit statistics

$$\text{AIC} = -2\ell + 2p, \qquad \text{BIC} = -2\ell + p \log N_{\text{obs}}$$

Two conventions worth stating because they are not universal:

- **$p$ counts free parameters only.** $\delta_1$ is fixed and excluded. `06`
  cross-checks the arithmetic count against Apollo's.
- **$N_{\text{obs}} = 3492$ rows, not 388 respondents.** BIC's penalty is
  therefore $\log 3492 = 8.16$ rather than $\log 388 = 5.96$. Using respondents
  would penalise parameters less and shift the selection toward larger models.
  Rows is the convention in the discrete choice literature and is used
  consistently across every model here, which is what makes the twelve-model
  BIC table comparable.

$\rho^2$ is reported against equal shares, $\ell_0 = -N_{\text{obs}} \log 2$,
which is the right null for two unlabelled alternatives with a 49.7/50.3 split.
