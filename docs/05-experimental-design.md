# 5. Experimental design

`09_idefix_design.R` answers the question the analysis provokes: given what
these travellers turned out to care about, what should the *next* survey ask?

---

## 5.1 What "efficient" means

A choice experiment's design determines how much information each respondent's
answers carry. Formally, the design $\mathbf{X}$ determines the Fisher
information matrix $\mathbf{I}(\beta, \mathbf{X})$, and the asymptotic
covariance of the estimator is its inverse. **D-efficiency** minimises

$$D\text{-error} = \det\!\big(\mathbf{I}(\beta, \mathbf{X})^{-1}\big)^{1/p}$$

the geometric mean of the eigenvalues of the covariance matrix — a scalar
summary of "how small are the confidence ellipsoids".

### The circularity, and the way out

$\mathbf{I}$ depends on $\beta$, the very thing the survey is meant to
estimate. A design optimised for the wrong $\beta$ can be **worse than
random**, so priors are not a technicality here; they are the whole method.

Two standard responses:

- **Locally optimal**: assume a single $\beta_0$. Efficient if you are right,
  fragile if you are not.
- **Bayesian D-optimal**: assume a *distribution* over $\beta$ and minimise the
  expected D-error,
  $$D_B\text{-error} = \int \det\!\big(\mathbf{I}(\beta, \mathbf{X})^{-1}\big)^{1/p} f(\beta)\,d\beta$$
  approximated by averaging over draws.

This project uses the Bayesian version, because it has something better than a
guess: a fitted model with a covariance matrix.

---

## 5.2 Propagating the priors exactly

idefix dummy-codes each attribute with the **first level as reference**, so the
design parameters are utility differences between level $j$ and level 1:

$$\beta^{\text{design}}_{a,j} = \beta_a \cdot (\ell_{aj} - \ell_{a1})$$

That is a **linear map** $\mathbf{A}$ of the four MNL coefficients, so both
moments follow exactly:

$$\boldsymbol{\mu} = \mathbf{A}\,\hat{\beta}, \qquad \boldsymbol{\Sigma} = \mathbf{A}\,\hat{\mathbf{V}}\,\mathbf{A}^{\!\top}$$

with $\hat{\mathbf{V}}$ the **robust** covariance matrix from `02`.

```r
A <- matrix(0, nrow = 4 * (N_LEVELS - 1L), ncol = 4)
row <- 1L
for (i in seq_along(ROUTE_ATTRS)) {
  lv <- levels_list[[ROUTE_ATTRS[i]]]
  for (j in 2:N_LEVELS) {
    A[row, i] <- lv[j] - lv[1]
    row <- row + 1L
  }
}
prior_mean <- as.vector(A %*% as.numeric(b_mnl))
prior_cov  <- A %*% V_mnl %*% t(A)
```

Using $\mathbf{A}\hat{\mathbf{V}}\mathbf{A}^\top$ rather than independent
guesses matters because it carries the **correlation** between the cost and
time coefficients — the correlation that drives every willingness-to-pay
ratio. A design robust to independent uncertainty in each coefficient is not
robust to the uncertainty that actually exists.

Draws come from a Cholesky factorisation:

```r
L_chol    <- t(chol(prior_cov + diag(1e-10, nrow(prior_cov))))
par_draws <- t(prior_mean + L_chol %*% matrix(rnorm(p * N_DRAWS), nrow = p))
```

The ridge is defensive: $\mathbf{A}$ has full row rank only if no two level
gaps coincide, which is not guaranteed, and `chol()` requires positive
definiteness. At $10^{-10}$ it keeps the factorisation defined without
materially moving the draws.

### Choosing the levels

Three levels per attribute, at the 10th, 50th and 90th percentiles of the
observed data:

| Attribute | Levels |
|---|---|
| Travel time | 10, 37, 118 min |
| Cost | 3, 11, 47 CHF |
| Headway | 15, 30, 60 min |
| Interchanges | 0, 1, 2 |

Wider levels buy statistical efficiency — the information matrix grows with
attribute spread — at the cost of asking respondents about journeys that do not
exist. Quantiles keep every level inside observed experience. Interchanges are
a small count, where quantiles would collapse onto repeated values, so the
observed integers are used directly.

---

## 5.3 The modified Fedorov algorithm

Finding the optimal design is a combinatorial problem: choose 24 profiles (12
sets × 2 alternatives) from the $3^4 = 81$ candidates. Exhaustive search is
out of reach.

**Modified Fedorov** is a greedy exchange heuristic:

```
start from a random design
repeat:
    for each position in the design:
        for each candidate profile:
            if swapping it in lowers the D-error: keep the swap
until no swap improves the design
```

It converges to a **local** optimum. That is the same problem the latent class
likelihood has, and it gets the same answer: restart from `n.start = 12`
random designs and keep the best.

```r
design <- Modfed(
  cand.set  = Profiles(lvls = rep(3, 4), coding = rep("D", 4)),
  n.sets    = 12, n.alts = 2,
  alt.cte   = rep(0, 2),        # unlabelled: no ASCs, matching the estimation model
  par.draws = par_draws,
  n.start   = 12,
  parallel  = FALSE
)
```

### `parallel = FALSE` is deliberate

`Modfed` parallelises across starts by default, and the workers draw their
random start designs from **RNG streams that `set.seed()` in the parent does
not reach**. The result is a design that cannot be reproduced from the seed.

For a script whose entire purpose is to replace an irreproducible artefact
with a reproducible one, that is disqualifying. Serial is slower (498 s here)
and correct.

---

## 5.4 Evaluating the result

A D-error means nothing on its own — only against an alternative. The
comparison is 200 randomly drawn designs of the same shape, which is what a
survey built without optimisation would give.

| Statistic | Value |
|---|---|
| Bayesian D-error, chosen design | **1.115** |
| Bayesian D-error, random median (200) | 15.380 |
| Bayesian D-error, random best (200) | 3.531 |
| Improvement vs median random | **92.75%** |
| A-error | 9.474 |
| Orthogonality | 0.267 |
| Attribute-level overlap rate | **0.000** |

The optimised design is roughly fourteen times more efficient than a typical
random one, and better than the best of 200 random draws by a factor of three.

```r
stopifnot(unname(best$DB.error) <= median(random_errors, na.rm = TRUE))
```

A design that cannot beat the median random draw is not worth the compute, and
the assertion stops the pipeline rather than writing out a table that looks
authoritative.

### Zero overlap

**Attribute-level overlap** is the share of (set, attribute) pairs where both
alternatives show the same level. An overlapping attribute contributes nothing
to identifying its coefficient in that task — the respondent is not being asked
to trade it off. Zero overlap across all 12 sets and 4 attributes means every
task asks about every attribute.

### Level balance

Near-uniform but not exact, which is correct: D-efficiency does not require
balance, and forcing it would cost information.

| Attribute | Level shares |
|---|---|
| Travel time | 0.333 / 0.417 / 0.250 |
| Cost | 0.333 / 0.417 / 0.250 |
| Headway | 0.292 / 0.375 / 0.333 |
| Interchanges | 0.292 / 0.375 / 0.333 |

### Decoding for humans

The optimiser works on the dummy matrix; nobody can read a survey off it.
`decode_row()` inverts the coding back to the levels a respondent would
actually see, and `09_design.csv` is the design as it would be fielded:

| Set | Alt | tt | tc | hw | ch |
|---|---|---|---|---|---|
| 1 | 1 | 37 | 3 | 15 | 1 |
| 1 | 2 | 10 | 11 | 30 | 2 |
| 2 | 1 | 118 | 11 | 60 | 0 |
| 2 | 2 | 37 | 47 | 30 | 2 |

An all-zero dummy row means the reference level, which is why `decode_row`
tests `sum(d) == 0` before looking for a `1`.

---

## 5.5 Provenance

`outputs/models/idefix_design_D.rds` **existed in this project before
`09_idefix_design.R` did.** It was produced in an ad-hoc session and no code
generated it — precisely the kind of artefact a replication must not contain,
because there is no way to check it, change it, or explain it.

This script regenerates it from a fixed seed and documented priors. The recipe
matches the orphaned object's structure — 12 choice sets, 2 alternatives, 4
attributes at 3 levels, dummy coded, no ASC, Bayesian D-optimality, 12 random
starts — but the numbers differ, because the original seed is unrecoverable.

That trade is deliberate. **A reproducible design that differs is worth more
than an irreproducible one that does not.**
