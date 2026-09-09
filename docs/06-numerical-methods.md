# 6. Numerical methods

Every bug in this project that survived code review was a numerical or shape
issue that produced **plausible-looking numbers rather than an error**. None
of them threw. Each is documented here with its symptom, because recognising
the symptom is the only reliable defence.

---

## 6.1 Probabilities that are exactly zero

### The problem

The obvious way to write a logit probability is

```r
p1 <- 1 / (1 + exp(-dv))
```

In IEEE double precision, `exp(709.8)` is the largest representable value and
`exp(710)` is `Inf`. So for `dv < -745` or so:

```r
1 / (1 + exp(745))   # -> 1 / Inf -> 0    exactly
log(0)               # -> -Inf
```

The likelihood becomes `-Inf` and the entire respondent is lost.

### Why this is not hypothetical here

The mixed logit draws coefficients from a **lognormal**, which has a heavy
right tail on $|\beta|$. With $\mu_{ch} = 1.17$ and $\sigma = 1.50$, the 95th
percentile of $|\beta_{ch}|$ is about 38; interchange differences reach 2, so
that term alone contributes $\pm 77$. The cost coefficient's log-scale standard
deviation is 1.89 and cost differences reach 44 CHF. Observed
$\max|\Delta_{nt}|$ across draws exceeded **60 000**.

Extreme utility differences are not a bug in the model. They are what a
lognormal tail *means*: a small number of respondents who are nearly
deterministic. The arithmetic has to survive them.

### The fix

`plogis(x, log.p = TRUE)` computes $\log \Lambda(x)$ in one step, without ever
forming the probability:

```r
1 / (1 + exp(-(-800)))      # 0
log(1 / (1 + exp(800)))     # -Inf
plogis(-800, log.p = TRUE)  # -800     exact
```

For large negative $x$, $\log \Lambda(x) \to x$, which is what the function
returns. Both scorers use it.

### The symptom to recognise

`NA` or `-Inf` in a likelihood column while the **hit rate looks fine**. Hit
rate depends only on the sign of $\Delta$, which is unaffected by overflow, so
a broken likelihood sits next to a plausible accuracy figure. That combination
should always be investigated.

---

## 6.2 `ifelse()` silently drops your draws

### The problem

This is the most dangerous line of code in the project's history:

```r
pch <- ifelse(db$choice == 1L, p1, 1 - p1)
```

`ifelse()` returns a result **shaped like its `test` argument**. `db$choice == 1L`
is a plain vector of length $N_{\text{obs}}$. When `p1` is an
$N_{\text{obs}} \times R$ matrix of draws, the result is a **vector of length
$N_{\text{obs}}$** containing only the first column.

No error. No warning. The mixed logit was being scored with **one draw out of
five hundred**.

```r
M    <- matrix(1:12, nrow = 4)      # 4 obs x 3 draws
test <- c(TRUE, FALSE, TRUE, FALSE)
dim(ifelse(test, M, -M))            # NULL   <- collapsed
length(ifelse(test, M, -M))         # 4      <- draw 1 only
dim(ifelse(test, 1, -1) * M)        # 4 3    <- correct
```

### The symptom to recognise

**A Monte Carlo estimate that does not change when you add draws.**

| Draws | Simulated LL |
|---|---|
| 200 | $-21211.53$ |
| 1000 | $-21211.53$ |
| 4000 | $-21211.53$ |

Byte-identical across a twentyfold increase in $R$ is impossible for a
simulation. Any Monte Carlo quantity should move — either converging or
drifting — and one that does not is not a Monte Carlo quantity.

The internal tell was `dim(by_resp)` returning `388 1` where `388 500` was
expected.

### The fix

Multiply by a sign, which preserves matrix shape, and assert the shape:

```r
chosen_sign <- function(db) ifelse(db$choice == 1L, 1, -1)   # vector, fine
log_pch <- plogis(chosen_sign(db) * dv_draws, log.p = TRUE)
stopifnot(identical(dim(log_pch), dim(dv_draws)))
```

`chosen_sign` still uses `ifelse`, safely: its arguments are scalars, so the
result is a vector of length $N_{\text{obs}}$, and `vector * matrix` recycles
down the rows — which is exactly right, because the sign varies by observation
and not by draw.

### Why the test suite missed it

The original test scored identical draws and checked the result matched the
closed form:

```r
flat <- panel_logit_score(matrix(rep(dv, 5), ncol = 5), db, "flat")
expect_equal(flat$LL, binary_logit_score(dv, db, "exact")$LL)
```

**It passed against the bug.** When every draw is identical, draw 1 *is* the
right answer. The test verified the log-sum-exp arithmetic and was blind to the
thing that was actually broken.

The replacements fail on it:

```r
# 1. An analytic two-draw answer the collapsed version cannot produce
dv <- cbind(rep(0, 6), ifelse(db$choice == 1L, 4, -4))
expected <- 2 * log((0.5^3 + plogis(4)^3) / 2)
expect_equal(panel_logit_score(dv, db, "two draws")$LL, expected)
expect_false(isTRUE(all.equal(s$LL, 2 * log(0.5^3))))   # the buggy answer

# 2. The likelihood must move when draws are added
expect_false(isTRUE(all.equal(ll_at_R5, ll_at_R200)))
```

**The lesson generalises:** a test whose fixture is degenerate in the same
dimension as the bug cannot catch the bug. Varying draws were essential, and
the original fixture deliberately removed that variation to make the expected
value easy to compute.

---

## 6.3 Log-sum-exp

### The problem

The panel likelihood needs

$$\hat\ell_n = \log\left(\frac{1}{R}\sum_r e^{s_{nr}}\right), \qquad s_{nr} = \sum_t \log P(y_{nt} \mid \beta^{(r)}_n)$$

$s_{nr}$ is a sum of nine log-probabilities, so it is routinely $-50$ or lower
and can reach $-700$. Exponentiating directly underflows every term to zero and
returns $\log(0) = -\infty$.

### The fix

Factor out the row maximum:

$$\log \sum_r e^{s_r} = m + \log \sum_r e^{s_r - m}, \qquad m = \max_r s_r$$

```r
by_resp <- rowsum(log_pch, group = db$ID, reorder = TRUE)
mx      <- apply(by_resp, 1, max)
ll_n    <- mx + log(rowMeans(exp(by_resp - mx)))
```

The largest term becomes exactly $e^0 = 1$, so nothing overflows, and terms far
below the maximum underflow to zero — which is harmless, since they contribute
nothing to the sum anyway.

Note `mx` has length $N_{\text{resp}}$ and `by_resp` is
$N_{\text{resp}} \times R$; R recycles the vector **down the columns**, which
is the correct orientation. This is the same recycling rule that makes
`chosen_sign(db) * dv_draws` work, and the same rule that makes `ifelse` fail —
worth internalising.

### Summing on the log scale is also the point

Note that `rowsum` operates on `log_pch`, not on probabilities. Computing
$\prod_t P$ first and then taking the log would underflow before the sum is
ever formed: nine probabilities of $10^{-300}$ multiply to $10^{-2700}$, which
is zero as a double. A test asserts a finite result for exactly that case.

---

## 6.4 List versus atomic subsetting

### The problem

```r
x <- list(a = 1);         x[["missing"]]   # NULL
x <- c(a = 1);            x[["missing"]]   # Error: subscript out of bounds
```

Apollo returns `model$estimate` as a **named numeric vector**. The helper

```r
lc_par_value <- function(est, name) {
  if (is.null(est[[name]])) 0 else unname(est[[name]])
}
```

worked on the synthetic list fixtures in the test suite and **threw** on real
Apollo output.

This is not an edge case in the covariate allocation model: class 1 has no
$\delta_1$ and no $\gamma_{c1}$ *by construction* — its whole allocation index
is normalised to zero — so looking up an absent parameter is the normal path,
not a mistake.

### The fix

Test membership rather than relying on the return value:

```r
lc_par_value <- function(est, name) {
  if (!name %in% names(est)) return(0)
  unname(est[[name]])
}
```

The regression test passes `unlist(fake_lc3(covariate_alloc = TRUE))` — an
atomic vector — and asserts `expect_type(vec, "double")` so the fixture cannot
silently revert to a list.

---

## 6.5 Cholesky parameterisation

A covariance matrix must be positive semi-definite. Estimating $\Sigma$
directly means the optimiser can propose invalid matrices, which then need
constraints or projection.

Estimating the **lower-triangular Cholesky factor** $\mathbf{L}$ removes the
problem: $\mathbf{L}\mathbf{L}^{\!\top}$ is positive semi-definite for *any*
real $\mathbf{L}$. There is no constraint to violate.

$$\log|\beta_n| = \boldsymbol\mu + \mathbf{L}\mathbf{z}_n, \qquad \mathbf{z}_n \sim \mathcal{N}(\mathbf{0}, \mathbf{I})$$

```r
randcoeff[["b_tt"]] <- -exp(mu_tt + L_tt_tt * d_tt)
randcoeff[["b_tc"]] <- -exp(mu_tc + L_tc_tt * d_tt + L_tc_tc * d_tc)
randcoeff[["b_hw"]] <- -exp(mu_hw + L_hw_tt * d_tt + L_hw_tc * d_tc + L_hw_hw * d_hw)
randcoeff[["b_ch"]] <- -exp(mu_ch + L_ch_tt * d_tt + L_ch_tc * d_tc + L_ch_hw * d_hw + L_ch_ch * d_ch)
```

Recovering the interpretable quantities:

$$\Sigma = \mathbf{L}\mathbf{L}^{\!\top}, \qquad \sigma_a = \sqrt{\Sigma_{aa}}, \qquad \rho_{ab} = \frac{\Sigma_{ab}}{\sigma_a \sigma_b}$$

**A sign caveat.** $\mathbf{L}$ is identified only up to the sign of each
column, so individual $L_{ij}$ are not directly interpretable — a column-sign
flip gives an identical $\Sigma$ and an identical likelihood. Report $\Sigma$,
$\sigma$ and $\rho$, never the raw factor entries.

### Lognormal moments

For $\beta = -e^{\mu + \sigma z}$:

$$\text{median} = -e^{\mu}, \quad \mathbb{E}[\beta] = -e^{\mu + \sigma^2/2}, \quad \text{SD} = e^{\mu + \sigma^2/2}\sqrt{e^{\sigma^2} - 1}$$

Quantiles invert directly, with the tails swapped by the negation:

$$q_{0.05}(\beta) = -e^{\mu + \sigma \Phi^{-1}(0.95)}$$

### Willingness to pay is lognormal too

$$\log\!\left(\frac{\text{VTT}}{60}\right) = \mu_{tt} - \mu_{tc} + (\mathbf{L}_{tt\cdot} - \mathbf{L}_{tc\cdot})\mathbf{z}$$

a normal random variable, so VTT is lognormal with

$$\sigma_{\text{VTT}} = \|\mathbf{L}_{tt\cdot} - \mathbf{L}_{tc\cdot}\|_2$$

and exact quantiles — no simulation needed. This is one of the concrete
benefits of choosing lognormal over normal coefficients.

---

## 6.6 Reproducibility of random numbers

| Issue | Where | Resolution |
|---|---|---|
| Parallel workers ignore `set.seed()` | `Modfed` in `09` | `parallel = FALSE` |
| Draw seeds must be stable across resumes | `08` | Seeds derived from `SEED + 1000K + seed`, never from loop position |
| Resampled respondents must be identical when a cell is re-entered | `08` | `set.seed(SEED + seed)` before sampling IDs |
| Scoring draws must not reuse estimation draws | `04` | Different seed (`SEED + 7`); reusing them would flatter the model |

The `09` case is worth dwelling on. `Modfed(parallel = TRUE)` spawns workers
whose RNG streams the parent's `set.seed()` never reaches, so the "best of 12
random starts" is different on every run. The script exists specifically to
replace an irreproducible artefact, so serial execution is not a preference
but a requirement.

---

## 6.7 Checklist

When a numerical result looks wrong, in order:

1. **Does a Monte Carlo estimate change with $R$?** If not, draws are being
   dropped.
2. **Does the scorer reproduce the model's own in-sample likelihood?** If not,
   the scorer is wrong, not the model.
3. **Are there `NA` or `-Inf` next to plausible hit rates?** Overflow.
4. **Do `dim()` calls return what you expect after `ifelse`, `apply` or
   `sapply`?** Shape collapse.
5. **Is a likelihood suspiciously worse than a coin flip while accuracy is
   good?** Directionally right, wildly overconfident — usually a scale or
   averaging error.
