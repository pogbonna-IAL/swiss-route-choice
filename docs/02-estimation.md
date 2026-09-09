# 2. Estimation

How each likelihood is actually maximised, and the machinery the awkward cases
need.

---

## 2.1 The optimiser

All models are estimated by Apollo using **BGW** (Bunch, Gay & Welsch), a
trust-region quasi-Newton method for nonlinear least squares and maximum
likelihood. Apollo supplies analytical gradients where it can and computes the
covariance matrix from a numerical Jacobian of the analytical gradient.

Two covariance matrices are reported:

- **Classical**, $-\mathbf{H}^{-1}$ from the Hessian at the optimum.
- **Robust** (sandwich), $\mathbf{H}^{-1} \mathbf{G} \mathbf{H}^{-1}$, valid
  when the likelihood is misspecified.

For the MNL the robust errors run 1.6–1.9× the classical ones. That ratio is
reported explicitly as `se_inflation` in `02_mnl_estimates.csv`, and it has a
clear interpretation: the classical errors treat a respondent's nine tasks as
nine independent observations, and the inflation factor is the price of that
pretence. Willingness-to-pay standard errors in `02` are computed on the
**robust** matrix for that reason.

---

## 2.2 The panel product

Within `apollo_probabilities`, the order of operations matters:

```r
P[["model"]] <- apollo_mnl(mnl_settings, functionality)   # per-task probability
P <- apollo_panelProd(P, apollo_inputs, functionality)    # product over tasks
P <- apollo_avgInterDraws(P, apollo_inputs, functionality) # average over draws (04 only)
P <- apollo_prepareProb(P, apollo_inputs, functionality)
```

For the MNL, `apollo_panelProd` does not change the estimates — the likelihood
factorises over rows either way. It is stated explicitly because it is the
structure the mixed logit depends on, and because getting it wrong there is
not a rounding difference but a different model
(see [§2.3](#23-simulated-maximum-likelihood)).

For latent class models the product is taken **inside** each class before the
classes are mixed:

```r
for (s in 1:K) {
  P[[paste0("Class_", s)]] <- apollo_mnl(mnl_settings, functionality)
  P[[paste0("Class_", s)]] <- apollo_panelProd(P[[paste0("Class_", s)]], ...)
}
P[["model"]] <- apollo_lc(lc_settings, apollo_inputs, functionality)
```

This encodes the substantive claim: a respondent belongs to **one** class for
all nine tasks. Mixing per task instead would say class membership is redrawn
every question, which is a different and far weaker model.

---

## 2.3 Simulated maximum likelihood

The mixed logit has no closed form. The respondent's sequence probability is
an integral over the taste distribution:

$$L_n = \int \left[\prod_{t=1}^{T_n} P(y_{nt} \mid \beta)\right] f(\beta \mid \theta)\, d\beta$$

approximated by averaging over $R$ draws:

$$\hat{L}_n = \frac{1}{R}\sum_{r=1}^{R} \prod_{t=1}^{T_n} P(y_{nt} \mid \beta_{n}^{(r)})$$

### The order of operations is the model

**Product over tasks first, then average over draws.** Reversing them,

$$\frac{1}{R}\sum_r \prod_t P \quad\text{vs}\quad \prod_t \frac{1}{R}\sum_r P$$

gives two different models. The first says each respondent has one taste
vector used for all nine answers; the second integrates each task
independently and throws away the panel structure entirely — it collapses to
something close to an MNL with a strange error term. `apollo_panelProd` must
precede `apollo_avgInterDraws`, and that is why the comment sits above it in
the code.

### Inter- versus intra-respondent draws

```r
apollo_draws <- list(
  interDrawsType = "mlhs", interNDraws = 200,
  interNormDraws = c("d_tt", "d_tc", "d_hw", "d_ch"),
  intraNDraws    = 0
)
```

`interNDraws = 200` with `intraNDraws = 0` says a respondent's taste vector is
drawn once and persists across their nine tasks. Intra-respondent draws would
redraw it every task — a claim about within-person instability rather than
between-person heterogeneity.

### Why MLHS rather than Halton

**MLHS** (Modified Latin Hypercube Sampling) stratifies each dimension into
$R$ equal-probability intervals, samples one point per interval, and shuffles
the dimensions independently. Halton sequences are cheaper but their
low-dimensional projections become visibly correlated as dimensions increase —
a well-documented failure that produces systematically biased estimates. With
four random parameters, MLHS is the safer default.

### Simulation bias

$\log(\hat{L}_n)$ is a **biased** estimator of $\log(L_n)$ even though
$\hat{L}_n$ is unbiased for $L_n$, because $\log$ is concave: by Jensen's
inequality the simulated log-likelihood is biased *downward*, and the bias
shrinks as $R$ grows.

This is measurable here. Scoring the fitted MXL-C against its own in-sample
data with plain i.i.d. draws:

| Draws | Simulated LL | Apollo's value |
|---|---|---|
| 200 | $-1421.20$ | $-1404.47$ |
| 1000 | $-1408.89$ | $-1404.47$ |
| 4000 | $-1409.38$ | $-1404.47$ |

Apollo reaches $-1404.47$ with only 200 draws because MLHS draws are
stratified; plain i.i.d. draws need roughly five times as many to get close.
`04` therefore estimates with 200 MLHS draws and *scores* with 2000 i.i.d.
draws, and the docstring records the calibration above so the choice is
auditable rather than arbitrary.

**This table is also a test.** A scoring routine that cannot reproduce the
model's own in-sample likelihood is broken, and that is exactly how the
single-draw bug in [Numerical methods §6.2](06-numerical-methods.md#62-ifelse-silently-drops-your-draws)
was found.

---

## 2.4 Multi-start search

The MNL likelihood is globally concave. **Latent class likelihoods are not.**
They are multi-modal, and a single set of starting values routinely converges
to a local optimum that looks perfectly respectable — it converges, it
produces standard errors, and it is simply not the best solution.

The only practical defence is to start many times and keep the best.

### The algorithm — `lc_estimate_best()` in `lc_helpers.R`

```
for each start s in 1..S:
    validate inputs
    estimate with hessianRoutine = "none"      # skip the covariance matrix
    record (converged, LL, error message)
    if LL > best so far: best <- this fit

re-validate inputs
re-estimate from best$estimate WITH the covariance matrix
record whether the covariance matrix is usable
```

Four design decisions worth stating:

1. **The search pass skips the Hessian.** Computing a covariance matrix at
   every one of fifty starts roughly doubles the cost for information that is
   discarded forty-nine times out of fifty. Only the winner pays for it.

2. **Errors are kept, not swallowed.** A start that fails because the
   specification is wrong must not be indistinguishable from ordinary
   non-convergence. `run_log` carries the error message, and both `05` and `06`
   write the whole per-start log to disk (`05_lccov_runs.csv`,
   `06_lc_runs.csv`).

3. **Inputs are re-validated before the final fit.** `05` used to re-estimate
   the winning start against whatever `apollo_inputs` happened to be left in
   the global environment by the last draw that validated. That was safe only
   by accident. One shared driver now makes it safe by construction.

4. **A failed covariance matrix does not lose the fit.** In small samples the
   Hessian often fails even when the point estimates are fine. The result
   records `hessian_ok = FALSE` and keeps the point estimates, because "the
   optimiser converged but the covariance matrix is unusable" is precisely the
   distinction `07_lc_stability.R` is built to measure.

### Two different starting-value schemes

`lc_helpers.R` provides both, and which one is right depends on the model.

**`lc_start_values()` — perturb the MNL solution.** Class $k$ starts at the
full-sample MNL vector scaled by a class-level factor and jittered per
coefficient:

$$\beta^{(0)}_{ak} = \beta^{\text{MNL}}_a \cdot s_k \cdot e^{\eta_{ak}}, \qquad s_k \sim U(0.4, 2.0),\ \eta \sim \mathcal{N}(0, 0.25^2)$$

The scale factor is strictly positive, so every start begins with all four
coefficients correctly signed. Used by `06` and `08`.

**`lc_start_uniform()` — naive dispersed starts**, $U(-0.1, 0.1)$ on every
parameter. Used by `05`.

The choice is empirical, and the code records the experiment: on LC2, uniform
starts reached the best optimum in **8 of 20** tries against **1 of 20** for
the perturbation scheme. Anchoring on the MNL biases the search toward the
MNL-like local optimum — which is exactly the one worth avoiding, since a
latent class model that collapses toward the MNL has found nothing.

The first draw of `lc_start_values` is deterministic (`seq(0.5, 1.8, length.out = K)`,
no jitter) so at least one start is reproducible independently of the RNG
stream.

---

## 2.5 Generating Apollo's model code

Apollo **parses the source** of `apollo_probabilities` and `apollo_lcPars` to
work out which parameters exist and how they are used. A $K$-class model
therefore cannot be assembled at run time with `get()` and `assign()` — the
static checks do not follow that indirection, and validation fails.

The solution is to write the functions out as literal source text for a given
$K$ and evaluate them into the global environment:

```r
lc_make_lcPars <- function(K) {
  code <- c(
    "apollo_lcPars <- function(apollo_beta, apollo_inputs) {",
    "  lcpars <- list()",
    sprintf('  lcpars[["b_%s"]] <- list(%s)', a, paste0("b_", a, "_", 1:K, collapse = ", ")),
    ...
  )
  eval(parse(text = paste(code, collapse = "\n")), envir = globalenv())
}
```

`lc_install(K)` installs the pair for the constant-only model;
`lc_install_cov(K)` installs the covariate-allocation variant. Both clear a
stale `apollo_lcPars` when $K = 1$, so a one-class model cannot inherit an
allocation function from a previous iteration of the loop.

This code lived in **two places** — `lc_helpers.R` and a verbatim copy inside
`06_lc_multiclass.R` — meaning two definitions of the function that writes the
likelihood could drift apart with nothing to catch it. There is now one copy.

### Guards on the generated code

- `lc_make_lcPars(K)` and `lc_make_lcPars_cov(K)` both require $K \ge 2$.
  With one class there is nothing to allocate, and the covariate version's
  `for (k in 2:K)` would run *backwards* at $K = 1$ (`2:1` is `c(2, 1)`),
  emitting code for a `class_2` that does not exist. Refusing is the only
  correct behaviour; the callers route $K = 1$ to the no-allocation path.
- `tests/testthat/test-lc-codegen.R` asserts the generated source contains
  `V[["class_1"]] <- 0`, every $\delta_k$ and every $\gamma_{ck}$ for
  $k \ge 2$, and that the one-class model contains no `apollo_lc(` call.

---

## 2.6 Restricted models by parameter fixing

`03` needs five nested models. Rather than write five utility functions — five
chances to introduce a discrepancy — it writes **one** and estimates the
restrictions by pinning parameters:

```r
estimate_variant <- function(fix, name) {
  beta <- apollo_beta
  beta[fix] <- 0
  assign("apollo_fixed", fix, envir = globalenv())
  ...
}
```

Every model in the comparison is then *provably* the same specification under a
constraint, and the LR degrees of freedom are just `length(fix)`. The
fully-restricted case is checked against `02`'s independently estimated
log-likelihood, so the claim is verified rather than assumed.

### Why the global environment

`apollo_estimate()` requires `apollo_beta`, `apollo_fixed`,
`apollo_probabilities` and `apollo_inputs` to be resolvable by name in the
calling frame. Helper functions therefore `assign(..., envir = globalenv())`
rather than passing locals. This is Apollo's design, not a choice — and it is
the reason `run_all.R` runs every script in its **own process**
(see [Architecture §7.3](07-architecture.md#73-process-isolation)).

---

## 2.7 Convergence and identification checks

Assertions that stop the pipeline rather than write a wrong number:

| Script | Check |
|---|---|
| `02` | $\text{LL}(0) = -N_{\text{obs}} \log 2$ exactly — the null of equal shares |
| `02` | `mlogit` reproduces every coefficient to 6 dp and the LL to $10^{-10}$ |
| `03` | Fixing all twenty gammas reproduces `02`'s LL to $10^{-6}$ |
| `06` | The arithmetic parameter count $4K + (K-1)$ equals what Apollo estimated |
| `09` | The efficient design beats the median of 200 random designs |

The `06` check deserves comment. `lc_n_par(K)` is an arithmetic definition used
throughout the comparison tables, while `length(model$estimate) - 1` is what
Apollo actually estimated. They are computed independently and must agree; if
someone changes the specification and updates only one, the pipeline stops
instead of silently producing AIC and BIC values computed on the wrong
parameter count.
