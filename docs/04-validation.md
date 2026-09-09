# 4. Validation

Every model here is compared to every other on likelihood, and likelihood
always improves with parameters. Four independent checks keep that honest.

---

## 4.1 The independent cross-check

Apollo can only tell you that *its own* optimiser converged. It cannot tell
you that the utility function says what you meant — a transposed alternative
or a dropped attribute gives a perfectly self-consistent Apollo run.

`02` therefore re-estimates the identical specification in **mlogit**: a
different package, by different authors, from a differently shaped dataset.

```r
mlogit_long <- bind_rows(
  database %>% transmute(obs = row_number(), alt = 1L, chosen = choice == 1L,
                         tt = tt1, tc = tc1, hw = hw1, ch = ch1),
  database %>% transmute(obs = row_number(), alt = 2L, chosen = choice == 2L,
                         tt = tt2, tc = tc2, hw = hw2, ch = ch2)
) %>% arrange(obs, alt)

model_mlogit <- mlogit(chosen ~ tt + tc + hw + ch | 0, data = mlogit_idx)
```

Apollo takes the data **wide** (one row per task, attributes suffixed by
alternative); mlogit takes it **long** (one row per alternative). Reshaping is
part of the test: a mistake in the wide layout will not survive being melted
and re-fitted. The `| 0` suppresses alternative-specific constants to match
`apollo_beta`.

Both maximise the same globally concave likelihood, so they must agree to
optimiser tolerance:

```r
stopifnot(ll_diff < 1e-4, max(crosscheck$rel_diff) < 1e-3)
```

Observed: coefficients identical to six decimal places, log-likelihood
differing by $8.07 \times 10^{-11}$.

Anything larger is a specification difference, not numerical noise, and it
stops the pipeline rather than being written to a table nobody reads.

---

## 4.2 Hold-out design

In-sample fit cannot distinguish a model that generalises from one that has
fitted noise. `02`, `03` and `04` are all scored on 194 respondents they never
saw.

### Split on respondent, never on task

```r
holdout_split <- function(db, frac = 0.5, seed = SEED) {
  set.seed(seed)
  ids       <- unique(db$ID)
  train_ids <- sample(ids, floor(frac * length(ids)))
  list(train = db %>% filter(ID %in% train_ids),
       test  = db %>% filter(!ID %in% train_ids), train_ids = train_ids)
}
```

**This is the single most important design decision in the validation.**
Putting some of a person's nine tasks in training and the rest in the hold-out
would let the model exploit that person's own revealed taste — and taste
heterogeneity is precisely what `03`, `04`, `05` and `06` are trying to
justify. A task-level split would hand every one of them a free win.

`tests/testthat/test-model-helpers.R` asserts the two sides share no
respondent, that respondent counts and row counts partition exactly, and that
every respondent keeps all nine tasks.

All three scripts use the **same** function with the **same** seed, so their
validation tables refer to the same 194 people and are comparable by
construction rather than by coincidence.

### Baselines must be fitted on the training half

An early version of `04` scored a training-half mixed logit against the
**full-sample** MNL — a baseline that had already seen every hold-out
respondent. Beating it would have proved nothing.

The fix uses `fit_binary_logit()`, an exact maximum-likelihood binary logit in
the attribute differences:

```r
nll <- function(b) {
  dv <- as.vector(X %*% b)
  -sum(y * plogis(dv, log.p = TRUE) + (1 - y) * plogis(-dv, log.p = TRUE))
}
gr <- function(b) -as.vector(crossprod(X, y - plogis(dv)))
optim(rep(0, 4), nll, gr, method = "BFGS", control = list(reltol = 1e-12))
```

Analytical gradient, globally concave, milliseconds to fit. It avoids swapping
Apollo's global model definition mid-script — exactly the state leak that
`run_all.R` uses separate processes to prevent — and a test asserts it
reproduces `02`'s Apollo estimates to $10^{-5}$.

---

## 4.3 Scoring

### Closed form (`02`, `03`)

For a binary logit the chosen alternative's probability has a closed form:

$$\log P(y_{nt}) = \log \Lambda(s_{nt}\, \Delta_{nt}), \qquad s_{nt} = \begin{cases} +1 & y_{nt} = 1 \\ -1 & y_{nt} = 2\end{cases}$$

```r
chosen_sign <- function(db) ifelse(db$choice == 1L, 1, -1)
log_pch <- plogis(chosen_sign(db) * dv, log.p = TRUE)
```

### Panel simulated (`04`)

The mixed logit likelihood does not factorise over tasks, so scoring must
mirror estimation: product within respondent, then average over draws.

$$\hat\ell_n = \log\left(\frac{1}{R}\sum_{r=1}^{R} \exp\Big(\sum_{t} \log P(y_{nt} \mid \beta_n^{(r)})\Big)\right)$$

Computed with the **log-sum-exp** trick, because a respondent's nine
probabilities multiply to a number that underflows to zero:

```r
by_resp <- rowsum(log_pch, group = db$ID, reorder = TRUE)   # respondents x draws
mx      <- apply(by_resp, 1, max)
ll_n    <- mx + log(rowMeans(exp(by_resp - mx)))
```

Subtracting the row maximum before exponentiating guarantees the largest term
is exactly 1, so nothing overflows and the smallest terms underflow harmlessly
to zero. See [Numerical methods](06-numerical-methods.md) for why this is not
optional.

### What is reported, and why more than one number

| Metric | What it catches |
|---|---|
| `LL_per_obs` | Overall fit, comparable across samples of different size |
| `rho2_vs_coin` | Same, scaled against the 50/50 null |
| `hit_rate` | Share of choices predicted correctly |
| `share_alt1_predicted` vs `observed` | Aggregate calibration |

These are reported **together** because they disagree in informative ways. A
model can recover the aggregate market share while getting every individual
task wrong, so share alone is worthless. And a model can have a good hit rate
while being wildly overconfident — high hit rate with a terrible likelihood is
the signature of predictions that are directionally right and far too certain.
That exact combination is what exposed the single-draw bug.

### Results

| Model | Sample | LL/obs | $\rho^2$ | Hit rate |
|---|---|---|---|---|
| MNL | training | $-0.4729$ | 0.318 | 0.793 |
| MNL | **hold-out** | $-0.4823$ | 0.304 | 0.791 |
| MNL + covariates | **hold-out** | $-0.4618$ | 0.334 | 0.803 |
| MXL-C | **hold-out** | $\mathbf{-0.4037}$ | **0.418** | 0.799 |

The MNL's training-to-hold-out gap is only 0.0094 LL per observation — four
parameters cannot overfit 1746 observations. Both richer models improve on it
out of sample, so their extra parameters are buying real structure rather than
noise. The mixed logit's margin is large.

This ordering matches the BIC ranking, which is reassuring but not
guaranteed — information criteria penalise parameters by a formula, the
hold-out penalises them by whether they help.

---

## 4.4 Assertions in the pipeline

| Script | Assertion | Catches |
|---|---|---|
| `01` | $N$ respondents $=$ `FULL_N` | Data changing under a hardcoded constant |
| `02` | $\ell(0) = -N \log 2$ | Misspecified utility |
| `02` | mlogit agreement | Coding error Apollo cannot see |
| `03` | Restricted model reproduces `02` | Interactions not entering as claimed |
| `06` | Parameter count matches Apollo's | Comparison table desynchronised from the models |
| `09` | Beats the median random design | Design search failing silently |

---

## 4.5 The small-sample experiment

`08_small_sample_experiment.R` is the largest computation in the project — 485
model fits — and answers the question the rest of the analysis provokes: **how
much of this survives a smaller study?**

### Design

- **Sample sizes** 250, 150, 100, 75, 50, 30 respondents; **20 seeds** each.
- **Models** MNL, LC2, LC3, LC4 in every replication.
- **Sampling is on ID**, so a drawn respondent keeps all nine tasks.
- **$N = 388$ is not a sampled cell.** Drawing 388 from 388 returns the whole
  panel for every seed. It is estimated once, with seed 0, as the benchmark
  that bias and coverage are measured against.
- **Every fit is canonicalised** before anything is recorded
  (see [Post-estimation §3.1](03-post-estimation.md#31-label-switching)).

### Resume protocol

Results are appended to CSV after each fit and completed `(n, seed, K)`
combinations are skipped on restart. Keying on the **triple** rather than the
cell matters: widening `K_SET` to add LC4 to an experiment that already ran
MNL/LC2/LC3 would otherwise mark every existing cell incomplete and
re-estimate hundreds of models that are already correct and on disk.

Resampled respondents are drawn from `SEED + seed`, so a cell re-entered to add
a model gets exactly the same people as its existing fits.

### What it finds

**Coefficient recovery (MNL, $N = 250$):** bias essentially zero, coverage
0.90–1.00. The MNL is well-behaved all the way down.

**Latent class stability:**

| $N$ | $K$ | Usable covariance | Median smallest class | Sign reversals /20 | Verdict |
|---|---|---|---|---|---|
| 250 | 2 | 100% | 43.1% | 0 | Yes |
| 250 | 4 | 90% | 16.9% | 0 | Mostly |
| 150 | 3 | 95% | 20.9% | 2 | Weak |
| 150 | 4 | 60% | 14.8% | 4 | No |
| 100 | 2 | 100% | 40.4% | 4 | No |
| 100 | 4 | 55% | 14.5% | 11 | No |
| 50 | 4 | 20% | 10.1% | 14 | meaningless |
| 30 | 4 | 15% | 13.2% | 18 | meaningless |

A **sign reversal** is a replication in which any class returns a positive
coefficient on time, cost, headway or interchanges — an impossibility that can
only mean the class is fitting noise. At $N = 100$, four of twenty two-class
fits do this. At $N = 30$ with four classes, eighteen of twenty do.

**Class selection collapses too.** Share of replications in which BIC picks
each $K$:

| $N$ | $K=1$ | $K=2$ | $K=3$ | $K=4$ |
|---|---|---|---|---|
| 250 | – | – | 0.05 | **0.95** |
| 150 | – | 0.05 | 0.35 | **0.60** |
| 100 | – | 0.15 | **0.80** | 0.05 |
| 75 | – | 0.45 | 0.45 | 0.10 |
| 50 | 0.10 | **0.50** | 0.35 | 0.05 |
| 30 | **0.45** | 0.35 | 0.20 | – |

BIC recovers the full-sample answer ($K = 4$) in 95% of replications at
$N = 250$ and in **none** at $N = 30$, where it picks a single class nearly
half the time. Small-sample BIC does not merely get noisier; it is
systematically biased toward too few classes.

### Two kinds of instability, never mixed

`07_lc_stability.R` assembles the verdict table, and its central design point
is that the $N = 388$ row measures something **different** from the others:

- **$N < 388$: sampling instability.** The same model refitted to 20 different
  draws of respondents. This is what a researcher with one dataset of that
  size would face without ever knowing it.
- **$N = 388$: starting-value instability.** The full panel has no sampling
  variation — it is every respondent there is. What varies is only the
  optimiser's starting point.

Reading a single "convergence rate" across both would be meaningless. The rows
are labelled by `source`, carry their own verdict scale ("20 starts enough" /
"50+ starts" / "100+ starts" rather than "Yes" / "No"), and columns that do not
apply print `--` instead of a number that invites a false comparison.

The table also flags $K = 5$ as having **no resampling counterpart**: `06`
sweeps five classes while `08` covers four, so that row rests on
starting-value evidence alone.

### Convergence is not usability

The experiment's most consequential finding is not a number, it is a category
error the code used to make. **Every single small-sample fit converges.**
Across 560 fits at n from 25 down to 5, not one returned a non-finite
log-likelihood. The optimiser reports success, hands back a log-likelihood,
class shares and a parameter vector, and the pipeline printed:

```
   K=4  converged LL =     -84.794  starts ok 8   21.4s
```

At n = 20 that word appeared over four-class fits whose Hessian was singular
in **every** replication. Convergence means the optimiser stopped moving. It
says nothing about whether two classes collapsed onto each other, whether a
class holds one respondent, whether a coefficient came back economically
impossible, or whether the covariance matrix exists at all.

`lc_fit_verdict()` is the single definition of usable, and 05, 06 and 08 all
route through it so no script can quietly disagree:

| Flag | Meaning |
|---|---|
| `no_fit` | No start produced a finite log-likelihood |
| `no_covariance` | Hessian singular — every s.e., t-ratio and interval is **undefined** |
| `wrong_sign` | A route attribute has a positive coefficient; the class is fitting noise |
| `degenerate_class` | A class holds under 5% of respondents and cannot be interpreted |
| `classes_collapsed` | Two classes agree within 10% on every attribute — one class fitted twice |

Class separation is measured **scale-free** — the largest per-attribute
difference divided by the larger magnitude — because `b_ch` is roughly thirty
times `b_hw`, and a raw distance would be driven entirely by interchanges.

The console output is now unmissable rather than a buffered `warning()` that
R prints at the end of a run under a "There were 50 or more warnings" line:

```
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !! UNUSABLE -- n=20 seed=3 K=4 -- NO COVARIANCE
     !!   the Hessian is singular -- every standard error, t-ratio and
     !!   confidence interval in this fit is UNDEFINED
     !! Any number derived from this fit is not a result.
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```

The marker is greppable on purpose: `grep UNUSABLE` over a run log finds every
one. The verdict is also carried into `08_fits.csv`, `06_lc_comparison.csv` and
`05_lccov_comparison.csv`, so a reader of the tables sees what a reader of the
console saw &mdash; a comparison table that ranks models on BIC without saying
which of them are estimable is a trap. `06` now says so explicitly when BIC's
pick is not usable.

### A note on what `converged` means

`08` records `ok = TRUE` whenever *any* start returns a finite likelihood, and
that is `1.000` in all 24 resampling cells. It is retained only as a check
that no cell failed outright — it is **not** a convergence rate and the table
says so. The honest measures are `mean_starts_ok` (how many of the eight
attempted starts converged) and `usable` (share of replications with a usable
covariance matrix), and the verdict rests on the latter.
