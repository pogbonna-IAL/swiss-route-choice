# 8. What broke

The latent class models were pushed from 250 respondents down to 5 — twenty
resampled panels at each size, four models in each, **1,044 fits in total**.
This is what did and did not reproduce.

The samples are **nested**: for a given seed the N = 5 panel is a subset of the
N = 30 panel is a subset of the N = 250 panel. Shrinking N removes respondents
from a fixed ordering rather than drawing an unrelated group, so every question
below is about the same people, progressively fewer.

---

## 8.1 The headline

> **1,044 fits. 1,044 converged. 440 usable.**

Not one fit failed to return a finite log-likelihood — not at 5 respondents,
not with 19 parameters fitted to 45 observations. Every single one produced a
log-likelihood, a parameter vector and a set of class shares. **The failure is
never loud.** In 17 of the 52 (N, K) cells, not one replication of twenty was
usable, and every fit in those cells converged.

How the 604 unusable fits failed:

| Verdict | Count |
|---|---|
| Singular Hessian — every standard error undefined | 438 |
| Wrong-signed coefficient | 166 |
| No fit at all | **0** |

This is why the pipeline no longer prints `converged`. It prints the verdict.

---

## 8.2 What reproduced

**The MNL, and almost nothing else.**

| N | obs/par | Usable | VTT median | VTT within 50% of benchmark |
|---|---|---|---|---|
| 250 | 562 | 100% | 27.9 | 100% |
| 150 | 338 | 100% | 27.5 | 95% |
| 100 | 225 | 100% | 28.0 | 95% |
| 75 | 169 | 100% | 28.8 | 90% |
| 50 | 113 | 95% | 29.0 | 85% |
| 30 | 68 | 90% | 24.4 | 65% |
| 15 | 34 | 75% | 15.2 | 45% |
| 5 | 11 | 60% | 12.3 | 30% |

The full-panel value of travel time is 27.21 CHF/h. The MNL's median tracks it
to within 6% down to **N = 30** and remains usable in 60% of replications at
**N = 5** — four parameters on 45 observations. Its 95% intervals cover the
benchmark 95% of the time at N = 250 and 150.

That is the entire list of things that reproduce. Note what it implies: the
degradation below is **not** simply "few observations". The MNL has the same
few observations and survives. It is a property of the latent class structure.

Two apparent successes that are not evidence:

- **The likelihood ordering LC4 > LC3 > LC2 > LC1 holds at every N.** It is
  mechanically guaranteed by nesting — more classes cannot fit worse. It says
  nothing about whether the classes mean anything.
- **Classes never collapsed onto each other.** The `classes_collapsed`
  diagnostic fires in **0%** of replications at every sample size. The classes
  stay numerically distinct all the way down while becoming meaningless. The
  obvious symptom never appears, which is precisely why it is not a sufficient
  check.

---

## 8.3 What did not reproduce, and in what order

The failures arrive in a definite sequence. The largest N at which each
pathology affects more than a fifth of replications:

| Failure | K = 2 | K = 3 | K = 4 |
|---|---|---|---|
| **Classes do not match the benchmark's** | N = 30 | **N = 250** | **N = 250** |
| Coefficient leaves the plausible range | N = 25 | N = 75 | N = 150 |
| Singular Hessian | N = 25 | N = 50 | N = 150 |
| Wrong-signed coefficient | N = 50 | N = 100 | N = 100 |
| Class too thin to interpret | — | N = 5 | N = 10 |

**Read the first row.** The first thing to break is not convergence, not the
standard errors, and not the signs. It is **class identity**, and it breaks at
N = 250 — a sample size nobody would think twice about, where every
conventional diagnostic is green.

At N = 250 the three-class model is usable in 100% of replications. It also
finds a materially different three-class partition from the full panel in 85%
of them. Comparing the benchmark against a typical N = 250 replication:

| | b_tt | b_tc | b_hw | b_ch |
|---|---|---|---|---|
| **Benchmark class 1** | −0.341 | −1.04 | −0.069 | −3.93 |
| **Benchmark class 2** | −0.227 | −1.83 | −0.047 | −1.72 |
| Replication class 1 | −0.206 | −1.75 | −0.049 | −1.82 |
| Replication class 2 | −0.131 | −0.31 | −0.054 | −2.48 |

The replication's class 1 is the benchmark's class 2. Its class 2 is not
clearly either. Canonical ordering by `b_tt` is doing its job — the classes are
sorted — but sorting cannot align two fits when the underlying partitions
differ, and at N = 250 they already differ.

### Coverage fails earliest and worst

Share of replications whose nominal **95%** interval contains the full-panel
value:

| N | MNL | LC2 | LC3 | LC4 |
|---|---|---|---|---|
| 250 | 0.95 | 0.54 | 0.38 | 0.60 |
| 150 | 0.95 | 0.59 | 0.42 | 0.34 |
| 100 | 0.91 | 0.56 | 0.46 | 0.33 |
| 50 | 0.76 | 0.49 | 0.37 | 0.12 |
| 30 | 0.88 | 0.42 | 0.34 | 0.08 |
| 15 | 0.71 | 0.36 | 0.10 | 0.00 |

At N = 250, where LC2 is usable in every replication and shows no other
pathology, its 95% intervals are really about 54% intervals. The MNL's are
honest at 95%. **Nothing in the ordinary output distinguishes the two.**

### Individual classes

Share of replications recovering each benchmark class — the closest fitted
class within 50% of that class's own magnitude, irrespective of labelling:

| N | LC3 class 1 | class 2 | class 3 |
|---|---|---|---|
| 250 | 40% | 60% | 80% |
| 150 | 60% | 55% | 80% |
| 100 | 45% | 50% | 70% |
| 50 | 35% | 25% | 35% |
| 30 | 25% | 10% | 30% |
| 10 | 15% | 10% | 5% |

No class of the three-class solution is recovered in more than 80% of
replications at **any** sample size tested, including 250.

---

## 8.4 A failure on the full panel

The small-sample sweep surfaced something that had been true all along and
invisible: **`LCcov3` and `LCcov4` have singular Hessians on the complete
388-respondent panel.** All 24 and all 34 of their standard errors are `NA`,
and `05_lccov_parameters.csv` had been carrying those NA columns since the
first run.

The only signal was an R `warning()`, which R buffers and prints at the end of
a session under a "There were 13 warnings" line.

Be precise about what this invalidates. The LR *statistic* uses only
log-likelihoods and is still computable. What fails is everything else: a
singular Hessian means the model is **not locally identified**, so the optimum
is a ridge rather than a point, the coefficients on it are arbitrary, the
degrees of freedom overstate the free parameters, and BIC is penalising a
parameter count the model does not have. Of the three covariate-allocation
models, only K = 2 is a result. The other two are a diagnosis.

---

## 8.5 Two of these metrics were wrong first

Worth recording, because both produced confident nonsense that looked like a
finding.

**Relative RMSE divided by each parameter's own benchmark value.** Class
coefficients pass near zero, so the denominator did, and the table reported
relative RMSEs in the hundreds — describing the denominator, not the estimate.
Now normalised by the attribute's scale across classes.

**Class recovery used an absolute distance threshold.** In scaled units the
near-zero class sits close to the origin, so *any* badly estimated class shrunk
toward zero landed beside it. The first version of the table reported that one
class reproduced in 100% of replications and another in 10% — a clean,
memorable, entirely manufactured asymmetry. With the threshold made relative to
each class's own magnitude, the split is 75%/100% at N = 250 and the real
finding is far duller: recovery declines for every class at much the same rate.

Both are the same error in different clothes: **a normalisation that is not
scale-fair across the things being compared will invent structure.** The
symptom in each case was a number too clean to be true.

---

## 8.6 What to take from this

1. **Convergence is not evidence.** Every fit converged, at every sample size,
   including four classes fitted to five people.
2. **The first thing to fail is the thing hardest to see.** Class identity goes
   at N = 250, long before signs, standard errors or convergence.
3. **Reported uncertainty is not the uncertainty.** LC2's 95% intervals cover
   at 54% where the MNL's cover at 95%.
4. **The obvious diagnostic never fired.** Classes never collapsed. Checking
   only for the failure you expect will find nothing.
5. **Fewer observations is not the explanation.** The MNL survives the same
   sample sizes that destroy the latent class models.

For practice: a latent class model on a few hundred respondents can be
estimated, but the *classes* should not be treated as reproducible objects, and
their standard errors should not be reported at face value without a resampling
check of the kind in `08_small_sample_experiment.R`.
