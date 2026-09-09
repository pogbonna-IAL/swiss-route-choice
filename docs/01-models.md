# 1. Models

Four specifications, each relaxing the same assumption in a different way.
All four are fitted to identical data, so the likelihoods are directly
comparable.

---

## 1.1 The data-generating story

Each respondent $n$ sees $T_n = 9$ tasks. Each task offers two **unlabelled**
public transport routes described by four attributes: travel time (min), cost
(CHF), headway (min) and number of interchanges. The respondent picks one.

Random utility theory says the respondent attaches a utility to each route and
picks the larger:

$$U_{ntj} = V_{ntj} + \varepsilon_{ntj}$$

with $V$ the systematic part and $\varepsilon$ an unobserved shock. Assuming
$\varepsilon$ is i.i.d. type-I extreme value gives the logit probability

$$P(y_{nt} = j) = \frac{e^{V_{ntj}}}{e^{V_{nt1}} + e^{V_{nt2}}}$$

Every model in this project keeps that shape and changes only what goes into
$V$.

### Why there is no alternative-specific constant

The alternatives are **unlabelled** — "Route A" and "Route B", not "bus" and
"train". There is nothing for a constant to represent except position bias,
and `01_data_audit.R` measures the choice split at 49.7 / 50.3. An ASC would
be fitting the absence of an effect. Dropping it also makes the two
alternatives exchangeable, which is what justifies **generic** coefficients:
one minute of travel time means the same thing in either position.

### Why only differences matter

With two always-available alternatives and generic coefficients, the logit
probability depends on $V_{nt1}$ and $V_{nt2}$ only through their difference:

$$P(y_{nt} = 1) = \frac{1}{1 + e^{-(V_{nt1} - V_{nt2})}} = \Lambda(\Delta_{nt})$$

where $\Lambda$ is the logistic CDF and

$$\Delta_{nt} = \sum_a \beta_a \,(x_{nt1a} - x_{nt2a})$$

This is not a simplification imposed for convenience; it is exactly what the
model implies. It is why `attr_diff_matrix()` in `model_helpers.R` returns the
$N_{\text{obs}} \times 4$ matrix of attribute differences and why the hold-out
scoring never needs to reconstruct both utilities.

### Why the likelihood is taken over respondents, not rows

The nine tasks from one person are not nine independent observations. Whatever
makes someone unusually time-sensitive is present in all nine. Every model
here therefore forms the **panel product**

$$L_n(\theta) = \prod_{t=1}^{T_n} P(y_{nt} \mid \theta)$$

and the sample log-likelihood is $\sum_n \log L_n$.

For the MNL this changes nothing — the product factorises and the maximum is
the same either way — but it is not decorative. It is the structure the mixed
logit integrates over, and taking it in the wrong order there silently
destroys the model (see [Estimation §2.3](02-estimation.md#23-simulated-maximum-likelihood)).
It also drives the robust standard errors in `02`, which come out 1.6–1.9×
the classical ones precisely because the classical ones pretend the nine tasks
are independent.

---

## 1.2 Baseline MNL — `02_mnl.R`

$$V_{ntj} = \beta_{tt} x^{tt}_{ntj} + \beta_{tc} x^{tc}_{ntj} + \beta_{hw} x^{hw}_{ntj} + \beta_{ch} x^{ch}_{ntj}$$

**4 parameters.** One taste vector for the entire population.

The log-likelihood is globally concave, so there is a single maximum and
starting values do not matter. Zero starts are used deliberately, and this is
the last model in the project for which that is safe.

**Identification.** No scale or location normalisation is needed: with no
constant and generic coefficients there is nothing to normalise. `apollo_fixed`
is empty.

**Results.** All four coefficients negative and strongly significant. Converted
to money by ratios against the cost coefficient:

| Measure | Estimate | 95% CI |
|---|---|---|
| Value of travel time | 27.21 CHF/h | 20.67 – 33.74 |
| Value of headway | 17.05 CHF/h | 10.95 – 23.15 |
| Cost of an interchange | 8.74 CHF | 5.73 – 11.75 |
| Interchange in time terms | 19.27 min | 15.21 – 23.34 |

$\text{LL} = -1665.69$, BIC $= 3364.01$.

---

## 1.3 Observed heterogeneity — `03_mnl_covariates.R`

The cheapest possible relaxation: each coefficient becomes a linear function
of things the survey recorded.

$$\beta_{a,n} = \beta_a + \gamma_{a,\text{inc}} \widetilde{\log y_n} + \gamma_{a,\text{car}} c_n + \gamma_{a,\text{com}} m_n + \gamma_{a,\text{shop}} s_n + \gamma_{a,\text{bus}} b_n$$

for each $a \in \{tt, tc, hw, ch\}$. **24 parameters** (4 base + 20 gammas).

**Centring.** $\widetilde{\log y_n} = \log y_n - \overline{\log y}$. Without
centring, $\beta_a$ would be the coefficient for a respondent earning 1 CHF —
a point far outside the data, where the linear extrapolation is meaningless.
Centred, $\beta_a$ is the coefficient for a respondent of average income.

**The reference purpose.** `01_data_audit.R` verifies that the four purpose
dummies sum to exactly 1 for all 388 respondents. Including all four alongside
a base coefficient would be exact collinearity, so *leisure* is the omitted
reference and every purpose gamma reads as a contrast against a leisure trip.

**Nesting.** The baseline is this model with all twenty gammas at zero. That
makes the comparison a clean likelihood-ratio test on 20 df — and it is
*enforced*, not asserted: `03` estimates the restricted model by fixing the
gammas with `apollo_fixed` and checks that its log-likelihood reproduces `02`'s
to within $10^{-6}$. If the interaction terms were not entering the utility the
way the script claims, that assertion fails and the pipeline stops.

**Per-attribute tests.** Each attribute's five interactions are dropped as a
block, giving a test of "does sensitivity to *this* attribute vary across
respondents at all":

| Dropped block | LR statistic | df | p |
|---|---|---|---|
| All twenty (vs baseline) | 173.32 | 20 | $2\times10^{-26}$ |
| Travel time interactions | 42.42 | 5 | $4.8\times10^{-8}$ |
| Cost interactions | 17.49 | 5 | 0.0037 |
| Headway interactions | 31.97 | 5 | $6.0\times10^{-6}$ |
| Interchange interactions | 9.31 | 5 | 0.097 |

Time, cost and headway sensitivity vary with who the respondent is.
Interchange aversion does not, at conventional levels — a substantive result,
and the sort of thing an overall LR test alone would hide.

$\text{LL} = -1579.03$, BIC $= 3353.85$. Note that BIC barely improves on the
baseline despite a decisive LR test: twenty parameters is expensive.

---

## 1.4 Continuous unobserved heterogeneity — `04_mixed_logit.R`

Observed characteristics explain only some taste variation. The mixed logit
lets each respondent draw a whole coefficient vector from a distribution and
estimates only that distribution's parameters.

$$\beta_{a,n} = -\exp\!\big(\mu_a + (\mathbf{L}\, \mathbf{z}_n)_a\big), \qquad \mathbf{z}_n \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_4)$$

### Why negative lognormal rather than normal

A normal distribution places mass on both sides of zero, so a share of
respondents would be estimated to *prefer* journeys that are longer, dearer,
and have more interchanges. That is not heterogeneity; it is a
misspecification that the distribution guarantees by construction. `06` already
treats a wrong-signed class coefficient as a symptom rather than a finding,
and it would be incoherent to accept by assumption what is flagged as a defect
elsewhere.

The lognormal has two further benefits:

- Willingness to pay becomes a ratio of two lognormals, which is itself
  lognormal, so its quantiles are exact rather than simulated.
- $\log|\beta_a|$ is the natural parameter, so $\mu_a = \log|\beta_a^{\text{MNL}}|$
  is an excellent starting value: it puts the median respondent exactly where
  the fixed-coefficient model put everybody.

Its cost is a heavy right tail on $|\beta|$, which is the direct cause of the
overflow problem documented in [Numerical methods §6.1](06-numerical-methods.md#61-probabilities-that-are-exactly-zero).

### Two variants

| | $\mathbf{L}$ | Parameters | LL | BIC |
|---|---|---|---|---|
| **MXL-I** | diagonal | 8 | $-1443.67$ | 2952.61 |
| **MXL-C** | lower triangular | 14 | $-1404.47$ | **2923.16** |

$\mathbf{L}$ is the Cholesky factor of the covariance of the *logged absolute*
coefficients. Parameterising the factor rather than the covariance matrix keeps
the estimate positive semi-definite by construction — there is no constraint
for the optimiser to violate and no need to project back onto the valid set.

### Why the MNL comparison is not a likelihood-ratio test

MXL-I nests the MNL only in the limit $\mathbf{L} \to \mathbf{0}$, which sits
on the **boundary** of the parameter space. The usual $\chi^2$ reference
distribution does not apply there; the correct asymptotic distribution is a
mixture of chi-squareds. Rather than get that subtly wrong, `04` reports the
MNL row for reference and rests the comparison on BIC and out-of-sample fit.

MXL-C nests MXL-I in the **interior** (six off-diagonals equal to zero, not on
a boundary), so *that* comparison is a valid LR test: $\text{LR} = 78.40$ on
6 df, $p \approx 0$.

### What the model actually says

Taste correlations on the log scale (MXL-C):

| | tt | tc | hw | ch |
|---|---|---|---|---|
| **tt** | 1.00 | 0.90 | 0.72 | 0.83 |
| **tc** | 0.90 | 1.00 | 0.65 | 0.74 |
| **hw** | 0.72 | 0.65 | 1.00 | 0.76 |
| **ch** | 0.83 | 0.74 | 0.76 | 1.00 |

Uniformly high and positive: someone sensitive to travel time is also
sensitive to cost, headway and interchanges. In a random utility model this
pattern is the signature of **scale heterogeneity** — respondents differing in
how deterministically they choose, rather than in what they want. Six extra
parameters is a large price for that, and it is worth reading the MXL-C result
with that interpretation in mind.

The headline is not a number but a spread:

| | Median | 10th pct | 90th pct |
|---|---|---|---|
| Value of travel time | 24.86 CHF/h | 7.28 | 84.82 |

The MNL's single 27.21 CHF/h sits near the median of a distribution spanning
more than a factor of ten. That dispersion is the finding.

---

## 1.5 Discrete unobserved heterogeneity — `06_lc_multiclass.R`, `05_lc_2class.R`

Instead of a continuous distribution, a small number of discrete types. Each
class $k$ has its own full taste vector, and each respondent belongs to one —
we just never observe which.

$$L_n = \sum_{k=1}^{K} \pi_{k} \prod_{t=1}^{T_n} P(y_{nt} \mid \beta_k)$$

The class probability is a softmax over allocation utilities:

$$\pi_{nk} = \frac{e^{W_{nk}}}{\sum_{l} e^{W_{nl}}}$$

### Constant-only allocation (`06`)

$W_{nk} = \delta_k$, identical for every respondent. Only differences between
classes are identified, so $\delta_1 \equiv 0$ (imposed with `apollo_fixed`).

**Parameters:** $4K + (K-1)$.

| Model | Par | LL | BIC | Class shares | Starts at best |
|---|---|---|---|---|---|
| LC1 | 4 | $-1665.69$ | 3364.01 | 1.000 | 1/1 |
| LC2 | 9 | $-1552.53$ | 3178.49 | 0.305 / 0.695 | 2/50 |
| LC3 | 14 | $-1490.22$ | 3094.65 | 0.359 / 0.183 / 0.458 | 8/50 |
| LC4 | 19 | $-1449.52$ | **3054.06** | 0.300 / 0.176 / 0.230 / 0.294 | 16/50 |
| LC5 | 24 | $-1431.01$ | 3057.81 | 0.246 / 0.180 / 0.237 / 0.184 / 0.153 | 1/50 |

LC1 reproduces the MNL exactly, as it must — a useful internal consistency
check that the latent class machinery reduces correctly.

The **`at_best` column is the most informative one in the table.** It counts
how many of the fifty random starts rediscovered the reported optimum. LC5's
1/50 means the reported log-likelihood was found once and could easily have
been missed; LC2's 2/50 is barely better. A latent class result without this
column is a result you cannot assess.

### Covariate allocation (`05`)

$$W_{n1} = 0, \qquad W_{nk} = \delta_k + \sum_c \gamma_{ck}\, z_{nc} \quad (k \ge 2)$$

Class 1's **entire allocation index** is normalised to zero, not just its
constant. **Parameters:** $4K + (K-1)(1 + C)$ with $C = 5$ covariates.

Each model nests its constant-only counterpart (set every $\gamma$ to zero),
giving a clean LR test on $(K-1) \times 5$ df:

| $K$ | Par | LL | BIC | LR | p | BIC better? | Usable? |
|---|---|---|---|---|---|---|---|
| 2 | 14 | $-1544.60$ | 3203.42 | 15.86 | 0.0072 | no | yes |
| 3 | 24 | $-1477.12$ | 3150.05 | 26.18 | 0.0035 | no | **NO** |
| 4 | 34 | $-1427.72$ | 3132.82 | 43.61 | 0.00013 | no | **NO** |

**Read that table carefully — twice.**

First: every LR test is significant and *every* BIC gets worse. The covariates
shift class membership detectably but do not pay for the parameters they cost.
Reporting only the p-values would tell a misleadingly positive story.

Second, and more seriously: **`LCcov3` and `LCcov4` have singular Hessians on
the full 388-respondent panel.** All 24 and all 34 of their standard errors are
`NA`. This was true from the first run and went unnoticed for a long time,
because the only signal was an R `warning()` that gets buffered to the end of
the session and printed under a "There were 13 warnings" line.

Be precise about what that invalidates. The LR *statistic* uses only
log-likelihoods, so it is still computable, and the two log-likelihoods are
real. What fails is everything else:

- a singular Hessian means the model is **not locally identified** — there is a
  direction in parameter space along which the likelihood does not change;
- so the reported optimum is a ridge rather than a point, and the coefficient
  values on it are arbitrary;
- the degrees of freedom overstate the number of parameters the model actually
  has, and BIC is penalising a count the model does not carry;
- no coefficient, odds ratio or membership profile from those two fits can be
  interpreted at all.

Only the $K = 2$ row is a result. The other two rows are a diagnosis.

This is the reason `lc_fit_verdict()` exists and why every table now carries a
`usable` column: a model that converges, reports a log-likelihood, and produces
a plausible-looking parameter table can still be telling you nothing, and
nothing about the ordinary output distinguishes it from one that is.

`05` also computes a `lr_trustworthy` flag. A likelihood-ratio test compares
two **global** maxima; if either model's optimum was found by only a handful of
starts, the test is being run against a number the search may simply have
failed to beat, and the p-value is optimistic. Rows where either side has
fewer than three starts at the optimum are marked, and here the $K=2$
comparison fails that check (2/50 on the constant-only side).

---

## 1.6 Where this leaves us

All twelve models on one BIC scale:

| Rank | Model | Family | Par | LL | BIC |
|---|---|---|---|---|---|
| 1 | **MXL-C** | continuous | 14 | $-1404.47$ | **2923.16** |
| 2 | MXL-I | continuous | 8 | $-1443.67$ | 2952.61 |
| 3 | LC4 | discrete | 19 | $-1449.52$ | 3054.06 |
| 4 | LC5 | discrete | 24 | $-1431.01$ | 3057.81 |
| 5 | LC3 | discrete | 14 | $-1490.22$ | 3094.65 |
| 6 | LCcov4 | discrete + observed | 34 | $-1427.72$ | 3132.82 |
| … | … | | | | |
| 11 | MNL | fixed | 4 | $-1665.69$ | 3364.01 |
| 12 | LC1 | discrete | 4 | $-1665.69$ | 3364.01 |

Continuous heterogeneity wins, and it wins on the hold-out too
(see [Validation](04-validation.md)), which BIC alone cannot establish.

But the ranking is not the most important output of this project. That is in
[the small-sample experiment](04-validation.md#45-the-small-sample-experiment):
at $N = 100$, four of twenty resampled LC2 fits return a wrong-signed
coefficient, and by $N = 30$ BIC picks a single class 45% of the time. Almost
none of this structure is recoverable from a study a quarter the size.
