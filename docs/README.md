# Documentation

Reference documentation for the Swiss route choice: what each
model says, how it is estimated, and why the code is arranged the way it is.

This is the canonical documentation and it lives with the code. If you change
a specification, change the file here that describes it.

## Reading order

| | Document | What it covers |
|---|---|---|
| 1 | [Models](01-models.md) | The four specifications, their identification, and what each one claims about travellers |
| 2 | [Estimation](02-estimation.md) | Likelihoods, the panel product, simulated maximum likelihood, multi-start search, Apollo code generation |
| 3 | [Post-estimation](03-post-estimation.md) | Label switching and canonical ordering, contrast standard errors, the delta method, class shares |
| 4 | [Validation](04-validation.md) | Hold-out design, scoring, the independent cross-check, the small-sample experiment |
| 5 | [Experimental design](05-experimental-design.md) | Bayesian D-efficiency, modified Fedorov, prior propagation |
| 6 | [Numerical methods](06-numerical-methods.md) | Overflow, underflow, log-sum-exp, and four traps that fail silently |
| 7 | [Architecture](07-architecture.md) | File map, dependency graph, caching, process isolation, the orchestration layer |
| 8 | [**What broke**](08-what-broke.md) | Pushing the latent class models to N = 5: what did and did not reproduce |

If you only read one thing, read [What broke](08-what-broke.md): 1,044 fits,
1,044 converged, 440 usable. Then read
[Numerical methods](06-numerical-methods.md).
Every bug in this project that survived code review was a numerical or shape
issue that produced plausible-looking numbers rather than an error.

## Notation

Used consistently across all eight documents.

| Symbol | Meaning |
|---|---|
| $n = 1 \dots N$ | respondent, $N = 388$ |
| $t = 1 \dots T_n$ | choice task within respondent, $T_n = 9$ for all $n$ |
| $j \in \{1, 2\}$ | alternative (unlabelled route) |
| $a \in \{tt, tc, hw, ch\}$ | attribute: travel time, cost, headway, interchanges |
| $x_{ntja}$ | value of attribute $a$ for alternative $j$ in task $t$ of respondent $n$ |
| $y_{nt} \in \{1,2\}$ | the alternative chosen |
| $\beta_a$ | taste coefficient on attribute $a$ |
| $z_n$ | respondent-level covariates (income, car access, trip purpose) |
| $k = 1 \dots K$ | latent class |
| $r = 1 \dots R$ | simulation draw |
| $V_{ntj}$ | systematic utility |
| $\Delta_{nt}$ | utility difference $V_{nt1} - V_{nt2}$ |

## The one-paragraph summary

Four routes to the same question — *travellers are not identical, so how
should a model say so?* A baseline MNL gives everyone the same four
coefficients. `03` lets them vary with observed characteristics, `04` lets
them vary continuously in ways nothing observed explains, and `05`/`06` let
them vary in discrete classes. On this data the continuous specification wins
decisively (BIC 2923 against the baseline's 3364, and it is the only
alternative that also wins out of sample).

The more important result is in `07`, `08` and `11`. Across 1,044 fits at
sample sizes from 250 respondents down to 5, **every single one converged and
fewer than half were usable.** The first thing to fail is not convergence,
not the standard errors and not the signs — it is class identity, and it goes
at $N = 250$, where every conventional diagnostic still reads green. At that
same size a two-class model's nominal 95% intervals contain the full-panel
value 54% of the time, where the MNL's contain it 95% of the time, and nothing
in the ordinary output distinguishes them. See
[What broke](08-what-broke.md).
