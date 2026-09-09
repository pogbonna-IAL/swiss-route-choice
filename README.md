# Swiss route choice

Four different ways of letting travellers differ from one another, fitted to
the same public route choice panel, plus an experiment on how much of the
resulting structure survives a smaller sample.

The data is `apollo_swissRouteChoiceData`, which ships with the
[Apollo](http://www.apollochoicemodelling.com/) package: 388 Swiss
respondents, 9 stated-preference tasks each, two unlabelled public transport
routes per task described by travel time, cost, headway and number of
interchanges. Nothing is read from disk — there is no `data/` directory.

## Documentation

Full reference documentation lives in [`docs/`](docs/README.md): what each
model says, how it is estimated, and why the code is arranged the way it is.

| | | |
|---|---|---|
| [Models](docs/01-models.md) | [Estimation](docs/02-estimation.md) | [Post-estimation](docs/03-post-estimation.md) |
| [Validation](docs/04-validation.md) | [Experimental design](docs/05-experimental-design.md) | [Numerical methods](docs/06-numerical-methods.md) |
| [Architecture](docs/07-architecture.md) | [**What broke**](docs/08-what-broke.md) | |

## The question

A multinomial logit gives every traveller the same four coefficients. That is
almost certainly wrong. There are three standard ways to relax it, and this
project fits all three on identical data so they can be read against each
other:

| | Approach | Script |
|---|---|---|
| **Observed** | Coefficients vary with income, car access and trip purpose | `03_mnl_covariates.R` |
| **Continuous unobserved** | Each respondent draws a taste vector from a lognormal | `04_mixed_logit.R` |
| **Discrete unobserved** | Respondents belong to one of K latent classes | `06_lc_multiclass.R`, `05_lc_2class.R` |

Then the question that matters more than any of them: **how much of this is
real?** `08_small_sample_experiment.R` resamples panels from 250 respondents
down to 5 and refits everything 20 times over. `07_lc_stability.R` turns the
plausible-study range into a verdict per sample size; `11_breakdown.R` pushes
past it and asks *how* latent class estimation fails, not just whether.

The samples are **nested**: for a given seed the N = 5 panel is a subset of the
N = 30 panel is a subset of the N = 250 panel. Shrinking N removes respondents
from a fixed ordering rather than drawing an unrelated group, so "where does it
break" is a question about the same people, progressively fewer.

## ⚠ The file numbering is not the run order

Two scripts read tables that a higher-numbered script writes:

- `05_lc_2class.R` reads `06_lc_comparison.csv` → **must run after 06**
- `07_lc_stability.R` reads `08_fits.csv`, `08_parameters.csv`,
  `06_lc_comparison.csv`, `06_lc_parameters.csv` → **must run after 06 and 08**

Running the scripts in filename order fails at `05`. The real order lives in
`R/run_all.R`, which checks each script's inputs before running it and stops
with a useful message rather than a missing-file error inside someone else's
code.

The numbering is kept rather than fixed because it is baked into every output
filename (`05_lccov_*.csv`, `06_lc_*.csv`, …) and into 1,044 latent class fits
already on disk. Renaming would orphan all of it to solve a problem that one
ordered list solves.

## Quickstart

```sh
git clone https://github.com/pogbonna-IAL/swiss-route-choice.git
cd swiss-route-choice

Rscript -e 'renv::restore()'      # 1. install pinned packages   (20-60 min, once)
Rscript tests/run_tests.R         # 2. check the environment     (~30 s)
Rscript R/run_all.R               # 3. run everything            (~4 min, see below)
```

**Do step 2 before step 3.** It runs the full test suite without estimating
anything, so it tells you in half a minute whether the environment is sound
rather than failing an hour into a latent class sweep.

There is no data to download. The estimation data is
`apollo_swissRouteChoiceData`, which ships with the Apollo package that `renv`
installs in step 1.

### What a run costs

The expensive fits are committed, so a fresh clone does **not** re-estimate
them:

| Scenario | Time | What happens |
|---|---|---|
| `Rscript R/run_all.R` on a fresh clone | **~4 min** | Every model and the design load from cache; only the cheap scripts recompute |
| After deleting `outputs/models/` | ~100 min | `04` (20), `06` (45), `05` (30), `09` (5) re-estimate; `08` still resumes from its committed CSVs |
| `REFIT=1 Rscript R/run_all.R` | **~6 hours** | Everything from scratch, including all 1,044 fits in `08` |

Only the last is a true from-nothing reproduction. It is worth doing once if
you are checking the work; it is not worth doing to read the results.

### Verifying it reproduced

Running to completion is not the same as reproducing. After a run:

```sh
Rscript tests/run_tests.R
```

`tests/testthat/test-reproduction.R` compares twenty headline quantities —
log-likelihoods, BIC values, the value of travel time, the design's D-error,
the small-sample usable-fit counts — against
`tests/reference/expected_results.csv`, each at a stated tolerance. A failure
names the quantity, what was expected, what was produced, and by how much it
missed. `run_all.R` also runs this check as its final step.

If you deliberately change what the pipeline should produce, regenerate the
reference with `Rscript R/build_reference.R` and commit the diff so the change
is reviewable.

## Running it

```sh
Rscript R/run_all.R              # everything, in dependency order
Rscript R/run_all.R --list       # the plan, then exit
Rscript R/run_all.R --resume     # restart at the first step that did not finish
Rscript R/run_all.R --from 06    # 06 and everything after it
Rscript R/run_all.R 02 06 05     # just those, still in dependency order
REFIT=1 Rscript R/run_all.R      # ignore every cache, re-estimate everything
Rscript tests/run_tests.R        # test suite + reproduction check
```

### What a run tells you

`run_all.R` opens with a **pre-flight**: which steps will hit cache, which will
re-estimate, and how long it will take &mdash; taken from what each step
actually took last time, not from a guess. Then, per step:

```
[05] 05_lc_2class.R -- covariate class allocation, K = 2..4
     done in 0.16 min   !! 3 critical
     -> unusable model fit (x2), log line 134
     -> model not estimable (x1), log line 336
```

The scripts have always printed those findings into their own logs. Nothing
read them back, so the summary line said `ok` for a step that had just declared
two models inestimable. The orchestrator now scans each log for known markers
and surfaces them with a line number.

**Critical findings do not fail the run.** `LCcov3` being inestimable is a
*result*; a pipeline that refused to finish over one would be unusable. The one
thing that does fail the run, with a non-zero exit code, is **reproduction
drift** &mdash; the pipeline no longer producing the numbers this repository
claims it produces.

Every run writes `outputs/runs/<timestamp>/` containing per-step logs (each
stamped with start, finish and exit code), a machine-readable `run.json`, and a
`summary.md`. The last ten are kept.

Notifications are opt-in and off by default:

```sh
NOTIFY_DESKTOP=1 Rscript R/run_all.R                    # toast when it ends
NOTIFY_WEBHOOK=https://hooks.slack.com/... Rscript R/run_all.R
```

Each script runs in its **own R process**. Apollo keeps its model definition in
the global environment — `apollo_probabilities`, `apollo_randCoeff`,
`apollo_lcPars` — and a leftover `apollo_randCoeff` from `04` would silently
turn `05`'s latent class model into a mixed logit. Process isolation is the
only reliable guard.

**Every expensive fit is cached and resumed.** `08` appends to CSV after each
cell and skips completed `(n, seed, K)` combinations; `04`, `05`, `06` and `09`
save their fitted objects and reload them on the next run. An interrupted run
therefore costs you the model it was working on, not the whole sweep — which
matters, because these scripts are long enough that something will interrupt
them.

Those caches are **committed** (711 KB), which is what makes a fresh clone take
four minutes instead of a hundred. The rule this repo applies is: track what is
expensive to reproduce and cheap to store, ignore what is cheap to reproduce.
See `.gitignore`, which lists what is kept and why.

## Troubleshooting

**`renv` says the library is out of sync.** Run `Rscript -e 'renv::restore()'`.
`run_all.R` checks this before it starts and warns rather than letting you
discover it an hour in.

**A different R version.** `renv.lock` pins R 4.4.1 and `00_setup.R` says so on
startup if you are on something else. Results should still reproduce; the risk
is that CRAN has no binaries for your version, in which case `renv::restore()`
builds from source and you will need a toolchain (Rtools on Windows).

**The process gets killed part-way through.** These scripts are memory-hungry;
`05` and `06` were each killed repeatedly during development on an 8 GB
machine. Nothing is lost — re-run the same command and it resumes from the
last completed model. If it keeps happening, run one step at a time:
`Rscript R/run_all.R 06`, then `Rscript R/run_all.R 05`.

**A reproduction check fails.** The message names the quantity and the size of
the drift. Small drift in a log-likelihood usually means a package update moved
an optimiser; a large or exact-count difference means the specification or the
experiment changed.

## Layout

```
R/
  00_setup.R                    paths, shared constants, IO helpers, provenance
  lc_helpers.R                  latent class code generation, canonical ordering,
                                  multi-start driver, delta-method valuations
  model_helpers.R               hold-out splitting and scoring, shared by 02/03/04
  01_data_audit.R               panel structure, missingness, dominance, non-traders
  02_mnl.R                      baseline MNL + mlogit cross-check + hold-out
  03_mnl_covariates.R           observed heterogeneity
  04_mixed_logit.R              continuous unobserved heterogeneity
  05_lc_2class.R                covariate class allocation, K = 2..4
  06_lc_multiclass.R            LC1–LC5, constant-only allocation
  07_lc_stability.R             stability verdicts
  08_small_sample_experiment.R  resampling experiment, N = 250 down to 5
  09_idefix_design.R            Bayesian D-efficient design for a follow-up
  11_breakdown.R                where latent class estimation stops working
  10_report.R                   assembles every table into outputs/report.md
  build_docs.R                  bundles docs/ into one markdown file
  run_all.R                     the real run order
tests/
  run_tests.R                   entry point
  testthat/                     assertions on the pure helpers
outputs/
  models/   figures/   tables/   logs/
  report.md                     the analysis report (generated)
  documentation.md              docs/ bundled into one file (generated)
  session_info.txt              package versions for the last full run
```

### Markdown deliverables

Every page published as an Artifact has a markdown source in this repo, and the
markdown is canonical &mdash; the page is a rendering of it, not the other way
round. `10_report.R` regenerates the two derived ones on every run:

| Markdown | Contents |
|---|---|
| `outputs/report.md` | The analysis report: every table and figure, all twelve models on one BIC scale, the hold-out comparison and the small-sample breakdown |
| `outputs/documentation.md` | The eight `docs/` files bundled into one portable document, with cross-links rewritten to in-document anchors |
| `docs/08-what-broke.md` | The small-sample breakdown write-up (hand-written, canonical) |

The two generated files are gitignored because they are rebuilt from tracked
sources; `docs/*.md` is tracked.

## How the models are kept honest

Estimation code that only checks itself is not checked. Each of these is
asserted in the script that produces it, and stops the pipeline when it fails:

- **`02`** — log-likelihood at zero coefficients equals `-N log 2`; an
  independent `mlogit` implementation reproduces the estimates to six decimal
  places and the log-likelihood to 1e-10.
- **`03`** — fixing all twenty interaction terms at zero must reproduce `02`'s
  log-likelihood exactly, or the interactions are not entering the utility the
  way the script claims.
- **`06`** — the arithmetic parameter count (`4K + K-1`) must match what Apollo
  actually estimated, so a respecification cannot desynchronise the comparison
  table from the models.
- **`09`** — the efficient design must beat the median of 200 random designs of
  the same shape.
- **`04`** — the simulated panel likelihood used for scoring is checked against
  Apollo's own in-sample number: 2000 draws reproduce it to within ordinary
  simulation bias. An earlier version silently scored with a single draw (see
  below), and this is the check that catches that class of error.
- **everything** — no model is compared on in-sample fit alone. `02`, `03` and
  `04` are scored on a hold-out of 194 respondents, split on **respondent, not
  on task**: putting some of a person's nine tasks in training and the rest in
  the hold-out would let the model exploit that person's own revealed taste,
  which is exactly what these models are trying to justify.

### Canonical class ordering

Latent class labels are identified only up to permutation — two runs can find
the same optimum with classes 1 and 2 swapped. Every latent class model is put
in **canonical order** (classes sorted by travel time coefficient, most
negative first) before anything is reported, so "class 1" means the same thing
in `05`, `06` and `08`.

Reordering changes the reference class, which turns every allocation parameter
into a *contrast* against the new reference. Their standard errors are computed
as such (`V_kk + V_jj - 2 V_kj`), not permuted — permuting them would misstate
the uncertainty by however much the two parameters covary.

### A trap worth knowing about

`ifelse()` returns a result shaped like its **test** argument. In

```r
ifelse(db$choice == 1L, p, 1 - p)
```

`db$choice == 1L` is a plain vector, so if `p` is an observations-by-draws
matrix the result is a **vector holding only the first draw** — no error, no
warning, and a mixed logit hold-out likelihood that does not change when you
add draws. The scoring code multiplies by a sign instead
(`plogis(chosen_sign(db) * dv, log.p = TRUE)`) and asserts the shape it got
back.

The same function computes log-probabilities directly rather than as `log(p)`:
`exp(745)` overflows to `Inf`, so `1/(1 + exp(-dv))` is exactly zero below
about `dv = -745`, and the lognormal coefficient tails in `04` reach well past
that.

## Notes on the design script

`outputs/models/idefix_design_D.rds` existed in this project before
`09_idefix_design.R` did: it was produced in an ad-hoc session and no code
generated it. `09` regenerates it from a fixed seed and documented priors
taken from `02`. The recipe matches the original object's structure, but the
numbers differ because the original seed is unrecoverable — a reproducible
design that differs is worth more than an irreproducible one that does not.
