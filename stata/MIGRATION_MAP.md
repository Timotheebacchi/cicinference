# Migration Map: R package `quantcdf.inference` to Stata/Mata

Status: Phase 1 inventory only. No Stata implementation has been started.

Repository inspected:

- Package metadata: `DESCRIPTION`, `NAMESPACE`
- R sources: `R/cic.R`, `R/data-generating.R`, `R/utils.R`, `R/print_summary.R`
- C++ sources: `src/Density_cic_boostrap_functions.cpp`
- Documentation and examples: `README.md`, `vignettes/introduction-cic.Rmd`, `man/*.Rd`, `inst/CITATION`
- Tests: `tests/testthat/test-cic.R`

## Exported R functions

Primary exported functions from `NAMESPACE`:

| R function | Source | Role | Proposed Stata/Mata target |
|---|---|---|---|
| `fit()` | `R/cic.R` | Main user-facing quantile estimator and CI constructor | `stata/ado/cicinference.ado`, `stata/mata/cicinference_main.mata` |
| `qY_dgp()` | `R/data-generating.R` | Quantile function for simulation DGP | Optional: `stata/mata/cicinference_dgp.mata` or test-only `.do`/Mata helper |
| `theta_true()` | `R/data-generating.R` | Closed-form true DGP parameter | Optional: `stata/mata/cicinference_dgp.mata` or reference-output scripts |
| `sim_dgp()` | `R/data-generating.R` | Simulation helper for examples/tests | `stata/tests/make_reference_outputs.R`; optional Stata example helper |

Registered S3 methods:

| R method | Source | Role | Proposed Stata/Mata target |
|---|---|---|---|
| `summary.quantcdf_fit()` | `R/print_summary.R` | Printed estimator summary | `stata/ado/cicinference.ado`, `stata/ado/cicinference.sthlp` |
| `coef.quantcdf_fit()` | `R/print_summary.R` | Extract point estimate | Stored result, likely `e(theta)` or `r(theta)` |
| `confint.quantcdf_fit()` | `R/print_summary.R` | Extract CI matrix | Stored results matrix, likely `e(ci)` or `r(ci)` |

Note: `NAMESPACE` exports class methods for `quantcdf_fit`; the R package stores the split/no-split bandwidth multiplier as `epsilon_n` and the Athey-Imbens bandwidth as `h_ai`.

## Internal R functions

| Internal function | Source | Role | Proposed Stata/Mata target |
|---|---|---|---|
| `.prepare_left_quantile()` | `R/utils.R` | Builds left-continuous empirical quantile evaluator | `stata/mata/cicinference_quantiles.mata` |
| `.make_density_estimator()` | `R/utils.R` | Cached wrapper around interval counts and density conversion | `stata/mata/cicinference_density.mata` |
| `.compute_eta_from_f()` | `R/utils.R` | Aggregates KDE variance term from density estimates | `stata/mata/cicinference_variance.mata` |
| `.fast_eta()` | `R/utils.R` | Fast double-integral/eta variance term | `stata/mata/cicinference_variance.mata` |
| `.panel_no_split_estimator()` | `R/utils.R` | Panel-data paired-sample no-split estimator | `stata/mata/cicinference_main.mata`, `stata/mata/cicinference_variance.mata` |
| `.panel_split_estimator()` | `R/utils.R` | Panel-data paired-sample split estimator | `stata/mata/cicinference_main.mata`, `stata/mata/cicinference_variance.mata` |
| `.quantcdf_table()` | `R/print_summary.R` | Formats standard errors, test statistics, p-values | `stata/ado/cicinference.ado` |
| `log_timing()` | nested in `fit()` | Optional timing messages | `stata/ado/cicinference.ado` |
| `.make_density_estimator()$estimate()` | closure | Cached density evaluation for one `(epsilon, pointwise)` key | `stata/mata/cicinference_density.mata` |
| `.make_density_estimator()$reset()` | closure | Clears R cache | Usually unnecessary in Stata; use local Mata structs/vectors |

## Rcpp/C++ components

All compiled functions are in `src/Density_cic_boostrap_functions.cpp` and are called from R through Rcpp.

| C++ function | R caller | Role | Complexity | Proposed Mata target |
|---|---|---|---|---|
| `f_y_hat_epnechikov()` | `fit(method = "ai")` | Epanechnikov density estimate using sorted centered Sample 1, prefix sums, and AI bandwidth | `O(n log n + m log n)` time, `O(n)` memory | `stata/mata/cicinference_density.mata` |
| `rect_counts_rcpp()` | `.make_density_estimator()` | Counts sorted observations in `[x_eval - h, x_eval + h]` via binary search | `O(m log n)` time | `stata/mata/cicinference_density.mata` |
| `counts_to_density()` | `.make_density_estimator()` | Converts rectangle counts to `count / (2*n*h)` | `O(m)` time | `stata/mata/cicinference_density.mata` |
| `boot_core()` | `fit(method = "bse"|"bpc")` | Nonparametric bootstrap for theta using multinomial counts for sorted Sample 1 and Sample 3, and resampled Sample 2 | `O(B*(n1+n3+n2*(log n1+log n3)))` time | `stata/mata/cicinference_bootstrap.mata` |

Important spelling note: the exported C++ function name is `f_y_hat_epnechikov`, matching the current R call, even though the conventional spelling is Epanechnikov.

## Function dependency graph

```text
fit()
  |-- stats::qnorm()
  |-- input validation, warning, timing setup
  |-- default epsilon_n: function(n) 1/log(n)
  |-- cross-sectional baseline, if panel_data == FALSE
  |     |-- stats::ecdf(Z)
  |     |-- FZ(X) -> Uhat
  |     |-- .prepare_left_quantile(Y)
  |     |-- qcdf_transform = FY_left_inverse(Uhat)
  |     |-- theta_hat = mean(qcdf_transform)
  |     |-- eps_hat = mean((theta_hat - qcdf_transform)^2)
  |     `-- h = epsilon_n * Uhat * (1 - Uhat)
  |
  |-- method "no-split", cross-sectional
  |     |-- sort(Y), diff()
  |     |-- FYhat = (1:(n1-1))/n1
  |     |-- .make_density_estimator(sort(Uhat), FYhat)
  |     |     |-- rect_counts_rcpp()
  |     |     `-- counts_to_density()
  |     |-- .fast_eta()
  |     `-- sigma_sq = lbda1_3 * eta + lbda2 * eps_hat
  |
  |-- method "split", cross-sectional
  |     |-- first/last half sample splits for Y, X, Z
  |     |-- stats::ecdf(Z_half_1), stats::ecdf(Z_half_2)
  |     |-- Uhat1 = FZ1(X_first_half)
  |     |-- Uhat2 = FZ2(X_last_half)
  |     |-- .make_density_estimator(sort(Uhat1), FYhat_split)
  |     |-- .make_density_estimator(sort(Uhat2), FYhat_split)
  |     |-- .fast_eta()
  |     `-- sigma_sq_split = lbda1_3 * eta_hat_split + lbda2 * eps_hat
  |
  |-- method "kde", cross-sectional
  |     |-- order(Uhat)
  |     |-- findInterval(grid - 1e-12, U_sort) + 1
  |     |-- f_y_hat_epnechikov(Y, qcdf_transform, h)
  |     |-- .compute_eta_from_f()
  |     `-- sigma_sq_kde = 2 * eta_ai + eps_hat
  |
  |-- method "bse" or "bpc", cross-sectional only
  |     |-- boot_core(sort(Y), X, sort(Z), B)
  |     |-- stats::sd() for "bse"
  |     `-- stats::quantile() for "bpc"
  |
  |-- method "no-split", panel_data == TRUE
  |     `-- .panel_no_split_estimator()
  |
  `-- method "split", panel_data == TRUE
        `-- .panel_split_estimator()

.panel_no_split_estimator()
  |-- split paired (Y,Z) into four blocks, X into two blocks
  |-- stats::ecdf() on Z blocks
  |-- .prepare_left_quantile() on Y blocks
  |-- theta_hat = mean(theta_1, theta_2)
  |-- .make_density_estimator() for U_1, U_2
  |-- .fast_eta()
  `-- panel residual variance term and standard error

.panel_split_estimator()
  |-- same core dependencies as .panel_no_split_estimator()
  `-- different asymptotic scaling with N = min(n1, n2)

summary.quantcdf_fit()
  |-- .quantcdf_table()
  |     |-- stats::qnorm()
  |     `-- stats::pnorm()
  `-- print()

confint.quantcdf_fit()
  `-- stored object$ci

coef.quantcdf_fit()
  `-- stored object$theta_hat

sim_dgp()
  |-- optionally set.seed()
  |-- stats RNG: runif(), rnorm(), rbeta()
  |-- stats::qnorm()
  `-- qY_dgp()

theta_true()
  `-- beta()
```

## Major function ledger

### `fit()`

- Inputs: numeric vectors `Y`, `X`, `Z`; `method`; bootstrap count `B`; `epsilon_n` scalar/function/NULL; CI `level`; `panel_data`; `timings`.
- Outputs: S3 list with `theta_hat`, `ci`, `level`, sample sizes `n`, requested `method`, `panel_data`, and `epsilon_n`.
- Side effects: warnings for coercion, NA/non-finite inputs, tied `Y`, small `B`; timing messages; bootstrap consumes R RNG state.
- Mathematical role: estimates `theta = E[F_Y^{-1}(F_Z(X))]` and constructs confidence intervals using no-split, split, KDE, bootstrap-SE, and bootstrap-percentile methods. In panel mode, uses paired `(Y,Z)` splitting logic and supports only no-split and split.
- Complexity: dominated by sorting, empirical CDF evaluation, density counts, and bootstrap. Cross-sectional no-split/split are roughly `O(n log n)`; KDE is `O(n log n)`; bootstrap is `O(B*n*log n)` with separate `n1`, `n2`, `n3` terms.
- Target: ado wrapper for syntax/options/display plus Mata implementation in `cicinference_main.mata`, with calls to quantile, density, variance, and bootstrap modules.

### `qY_dgp()`

- Inputs: probabilities `t`, tail parameters `d1`, `d2`.
- Outputs: numeric vector `-t^(-d1) + (1-t)^(-d2)` with special cases `d1 == 0` and `d2 == 0`.
- Side effects: none.
- Mathematical role: DGP quantile function for simulation.
- Complexity: `O(n)` for vector length `n`.
- Target: reference/testing helper, not required for the main ado command.

### `theta_true()`

- Inputs: boundary parameters `b1`, `b2`; tail parameters `d1`, `d2`.
- Outputs: scalar true DGP theta using beta functions.
- Side effects: none.
- Mathematical role: analytic benchmark for simulations.
- Complexity: `O(1)`.
- Target: reference/testing helper.

### `sim_dgp()`

- Inputs: `n`, `b1`, `b2`, `d1`, `d2`, optional `seed`, `panel_data`.
- Outputs: list with `Y`, `X`, `Z`.
- Side effects: calls `set.seed(seed)` when supplied; consumes RNG state.
- Mathematical role: Monte Carlo DGP: `W ~ U(0,1)`, `Y = qY_dgp(W)`, `V ~ Beta(1-b1, 1-b2)`, `X = qnorm(V)`, and `Z = rnorm(n)` unless `panel_data`, where `Z = qnorm(W)`.
- Complexity: `O(n)`.
- Target: reference-output generation first; optional Stata example support later.

### `.prepare_left_quantile()`

- Inputs: numeric vector `x`.
- Outputs: closure evaluating `xs[pmax(ceiling(p*n), 1)]` after radix-sorting `x`.
- Side effects: none after closure creation.
- Mathematical role: left-continuous empirical quantile with convention `F_Y^{-1}(0) = Y_(1)`.
- Complexity: setup `O(n log n)`, evaluation `O(m)` for `m` probabilities.
- Target: `cicinference_quantiles.mata`.

### `.make_density_estimator()`

- Inputs: sorted rank vector `U_sorted`; empirical grid `FYhat_split`.
- Outputs: closure with `estimate(eps, pointwise)` returning rectangle-count density estimates.
- Side effects: mutates an internal R cache by `(eps, pointwise)` key.
- Mathematical role: estimates density of `U = F_Z(X)` on the `FYhat` grid using variable bandwidth `h = eps * (F*(1-F))^pointwise`.
- Complexity: each uncached estimate is `O(m log n)` through binary-search counts plus `O(m)` conversion.
- Target: `cicinference_density.mata`; cache can be replaced by local reuse of Mata vectors.

### `.fast_eta()`

- Inputs: `Ydiff1`, `fUhat1`, `Ydiff2`, `fUhat2`, `FYhat`.
- Outputs: scalar eta term.
- Side effects: none.
- Mathematical role: computes a double-integral variance component without forming an `n x n` matrix. Uses `u = Ydiff1*fUhat1`, `v = Ydiff2*fUhat2`, reverse cumulative sums, and grid increments.
- Complexity: `O(m)` time and memory for `m = length(FYhat)`.
- Target: `cicinference_variance.mata`.

### `.compute_eta_from_f()`

- Inputs: density values `f_vals`; unsorted `Uhat`; sort index `idx_sort`; interval positions `k`; logical mask `ok`; sample size `n`.
- Outputs: scalar eta term for KDE method.
- Side effects: none.
- Mathematical role: computes `mean((T1 - C2)^2)` where `C2 = mean(Uhat/f_vals)` and `T1` is based on reverse cumulative sums of sorted inverse densities.
- Complexity: `O(n)` if sort/order and interval positions are precomputed.
- Target: `cicinference_variance.mata`.

### `.panel_no_split_estimator()`

- Inputs: paired `Y`, `Z`, independent `X`, scalar `epsilon_n`.
- Outputs: list with `theta_hat`, `sigma_sq`, `se`, `q`, `r`.
- Side effects: errors if `length(Y) != length(Z)` or samples are too small.
- Mathematical role: paired panel estimator using four blocks of `(Y,Z)` and two blocks of `X`; cross-combines empirical CDFs and quantiles to account for dependence between `Y` and `Z`.
- Complexity: roughly `O(n log n)` from sorting/CDF/density operations.
- Target: `cicinference_main.mata` plus shared quantile/density/variance helpers.

### `.panel_split_estimator()`

- Inputs: paired `Y`, `Z`, independent `X`, scalar `epsilon_n`.
- Outputs: list with `theta_hat`, `sigma_sq`, `se`, `q`, `r`.
- Side effects: same validation errors as panel no-split.
- Mathematical role: panel split estimator using the same block construction as panel no-split but asymptotic scaling with `N = min(n1,n2)`.
- Complexity: roughly `O(n log n)`.
- Target: `cicinference_main.mata` plus shared quantile/density/variance helpers.

### `.quantcdf_table()`

- Inputs: `quantcdf_fit` object.
- Outputs: data frame with method, standard error, t statistic, p-value, CI bounds, and length.
- Side effects: none.
- Mathematical role: display-only post-processing; derives standard errors from CI length except for percentile bootstrap.
- Complexity: `O(k)` for number of methods.
- Target: ado display and stored-results formatting.

### `summary.quantcdf_fit()`

- Inputs: fitted object, `digits`, `...`.
- Outputs: invisibly returns object.
- Side effects: prints summary to console.
- Mathematical role: display only.
- Complexity: `O(k)`.
- Target: `cicinference.ado` display block and `.sthlp` examples.

### `coef.quantcdf_fit()`

- Inputs: fitted object.
- Outputs: scalar `theta_hat`.
- Side effects: none.
- Mathematical role: extraction only.
- Complexity: `O(1)`.
- Target: stored result scalar.

### `confint.quantcdf_fit()`

- Inputs: fitted object; optional `level`.
- Outputs: matrix of lower/upper CI bounds with method row names.
- Side effects: warning if requested level differs from stored object level.
- Mathematical role: extraction only.
- Complexity: `O(k)`.
- Target: stored result matrix.

## Algorithmic components

1. Input validation and sanitization
   - Validate numeric vectors, minimum lengths, methods, CI level, bootstrap `B`, `epsilon_n`, and panel method restrictions.
   - Current R code warns about NA/non-finite values but does not explicitly remove them before estimation.

2. Empirical CDF and rank transformation
   - Compute `Uhat = F_Z(X)` using R's `stats::ecdf`.
   - Stata/Mata must reproduce R `ecdf` behavior at ties and boundaries.

3. Left-continuous empirical quantile
   - Compute `F_Y^{-1}(u) = Y_(ceil(u*n))`, with `u = 0` mapped to the first order statistic.
   - This convention is central to theta and bootstrap equivalence.

4. Point estimate
   - Cross-sectional: `theta_hat = mean(F_Y^{-1}(F_Z(X)))`.
   - Panel: cross-fit two block-specific estimates and average them.

5. Residual/moment term
   - Cross-sectional `eps_hat = mean((theta_hat - qcdf_transform)^2)`.
   - Panel `eps_panel = mean((eps_1^2 + eps_2^2)/2)`.

6. Rectangle-count density estimator
   - Count observations of sorted `Uhat` in variable intervals around `FYhat`.
   - Convert to density via `count/(2*n*h)`.
   - This is used by no-split, split, and panel variance estimators.

7. Fast eta aggregation
   - Uses sorted `Y` gaps, density estimates, reverse cumulative sums, and grid increments.
   - Avoids the old matrix construction and must be implemented in Mata.

8. KDE localized density
   - Epanechnikov estimate over `Y` evaluated at `qcdf_transform`, with bandwidth `h = epsilon_n*Uhat*(1-Uhat)`.
   - C++ centers `Y` before moment calculations to reduce cancellation.

9. KDE eta aggregation
   - Sort `Uhat`, calculate interval indices using `findInterval(grid - 1e-12, U_sort) + 1`, then compute inverse-density cumulative terms.

10. Bootstrap
    - Resamples `Y` and `Z` via multinomial count arrays over sorted values.
    - Resamples `X` observation-by-observation.
    - Computes bootstrap theta values, then either percentile intervals (`bpc`) or SE intervals (`bse`).

11. Panel splitting
    - Requires `length(Y) == length(Z)`.
    - Uses `q = floor(n1/4)`, `r = floor(n2/2)`.
    - Blocks: `z_1`, `z_2`, `y_1`, `y_2`, `x_1`, `x_2`.
    - Supports only `no-split` and `split`.

12. Display and stored results
    - R returns an S3 object; Stata should decide between `e()` and `r()`.
    - Since the command reports an estimate, standard errors, confidence intervals, and p-values, `e()` is likely natural, but this should be finalized before implementation.

## Interface, display, and documentation components

| Component | Current R behavior | Stata target |
|---|---|---|
| Main command | `fit(sample1, sample2, sample3, method=..., ...)` | `cicinference y x z, method(...) ...` |
| Multiple methods | R accepts a character vector and returns one CI row per method | Stata option may accept one or more methods, or repeatable/multi-token `method()` syntax |
| Bootstrap options | `B`, with minimum coerced to 200 after warning | `bootstrap(#)` or `breps(#)`; preserve minimum/warning behavior |
| Bandwidth | `epsilon_n = NULL` means `1/log(n2)`; function values allowed in R | Stata should support scalar `epsilon(#)` and default `1/log(n2)`; function-valued option likely not portable |
| Panel option | `panel_data = TRUE` only for no-split/split | `panel` option or `paneldata` option |
| Timing | `timings = TRUE` emits messages | `timings` option in ado |
| Summary | S3 print table | Native Stata results table |
| Extractors | `coef()`, `confint()` | Stored scalars/matrices |
| Help/docs | roxygen Rd, README, vignette | `stata/ado/cicinference.sthlp`, `stata/README-stata.md` |

## Proposed Stata/Mata target files

Initial target layout from `AGENTS.md`:

```text
stata/
  ado/
    cicinference.ado
    cicinference.sthlp

  mata/
    cicinference_utils.mata
    cicinference_quantiles.mata
    cicinference_density.mata
    cicinference_variance.mata
    cicinference_bootstrap.mata
    cicinference_main.mata

  tests/
    make_reference_outputs.R
    test_nosplit.do
    test_split.do
    test_kde.do
    test_bootstrap.do
    benchmark.do

  examples/
    basic_usage.do
    replication_example.do

  README-stata.md
```

Suggested ownership:

| File | Responsibilities |
|---|---|
| `stata/ado/cicinference.ado` | Syntax parsing, `if`/`in`, missing handling, option validation, Mata calls, display, stored results |
| `stata/ado/cicinference.sthlp` | Stata help, examples, stored results, equivalence notes |
| `stata/mata/cicinference_utils.mata` | Common checks, index/block helpers, timing hooks, safe numeric helpers |
| `stata/mata/cicinference_quantiles.mata` | Sorting, empirical CDF evaluation, left-continuous quantiles, tie handling |
| `stata/mata/cicinference_density.mata` | Rectangle counts, density conversion, Epanechnikov KDE |
| `stata/mata/cicinference_variance.mata` | `eps_hat`, `.fast_eta`, KDE eta, asymptotic variance scalings |
| `stata/mata/cicinference_bootstrap.mata` | Bootstrap sampling, sorted-count CDF logic, theta bootstrap values |
| `stata/mata/cicinference_main.mata` | Cross-sectional and panel estimator orchestration |
| `stata/tests/make_reference_outputs.R` | Future deterministic R reference outputs |
| `stata/tests/*.do` | Future Stata numerical equivalence tests |

## Difficulty ranking for migration

| Rank | Component | Difficulty | Main reason |
|---:|---|---|---|
| 1 | Display/extractors (`summary`, `coef`, `confint`) | Low | Mostly formatting and stored-result design |
| 2 | DGP helpers (`qY_dgp`, `theta_true`, `sim_dgp`) | Low to Medium | Formulae are simple; exact RNG equivalence is not required unless comparing simulated draws |
| 3 | Input parsing and validation | Medium | Need natural Stata syntax while preserving R defaults and warnings |
| 4 | Empirical CDF and left quantile primitives | Medium | Tie, boundary, missing, and `ceil(u*n)` behavior must match R |
| 5 | Cross-sectional point estimate | Medium | Straightforward once CDF/quantile primitives are correct |
| 6 | `.fast_eta` variance term | Medium | Algorithm is compact but sensitive to indexing/grid conventions |
| 7 | No-split and split variance estimators | Medium to High | Must reproduce sample splitting, scaling, and density grid exactly |
| 8 | Panel no-split/split estimators | High | Block indexing and two different variance scalings are easy to misalign |
| 9 | Rectangle-count density estimator | High | Requires efficient Mata binary search/two-pointer logic and exact interval inclusivity |
| 10 | KDE method | High | Must reproduce centered Epanechnikov prefix-sum implementation and `findInterval(... - 1e-12)` behavior |
| 11 | Bootstrap (`boot_core`) | Very High | R RNG, multinomial count mechanics, sorted-count quantiles, percentile type, and Stata RNG differences make exact replication hard |
| 12 | End-to-end numerical equivalence tests | Very High | Existing R tests are partly stale; reference outputs must be generated from the current implementation before porting |

## First implementation milestone

Do not implement this in the current step. The first future implementation milestone should be:

1. Create `stata/tests/make_reference_outputs.R`.
2. Generate deterministic reference outputs from the current R package for:
   - small, medium, and moderately large cross-sectional datasets;
   - no-split, split, KDE, BSE, and BPC methods;
   - panel no-split and panel split;
   - tied values, boundary ranks, odd sample sizes, and missing-value probes.
3. Implement only the foundational Mata utilities needed to reproduce the cross-sectional point estimate and no-split variance:
   - sorting;
   - R-compatible empirical CDF evaluation;
   - left-continuous quantile evaluation;
   - rectangle-count density;
   - `.fast_eta`.
4. Validate a first `method(no-split)` path against the R reference outputs before adding split, KDE, panel, or bootstrap.

Success criterion for milestone 1: Stata/Mata reproduces R `theta_hat`, no-split CI bounds, standard error, and intermediate diagnostics on fixed reference datasets within documented numerical tolerances.

## Migration risks

1. R/Stata empirical CDF differences
   - R's `ecdf` and Stata/Mata custom code must agree on ties and endpoints. This affects every estimator.

2. Left-quantile convention
   - The R package uses `ceil(p*n)` and maps `p=0` to the first order statistic. Any shift to Stata's built-in percentile definitions would change results.

3. Missing values
   - Current R code warns but does not clearly drop missing/non-finite values. Stata's default missing ordering is different from R's NA behavior, so missing handling must be explicitly specified and tested.

4. Bootstrap RNG equivalence
   - Stata will not naturally reproduce R's `unif_rand()` stream. Validation should compare saved R reference outputs or use fixed bootstrap index draws if exact bootstrap replication is required.

5. Density interval inclusivity
   - C++ uses `lower_bound(x-h)` and `upper_bound(x+h)`, i.e. closed intervals at both ends for values equal to the upper bound. Mata must match this exactly.

6. KDE numerical details
   - The C++ implementation centers `Y`, uses prefix sums, and applies `grid - 1e-12` before `findInterval`. These small numerical choices affect edge cases.

7. Panel variance scaling
   - `.panel_no_split_estimator()` and `.panel_split_estimator()` are structurally similar but use different scaling formulas. This should be tested with intermediate outputs, not only final CIs.

8. Documentation/test drift
   - Tests and docs still refer to fields/options/classes such as `h` and `"cic"` that do not match the current code. Reference outputs should be based on the installed/current implementation, and these discrepancies should be resolved or documented before using the tests as a contract.
