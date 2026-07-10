# cicinference

Variance estimators for the Changes-in-Changes (CiC) design.

## Requirements

R >= 4.0, a C++ compiler (variance estimators are implemented via Rcpp).

## Installation

Not available on CRAN:

```r
devtools::install_github("Timotheebacchi/cicinference")
```

## Quick Start

```r
library(cicinference)

d <- sim_dgp(n = 5000, b1 = 0.05, b2 = 0.05, d1 = 0.05, d2 = 0.05, seed = 2026)
fit <- cic_inference(d$Y, d$X, d$Z, method = "no-split")
summary(fit)
```

## Usage

### `cic_inference()`

```r
cic_inference(
  Y,          # numeric vector: outcome
  X,          # numeric vector: treatment/endogenous variable
  Z,          # numeric vector: instrument/exogenous variable
  method,     # character vector, one or more of:
              #   "no-split" - nonparametric, full sample
              #   "split"    - sample-splitting variance estimator
              #   "kde"      - Epanechnikov KDE variance estimator (Athey-Imbens)
              #   "bse"      - bootstrap standard-error
              #   "bpc"      - bootstrap percentile
  B = 1000,       # bootstrap replications; used only for "bse"/"bpc" (B >= 200 recommended)
  epsilon_n,      # KDE bandwidth multiplier: h_{n2,u} = epsilon_n * u * (1-u).
                  # Default 1/log(n2), per the manuscript's recommendation.
  level = 0.95,   # confidence level
  panel_data = FALSE,  # TRUE for the paired-sample (Y, Z) panel estimator
  timings = FALSE      # print elapsed time per computation block
)
```

**Value**: an object of class `cic_fit` with, for each requested `method`, a point
estimate `theta_hat`, standard error `se`, and a `level`-confidence interval.
`summary()` prints these in a table; individual fields are accessible via
`fit$<method>$theta_hat`, `fit$<method>$se`, etc.

### `sim_dgp()`

Simulates data from the Monte Carlo design used in the manuscript.

```r
sim_dgp(
  n,                 # sample size
  b1 = 0, b2 = 0.05, # boundary parameters, must satisfy b1, b2 < 1
  d1 = 0, d2 = 0.05, # tail parameters, must satisfy d1 < 1 - b1, d2 < 1 - b2
  seed = NULL,
  panel_data = FALSE # if TRUE, Y and Z are paired (for the panel workflow)
)
```

`theta_true(b1, b2, d1, d2)` returns the true parameter for this DGP, useful as
a benchmark in Monte Carlo checks.

## Full Example

```r
library(cicinference)

b1 <- 0; b2 <- 0.05; d1 <- 0; d2 <- 0.05
seed <- 2026

theta0 <- theta_true(b1, b2, d1, d2)
cat("theta_true =", theta0, "\n")

## Large sample
data1 <- sim_dgp(1e6, b1, b2, d1, d2, seed)
fit <- cic_inference(data1$Y, data1$X, data1$Z,
                      method = c("no-split", "split", "kde"), timings = TRUE)
summary(fit)

## Smaller sample
data2 <- sim_dgp(1e4, b1, b2, d1, d2, seed)
fit1 <- cic_inference(data2$Y, data2$X, data2$Z,
                       method = c("no-split", "split", "kde", "bse", "bpc"),
                       timings = TRUE)
summary(fit1)

## Panel data: fit_panel accounts for the Y-Z dependence, fit_nopanel ignores
## it and is shown here only for comparison (it understates precision).
data3 <- sim_dgp(1e6, b1, b2, d1, d2, seed, panel_data = TRUE)
fit_panel   <- cic_inference(data3$Y, data3$X, data3$Z,
                              method = c("no-split", "split"), panel_data = TRUE,
                              timings = TRUE)
fit_nopanel <- cic_inference(data3$Y, data3$X, data3$Z,
                              method = c("no-split", "split"), timings = TRUE)
summary(fit_panel)
summary(fit_nopanel)
```

## Assumptions

The package targets the CiC setup of Chhor et al. (2026). Before interpreting
output, check the input data against the manuscript's assumptions, in
particular Assumption 3 for bandwidth choice. For `sim_dgp()`, parameters must
satisfy `b1, b2 < 1` and `d1 < 1 - b1`, `d2 < 1 - b2`, or `theta_true()` is not
well-defined.

## Notes

- `"no-split"`, `"split"`, and `"kde"` scale to samples of size 10^6+; each
  should complete in well under a minute. If not, something is wrong — please
  open an issue.
- The default KDE bandwidth `epsilon_n = 1/log(n2)` follows the manuscript's
  recommendation. Custom values must still satisfy Assumption 3.
- Bootstrap methods (`"bse"`, `"bpc"`) scale to millions of observations but
  take noticeably longer; prefer them on moderate sample sizes.

## References

Chhor, J., D'Haultfoeuille, X., L'Hour, J., & Mugnier, M. (2026). *Asymptotic
Properties of Empirical Quantile-Based Estimators*. [arXiv:2607.00219](https://arxiv.org/abs/2607.00219)

Athey, S. & Imbens, G. W. (2006). Identification and Inference in Nonlinear
Difference-in-Differences Models. *Econometrica*, 74(2), 431-497.
[doi:10.1111/j.1468-0262.2006.00668.x](https://doi.org/10.1111/j.1468-0262.2006.00668.x)

## Stata implementation

A native Stata implementation of `cicinference` is currently under development.

The implementation and migration tools are located in the `stata/` directory.

For development instructions, testing procedures, and release packaging, see:

- `stata/README-stata.md`


## License

MIT + file LICENSE

## Authors

Julien Chhor, Xavier D'Haultfoeuille, Jeremy L'Hour, Martin Mugnier,
Timothée Bacchi (Research Assistant)