# quantcdf.inference

Inference for empirical quantile-based estimators of
`E[F_Y^{-1}(F_Z(X))]`, including Changes-in-Changes applications.

## Requirements

R >= 4.0 and a C++ compiler. Computational helpers are implemented with Rcpp.

## Installation

The package is not on CRAN yet. Until it is, install it with `devtools`:

```r
install.packages("devtools")
devtools::install_github("Timotheebacchi/cicinference")
```

From a local checkout, use `devtools::install_local(".")`.

## Quick Start

```r
library(quantcdf.inference)

d <- sim_dgp(n = 5000, b1 = 0.05, b2 = 0.05, d1 = 0.05, d2 = 0.05, seed = 2026)
out <- fit(d$Y, d$X, d$Z, method = "no-split")
summary(out)
```

## Usage

### `fit()`

```r
fit(
  sample1,    # Sample 1, corresponding to Y
  sample2,    # Sample 2, corresponding to X
  sample3,    # Sample 3, corresponding to Z
  method,     # one or more of "no-split", "split", "ai", "bse", "bpc"
  B = 1000,
  epsilon_n = NULL,
  level = 0.95,
  panel_data = FALSE,
  timings = FALSE
)
```

`fit()` estimates `theta = E[F_Y^{-1}(F_Z(X))]` with empirical distribution
and left-continuous empirical quantile functions. The convention at probability
zero is `F_Y^{-1}(0) = Y_(1)`.

The split and no-split density estimators use
`h_n2(u) = epsilon_n2 * u * (1 - u)` with default
`epsilon_n2 = 1 / log(n2)`. The Athey-Imbens method uses an Epanechnikov
density estimate of Sample 1 with `h_AI = 1.06 * n1^(-1/5) / sd(Sample 1)`.

### `cic_fit()`

```r
cic_fit(
  Y11,
  Y10,
  Y01,
  Y00,
  method,
  B = 1000,
  epsilon_n = NULL,
  level = 0.95,
  timings = FALSE
)
```

`cic_fit()` estimates `mean(Y11) - fit(Y01, Y10, Y00)$theta_hat`. For analytic
methods, it adds the independent repeated-cross-section contribution of `Y11`
on the same `sqrt(N)` variance scale. For bootstrap methods, each replication
resamples all four samples independently and recomputes the full contrast.

### `sim_dgp()`

```r
sim_dgp(
  n,
  b1 = 0, b2 = 0.05,
  d1 = 0, d2 = 0.05,
  seed = NULL,
  panel_data = FALSE
)
```

`theta_true(b1, b2, d1, d2)` returns the true parameter for this DGP.

## Example

```r
library(quantcdf.inference)

b1 <- 0
b2 <- 0.05
d1 <- 0
d2 <- 0.05

d <- sim_dgp(10000, b1, b2, d1, d2, seed = 2026)
theta0 <- theta_true(b1, b2, d1, d2)

out <- fit(d$Y, d$X, d$Z, method = c("no-split", "split", "ai", "bse", "bpc"))
summary(out)

theta0
coef(out)
confint(out)
```

```r
Y11 <- d$Y + 0.5
cic <- cic_fit(Y11, d$X, d$Y, d$Z, method = c("no-split", "ai", "bse"))
summary(cic)
```

## References

Chhor, J., D'Haultfoeuille, X., L'Hour, J., & Mugnier, M. (2026).
*Asymptotic Properties of Empirical Quantile-Based Estimators*.
[arXiv:2607.00219](https://arxiv.org/abs/2607.00219)

Athey, S. & Imbens, G. W. (2006). Identification and Inference in Nonlinear
Difference-in-Differences Models. *Econometrica*, 74(2), 431-497.
[doi:10.1111/j.1468-0262.2006.00668.x](https://doi.org/10.1111/j.1468-0262.2006.00668.x)

## Stata Implementation

A native Stata implementation is under development in `stata/`. Its command
name and files are intentionally preserved during the R package rename.

## Authors And Affiliations

Paper authors:

- Julien Chhor: Toulouse School of Economics, University of Toulouse Capitole, Toulouse, France
- Xavier D'Haultfoeuille: CREST-ENSAE
- Jérémy L'Hour: CFM & CREST-ENSAE
- Martin Mugnier: Paris School of Economics

Research assistance:

- Timothée Bacchi: research assistant of Martin Mugnier
