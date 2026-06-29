  # cic: Changes-in-Changes Estimator

  `cic` is an R package for the Changes-in-Changes estimator and asymptotic inference for empirical quantile-based estimators. It computes the plug-in estimate, provides several confidence interval methods, and includes a diagnostic helper to check whether the input data look compatible with the model assumptions.

  ## Features

  - Point estimation of the CiC parameter from `Y`, `X`, and `Z`
  - Confidence intervals with five methods:
    - `"no-split"`
    - `"split"`
    - `"kde"`
    - `"bse"`
    - `"bpc"`
  - Diagnostic checks via `check_cic_assumptions()`
  - Rcpp-backed computation for the core routines, with a pure R fallback where available
  - Simulation helpers `sim_dgp()`, `qY_dgp()`, and `theta_true()`

  ## Installation

  ```r
  devtools::install_github("Timotheebacchi/cic_package")
  ```

  ## Quick Start

  ```r
  library(cic)

  set.seed(42)
  d <- sim_dgp(500)

  diag <- check_cic_assumptions(d$Y, d$X, d$Z)
  diag$pass_all

  fit <- cic(d$Y, d$X, d$Z, method = c("no-split", "split", "kde", "bse" , "bpc")

  fit
  fit$ci
  summary(fit)
  ```

  ## Diagnostics

  `check_cic_assumptions()` returns a list with:

  - `pass_all`: overall logical result
  - `metrics`: sample ratios, tail indices, boundary estimates, and related checks
  - `messages`: warnings or success messages

  It is designed as a quick pre-check before running `cic()` on empirical data.

  ## Methods

  | Method | Description |
  |---|---|
  | `"no-split"` | Full-sample nonparametric variance estimator |
  | `"split"` | Sample-splitting variance estimator |
  | `"kde"` | Epanechnikov KDE-based variance estimator |
  | `"bse"` | Bootstrap standard-error interval |
  | `"bpc"` | Bootstrap percentile interval |

  You can supply one method or several methods at once; `cic()` returns the intervals in the order requested.

  ## Package Contents

  - `cic()` for estimation and confidence intervals
  - `check_cic_assumptions()` for diagnostics
  - `sim_dgp()` for a reproducible data-generating process
  - `qY_dgp()` and `theta_true()` for simulation support

  ## Reference

  Chhor, J., D'Haultfœuille, X., L'Hour, J., & Mugnier, M. (2026). Asymptotic Properties of Empirical Quantile-Based Estimators. Manuscript.

  ## License

  MIT + file LICENSE

  ## Author

  Timothée Bacchi
