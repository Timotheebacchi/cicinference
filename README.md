  # cic: Changes-in-Changes Estimator

  An R package implementing the plug-in **Changes-in-Changes (CiC) estimator** from Mugnier & d'Haultefoeuille & Chhor & l'Hour with asymptotic confidence intervals and bootstrap alternatives. 
  It also constructs some tests to look if the estimator is usable with the data given and respects assumptions given in Asymptotic Properties of Empirical Quantile-Based Estimators (2026)

  ## Features

  - **Core estimator**: Average of quantile-transformed ranks
  - **Confidence intervals via three methods**:
    - No-split / Full-sample (nonparametric, uses all data)
    - Bootstrap standard-error (bse)
    - Bootstrap percentile (bpc)
  - **Performance**: C++ backend via Rcpp
  - **Validated**: Comprehensive Monte Carlo tests with known true parameter

  ## Installation

  ```r
  # Install from GitHub
  devtools::install_github("Timotheebacchi/cic_package")
  ```

  ## Quick Start

  ```r
  library(cic)

  # Simulate from the DGP (Athey & Imbens 2006)
  set.seed(42)
  n <- 500
  W <- runif(n)
  Y <- -W^(-0.2) + (1-W)^(-0.05)  # Quantile function
  X <- qnorm(rbeta(n, 1, 1.05))   # Covariate
  Z <- rnorm(n)                     # Instrument

  # Check data consistency with CiC assumptions
  diag <- check_cic_assumptions(Y, X, Z)
  if (diag$pass_all) {
    cat("Data passes all diagnostic checks.\n")
  }

  # Estimate with no-split CI
  fit <- cic(Y, X, Z, method = "no-split")

  # View results
  print(fit)

  # Get confidence interval
  fit$ci

  #Get everything
  summary(fit)
  ```

  ## Diagnostic Tool

  Before estimating the CiC model, use `check_cic_assumptions()` to validate that your data satisfy the theoretical requirements:

  ```r
  # Evaluate data against CiC assumptions
  diag <- check_cic_assumptions(Y, X, Z)

  # Returns a list with:
  # - $pass_all: Boolean indicating if all checks passed
  # - $metrics: Calculated diagnostic values (sample ratios, tail indices, etc.)
  # - $messages: Warnings or success messages
  ```

  The diagnostic evaluates:

  1. **Sampling balance** (Assumption 1): Checks sample size ratios and autocorrelation via Ljung-Box test
  2. **Continuity** (Assumption 2): Tests for absolute continuity via uniqueness ratios
  3. **Tail behavior & convergence** (Assumption 2): Estimates tail indices and boundary densities; flags severe convergence issues if $b + d \geq 0.5$

  ## Methods

  The `cic()` function supports:

  | Method | Speed | Notes |
  |--------|-------|-------|
  | `"no-split"` | Fast | Nonparametric, uses full sample |
  | `"bse"` | Medium | Bootstrap standard-error |
  | `"bpc"` | Medium | Bootstrap percentile |

  ## Code Quality & Improvements

  This package has been audited and enhanced with:

  - **Robust input validation**: B parameter sanitization with explicit NA/NULL checks
  - **Graceful error handling**: Informative messages when Rcpp compiled code is unavailable
  - **Bootstrap implementation**: Pure R fallback via `.boot_core()` for maximum compatibility
  - **Diagnostic capabilities**: Comprehensive `check_cic_assumptions()` function for empirical validation

  ## References

  - Athey, S., & Imbens, G. W. (2006). Identification and inference in nonlinear difference-in-differences models. Econometrica, 74(2), 431–497.

  ## License

  MIT + file LICENSE

  ## Author

  Martin Mugnier and Timothée Bacchi
