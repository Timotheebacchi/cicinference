  # cic.newassumptions.newvarianceestimator: Changes-in-Changes Estimator

  `cic.newassumptions.newvarianceestimator` is an R package for the Changes-in-Changes estimator and asymptotic inference for empirical quantile-based estimators. It computes the plug-in estimate, provides several confidence interval methods. This package is based on the inference methods proposed in https://arxiv.org/abs/2607.00219

## Features

* **Point Estimation:** Highly optimized calculation of the CiC parameter from outcome `Y`, endogenous treatment `X`, and instrument/exogenous `Z`.
* **Flexible Inference:** Support for 5 different confidence interval estimation methods:
  * `"no-split"` (Full sample nonparametric method — **fastest**)
  * `"split"` (Sample-splitting variance estimator)
  * `"kde"` (Epanechnikov Kernel Density Estimation variance estimator)
  * `"bse"` (Bootstrap standard-error method)
  * `"bpc"` (Bootstrap percentile method)
* **Panel Data Support:** Optional `panel_data = TRUE` workflow tailored for paired `(Y, Z)` samples (supported for `"no-split"` and `"split"` methods).
* **Execution Profiling:** Built-in `timings = TRUE` flag to log elapsed time across major computation blocks.
* **C++ Acceleration:** Powered by `Rcpp` for core routines, ensuring peak execution speed.

The `cic()` function returns a structured S3 object of class `"cic"`. Under the hood, it is a named list, allowing you to easily extract results for custom plots, tables, or Monte Carlo simulations.

### Object Components

| Element | Type | Description |
| :--- | :--- | :--- |
| `theta_hat` | `numeric` | The point estimate of the Changes-in-Changes treatment effect parameter. |
| `ci` | `data.frame` | A clean data frame containing the calculated confidence intervals with columns: `method`, `lower`, `upper`, and `length`. |
| `n` | `named numeric` | A vector tracking the sample sizes for each input: `n1` ($Y$), `n2` ($X$), and `n3` ($Z$). |
| `method` | `character` | A vector containing the names of all the estimation methods evaluated in this run. |
| `h` | `numeric` | The exact bandwidth value utilized for the density estimation (returns `NA` if only bootstrap methods were used). |
| `level` | `numeric` | The confidence level specified for the intervals (e.g., `0.95`). |
| `panel_data` | `logical` | A boolean flag indicating whether the paired panel-data workflow was used (`TRUE`) or not (`FALSE`). |

  ## Installation

  ```r
  devtools::install_github("Timotheebacchi/cic_package")
  ```

  ## Quick Start

  ```r
  library(cic.newassumptions.newvarianceestimator)

  set.seed(2026) #To match the code of the manuscript
  d1 <- sim_dgp(1000000)
  
  #For Big datasets
  fit <- cic(
      d1$Y, d1$X, d1$Z,
      method = c("no-split", "split", "kde"),
      timings = TRUE
    )
    summary(fit)

  d2 <- sim_dgp(2000)
  #For smaller datasets
  fit1 <- cic(
    d2$Y, d2$X, d2$Z,
    method = c("no-split", "split", "kde", "bse", "bpc"),
    timings = TRUE
  )
  summary(fit1)
  d3 <- sim_dgp(100000, panel_data = TRUE)
  #With Panel data
  fit_panel <- cic(
      d3$Y, d3$X, d3$Z,
      method = c("no-split","split"),
      panel_data = TRUE,
      timings = TRUE
    )
    
  fit_nopanel  <- cic(
      d3$Y, d3$X, d3$Z,
      method = c("no-split","split"),
      timings = TRUE
    )
    summary(fit_panel)
    summary(fit_nopanel)
  ```

  ## Validation and warnings

  The package performs lightweight input validation and emits clear warnings for common
  issues (e.g., non-numeric inputs, mismatched lengths, or invalid bandwidth
  values). Use `sim_dgp()` and the estimation `cic()` for simulation and
  inference; inspect warnings to help debug input problems. Use `sim_dgp_panel()` tu understand the differences between panel and npn_panel data samples.

  ## Assumptions Reminder

  The package is built for the CiC setup described in [the manuscript](https://arxiv.org/abs/2607.00219). Before
  interpreting the output, the input data should be checked against the main
  assumptions used by the estimator:

  - the observations should look approximately i.i.d. and continuous enough for
    the rank transformation to make sense;
  - the empirical quantile function of $Y$ should stay below a boundary envelope
    of the form $C_Y t^{-d_1}(1-t)^{-d_2}$ on the interior of $(0,1)$;
  - the transformed covariate ranks $U = F_Z(X)$ should admit a smooth density
    that can be screened against $C_U u^{-b_1}(1-u)^{-b_2}$;
  - the outcome distribution should have tails that are not too heavy;
  - the combined tail and boundary behavior should stay within the rate
    conditions required by the asymptotic theory;
  - the smoothing bandwidth should be reasonable for the sample size.

  - The package provides only lightweight validation warnings; inspect warnings produced by `cic()` when running your data.
  - Use `summary(fit)` to view the estimation and inference table.
  - In practice, `cic()` can receive one method of inference or several methods at once, and it returns the intervals in the requested order. 

  
  ## Indication to use the package

  - The methods bse and bpc are betterto be used on relatively small samples to be computed fast. It is expected to wait for more than a minute if you wish to use these methods for samples of more than 10^5 elements.
  - Split,no-split and kde methods are completely usable with gigantic samples. If it runs for more than a minute, even for n = O(10^6) there must be a problem
  - The default bandwith for kde is h = 1/log(n2) as adviced in the manuscript but it is possible to change it. However take care to have a bandwidth which respects assumptions of the paper (assumptions 3)
  

  ## Package Contents

  - `cic()` for estimation and confidence intervals
  - `sim_dgp()` for a reproducible data-generating process
  - `qY_dgp()` and `theta_true()` for simulation support
  

  ## Reference

  Chhor, J., D'Haultfoeuille, X., L'Hour, J., & Mugnier, M. (2026). Asymptotic Properties of Empirical Quantile-Based Estimators. Manuscript : https://arxiv.org/abs/2607.00219 , DOI : 

  ## License

  MIT + file LICENSE

  ## Author

  Timothée Bacchi
