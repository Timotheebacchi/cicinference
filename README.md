  # cic_inference : Inference Changes-in-Changes Estimator 

Library of variance estimators in Changes-in-Changes (CiC) design

 ## Setup 

  While it is not available on CRAN : 
  ```r
  devtools::install_github("Timotheebacchi/cic_inference")
  ```



## Arguments 
```r
cic_inference(
Y Numeric vector, 
X Numeric vector,
Z Numeric vector, 
 method :  Character vector,
#Options are:
#        "no-split" : Nonparametric method using full sample
#       "split" : Sample-splitting variance estimator
#      "kde" : Epanechnikov KDE variance estimator
#     "bse" : Bootstrap standard-error method
#    "bpc" : Bootstrap percentile method
B Integer,
# Number of bootstrap replications (default: 1000, take at least B>=200). Only used for "bse" and "bpc".
epsilon_n Numeric,
#Bandwidth multiplier epsilon_n used in h_{n_2,u} = epsilon_n u(1-u). The default is epsilon_n = 1/log(n_2).
level Numeric,
# Confidence level for intervals (default: 0.95)
panel_data Logical,
# if TRUE, use the panel-data estimator based on a paired (Y, Z) sample.
timings Logical,
# if TRUE, print elapsed time after each major computation block.
)
```

```r
sim_dgp
(
n integer,
panel_data Logical 
#if True the data created is a panel_data sample
)
```
## Example

  ```r
  library(cic.newassumptions.newvarianceestimator)
  set.seed(2026) 
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
  fit1 <- cic_inference(
    d2$Y, d2$X, d2$Z,
    method = c("no-split", "split", "kde", "bse", "bpc"),
    timings = TRUE
  )
  summary(fit1)
  d3 <- sim_dgp(100000, panel_data = TRUE)
  #With Panel data
  fit_panel <- cic_inference(
      d3$Y, d3$X, d3$Z,
      method = c("no-split","split"),
      panel_data = TRUE,
      timings = TRUE
    )
    
  fit_nopanel  <- cic_inference(
      d3$Y, d3$X, d3$Z,
      method = c("no-split","split"),
      timings = TRUE
    )
    summary(fit_panel)
    summary(fit_nopanel)
  ```

  ## Assumptions Reminder

  The package is built for the CiC setup described in [the manuscript](https://arxiv.org/abs/2607.00219). Before interpreting the output, the input data should be checked against the main assumptions used by the estimator

  ## Indication to use the package

  - Split,no-split and kde methods are completely usable with gigantic samples. If it runs for more than a minute, even for n = O(10^6) there must be a problem
  - The default bandwith for kde is h = 1/log(n2) as adviced in the manuscript but it is possible to change it. However take care to have a bandwidth which respects assumptions of the paper (assumptions 3)
  - Bootstrap methods are usable on samples containing millions of elements but is expected to compute for a longer time in this case
  
  ## Reference

  Chhor, J., D'Haultfoeuille, X., L'Hour, J., & Mugnier, M. (2026). Asymptotic Properties of Empirical Quantile-Based Estimators. Manuscript : https://arxiv.org/abs/2607.00219 , DOI : 2607.00219

  ## License

  MIT + file LICENSE

  ## Author

    Julien Chhor,Xavier D'Haultfoeuille, Jeremy L'Hour, Martin Mugnier, Timothée Bacchi (Research Assistant)
