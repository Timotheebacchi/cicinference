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

  fit <- cic(d$Y, d$X, d$Z, method = c("no-split", "split", "kde", "bse", "bpc"))

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

  ## Rappel mathématique des méthodes

  Soit $F_Z$ la fonction de répartition empirique de $Z$ et $Q_Y(u)=F_Y^{-1}(u)$ la quantile left-continuous de $Y$. Le point estimateur commun à toutes les méthodes est

  $$
  \hat\theta = \frac{1}{n}\sum_{i=1}^n Q_Y(F_Z(X_i)).
  $$

  L’intervalle de confiance est ensuite construit autour de $\hat\theta$ avec des estimateurs de variance ou par bootstrap. En notation compacte, les méthodes asymptotiques utilisent une forme du type

  $$
  \hat\theta \pm z_{1-\alpha}\,\widehat{\mathrm{se}},
  $$

  où $\alpha = (1-\text{level})/2$ et $z_{1-\alpha}$ est le quantile normal. Le code sépare aussi les contributions de variance via

  $$
  \hat\sigma^2 = \lambda_{1,3}\,\hat\eta + \lambda_2\,\hat\varepsilon,
  $$

  avec $\hat\varepsilon = n^{-1}\sum_i (\hat\theta - Q_Y(F_Z(X_i)))^2$, $\lambda_{1,3} = N(n_1+n_3)/(n_1n_3)$, $\lambda_2 = N/n_2$ et $N = \min(n_1,n_2,n_3)$.

  | Méthode | Ce qu’elle fait | Différence principale |
  |---|---|---|
  | "no-split" | Estime la partie non paramétrique de variance sur l’échantillon complet, avec un lissage piloté par `h`. | C’est la méthode asymptotique la plus directe et la plus efficace en termes d’utilisation des données. |
  | "split" | Coupe les échantillons en deux moitiés, ré-estime les composantes de densité sur chaque moitié, puis forme l’intervalle avec une variance de sample splitting. | Elle réduit le biais de réutilisation des données, au prix d’une perte d’information et d’une variance souvent plus grande. |
  | "kde" | Approxime la densité à l’aide d’un noyau d’Epanechnikov appliqué aux scores $Q_Y(F_Z(X_i))$. | Elle remplace l’estimation par comptage local par une estimation de densité au noyau, utile comme alternative lissée. |
  | "bse" | Rééchantillonne les triplets $(Y,X,Z)$, calcule $\hat\theta^*$ à chaque réplication, puis prend l’écart-type bootstrap comme erreur-type. | L’intervalle est centré sur $\hat\theta$ et hérite de la dispersion empirique bootstrap. |
  | "bpc" | Rééchantillonne comme "bse", mais utilise directement les quantiles empiriques de $\hat\theta^*$ pour construire l’intervalle. | L’intervalle n’est pas forcé d’être symétrique autour de $\hat\theta$ et suit mieux l’asymétrie finie-échantillon. |

  En pratique, `cic()` peut recevoir une seule méthode ou plusieurs méthodes à la fois, et renvoie les intervalles dans l’ordre demandé. Pour lire rapidement les résultats: `no-split`, `split` et `kde` sont les trois variantes asymptotiques principales, tandis que `bse` et `bpc` servent surtout de contrepoints bootstrap pour comparer la robustesse en petit échantillon.

  ## Package Contents

  - `cic()` for estimation and confidence intervals
  - `check_cic_assumptions()` for diagnostics
  - `sim_dgp()` for a reproducible data-generating process
  - `qY_dgp()` and `theta_true()` for simulation support

  ## Reference

  Chhor, J., D'Haultfoeuille, X., L'Hour, J., & Mugnier, M. (2026). Asymptotic Properties of Empirical Quantile-Based Estimators. Manuscript.

  ## License

  MIT + file LICENSE

  ## Author

  Timothée Bacchi
