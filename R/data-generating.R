#' @title Quantile Function for Athey & Imbens (2006) DGP
#' @description Quantile function from the Monte Carlo section of Athey & Imbens (2006):
#'   \deqn{F_Y^{-1}(t) = -t^{-d1} + (1-t)^{-d2}}
#' with the convention that d1 = 0 => t^{-d1} = 1, d2 = 0 => (1-t)^{-d2} = 1.
#' @param t Vector of probabilities in [0,1]
#' @param d1 Left tail parameter (default: 0)
#' @param d2 Right tail parameter (default: 0.05)
#' @return Numeric vector of quantiles
#' @export
#' @examples
#' qY_dgp(0.5, d1 = 0, d2 = 0.05)
qY_dgp <- function(t, d1 = 0, d2 = 0.05) {
  term1 <- if (d1 == 0) 1 else t^(-d1)
  term2 <- if (d2 == 0) 1 else (1 - t)^(-d2)
  -term1 + term2
}

#' @title True CiC Parameter for Athey & Imbens (2006) DGP
#' @description Computes the true theta_0 parameter for the DGP:
#'   \deqn{theta_0 = [B(1-b1, 1-b2-d2) - B(1-b1-d1, 1-b2)] / B(1-b1, 1-b2)}
#' where B is the beta function.
#' @param b1 Left boundary parameter (default: 0)
#' @param b2 Right boundary parameter (default: 0.05)
#' @param d1 Left tail parameter (default: 0)
#' @param d2 Right tail parameter (default: 0.05)
#' @return Numeric scalar: true theta_0 value
#' @export
#' @examples
#' theta_true(b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05)
theta_true <- function(b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05) {
  (beta(1 - b1, 1 - b2 - d2) - beta(1 - b1 - d1, 1 - b2)) /
    beta(1 - b1, 1 - b2)
}

#' @title Simulate Data from Athey & Imbens (2006) DGP
#' @description Simulates a dataset from the Data Generating Process (DGP)
#' described in Athey & Imbens (2006).
#' @param n Sample size
#' @param b1 Left boundary parameter (default: 0)
#' @param b2 Right boundary parameter (default: 0.05)
#' @param d1 Left tail parameter (default: 0)
#' @param d2 Right tail parameter (default: 0.05)
#' @param seed Random seed (default: NULL)
#' @return A list with elements:
#' \item{Y}{Outcome variable}
#' \item{X}{Treatment/endogenous variable}
#' \item{Z}{Instrument/exogenous variable}
#' @export
#' @examples
#' set.seed(42)
#' d <- sim_dgp(500)
#' head(d)
sim_dgp <- function(n, b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  W <- runif(n)
  Y <- qY_dgp(W, d1, d2)
  Z <- rnorm(n)
  V <- rbeta(n, 1 - b1, 1 - b2)
  X <- qnorm(V)
  list(Y = Y, X = X, Z = Z)
}