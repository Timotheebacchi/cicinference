#' Changes-in-Changes estimator
#'
#' Computes the Changes-in-Changes (CiC) estimator
#' \deqn{\hat\theta = \frac{1}{n_2}\sum_{j=1}^{n_2} \hat F_Y^{-1}(\hat F_Z(X_j))}
#' and returns asymptotic confidence intervals based on the chosen variance
#' estimator(s).
#'
#' @param Y Numeric vector. Outcome sample of length \eqn{n_1}.
#' @param X Numeric vector. Covariate sample of length \eqn{n_2}.
#' @param Z Numeric vector. Instrument sample of length \eqn{n_3}.
#' @param method Character vector. One or more of `"split"`, `"bse"`, `"bpc"`.
#'   Defaults to all three.
#'   * `"split"` : asymptotic CI using the split-sample variance estimator.
#'   * `"bse"`   : bootstrap CI using the standard-error bootstrap.
#'   * `"bpc"`   : bootstrap CI using the percentile bootstrap.
#' @param B Integer. Number of bootstrap replications. Ignored if neither
#'   `"bse"` nor `"bpc"` is in `method`. Must be >= 200. Default: 999.
#' @param h Numeric or `NULL`. Bandwidth for the Epanechnikov KDE used in the
#'   `"split"` variance estimator. If `NULL` (default), Silverman's
#'   rule-of-thumb is applied automatically.
#' @param level Numeric in (0, 1). Confidence level. Default: 0.95.
#'
#' @return An object of class `"cic"`, a list containing:
#'   * `theta_hat`  : point estimate \eqn{\hat\theta}.
#'   * `ci`         : data frame with columns `method`, `lower`, `upper`,
#'     `length` for each requested method.
#'   * `level`      : confidence level used.
#'   * `n`          : named integer vector `c(n1, n2, n3)`.
#'   * `method`     : methods requested.
#'   * `h`          : bandwidth used (NA if only bootstrap methods requested).
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' Y <- rnorm(n)
#' Z <- rnorm(n)
#' X <- qnorm(pbeta(pnorm(Z), 0.95, 0.95))
#' fit <- cic(Y, X, Z)
#' print(fit)
#' summary(fit)
#'
#' @export
cic <- function(Y, X, Z,
                method = c("split", "bse", "bpc"),
                B      = 999L,
                h      = NULL,
                level  = 0.95) {

  # ── Input checks ─────────────────────────────────────────────────────────────
  method <- match.arg(method, several.ok = TRUE)
  stopifnot(
    is.numeric(Y), is.numeric(X), is.numeric(Z),
    length(Y) >= 4, length(X) >= 4, length(Z) >= 4,
    length(Y) %% 2 == 0, length(X) %% 2 == 0, length(Z) %% 2 == 0,
    is.numeric(level), level > 0, level < 1
  )
  if (any(c("bse", "bpc") %in% method)) {
    B <- as.integer(B)
    if (B < 200L) {
      warning("B < 200: the bootstrap standard error may be unstable. ",
              "Using B = 200.")
      B <- 200L
    }
  }

  n1 <- length(Y); n2 <- length(X); n3 <- length(Z)
  alpha <- (1 - level) / 2
  z_a   <- stats::qnorm(1 - alpha)

  # ── Point estimate ────────────────────────────────────────────────────────────
  FZ             <- stats::ecdf(Z)
  Uhat           <- FZ(X)
  qcdf_transform <- .prepare_left_quantile(Y)(Uhat)
  theta_hat      <- mean(qcdf_transform)
  eps_hat        <- mean((theta_hat - qcdf_transform)^2)

  # ── Bandwidth ─────────────────────────────────────────────────────────────────
  if (is.null(h)) h <- .default_bandwidth(Y)

  # ── Variance estimation ───────────────────────────────────────────────────────
  ci_rows <- list()
  N        <- min(n1, n2, n3)
  lbda1_3  <- N * (n1 + n3) / (n1 * n3)
  lbda2    <- N / n2

  if ("split" %in% method) {
    eps_bw <- 1 / log(n2)

    # Split Z and X into two halves
    FZ1   <- stats::ecdf(Z[1:(n3 / 2)])
    FZ2   <- stats::ecdf(Z[(n3 / 2 + 1):n3])
    Uhat1 <- FZ1(X[1:(n2 / 2)])
    Uhat2 <- FZ2(X[(n2 / 2 + 1):n2])

    Ysort1diff  <- diff(sort(Y[1:(n1 / 2)]))
    Ysort2diff  <- diff(sort(Y[(n1 / 2 + 1):n1]))
    FYhat_split <- seq_len(n1 / 2 - 1) / (n1 / 2)

    est1   <- .make_density_estimator(sort(Uhat1), FYhat_split)
    est2   <- .make_density_estimator(sort(Uhat2), FYhat_split)
    fUhat1 <- est1$estimate(eps_bw, pointwise = 1)
    fUhat2 <- est2$estimate(eps_bw, pointwise = 1)

    eta_split <- .fast_eta(Ysort1diff, fUhat1, Ysort2diff, fUhat2, FYhat_split)
    sigma_sq  <- lbda1_3 * eta_split + lbda2 * eps_hat
    se        <- sqrt(sigma_sq / N)

    ci_rows[["split"]] <- data.frame(
      method = "split",
      lower  = theta_hat - z_a * se,
      upper  = theta_hat + z_a * se,
      length = 2 * z_a * se
    )
    rm(Ysort1diff, Ysort2diff, FYhat_split, fUhat1, fUhat2, est1, est2)
  }

  if (any(c("bse", "bpc") %in% method)) {
    boot_vals <- boot_core(sort(Y), X, sort(Z), B = B)

    if ("bpc" %in% method) {
      q_lo <- stats::quantile(boot_vals, probs = alpha,     names = FALSE)
      q_hi <- stats::quantile(boot_vals, probs = 1 - alpha, names = FALSE)
      ci_rows[["bpc"]] <- data.frame(
        method = "bpc",
        lower  = q_lo,
        upper  = q_hi,
        length = q_hi - q_lo
      )
    }

    if ("bse" %in% method) {
      se_boot <- stats::sd(boot_vals)
      ci_rows[["bse"]] <- data.frame(
        method = "bse",
        lower  = theta_hat - z_a * se_boot,
        upper  = theta_hat + z_a * se_boot,
        length = 2 * z_a * se_boot
      )
    }
  }

  ci <- do.call(rbind, ci_rows[method])  # preserve requested order
  rownames(ci) <- NULL

  # ── Output ───────────────────────────────────────────────────────────────────
  structure(
    list(
      theta_hat = theta_hat,
      ci        = ci,
      level     = level,
      n         = c(n1 = n1, n2 = n2, n3 = n3),
      method    = method,
      h         = if ("split" %in% method) h else NA_real_
    ),
    class = "cic"
  )
}
