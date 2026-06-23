#' @title Changes-in-Changes Estimator
#' @description Computes the CiC estimator for nonlinear difference-in-differences models
#' as described in Athey & Imbens (2006).
#' @useDynLib cic, .registration = TRUE
#' @importFrom Rcpp evalCpp         
#' @import stats
#' @param Y Numeric vector: Outcome variable
#' @param X Numeric vector: Treatment/endogenous variable
#' @param Z Numeric vector: Instrument/exogenous variable
#' @param method Character vector: Estimation method(s). Options are:
#'   \itemize{
#'     \item{"no-split"}{Nonparametric method using full sample (fastest)}
#'     \item{"bse"}{Bootstrap standard-error method}
#'     \item{"bpc"}{Bootstrap percentile method}
#'   }
#' @param B Integer: Number of bootstrap replications (default: 200). Only used for "bse" and "bpc".
#' @param h Numeric: Bandwidth parameter for no-split method (default: NULL = auto)
#' @param level Numeric: Confidence level for intervals (default: 0.95)
#' @return An S3 object of class 'cic' with elements:
#'   \item{theta_hat}{Estimated CiC parameter}
#'   \item{ci}{Data frame with confidence intervals (columns: lower, upper, length, method)}
#'   \item{n}{Sample sizes (Y, X, Z)}
#'   \item{method}{Estimation method(s) used}
#'   \item{h}{Bandwidth used (for no-split)}
#'   \item{level}{Confidence level}
#' @examples
#' set.seed(42)
#' d <- sim_dgp(500)
#' fit <- cic(d$Y, d$X, d$Z, method = "no-split")
#' summary(fit)
#' @seealso \code{\link{check_cic_assumptions}}, \code{\link{sim_dgp}}
#' @references Athey, S., & Imbens, G. W. (2006). Identification and inference in
#'   nonlinear difference-in-differences models. Econometrica, 74(2), 431-497.
#' @export
cic <- function(Y, X, Z, method = c("no-split", "bse", "bpc"), B = 200, h = NULL, level = 0.95) {

  # ── Input checks ─────────────────────────────────────────────────────────────
  method <- match.arg(method, several.ok = TRUE)
  stopifnot(
    is.numeric(Y), is.numeric(X), is.numeric(Z),
    length(Y) >= 4, length(X) >= 4, length(Z) >= 4,
    is.null(h) || (is.numeric(h) && length(h) == 1 && h > 0),
    is.numeric(level), level > 0, level < 1
  )
  if (any(c("bse", "bpc") %in% method)) {
    # Sanitize B: check for NA and NULL before conversion
    if (is.null(B) || (is.na(B))) {
      stop("B must not be NULL or NA.")
    }
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

  if ("no-split" %in% method) {
    # Use the user-provided or default Silverman bandwidth in the no-split
    # adaptive density estimator. This makes fit$h the actual smoothing
    # parameter used by the full-sample variance estimator.
    eps_bw <- h

    # No-split: use full samples, not divided into halves
    Ysortdiff <- diff(sort(Y))
    FYhat     <- seq_len(n1 - 1) / n1
    
    est_full   <- .make_density_estimator(sort(Uhat), FYhat)
    fUhat      <- est_full$estimate(eps_bw, pointwise = 1)
    fUhat_unif <- est_full$estimate(eps_bw, pointwise = 0)
    fUhat_t2   <- est_full$estimate(2 * eps_bw, pointwise = 1)
    fUhat_d2   <- est_full$estimate(eps_bw / 2, pointwise = 1)

    eta_hat    <- .fast_eta(Ysortdiff, fUhat,      Ysortdiff, fUhat,      FYhat)
    eta_unif   <- .fast_eta(Ysortdiff, fUhat_unif, Ysortdiff, fUhat_unif, FYhat)
    eta_t2     <- .fast_eta(Ysortdiff, fUhat_t2,   Ysortdiff, fUhat_t2,   FYhat)
    eta_d2     <- .fast_eta(Ysortdiff, fUhat_d2,   Ysortdiff, fUhat_d2,   FYhat)

    eta_nosplit <- (eta_hat + eta_unif + eta_t2 + eta_d2) / 4
    sigma_sq    <- lbda1_3 * eta_nosplit + lbda2 * eps_hat
    se          <- sqrt(sigma_sq / N)

    ci_rows[["no-split"]] <- data.frame(
      method = "no-split",
      lower  = theta_hat - z_a * se,
      upper  = theta_hat + z_a * se,
      length = 2 * z_a * se
    )
    rm(Ysortdiff, FYhat, fUhat, fUhat_unif, fUhat_t2, fUhat_d2, est_full)
  }

  if (any(c("bse", "bpc") %in% method)) {
    # Use internal R implementation of bootstrap core
    # (fallback to pure R if Rcpp version unavailable)
    tryCatch({
      boot_vals <-  boot_core(sort(Y), X, sort(Z), B = B)
    }, error = function(e) {
      stop(
        "Bootstrap computation failed. This may occur if the compiled C++ code ",
        "is unavailable or corrupted. Please ensure the package is properly ",
        "installed and compiled: install.packages('cic', type='source'). ",
        "Original error: ", conditionMessage(e)
      )
    })

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
      h         = if ("no-split" %in% method) h else NA_real_
    ),
    class = "cic"
  )
}
