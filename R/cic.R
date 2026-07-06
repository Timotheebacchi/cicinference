#' @title Changes-in-Changes Estimator
#' @description Computes the CiC estimator for nonlinear difference-in-differences models
#' as described in the manuscript by Chhor, D'Haultfoeuille, L'Hour, and Mugnier.
#' @useDynLib cic.newassumptions.newvarianceestimator, .registration = TRUE
#' @importFrom Rcpp evalCpp         
#' @import stats
#' @param Y Numeric vector: Outcome variable
#' @param X Numeric vector: Treatment/endogenous variable
#' @param Z Numeric vector: Instrument/exogenous variable
#' @param method Character vector: Estimation method(s). Options are:
#'   \itemize{
#'     \item{"no-split"}{Nonparametric method using full sample (fastest)}
#'     \item{"split"}{Sample-splitting variance estimator}
#'     \item{"kde"}{Epanechnikov KDE variance estimator}
#'     \item{"bse"}{Bootstrap standard-error method}
#'     \item{"bpc"}{Bootstrap percentile method}
#'   }
#' @param B Integer: Number of bootstrap replications (default: 200). Only used for "bse" and "bpc".
#' @param h Numeric: Bandwidth multiplier epsilon_n used in h_{n_2,u} = epsilon_n u(1-u). The default is 1/log(n_2).
#' @param level Numeric: Confidence level for intervals (default: 0.95)
#' @param panel_data Logical: if TRUE, use the panel-data estimator based on a paired (Y, Z) sample.
#' @param timings Logical: if TRUE, print elapsed time after each major computation block.
#' @return An S3 object of class 'cic' with elements:
#'   \item{theta_hat}{Estimated CiC parameter}
#'   \item{ci}{Data frame with confidence intervals (columns: lower, upper, length, method)}
#'   \item{n}{Sample sizes (Y, X, Z)}
#'   \item{method}{Estimation method(s) used}
#'   \item{h}{Bandwidth used for the nonparametric density estimators}
#'   \item{level}{Confidence level}
#'   \item{panel_data}{Logical: TRUE when the panel-data estimator is used}
#' @examples
#' set.seed(42)
#' d <- sim_dgp(500)
#' fit <- cic(d$Y, d$X, d$Z, method = "no-split")
#' summary(fit)
#' @seealso \code{\link{sim_dgp}}
#' @references Chhor, J., D'Haultfoeuille, X., L'Hour, J., & Mugnier, M. (2026).
#'   Asymptotic Properties of Empirical Quantile-Based Estimators. Manuscript.
#' @export
cic <- function(Y, X, Z, method = c("no-split", "split", "kde", "bse", "bpc"), B = 200, h = NULL, level = 0.95, panel_data = FALSE, timings = FALSE) {

  # ── Input checks ─────────────────────────────────────────────────────────────
  method <- match.arg(method, several.ok = TRUE)
  stopifnot(
    is.numeric(Y), is.numeric(X), is.numeric(Z),
    length(Y) >= 4, length(X) >= 4, length(Z) >= 4,
    is.null(h) || (is.numeric(h) && length(h) == 1 && h > 0),
    is.numeric(level), level > 0, level < 1
  )

  # Lightweight input warnings and coercions
  if (!is.atomic(Y) || is.matrix(Y) || is.data.frame(Y)) {
    warning("Y should be a numeric vector; coercing to a numeric vector.")
    Y <- as.numeric(Y)
  }
  if (!is.atomic(X) || is.matrix(X) || is.data.frame(X)) {
    warning("X should be a numeric vector; coercing to a numeric vector.")
    X <- as.numeric(X)
  }
  if (!is.atomic(Z) || is.matrix(Z) || is.data.frame(Z)) {
    warning("Z should be a numeric vector; coercing to a numeric vector.")
    Z <- as.numeric(Z)
  }

  if (anyNA(Y) || anyNA(X) || anyNA(Z)) {
    warning("Input contains NA or non-finite values; these will be dropped or may cause errors.")
  }

  uniqY <- length(unique(Y)) / max(1L, length(Y))
  if (uniqY < 0.98) {
    warning("Y contains many tied values; results may be unreliable (uniqueness < 0.98).")
  }

  # Fixed logical check for panel methods compatibility
  if (isTRUE(panel_data) && !all(method %in% c("no-split", "split"))) {
    stop("panel_data = TRUE is currently supported only for method = 'no-split' or 'split'. Please set method accordingly.")
  }

  if (any(c("bse", "bpc") %in% method)) {
    if (is.null(B) || is.na(B)) {
      stop("B must not be NULL or NA.")
    }
    B <- as.integer(B)
    if (B < 200L) {
      warning("B < 200: the bootstrap standard error may be unstable. Using B = 200.")
      B <- 200L
    }
  }

  timings <- isTRUE(timings)
  timing_start <- proc.time()[["elapsed"]]
  timing_last <- timing_start
  log_timing <- function(stage) {
    if (!timings) return(invisible(NULL))
    timing_now <- proc.time()[["elapsed"]]
    message(sprintf(
      "cic timing [%s]: %.3f s (stage), %.3f s (total)",
      stage, timing_now - timing_last, timing_now - timing_start
    ))
    timing_last <<- timing_now
    invisible(NULL)
  }

  n1 <- as.double(length(Y)); n2 <- as.double(length(X)); n3 <- as.double(length(Z))
  alpha <- (1 - level) / 2
  z_a   <- stats::qnorm(1 - alpha)

  # ── Bandwidth initialization ──────────────────────────────────────────────────
  if (is.null(h)) h <- .default_bandwidth(n2)

  # ── Baseline setup (Conditional lazy execution) ───────────────────────────────
  ci_rows   <- list()
  theta_hat <- NULL 
  N         <- min(n1, n2, n3)
  lbda1_3   <- N * (n1 + n3) / (n1 * n3)
  lbda2     <- N / n2

  # We only compute cross-sectional matrices if we are NOT in panel data mode
  if (!isTRUE(panel_data)) {
    FZ             <- stats::ecdf(Z)
    Uhat           <- FZ(X)
    qcdf_transform <- .prepare_left_quantile(Y)(Uhat)
    theta_hat      <- mean(qcdf_transform)
    eps_hat        <- mean((theta_hat - qcdf_transform)^2)
  }

  # ── Variance & Estimation methods ────────────────────────────────────────────
  
  if ("no-split" %in% method) {
    if (isTRUE(panel_data)) {
      panel_fit <- .panel_no_split_estimator(Y, X, Z, h)
      theta_hat <- panel_fit$theta_hat
      se        <- panel_fit$se
      log_timing("panel no-split estimator")
    } else {
      Ysortdiff   <- diff(sort(Y))
      FYhat       <- seq_len(n1 - 1) / n1
      est_full    <- .make_density_estimator(sort(Uhat), FYhat)
      fUhat       <- est_full$estimate(h, pointwise = 1)
      eta_nosplit <- .fast_eta(Ysortdiff, fUhat, Ysortdiff, fUhat, FYhat)
      sigma_sq    <- lbda1_3 * eta_nosplit + lbda2 * eps_hat
      se          <- sqrt(sigma_sq / N)
      rm(Ysortdiff, FYhat, fUhat, est_full)
    }

    ci_rows[["no-split"]] <- data.frame(
      method = "no-split",
      lower  = theta_hat - z_a * se,
      upper  = theta_hat + z_a * se,
      length = 2 * z_a * se
    )
    log_timing("variance no-split")
  }

  if ("split" %in% method) {
    if (isTRUE(panel_data)) {
      panel_split_fit <- .panel_estimator(Y, X, Z, h)
      theta_hat       <- panel_split_fit$theta_hat
      se_split        <- panel_split_fit$se
      log_timing("panel split estimator")
    } else {
      n_half <- min(floor(n1 / 2), floor(n2 / 2), floor(n3 / 2))
      if (n_half < 2) {
        stop("split method requires at least two observations in the first half-sample.")
      }

      FZ1   <- stats::ecdf(Z[seq_len(n_half)])
      FZ2   <- stats::ecdf(Z[seq.int(n3 - n_half + 1L, n3)])
      Uhat1 <- FZ1(X[seq_len(n_half)])
      Uhat2 <- FZ2(X[seq.int(n2 - n_half + 1L, n2)])

      Ysort1diff  <- diff(sort(Y[seq_len(n_half)]))
      Ysort2diff  <- diff(sort(Y[seq.int(n1 - n_half + 1L, n1)]))
      FYhat_split <- seq_len(n_half - 1L) / n_half

      est1   <- .make_density_estimator(sort(Uhat1), FYhat_split)
      est2   <- .make_density_estimator(sort(Uhat2), FYhat_split)
      fUhat1 <- est1$estimate(h, pointwise = 1)
      fUhat2 <- est2$estimate(h, pointwise = 1)

      eta_hat_split  <- .fast_eta(Ysort1diff, fUhat1, Ysort2diff, fUhat2, FYhat_split)
      sigma_sq_split <- lbda1_3 * eta_hat_split + lbda2 * eps_hat
      se_split       <- sqrt(max(sigma_sq_split, 0) / N)
    }

    ci_rows[["split"]] <- data.frame(
      method = "split",
      lower  = theta_hat - z_a * se_split,
      upper  = theta_hat + z_a * se_split,
      length = 2 * z_a * se_split
    )
    log_timing("variance split")
  }

  if ("kde" %in% method) {
    n_kde    <- length(Uhat)
    idx_sort <- order(Uhat)
    U_sort   <- Uhat[idx_sort]
    grid     <- seq_len(n_kde) / n_kde
    k        <- findInterval(grid - 1e-12, U_sort) + 1L
    ok       <- k <= n_kde

    f_one  <- f_y_hat_epnechikov(Y, qcdf_transform, h)
    eta_ai <- .compute_eta_from_f(f_one, Uhat, idx_sort, k, ok, n_kde)

    sigma_sq_kde <- 2 * eta_ai + eps_hat
    se_kde       <- sqrt(max(sigma_sq_kde, 0) / N)

    ci_rows[["kde"]] <- data.frame(
      method = "kde",
      lower  = theta_hat - z_a * se_kde,
      upper  = theta_hat + z_a * se_kde,
      length = 2 * z_a * se_kde
    )
    log_timing("variance kde")
  }

  if (any(c("bse", "bpc") %in% method)) {
    tryCatch({
      boot_vals <- boot_core(sort(Y), X, sort(Z), B = B)
    }, error = function(e) {
      stop(
        "Bootstrap computation failed. This may occur if the compiled C++ code is unavailable. ",
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
    log_timing("bootstrap")
  }

  ci <- do.call(rbind, ci_rows[method])
  rownames(ci) <- NULL

  # ── Output ───────────────────────────────────────────────────────────────────
  structure(
    list(
      theta_hat  = theta_hat,
      ci         = ci,
      level      = level,
      n          = c(n1 = n1, n2 = n2, n3 = n3),
      method     = method,
      panel_data = isTRUE(panel_data),
      h          = if (any(c("no-split", "kde") %in% method)) h else NA_real_
    ),
    class = "cic"
  )
}