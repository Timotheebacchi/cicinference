#' Fit Empirical Quantile-Based Inference
#'
#' @description
#' Estimates \eqn{\theta = E[F_Y^{-1}(F_Z(X))]} using empirical distribution
#' and left-continuous empirical quantile functions, and reports confidence
#' intervals from the requested inference methods.
#'
#' @useDynLib quantcdf.inference, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import stats
#'
#' @param sample1 Numeric vector: Sample 1, corresponding to \eqn{Y}.
#' @param sample2 Numeric vector: Sample 2, corresponding to \eqn{X}.
#' @param sample3 Numeric vector: Sample 3, corresponding to \eqn{Z}.
#' @param method Character vector. One or more of:
#'   \itemize{
#'     \item{"no-split"}{Full-sample plug-in asymptotic variance estimator}
#'     \item{"split"}{Sample-splitting asymptotic variance estimator}
#'     \item{"ai"}{Athey-Imbens asymptotic variance estimator}
#'     \item{"bse"}{Bootstrap standard-error method}
#'     \item{"bpc"}{Bootstrap percentile method}
#'   }
#' @param ... Reserved for future extensions. Unused arguments are rejected.
#' @param B Integer: number of bootstrap replications. Only used for
#'   \code{"bse"} and \code{"bpc"}.
#' @param epsilon_n Numeric scalar or function. Bandwidth multiplier for the
#'   split and no-split density estimators, where
#'   \eqn{h_{n_2}(u) = \epsilon_{n_2} u(1-u)}. The default is
#'   \eqn{\epsilon_{n_2} = 1 / \log(n_2)}.
#' @param level Numeric confidence level.
#' @param panel_data Logical. If \code{TRUE}, use the paired-sample panel
#'   estimator for \code{"no-split"} and \code{"split"}.
#' @param timings Logical. If \code{TRUE}, print elapsed time after each major
#'   computation block.
#'
#' @return An S3 object of class \code{"quantcdf_fit"}.
#'
#' @examples
#' set.seed(42)
#' d <- sim_dgp(500)
#' out <- fit(d$Y, d$X, d$Z, method = "no-split")
#' summary(out)
#'
#' @seealso \code{\link{sim_dgp}}, \code{\link{cic_fit}}
#' @references Chhor, J., D'Haultfoeuille, X., L'Hour, J., & Mugnier, M. (2026).
#'   Asymptotic Properties of Empirical Quantile-Based Estimators. Manuscript.
#' @export
fit <- function(sample1,
                sample2,
                sample3,
                method = c("no-split", "split", "ai", "bse", "bpc"),
                ...,
                B = 1000,
                epsilon_n = NULL,
                level = 0.95,
                panel_data = FALSE,
                timings = FALSE) {
  .check_unused_dots(...)

  method <- match.arg(method, several.ok = TRUE)
  sample1 <- .validate_sample(sample1, "Sample 1")
  sample2 <- .validate_sample(sample2, "Sample 2")
  sample3 <- .validate_sample(sample3, "Sample 3")

  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("level must be a numeric scalar between 0 and 1.", call. = FALSE)
  }

  if (isTRUE(panel_data) && !all(method %in% c("no-split", "split"))) {
    stop("panel_data = TRUE is currently supported only for method = 'no-split' or 'split'.",
         call. = FALSE)
  }

  uniq1 <- length(unique(sample1)) / max(1L, length(sample1))
  if (uniq1 < 0.98) {
    warning("Sample 1 contains many tied values; results may be unreliable (uniqueness < 0.98).")
  }

  B <- .resolve_bootstrap_B(method, B)
  n1 <- as.double(length(sample1))
  n2 <- as.double(length(sample2))
  n3 <- as.double(length(sample3))
  N <- min(n1, n2, n3)
  epsilon_n <- .resolve_epsilon_n(epsilon_n, n2)

  alpha <- (1 - level) / 2
  z_a <- stats::qnorm(1 - alpha)

  timings <- isTRUE(timings)
  timing_start <- proc.time()[["elapsed"]]
  timing_last <- timing_start
  log_timing <- function(stage) {
    if (!timings) {
      return(invisible(NULL))
    }
    timing_now <- proc.time()[["elapsed"]]
    message(sprintf(
      "quantcdf timing [%s]: %.3f s (stage), %.3f s (total)",
      stage, timing_now - timing_last, timing_now - timing_start
    ))
    timing_last <<- timing_now
    invisible(NULL)
  }

  ci_rows <- list()
  h_ai <- NA_real_
  theta_hat <- NULL

  if (!isTRUE(panel_data)) {
    theta_parts <- .theta_components(sample1, sample2, sample3)
    theta_hat <- theta_parts$theta_hat
    Uhat <- theta_parts$Uhat
    qcdf_transform <- theta_parts$qcdf_transform
    eps_hat <- theta_parts$eps_hat
    lbda1_3 <- N * (n1 + n3) / (n1 * n3)
    lbda2 <- N / n2
  }

  if ("no-split" %in% method) {
    if (isTRUE(panel_data)) {
      panel_fit <- .panel_no_split_estimator(sample1, sample2, sample3, epsilon_n)
      theta_hat <- panel_fit$theta_hat
      sigma_sq <- panel_fit$sigma_sq
      se <- panel_fit$se
      log_timing("panel no-split estimator")
    } else {
      Ysortdiff <- diff(sort(sample1))
      FYhat <- seq_len(n1 - 1) / n1
      est_full <- .make_density_estimator(sort(Uhat), FYhat)
      fUhat <- est_full$estimate(epsilon_n, pointwise = 1)
      eta_nosplit <- .fast_eta(Ysortdiff, fUhat, Ysortdiff, fUhat, FYhat)
      sigma_sq <- lbda1_3 * eta_nosplit + lbda2 * eps_hat
      se <- sqrt(max(sigma_sq, 0) / N)
    }

    ci_rows[["no-split"]] <- .ci_row("no-split", theta_hat, se, z_a, sigma_sq)
    log_timing("variance no-split")
  }

  if ("split" %in% method) {
    if (isTRUE(panel_data)) {
      panel_split_fit <- .panel_split_estimator(sample1, sample2, sample3, epsilon_n)
      theta_hat <- panel_split_fit$theta_hat
      sigma_sq_split <- panel_split_fit$sigma_sq
      se_split <- panel_split_fit$se
      log_timing("panel split estimator")
    } else {
      n_half <- min(floor(n1 / 2), floor(n2 / 2), floor(n3 / 2))
      if (n_half < 2) {
        stop("split method requires at least two observations in the first half-sample.",
             call. = FALSE)
      }

      F3_1 <- stats::ecdf(sample3[seq_len(n_half)])
      F3_2 <- stats::ecdf(sample3[seq.int(n3 - n_half + 1L, n3)])
      Uhat1 <- F3_1(sample2[seq_len(n_half)])
      Uhat2 <- F3_2(sample2[seq.int(n2 - n_half + 1L, n2)])

      Ysort1diff <- diff(sort(sample1[seq_len(n_half)]))
      Ysort2diff <- diff(sort(sample1[seq.int(n1 - n_half + 1L, n1)]))
      FYhat_split <- seq_len(n_half - 1L) / n_half
      est1 <- .make_density_estimator(sort(Uhat1), FYhat_split)
      est2 <- .make_density_estimator(sort(Uhat2), FYhat_split)
      fUhat1 <- est1$estimate(epsilon_n, pointwise = 1)
      fUhat2 <- est2$estimate(epsilon_n, pointwise = 1)

      eta_hat_split <- .fast_eta(Ysort1diff, fUhat1, Ysort2diff, fUhat2, FYhat_split)
      sigma_sq_split <- lbda1_3 * eta_hat_split + lbda2 * eps_hat
      se_split <- sqrt(max(sigma_sq_split, 0) / N)
    }

    ci_rows[["split"]] <- .ci_row("split", theta_hat, se_split, z_a, sigma_sq_split)
    log_timing("variance split")
  }

  if ("ai" %in% method) {
    sd1 <- stats::sd(sample1)
    if (!is.finite(sd1) || sd1 <= 0) {
      stop("The AI method requires a positive finite empirical standard deviation for Sample 1.",
           call. = FALSE)
    }
    h_ai <- 1.06 * n1^(-1 / 5) / sd1
    if (!is.finite(h_ai) || h_ai <= 0) {
      stop("The AI bandwidth must be positive and finite.", call. = FALSE)
    }

    n_ai <- length(Uhat)
    idx_sort <- order(Uhat)
    U_sort <- Uhat[idx_sort]
    grid <- seq_len(n_ai) / n_ai
    k <- findInterval(grid - 1e-12, U_sort) + 1L
    ok <- k <= n_ai

    f_one <- f_y_hat_epnechikov(sample1, qcdf_transform, rep(h_ai, length(qcdf_transform)))
    eta_ai <- .compute_eta_from_f(f_one, Uhat, idx_sort, k, ok, n_ai)

    sigma_sq_ai <- 2 * eta_ai + eps_hat
    se_ai <- sqrt(max(sigma_sq_ai, 0) / N)

    ci_rows[["ai"]] <- .ci_row("ai", theta_hat, se_ai, z_a, sigma_sq_ai)
    log_timing("variance ai")
  }

  if (any(c("bse", "bpc") %in% method)) {
    tryCatch({
      boot_vals <- boot_core(sort(sample1), sample2, sort(sample3), B = B)
    }, error = function(e) {
      stop(
        "Bootstrap computation failed. This may occur if the compiled C++ code is unavailable. ",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    })

    if ("bpc" %in% method) {
      ci_rows[["bpc"]] <- .percentile_ci_row("bpc", boot_vals, alpha)
    }

    if ("bse" %in% method) {
      se_boot <- stats::sd(boot_vals)
      ci_rows[["bse"]] <- .ci_row("bse", theta_hat, se_boot, z_a, NA_real_)
    }
    log_timing("bootstrap")
  }

  ci <- do.call(rbind, ci_rows[method])
  rownames(ci) <- NULL

  structure(
    list(
      theta_hat = theta_hat,
      ci = ci,
      level = level,
      n = c(n1 = n1, n2 = n2, n3 = n3),
      method = method,
      panel_data = isTRUE(panel_data),
      epsilon_n = if (any(c("no-split", "split") %in% method)) epsilon_n else NA_real_,
      h_ai = h_ai
    ),
    class = "quantcdf_fit"
  )
}

#' Fit Changes-in-Changes Inference
#'
#' @description
#' Estimates the Changes-in-Changes contrast
#' \eqn{E[Y_{11}] - E[F_{Y_{01}}^{-1}(F_{Y_{00}}(Y_{10}))]}.
#' The counterfactual component is obtained from \code{\link{fit}} rather than
#' duplicated.
#'
#' @param Y11 Numeric vector: treated post-period sample.
#' @param Y10 Numeric vector: treated pre-period sample.
#' @param Y01 Numeric vector: control post-period sample.
#' @param Y00 Numeric vector: control pre-period sample.
#' @param method Character vector. One or more of \code{"no-split"},
#'   \code{"split"}, \code{"ai"}, \code{"bse"}, or \code{"bpc"}.
#' @param ... Reserved for future extensions. Unused arguments are rejected.
#' @inheritParams fit
#'
#' @return An S3 object of class \code{c("quantcdf_cic_fit", "quantcdf_fit")}.
#'
#' @examples
#' set.seed(42)
#' Y11 <- rnorm(200, 1)
#' Y10 <- rnorm(200)
#' Y01 <- rnorm(200)
#' Y00 <- rnorm(200)
#' cic_fit(Y11, Y10, Y01, Y00, method = "no-split")
#'
#' @export
cic_fit <- function(Y11,
                    Y10,
                    Y01,
                    Y00,
                    method = c("no-split", "split", "ai", "bse", "bpc"),
                    ...,
                    B = 1000,
                    epsilon_n = NULL,
                    level = 0.95,
                    timings = FALSE) {
  .check_unused_dots(...)

  method <- match.arg(method, several.ok = TRUE)
  Y11 <- .validate_sample(Y11, "Y11")
  Y10 <- .validate_sample(Y10, "Y10")
  Y01 <- .validate_sample(Y01, "Y01")
  Y00 <- .validate_sample(Y00, "Y00")

  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("level must be a numeric scalar between 0 and 1.", call. = FALSE)
  }

  B <- .resolve_bootstrap_B(method, B)
  alpha <- (1 - level) / 2
  z_a <- stats::qnorm(1 - alpha)

  analytic_methods <- intersect(method, c("no-split", "split", "ai"))
  theta_methods <- analytic_methods
  if (any(c("bse", "bpc") %in% method)) {
    theta_methods <- unique(c(theta_methods, "no-split"))
  }

  theta_fit <- fit(
    sample1 = Y01,
    sample2 = Y10,
    sample3 = Y00,
    method = theta_methods,
    B = B,
    epsilon_n = epsilon_n,
    level = level,
    timings = timings
  )

  cic_hat <- mean(Y11) - theta_fit$theta_hat
  n11 <- as.double(length(Y11))
  n10 <- as.double(length(Y10))
  n01 <- as.double(length(Y01))
  n00 <- as.double(length(Y00))
  N <- min(n11, n10, n01, n00)
  V11 <- N / n11 * (mean(Y11^2) - mean(Y11)^2)

  ci_rows <- list()

  for (analytic_method in analytic_methods) {
    theta_row <- theta_fit$ci[theta_fit$ci$method == analytic_method, , drop = FALSE]
    sigma_theta_sq <- theta_row$sigma_sq[[1L]]
    sigma_cic_sq <- sigma_theta_sq + V11
    se_cic <- sqrt(max(sigma_cic_sq, 0) / N)
    ci_rows[[analytic_method]] <- .ci_row(analytic_method, cic_hat, se_cic, z_a, sigma_cic_sq)
  }

  if (any(c("bse", "bpc") %in% method)) {
    boot_vals <- .cic_bootstrap_values(Y11, Y10, Y01, Y00, B)

    if ("bpc" %in% method) {
      ci_rows[["bpc"]] <- .percentile_ci_row("bpc", boot_vals, alpha)
    }

    if ("bse" %in% method) {
      se_boot <- stats::sd(boot_vals)
      ci_rows[["bse"]] <- .ci_row("bse", cic_hat, se_boot, z_a, NA_real_)
    }
  }

  ci <- do.call(rbind, ci_rows[method])
  rownames(ci) <- NULL

  structure(
    list(
      theta_hat = cic_hat,
      ci = ci,
      level = level,
      n = c(n11 = n11, n10 = n10, n01 = n01, n00 = n00),
      method = method,
      panel_data = FALSE,
      epsilon_n = theta_fit$epsilon_n,
      h_ai = theta_fit$h_ai,
      components = list(
        mean_Y11 = mean(Y11),
        counterfactual_theta_hat = theta_fit$theta_hat,
        counterfactual_fit = theta_fit,
        V11 = V11,
        N = N
      )
    ),
    class = c("quantcdf_cic_fit", "quantcdf_fit")
  )
}
