# Internal utilities.

.check_unused_dots <- function(...) {
  dots <- list(...)
  if (length(dots) == 0L) {
    return(invisible(NULL))
  }
  dot_names <- names(dots)
  dot_names[dot_names == ""] <- "<unnamed>"
  stop("Unused argument(s): ", paste(dot_names, collapse = ", "), call. = FALSE)
}

.validate_sample <- function(x, label) {
  if (!is.atomic(x) || is.matrix(x) || is.data.frame(x)) {
    warning(label, " should be a numeric vector; coercing to a numeric vector.")
    x <- as.numeric(x)
  }
  if (!is.numeric(x)) {
    stop(label, " must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    warning(label, " contains NA or non-finite values; results may be unreliable.")
  }
  if (length(x) < 4L) {
    stop(label, " must contain at least 4 observations.", call. = FALSE)
  }
  x
}

.resolve_epsilon_n <- function(epsilon_n, n2) {
  if (is.null(epsilon_n)) {
    epsilon_n <- function(n) 1 / log(n)
  }
  if (is.function(epsilon_n)) {
    epsilon_n <- epsilon_n(n2)
  }
  if (!is.numeric(epsilon_n) ||
      length(epsilon_n) != 1L ||
      !is.finite(epsilon_n) ||
      epsilon_n <= 0) {
    stop("epsilon_n must be a positive finite numeric scalar or a function returning one.",
         call. = FALSE)
  }
  epsilon_n
}

.resolve_bootstrap_B <- function(method, B) {
  if (!any(c("bse", "bpc") %in% method)) {
    return(as.integer(B))
  }
  if (is.null(B) || length(B) != 1L || is.na(B)) {
    stop("B must not be NULL or NA.", call. = FALSE)
  }
  B <- as.integer(B)
  if (B < 200L) {
    warning("B < 200: the bootstrap standard error may be unstable.")
    B <- 200L
  }
  B
}

# Left-continuous empirical quantile function with F_Y^{-1}(0) = Y_(1).
.prepare_left_quantile <- function(x) {
  xs <- sort(x, method = "radix")
  n <- length(xs)
  function(p) xs[pmax(ceiling(p * n), 1L)]
}

.theta_components <- function(sample1, sample2, sample3) {
  F3 <- stats::ecdf(sample3)
  Uhat <- F3(sample2)
  qcdf_transform <- .prepare_left_quantile(sample1)(Uhat)
  theta_hat <- mean(qcdf_transform)
  list(
    theta_hat = theta_hat,
    Uhat = Uhat,
    qcdf_transform = qcdf_transform,
    eps_hat = mean((theta_hat - qcdf_transform)^2)
  )
}

.ci_row <- function(method, estimate, se, z_a, sigma_sq = NA_real_) {
  data.frame(
    method = method,
    lower = estimate - z_a * se,
    upper = estimate + z_a * se,
    length = 2 * z_a * se,
    se = se,
    sigma_sq = sigma_sq
  )
}

.percentile_ci_row <- function(method, values, alpha) {
  q_lo <- stats::quantile(values, probs = alpha, names = FALSE)
  q_hi <- stats::quantile(values, probs = 1 - alpha, names = FALSE)
  data.frame(
    method = method,
    lower = q_lo,
    upper = q_hi,
    length = q_hi - q_lo,
    se = NA_real_,
    sigma_sq = NA_real_
  )
}

# Density estimator with cache. The no-split and split estimators use
# h_n2(u) = epsilon_n2 * u * (1 - u), represented on the FY grid here.
.make_density_estimator <- function(U_sorted, FYhat_split) {
  cache <- list()
  n <- length(U_sorted)
  list(
    estimate = function(eps, pointwise = 1) {
      key <- paste(eps, pointwise, sep = "_")
      if (is.null(cache[[key]])) {
        h_vals <- eps * (FYhat_split * (1 - FYhat_split))^pointwise
        cache[[key]] <<- list(
          counts = rect_counts_rcpp(U_sorted, FYhat_split, h_vals),
          h_vals = h_vals
        )
      }
      counts_to_density(cache[[key]]$counts, cache[[key]]$h_vals, n)
    },
    reset = function() {
      cache <<- list()
    }
  )
}

.compute_eta_from_f <- function(f_vals, Uhat, idx_sort, k, ok, n) {
  if (any(!is.finite(f_vals) | f_vals <= 0)) {
    stop("AI density estimates must be positive and finite.", call. = FALSE)
  }
  C2 <- mean(Uhat / f_vals)
  sinvf <- rev(cumsum(rev(1 / f_vals[idx_sort]))) / n
  T1 <- numeric(n)
  T1[ok] <- sinvf[k[ok]]
  mean((T1 - C2)^2)
}

.fast_eta <- function(Ydiff1, fUhat1, Ydiff2, fUhat2, FYhat) {
  u <- as.numeric(Ydiff1 * fUhat1)
  v <- as.numeric(Ydiff2 * fUhat2)
  term2 <- sum(u * FYhat) * sum(v * FYhat)
  Su <- rev(cumsum(rev(u)))
  Sv <- rev(cumsum(rev(v)))
  delta_V <- c(FYhat[1], diff(FYhat))
  sum(delta_V * Su * Sv) - term2
}

.cic_bootstrap_values <- function(Y11, Y10, Y01, Y00, B) {
  n11 <- length(Y11)
  n10 <- length(Y10)
  n01 <- length(Y01)
  n00 <- length(Y00)
  values <- numeric(B)

  for (b in seq_len(B)) {
    Y11_star <- Y11[sample.int(n11, n11, replace = TRUE)]
    Y10_star <- Y10[sample.int(n10, n10, replace = TRUE)]
    Y01_star <- Y01[sample.int(n01, n01, replace = TRUE)]
    Y00_star <- Y00[sample.int(n00, n00, replace = TRUE)]
    values[[b]] <- mean(Y11_star) -
      .theta_components(Y01_star, Y10_star, Y00_star)$theta_hat
  }

  values
}

.panel_no_split_estimator <- function(Y, X, Z, epsilon_n) {
  n1 <- length(Y)
  n2 <- length(X)

  if (n1 != length(Z)) {
    stop("panel_data = TRUE requires length(sample1) to equal length(sample3).")
  }

  q <- floor(n1 / 4)
  r <- floor(n2 / 2)
  if (q < 2 || r < 2) {
    stop("panel_data = TRUE requires at least 8 paired observations and 4 sample2 observations.")
  }

  z_1 <- seq_len(q)
  z_2 <- seq.int(q + 1L, 2L * q)
  y_1 <- seq.int(2L * q + 1L, 3L * q)
  y_2 <- seq.int(3L * q + 1L, 4L * q)
  x_1 <- seq_len(r)
  x_2 <- seq.int(r + 1L, 2L * r)

  FZ1 <- stats::ecdf(Z[z_1])
  FZ2 <- stats::ecdf(Z[z_2])
  FY1_q <- .prepare_left_quantile(Y[y_1])
  FY2_q <- .prepare_left_quantile(Y[y_2])

  theta_1 <- mean(FY1_q(FZ2(X[x_1])))
  theta_2 <- mean(FY2_q(FZ1(X[x_2])))
  theta_hat <- mean(c(theta_1, theta_2))

  U_1 <- FZ1(X[x_1])
  U_2 <- FZ2(X[x_2])
  FYhat <- seq_len(q - 1L) / q

  est_1 <- .make_density_estimator(sort(U_1), FYhat)
  est_2 <- .make_density_estimator(sort(U_2), FYhat)
  fU_1 <- est_1$estimate(epsilon_n, pointwise = 1)
  fU_2 <- est_2$estimate(epsilon_n, pointwise = 1)

  eta_panel <- .fast_eta(diff(sort(Y[y_1])), fU_1, diff(sort(Y[y_2])), fU_2, FYhat)

  eps_1 <- theta_hat - FY1_q(FZ2(X))
  eps_2 <- theta_hat - FY2_q(FZ1(X))
  eps_panel <- mean((eps_1^2 + eps_2^2) / 2)

  sigma_sq <- (2 * q / n1) * eta_panel + (q / n2) * eps_panel

  list(
    theta_hat = theta_hat,
    sigma_sq = sigma_sq,
    se = sqrt(max(sigma_sq, 0) / min(n1, n2)),
    q = q,
    r = r
  )
}

.panel_split_estimator <- function(Y, X, Z, epsilon_n) {
  n1 <- length(Y)
  n2 <- length(X)

  if (n1 != length(Z)) {
    stop("panel_data = TRUE requires length(sample1) to equal length(sample3).")
  }

  q <- floor(n1 / 4)
  r <- floor(n2 / 2)
  if (q < 2 || r < 2) {
    stop("panel_data = TRUE requires at least 8 paired observations and 4 sample2 observations.")
  }

  z_1 <- seq_len(q)
  z_2 <- seq.int(q + 1L, 2L * q)
  y_1 <- seq.int(2L * q + 1L, 3L * q)
  y_2 <- seq.int(3L * q + 1L, 4L * q)
  x_1 <- seq_len(r)
  x_2 <- seq.int(r + 1L, 2L * r)

  FZ1 <- stats::ecdf(Z[z_1])
  FZ2 <- stats::ecdf(Z[z_2])
  FY1_q <- .prepare_left_quantile(Y[y_1])
  FY2_q <- .prepare_left_quantile(Y[y_2])

  theta_1 <- mean(FY1_q(FZ2(X[x_1])))
  theta_2 <- mean(FY2_q(FZ1(X[x_2])))
  theta_hat <- mean(c(theta_1, theta_2))

  U_1 <- FZ1(X[x_1])
  U_2 <- FZ2(X[x_2])
  FYhat <- seq_len(q - 1L) / q

  est_1 <- .make_density_estimator(sort(U_1), FYhat)
  est_2 <- .make_density_estimator(sort(U_2), FYhat)
  fU_1 <- est_1$estimate(epsilon_n, pointwise = 1)
  fU_2 <- est_2$estimate(epsilon_n, pointwise = 1)

  eta_panel <- .fast_eta(diff(sort(Y[y_1])), fU_1, diff(sort(Y[y_2])), fU_2, FYhat)

  eps_1 <- theta_hat - FY1_q(FZ2(X))
  eps_2 <- theta_hat - FY2_q(FZ1(X))
  eps_panel <- mean((eps_1^2 + eps_2^2) / 2)

  N <- min(n1, n2)
  sigma_sq <- (2 * N / n1) * eta_panel + (N / n2) * eps_panel

  list(
    theta_hat = theta_hat,
    sigma_sq = sigma_sq,
    se = sqrt(max(sigma_sq, 0) / N),
    q = q,
    r = r
  )
}
