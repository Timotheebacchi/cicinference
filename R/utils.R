# ── Internal utilities ─────────────────────────────────────────────────────────
# These functions are not exported. They are called by cic() internally.

# Left-continuous quantile function (i.e. F_Y^{-1}(p) = Y_{(ceil(pn))})
# with the convention F_Y^{-1}(0) = Y_{(1)}.
# Returns a *function* that evaluates the quantile at any vector of probabilities.
.prepare_left_quantile <- function(x) {
  xs <- sort(x, method = "radix")
  n  <- length(xs)
  function(p) xs[pmax(ceiling(p * n), 1L)]
}

# Density estimator with cache.
# Avoids recomputing rect_counts_rcpp when only the bandwidth scalar (eps)
# or the pointwise exponent changes — the counts depend only on the shape
# h_vals = eps * (F * (1 - F))^pointwise, but the *binary-search count step*
# depends on h_vals, so we cache (counts, h_vals) keyed by (eps, pointwise).
#
# Arguments:
#   U_sorted : sorted vector of estimated ranks (Uhat)
#   FYhat    : empirical CDF grid, typically (1:(n-1)) / n
#
# Returns a list with:
#   $estimate(eps, pointwise) : evaluates the density at FYhat
#   $reset()                  : clears the cache
#
# NOTE: This function requires the Rcpp-compiled rect_counts_rcpp and
# counts_to_density functions. If these are unavailable, an error is raised
# with instructions for proper installation.
.make_density_estimator <- function(U_sorted, FYhat) {
  cache <- list()
  n     <- length(U_sorted)
  list(
    estimate = function(eps, pointwise = 1) {
      key <- paste(eps, pointwise, sep = "_")
      if (is.null(cache[[key]])) {
        h_vals <- eps * (FYhat * (1 - FYhat))^pointwise
        
        # Wrap Rcpp calls with error handling
        tryCatch({
          counts <- rect_counts_rcpp(U_sorted, FYhat, h_vals)
          density <- counts_to_density(counts, h_vals, n)
        }, error = function(e) {
          stop(
            "Density estimation failed due to unavailable compiled C++ code. ",
            "This function requires rect_counts_rcpp and counts_to_density from ",
            "the Rcpp binding. Please ensure the package is installed with ",
            "compilation: install.packages('cic.newassumptions.newvarianceestimator', type='source'). ",
            "Original error: ", conditionMessage(e)
          )
        })
        
        cache[[key]] <<- list(
          counts = counts,
          h_vals = h_vals
        )
      }
      counts_to_density(cache[[key]]$counts, cache[[key]]$h_vals, n)
    },
    reset = function() { cache <<- list() }
  )
}

.compute_eta_from_f <- function(f_vals, Uhat, idx_sort, k, ok, n) {
  C2     <- mean(Uhat / f_vals)
  sinvf  <- rev(cumsum(rev(1 / f_vals[idx_sort]))) / n
  T1     <- numeric(n)
  T1[ok] <- sinvf[k[ok]]
  mean((T1 - C2)^2)
}

# Computes the fast eta variance term without building an n x n matrix.
# This routine replaced the earlier matrix-based formulation and is the
# active implementation used by cic().
#
# Arguments:
#   Ydiff1, Ydiff2 : diff(sort(Y)) on each half-sample (or full sample)
#   fUhat1, fUhat2 : density estimates at the CDF grid points
#   FYhat          : CDF grid (1:(n-1)) / n
.fast_eta <- function(Ydiff1, fUhat1, Ydiff2, fUhat2, FYhat) {
  u       <- as.numeric(Ydiff1 * fUhat1)
  v       <- as.numeric(Ydiff2 * fUhat2)
  term2   <- sum(u * FYhat) * sum(v * FYhat)
  Su      <- rev(cumsum(rev(u)))
  Sv      <- rev(cumsum(rev(v)))
  delta_V <- c(FYhat[1], diff(FYhat))
  sum(delta_V * Su * Sv) - term2
}

# Default bandwidth multiplier for the local density estimators.
# The local bandwidth used by the estimator is h_{n_2,u} = eps_n * u * (1-u).
# The default multiplier follows the theoretical suggestion 1 / log(n_2).
.default_bandwidth <- function(n2) {
  1 / log(n2)
}

# Bootstrap core: resamples and recomputes point estimates
#
# Arguments:
#   Y_sorted : sorted outcome vector (length n1)
#   X        : covariate vector (length n2)
#   Z_sorted : sorted instrument vector (length n3)
#   B        : number of bootstrap replications
#
# Returns:
#   Numeric vector of length B containing bootstrap replications of theta_hat
.boot_core <- function(Y_sorted, X, Z_sorted, B) {
  n1 <- length(Y_sorted)
  n2 <- length(X)
  n3 <- length(Z_sorted)
  
  # Prepare quantile function
  qY <- .prepare_left_quantile(Y_sorted)
  
  # Bootstrap replications
  boot_vals <- numeric(B)
  
  for (b in seq_len(B)) {
    # Resample with replacement
    Y_b <- Y_sorted[sample.int(n1, n1, replace = TRUE)]
    X_b <- X[sample.int(n2, n2, replace = TRUE)]
    Z_b <- Z_sorted[sample.int(n3, n3, replace = TRUE)]
    
    # Recompute point estimate
    qY_b <- .prepare_left_quantile(sort(Y_b))
    FZ_b <- stats::ecdf(Z_b)
    Uhat_b <- FZ_b(X_b)
    boot_vals[b] <- mean(qY_b(Uhat_b))
  }
  
  boot_vals
}

# Estimate tail index using Hill-style log-log regression on quantile tail
#
# Arguments:
#   x        : numeric vector
#   tail     : "upper" or "lower"; which tail to estimate
#   quantile_cutoff : proportion for top/bottom quantile (default 0.10 for 10%)
#
# Returns:
#   Estimated tail index (bounded below at 0 for light tails)
.estimate_tail_index <- function(x, tail = "upper", quantile_cutoff = 0.10) {
  x <- na.omit(as.numeric(x))
  n <- length(x)
  xs <- sort(x)
  
  if (tail == "upper") {
    k <- max(1, ceiling(n * quantile_cutoff))
    if (k >= n) k <- n - 1
    # Upper k order statistics relative to k+1-th value
    tail_vals <- xs[(n - k + 1):n]
    reference <- xs[n - k]
    # Shift to positive domain for log
    tail_vals <- tail_vals - reference + 1e-10
  } else if (tail == "lower") {
    k <- max(1, ceiling(n * quantile_cutoff))
    if (k >= n) k <- n - 1
    # Lower k order statistics relative to k+1-th value (reversed)
    tail_vals <- xs[1:k]
    reference <- xs[k + 1]
    # Flip to positive domain for log
    tail_vals <- reference - tail_vals + 1e-10
  } else {
    stop("tail must be 'upper' or 'lower'")
  }
  
  if (length(tail_vals) < 2 || any(tail_vals <= 0)) {
    return(0)
  }
  
  # Hill estimator: α = 1 / mean(log(tail_vals / reference))
  # For power-law tails: d = α (tail index parameter)
  log_ratio <- log(tail_vals)
  hill_est <- mean(log_ratio)
  
  if (hill_est <= 0) {
    return(0)  # Light tails (Gaussian-like)
  }
  
  # d_est = 1 / hill_est
  d_est <- 1 / hill_est
  
  # Reasonable bounds for tail indices
  max(0, min(d_est, 5))
}

# Estimate boundary density parameters (power decay at boundaries)
#
# Arguments:
#   U : sorted U = F_Z(X) vector (on [0,1])
#
# Returns:
#   List with b1 (left boundary) and b2 (right boundary) estimates
.estimate_boundary_density <- function(U) {
  U <- na.omit(as.numeric(U))
  n <- length(U)
  
  # Left boundary (near 0): check density within first 5%
  left_threshold <- 0.05
  n_left <- max(1, ceiling(n * left_threshold))
  U_sorted <- sort(U)
  left_vals <- U_sorted[1:n_left]
  
  # Estimate power decay: count how many points are in first 1%, 2%, etc.
  # Power law: density ~ u^{-b}, so cumulative count ~ u^{1-b}
  # Use finite differences to estimate
  counts_1pct <- sum(U <= 0.01)
  counts_5pct <- sum(U <= 0.05)
  
  if (counts_5pct > 0 && counts_1pct > 0) {
    # Ratio: counts_1pct / counts_5pct ~ (0.01 / 0.05)^{1-b}
    # (1-b) * log(0.01/0.05) = log(counts_1pct / counts_5pct)
    ratio <- (counts_1pct + 1e-8) / (counts_5pct + 1e-8)
    log_ratio <- log(ratio)
    log_scale <- log(0.01 / 0.05)
    b1_est <- 1 - log_ratio / log_scale
    b1_est <- max(0, b1_est)  # Non-negative
  } else {
    b1_est <- 0
  }
  
  # Right boundary (near 1): check density within last 5%
  n_right <- max(1, ceiling(n * left_threshold))
  right_vals <- U_sorted[(n - n_right + 1):n]
  
  counts_99pct <- sum(U >= 0.99)
  counts_95pct <- sum(U >= 0.95)
  
  if (counts_95pct > 0 && counts_99pct > 0) {
    ratio <- (counts_99pct + 1e-8) / (counts_95pct + 1e-8)
    log_ratio <- log(ratio)
    log_scale <- log((1 - 0.99) / (1 - 0.95))
    b2_est <- 1 - log_ratio / log_scale
    b2_est <- max(0, b2_est)  # Non-negative
  } else {
    b2_est <- 0
  }
  
  list(b1 = b1_est, b2 = b2_est)
}

.screen_quantile_envelope <- function(Y, d1, d2, grid = seq(0.05, 0.95, length.out = 19)) {
  qhat <- .prepare_left_quantile(Y)(grid)
  envelope <- grid^(-d1) * (1 - grid)^(-d2)
  ratio <- abs(qhat) / pmax(envelope, .Machine$double.eps)
  reference <- stats::median(ratio[is.finite(ratio)])
  if (!is.finite(reference) || reference <= 0) reference <- 1

  list(
    grid = grid,
    qhat = qhat,
    envelope = envelope,
    ratio = ratio,
    reference = reference,
    max_to_median = max(ratio, na.rm = TRUE) / reference,
    exceed_rate = mean(ratio > 2 * reference, na.rm = TRUE)
  )
}

.panel_no_split_estimator <- function(Y, X, Z, h) {
  n1 <- length(Y)
  n2 <- length(X)

  if (n1 != length(Z)) {
    stop("panel_data = TRUE requires length(Y) to equal length(Z).")
  }

  q <- floor(n1 / 4)
  r <- floor(n2 / 2)
  if (q < 2 || r < 2) {
    stop("panel_data = TRUE requires at least 8 paired observations and 4 X observations.")
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
  fU_1 <- est_1$estimate(h, pointwise = 1)
  fU_2 <- est_2$estimate(h, pointwise = 1)

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

.screen_density_envelope <- function(U, b1, b2, h = NULL) {
  U <- na.omit(as.numeric(U))
  n <- length(U)
  if (n < 4) {
    return(list(
      grid = numeric(0),
      fhat = numeric(0),
      envelope = numeric(0),
      ratio = numeric(0),
      reference = NA_real_,
      max_to_median = NA_real_,
      exceed_rate = NA_real_
    ))
  }

  U_sorted <- sort(U, method = "radix")
  FYhat <- seq_len(n - 1L) / n
  h_diag <- h
  if (is.null(h_diag) || !is.finite(h_diag) || h_diag <= 0) {
    h_diag <- .default_bandwidth(n)
  }

  fhat <- tryCatch({
    .make_density_estimator(U_sorted, FYhat)$estimate(h_diag, pointwise = 1)
  }, error = function(e) {
    stats::density(U_sorted, from = 0, to = 1, n = length(FYhat), cut = 0)$y
  })

  envelope <- FYhat^(-b1) * (1 - FYhat)^(-b2)
  ratio <- fhat / pmax(envelope, .Machine$double.eps)
  reference <- stats::median(ratio[is.finite(ratio)])
  if (!is.finite(reference) || reference <= 0) reference <- 1

  list(
    grid = FYhat,
    fhat = fhat,
    envelope = envelope,
    ratio = ratio,
    reference = reference,
    max_to_median = max(ratio, na.rm = TRUE) / reference,
    exceed_rate = mean(ratio > 2 * reference, na.rm = TRUE)
  )
}

# Diagnostic helper removed by request: all functionality solely related to
# `check_cic_assumptions()` has been deleted from this file. The package now
# focuses on core estimation (`cic`) and basic input validation. If you need
# selective diagnostics in the future, reintroduce lightweight checks here.
#' @examples
#' set.seed(42)
#' d <- sim_dgp(500)
#' diag <- check_cic_assumptions(d$Y, d$X, d$Z)
#' print(diag)
#' @seealso \code{\link{cic}}, \code{\link{sim_dgp}}



