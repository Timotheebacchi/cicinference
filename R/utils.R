# ── Internal utilities ─────────────────────────────────────────────────────────
# These functions are not exported. They are called by cic_inference() internally.

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
            "compilation: install.packages('cicinference', type='source'). ",
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
# active implementation used by cic_inference().
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


.panel_no_split_estimator <- function(Y, X, Z, epsilon_n) {
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
  h_1 <- epsilon_n * U_1 * (1 - U_1)
  h_2 <- epsilon_n * U_2 * (1 - U_2)
  FYhat <- seq_len(q - 1L) / q

  est_1 <- .make_density_estimator(sort(U_1), FYhat)
  est_2 <- .make_density_estimator(sort(U_2), FYhat)
  fU_1 <- est_1$estimate(h_1, pointwise = 1)
  fU_2 <- est_2$estimate(h_2, pointwise = 1)

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

# Returns a list with point estimate, asymptotic variance, and standard error.
.panel_split_estimator <- function(Y, X, Z, epsilon_n) {
  n1 <- length(Y)
  n2 <- length(X)

  if (n1 != length(Z)) {
    stop("panel_data = TRUE requires length(Y) to equal length(Z).")
  }

  # Split n1 into 4 parts and n2 into 2 parts as defined in the LaTeX theorem
  q <- floor(n1 / 4)
  r <- floor(n2 / 2)
  if (q < 2 || r < 2) {
    stop("panel_data = TRUE requires at least 8 paired observations and 4 X observations.")
  }

  # 4 distinct splits for (Y, Z) pairs
  z_1 <- seq_len(q)
  z_2 <- seq.int(q + 1L, 2L * q)
  y_1 <- seq.int(2L * q + 1L, 3L * q)
  y_2 <- seq.int(3L * q + 1L, 4L * q)
  
  # 2 distinct splits for X
  x_1 <- seq_len(r)
  x_2 <- seq.int(r + 1L, 2L * r)

  # Compute empirical CDFs and left-quantiles
  FZ1 <- stats::ecdf(Z[z_1])
  FZ2 <- stats::ecdf(Z[z_2])
  FY1_q <- .prepare_left_quantile(Y[y_1])
  FY2_q <- .prepare_left_quantile(Y[y_2])

  # Step 2: Point estimation (\tilde{\theta} calculation)
  theta_1 <- mean(FY1_q(FZ2(X[x_1])))
  theta_2 <- mean(FY2_q(FZ1(X[x_2])))
  theta_hat <- mean(c(theta_1, theta_2))

  #  Ranks and density estimation (\check{f}_U^{(1)} and \check{f}_U^{(2)})
  U_1 <- FZ1(X[x_1])
  U_2 <- FZ2(X[x_2])
  h_1 <- epsilon_n * U_1 * (1 - U_1)
  h_2 <- epsilon_n * U_2 * (1 - U_2)
  FYhat <- seq_len(q - 1L) / q

  est_1 <- .make_density_estimator(sort(U_1), FYhat)
  est_2 <- .make_density_estimator(sort(U_2), FYhat)
  fU_1 <- est_1$estimate(h_1, pointwise = 1)
  fU_2 <- est_2$estimate(h_2, pointwise = 1)

  #  Double integral via the fast eta routine (\hat{E}[\eta^2])
  eta_panel <- .fast_eta(diff(sort(Y[y_1])), fU_1, diff(sort(Y[y_2])), fU_2, FYhat)

  # Compute residuals (\check{\eps}_i^{(j)}) over all X observations
  eps_1 <- theta_hat - FY1_q(FZ2(X))
  eps_2 <- theta_hat - FY2_q(FZ1(X))
  eps_panel <- mean((eps_1^2 + eps_2^2) / 2)

  # Asymptotic variance scaling matching the \widehat{\widetilde{\sigma}}^2 formula
  # We use N = min(n1, n2) as the standard package reference sample size.
  N <- min(n1, n2)
  
  # First term:  (2 * N / n1) * Integral
  # Second term: (N / n2^2) * sum(eps^2 / 2) = (N / n2) * mean(eps^2 / 2)
  sigma_sq <- (2 * N / n1) * eta_panel + (N / n2) * eps_panel

  list(
    theta_hat = theta_hat,
    sigma_sq = sigma_sq,
    se = sqrt(max(sigma_sq, 0) / N),
    q = q,
    r = r
  )
}







