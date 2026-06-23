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
            "compilation: install.packages('cic', type='source'). ",
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

# Default plug-in bandwidth for the Epanechnikov KDE.
# Rule: Silverman's rule-of-thumb adapted for the quantile density context.
# Can be overridden by the user via the h argument of cic().
.default_bandwidth <- function(Y) {
  n  <- length(Y)
  sd_Y <- stats::sd(Y)
  # Silverman (1986): h = 1.06 * sigma * n^{-1/5}
  1.06 * sd_Y * n^(-1/5)
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

#' @title Check CiC Model Assumptions
#' @description Validates that data satisfy the theoretical requirements for the CiC estimator.
#' @param Y Numeric vector: Outcome variable
#' @param X Numeric vector: Treatment/endogenous variable
#' @param Z Numeric vector: Instrument/exogenous variable
#' @return A list with elements:
#'   \item{pass_all}{Logical: TRUE if all checks passed}
#'   \item{metrics}{Named list of diagnostic metrics:
#'     \itemize{
#'       \item{lambda1, lambda2, lambda3}{Sample size ratios}
#'       \item{d1_left_tail, d2_right_tail}{Tail indices}
#'       \item{b1_left_boundary, b2_right_boundary}{Boundary density estimates}
#'       \item{autocorrelation_p}{Ljung-Box test p-value}
#'       \item{uniqueness_ratio}{Ratio of unique values}
#'     }
#'   }
#'   \item{messages}{Character vector: Warnings or success messages}
#' #' @export
check_cic_assumptions <- function(Y, X, Z) {
  # Input validation
  stopifnot(
    is.numeric(Y) && length(Y) >= 4,
    is.numeric(X) && length(X) >= 4,
    is.numeric(Z) && length(Z) >= 4
  )
  
  n1 <- length(Y)
  n2 <- length(X)
  n3 <- length(Z)
  N <- n1 + n2 + n3
  
  messages <- character(0)
  metrics <- list()
  
  # ─ Assumption 1: Sampling Evaluation ──────────────────────────────────────
  
  # Ratio test
  lambda1 <- n1 / N
  lambda2 <- n2 / N
  lambda3 <- n3 / N
  
  metrics$lambda1 <- lambda1
  metrics$lambda2 <- lambda2
  metrics$lambda3 <- lambda3
  
  if (any(c(lambda1, lambda2, lambda3) < 0.05)) {
    messages <- c(
      messages,
      paste0(
        "Warning [Assumption 1(iii)]: Sample size ratios are highly imbalanced. ",
        "Asymptotic convergence may be unstable. ",
        sprintf("(λ₁=%.3f, λ₂=%.3f, λ₃=%.3f)", lambda1, lambda2, lambda3)
      )
    )
  }
  
  # Autocorrelation test (Ljung-Box)
  lb_Y <- stats::Box.test(Y, lag = 1, type = "Ljung-Box")
  lb_X <- stats::Box.test(X, lag = 1, type = "Ljung-Box")
  lb_Z <- stats::Box.test(Z, lag = 1, type = "Ljung-Box")
  
  metrics$ljung_box_Y_pval <- lb_Y$p.value
  metrics$ljung_box_X_pval <- lb_X$p.value
  metrics$ljung_box_Z_pval <- lb_Z$p.value
  
  if (any(lb_Y$p.value < 0.05, lb_X$p.value < 0.05, lb_Z$p.value < 0.05)) {
    messages <- c(
      messages,
      paste0(
        "Warning [Assumption 1(i)]: Serial correlation detected via Ljung-Box test. ",
        "The i.i.d. assumption may be violated. ",
        sprintf(
          "(p-values: Y=%.4f, X=%.4f, Z=%.4f)",
          lb_Y$p.value, lb_X$p.value, lb_Z$p.value
        )
      )
    )
  }
  
  # ─ Assumption 2(i)/(ii): Continuity Evaluation ────────────────────────────
  
  uniqueness_Y <- length(unique(Y)) / length(Y)
  uniqueness_Z <- length(unique(Z)) / length(Z)
  
  metrics$uniqueness_Y <- uniqueness_Y
  metrics$uniqueness_Z <- uniqueness_Z
  
  if (uniqueness_Y < 0.98 || uniqueness_Z < 0.98) {
    messages <- c(
      messages,
      paste0(
        "Warning [Assumption 2(i)/(ii)]: High frequency of tied values detected. ",
        "The distributions may not be absolutely continuous, which violates ",
        "model smoothness rules. ",
        sprintf("(Y uniqueness: %.3f, Z uniqueness: %.3f)", uniqueness_Y, uniqueness_Z)
      )
    )
  }
  
  # ─ Assumption 2(ii-iv): Tail Heaviness & Convergence Rate ────────────────
  
  # Step A: Construct U-hat
  FZ <- stats::ecdf(Z)
  Uhat <- FZ(X)
  
  # Step B: Estimate tail indices
  d1 <- .estimate_tail_index(Y, tail = "lower", quantile_cutoff = 0.10)
  d2 <- .estimate_tail_index(Y, tail = "upper", quantile_cutoff = 0.10)
  
  metrics$d1_left_tail <- d1
  metrics$d2_right_tail <- d2
  
  # Step C: Estimate boundary density parameters
  boundary_est <- .estimate_boundary_density(Uhat)
  b1 <- boundary_est$b1
  b2 <- boundary_est$b2
  
  metrics$b1_left_boundary <- b1
  metrics$b2_right_boundary <- b2
  
  # Step D: Convergence check
  convergence_left <- b1 + d1
  convergence_right <- b2 + d2
  
  metrics$convergence_left <- convergence_left
  metrics$convergence_right <- convergence_right
  
  pass_all <- TRUE
  
  if (convergence_left >= 0.5 || convergence_right >= 0.5) {
    messages <- c(
      messages,
      paste0(
        "Convergence Alert [Assumption 2(iv)]: The combined tail heaviness of ",
        "your data (b + d ≥ 0.5) exceeds the theoretical threshold. ",
        "The estimator will experience severe convergence slowdowns or ",
        "non-normal asymptotics. ",
        sprintf(
          "(Left: b₁+d₁=%.3f, Right: b₂+d₂=%.3f)",
          convergence_left, convergence_right
        )
      )
    )
    pass_all <- FALSE
  }
  
  # If no messages, add success message
  if (length(messages) == 0) {
    messages <- "All checks passed. Data appear consistent with CiC assumptions."
  }
  
  # Return structured output
  structure(
    list(
      pass_all = pass_all,
      metrics = metrics,
      messages = messages
    ),
    class = "cic_diagnostics"
  )
}
#' @examples
#' set.seed(42)
#' d <- sim_dgp(500)
#' diag <- check_cic_assumptions(d$Y, d$X, d$Z)
#' print(diag)
#' @seealso \code{\link{cic}}, \code{\link{sim_dgp}}



