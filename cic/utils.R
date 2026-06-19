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
#   U_sorted    : sorted vector of estimated ranks (Uhat)
#   FYhat_split : empirical CDF grid, typically (1:(n-1)) / n
#
# Returns a list with:
#   $estimate(eps, pointwise) : evaluates the density at FYhat_split
#   $reset()                  : clears the cache
.make_density_estimator <- function(U_sorted, FYhat_split) {
  cache <- list()
  n     <- length(U_sorted)
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
    reset = function() { cache <<- list() }
  )
}

# Computes the eta term from a pre-computed density vector f_vals.
# Called three times in cic() with h/2, h, 2h to estimate eta_ai.
#
# Arguments:
#   f_vals   : density evaluated at qcdf_transform points  (length n)
#   Uhat     : estimated ranks FZ(X)                       (length n)
#   idx_sort : order(Uhat)
#   k        : findInterval indices
#   ok       : logical mask (k <= n)
#   n        : sample size
.compute_eta_from_f <- function(f_vals, Uhat, idx_sort, k, ok, n) {
  C2    <- mean(Uhat / f_vals)
  sinvf <- rev(cumsum(rev(1 / f_vals[idx_sort]))) / n
  T1    <- numeric(n)
  T1[ok] <- sinvf[k[ok]]
  mean((T1 - C2)^2)
}

# Computes the fast eta variance term without building an n x n matrix.
# See fast_eta() in the original file for the mathematical derivation.
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
