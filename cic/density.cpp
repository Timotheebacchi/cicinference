#include <Rcpp.h>
#include <algorithm>  // std::lower_bound, std::upper_bound
#include <cmath>      // std::sqrt, std::abs, std::ceil
using namespace Rcpp;

// ── Epanechnikov kernel density estimator ─────────────────────────────────────
// Evaluates the KDE with Epanechnikov kernel at each point in y,
// using the sample Y and bandwidth h.
// Used to estimate f_Y at the transformed quantile points.

// [[Rcpp::export]]
NumericVector f_y_hat_epanechnikov(NumericVector Y, NumericVector y, double h) {
  int n = Y.size();
  int m = y.size();
  NumericVector res(m);
  double inv_h  = 1.0 / h;
  double sqrt5  = std::sqrt(5.0);

  for (int i = 0; i < m; i++) {
    double s = 0.0;
    for (int j = 0; j < n; j++) {
      double u = (Y[j] - y[i]) * inv_h;
      if (std::abs(u) < sqrt5) {
        s += (1.0 - (u * u) / 5.0) * 3.0 / (4.0 * sqrt5);
      }
    }
    res[i] = s * inv_h / n;
  }
  return res;
}

// ── Rectangular counts (step 1/2 of adaptive rect. density) ──────────────────
// For each evaluation point x_eval[j] with half-bandwidth h_vals[j],
// counts how many points in X_sorted fall in (x - h, x + h).
// Uses binary search => O(m * log(n)) instead of O(m * n).
// X_sorted MUST be sorted in ascending order before calling.

// [[Rcpp::export]]
IntegerVector rect_counts_rcpp(NumericVector X_sorted, NumericVector x_eval,
                                NumericVector h_vals) {
  const int m = x_eval.size();
  IntegerVector counts(m);
  const double* xs   = X_sorted.begin();
  const double* xend = X_sorted.end();

  for (int j = 0; j < m; j++) {
    const double x = x_eval[j];
    const double h = h_vals[j];
    const double* lo = std::lower_bound(xs, xend, x - h);
    const double* hi = std::upper_bound(lo,  xend, x + h);
    counts[j] = (int)(hi - lo);
  }
  return counts;
}

// ── Rectangular density (step 2/2) ───────────────────────────────────────────
// Converts counts from rect_counts_rcpp into a density estimate:
//   f_hat[j] = counts[j] / (2 * n * h_vals[j])
// Separated from the count step so that rect_counts_rcpp can be cached
// and reused with different bandwidth scalings (see make_density_estimator in R).

// [[Rcpp::export]]
NumericVector counts_to_density(IntegerVector counts, NumericVector h_vals,
                                 int n) {
  const int m = counts.size();
  NumericVector f_hat(m);
  const double inv2n = 1.0 / (2.0 * n);
  for (int j = 0; j < m; j++)
    f_hat[j] = counts[j] * inv2n / h_vals[j];
  return f_hat;
}
