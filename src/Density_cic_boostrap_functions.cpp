#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector f_y_hat_epnechikov(NumericVector Y, NumericVector y, double h) {
  int n = Y.size();
  int m = y.size();
  NumericVector res(m);
  double inv_h = 1.0 / h;
  double sqrt5 = std::sqrt(5.0);

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

// [[Rcpp::export]]
IntegerVector rect_counts_rcpp(NumericVector X_sorted, NumericVector x_eval,
                               NumericVector h_vals) {
  const int m = x_eval.size();
  IntegerVector counts(m);
  const double* xs = X_sorted.begin();
  const double* xend = X_sorted.end();

  for (int j = 0; j < m; j++) {
    const double x = x_eval[j];
    const double h = h_vals[j];
    const double* lo = std::lower_bound(xs, xend, x - h);
    const double* hi = std::upper_bound(lo, xend, x + h);
    counts[j] = (int)(hi - lo);
  }
  return counts;
}

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

// [[Rcpp::export]]
NumericVector boot_core(NumericVector Ys, NumericVector Xs,
                        NumericVector Zs, int B) {
  int n1 = Ys.size(), n2 = Xs.size(), n3 = Zs.size();
  NumericVector results(B);
  std::vector<int> counts_y(n1), counts_z(n3), cdf_y_counts(n1);
  std::vector<double> cdf_z(n3);

  GetRNGstate();
  for (int b = 0; b < B; b++) {
    std::fill(counts_y.begin(), counts_y.end(), 0);
    std::fill(counts_z.begin(), counts_z.end(), 0);
    for (int i = 0; i < n1; i++) counts_y[(int)(unif_rand() * n1)]++;
    for (int i = 0; i < n3; i++) counts_z[(int)(unif_rand() * n3)]++;

    int cumul = 0;
    for (int j = 0; j < n1; j++) {
      cumul += counts_y[j];
      cdf_y_counts[j] = cumul;
    }

    cumul = 0;
    for (int j = 0; j < n3; j++) {
      cumul += counts_z[j];
      cdf_z[j] = (double)cumul / n3;
    }

    double s = 0.0;
    for (int i = 0; i < n2; i++) {
      double xb = Xs[(int)(unif_rand() * n2)];
      int pos = (int)(std::upper_bound(Zs.begin(), Zs.end(), xb) - Zs.begin());
      double u = (pos == 0) ? 0.0 : cdf_z[pos - 1];
      int rank_y = (int)std::ceil(u * n1);
      if (rank_y < 1) rank_y = 1;
      if (rank_y > n1) rank_y = n1;
      int idx_y = (int)(std::lower_bound(cdf_y_counts.begin(), cdf_y_counts.end(), rank_y) - cdf_y_counts.begin());
      s += Ys[idx_y];
    }
    results[b] = s / n2;
  }
  PutRNGstate();
  return results;
}