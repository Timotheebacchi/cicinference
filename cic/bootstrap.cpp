#include <Rcpp.h>
#include <algorithm>  // std::lower_bound, std::fill
#include <cmath>      // std::ceil
#include <vector>
using namespace Rcpp;

// ── Bootstrap core ────────────────────────────────────────────────────────────
// Computes B bootstrap replications of theta_hat.
// Each replication:
//   1. Resamples X with replacement => Xb
//   2. Resamples Z with replacement via multinomial counts => bootstrapped ecdf of Z
//   3. Evaluates theta_hat on (Ys, Xb, bootstrapped Z ecdf)
//
// Inputs:
//   Ys : sorted Y sample  (length n1)
//   Xs : X sample         (length n2, need not be sorted)
//   Zs : sorted Z sample  (length n3)
//   B  : number of bootstrap replications (>= 200 recommended)
//
// Returns a NumericVector of length B with the bootstrap estimates.

// [[Rcpp::export]]
NumericVector boot_core(NumericVector Ys, NumericVector Xs,
                        NumericVector Zs, int B) {
  int n1 = Ys.size(), n2 = Xs.size(), n3 = Zs.size();
  NumericVector results(B);
  std::vector<int>    counts(n3);
  std::vector<double> Xb(n2);

  GetRNGstate();
  for (int b = 0; b < B; b++) {

    // Resample X
    for (int i = 0; i < n2; i++)
      Xb[i] = Xs[(int)(unif_rand() * n2)];

    // Multinomial resample of Z => bootstrapped ecdf
    std::fill(counts.begin(), counts.end(), 0);
    for (int i = 0; i < n3; i++) counts[(int)(unif_rand() * n3)]++;
    std::vector<double> cdf_z(n3);
    int cumul = 0;
    for (int j = 0; j < n3; j++) {
      cumul   += counts[j];
      cdf_z[j] = (double)cumul / n3;
    }

    // Compute theta_hat on bootstrap sample
    double s = 0.0;
    for (int i = 0; i < n2; i++) {
      int pos = (int)(std::lower_bound(Zs.begin(), Zs.end(), Xb[i]) - Zs.begin());
      double u = (pos == 0) ? 0.0 : (pos >= n3) ? 1.0 : cdf_z[pos - 1];
      int idx_q = (int)std::ceil(u * n1);
      if (idx_q < 1)  idx_q = 1;
      if (idx_q > n1) idx_q = n1;
      s += Ys[idx_q - 1];
    }
    results[b] = s / n2;
  }
  PutRNGstate();
  return results;
}
