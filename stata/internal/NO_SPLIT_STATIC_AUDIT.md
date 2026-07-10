# No-Split Static Audit

Status: static audit only. Stata execution is pending.

## Files Audited

- `stata/mata/cicinference_quantiles.mata`
- `stata/mata/cicinference_density.mata`
- `stata/mata/cicinference_variance.mata`
- `stata/mata/cicinference_main.mata`
- `stata/ado/cicinference_nosplit.ado`
- `stata/tests/test_nosplit.do`

## Findings

1. Mata syntax
   - Functions are defined inside `mata:`/`end` blocks with `matastrict on`.
   - The no-split core returns a fixed 12-column row vector consumed by the
     temporary ado wrapper.
   - Actual Stata parsing has not been executed.

2. Function signatures
   - `cic_nosplit_fit(y, x, z, level)` accepts three column vectors and a
     confidence level on the R scale, e.g. `0.95`.
   - It intentionally implements only cross-sectional no-split.

3. Matrix dimensions
   - `cic_rank_quantile_transform()` returns an `n_x x 2` matrix.
   - No-split density uses `FYhat` and `diff(sort(Y))`, both length `n_y - 1`.
   - Tests compare scalar outputs and key intermediates against R references.

4. Indexing
   - Empirical CDF uses binary-search upper-bound counts.
   - Left quantile uses `ceil(p*n)` with `p=0` mapped to index 1.
   - Rectangle counts use lower-bound left endpoint and upper-bound right
     endpoint, matching the Rcpp closed-interval count.

5. Missing-value handling
   - The temporary ado wrapper uses `marksample` and `markout`, so it works on
     complete rows for the three variables.
   - This differs from the R function's separate-vector API and must be
     revisited for the final public command.

6. Empirical-CDF convention
   - `cic_ecdf_at()` implements `count(z <= x) / n`, matching R `stats::ecdf()`
     on finite inputs.

7. Left-quantile convention
   - `cic_left_quantile()` preserves the package convention
     `Y_(max(ceil(p*n), 1))`.

8. Returned result names
   - The temporary no-split wrapper stores `r(theta_hat)`, `r(se)`,
     `r(ci_lower)`, `r(ci_upper)`, `r(ci_length)`, `r(eps_hat)`, `r(eta)`,
     `r(sigma_sq)`, `r(epsilon_n)`, `r(n_y)`, `r(n_x)`, and `r(n_z)`.

9. Temporary-variable handling
   - The temporary wrapper uses a `tempname` matrix for Mata output.
   - It does not create temporary Stata variables.

## Pending

Run:

```bash
stata-mp -b do stata/tests/test_phase3_all.do
stata-mp -b do stata/tests/test_nosplit.do
```

Do not claim no-split passes until these commands complete successfully in
Stata.
