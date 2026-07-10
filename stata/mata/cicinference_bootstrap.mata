version 16.0

mata:
mata set matastrict on

/*
Inputs:
  sorted_values: ascending real column vector.
  probability: scalar probability in [0, 1].
Outputs:
  R default quantile type 7 at probability.
Assumptions:
  sorted_values is nonempty and sorted ascending.
*/
real scalar cic_quantile_type7_sorted(real colvector sorted_values,
                                      real scalar probability)
{
    real scalar n
    real scalar h
    real scalar lower_index
    real scalar upper_index
    real scalar weight

    n = rows(sorted_values)
    if (n < 1 | probability < 0 | probability > 1) {
        _error(3498)
    }
    if (n == 1) {
        return(sorted_values[1])
    }

    h = (n - 1) * probability + 1
    lower_index = floor(h)
    upper_index = ceil(h)
    weight = h - lower_index

    if (lower_index < 1) {
        return(sorted_values[1])
    }
    if (upper_index > n) {
        return(sorted_values[n])
    }

    return((1 - weight) * sorted_values[lower_index] +
           weight * sorted_values[upper_index])
}

/*
Inputs:
  values: real column vector.
  probability: scalar probability in [0, 1].
Outputs:
  R default quantile type 7 at probability.
Assumptions:
  values contains finite numbers.
*/
real scalar cic_quantile_type7(real colvector values,
                               real scalar probability)
{
    return(cic_quantile_type7_sorted(cic_sort_ascending(values), probability))
}

/*
Inputs:
  sorted_y: sorted outcome support from the original sample.
  x: original X sample in observation order.
  sorted_z: sorted Z support from the original sample.
  y_indices, x_indices, z_indices: B x n index matrices, 1-based.
Outputs:
  B x 1 vector of bootstrap theta values.
Assumptions:
  Matches boot_core() conditional on common resampling indices:
  Y and Z are represented by bootstrap counts over sorted supports; X is
  resampled as observations from the original X vector.
*/
real colvector cic_bootstrap_values_with_indices(real colvector sorted_y,
                                                 real colvector x,
                                                 real colvector sorted_z,
                                                 real matrix y_indices,
                                                 real matrix x_indices,
                                                 real matrix z_indices)
{
    real scalar n_y
    real scalar n_x
    real scalar n_z
    real scalar B
    real scalar b
    real scalar i
    real scalar draw_index
    real scalar pos
    real scalar u
    real scalar rank_y
    real scalar y_position
    real scalar running
    real colvector theta_boot
    real colvector counts_y
    real colvector counts_z
    real colvector cdf_y_counts
    real colvector cdf_z

    n_y = rows(sorted_y)
    n_x = rows(x)
    n_z = rows(sorted_z)
    B = rows(y_indices)

    if (cols(y_indices) != n_y | cols(x_indices) != n_x |
        cols(z_indices) != n_z | rows(x_indices) != B |
        rows(z_indices) != B) {
        _error(3200)
    }

    theta_boot = J(B, 1, .)

    for (b = 1; b <= B; b++) {
        counts_y = J(n_y, 1, 0)
        counts_z = J(n_z, 1, 0)

        for (i = 1; i <= n_y; i++) {
            draw_index = y_indices[b, i]
            if (draw_index < 1 | draw_index > n_y | draw_index >= .) {
                _error(3300)
            }
            counts_y[draw_index] = counts_y[draw_index] + 1
        }

        for (i = 1; i <= n_z; i++) {
            draw_index = z_indices[b, i]
            if (draw_index < 1 | draw_index > n_z | draw_index >= .) {
                _error(3300)
            }
            counts_z[draw_index] = counts_z[draw_index] + 1
        }

        cdf_y_counts = J(n_y, 1, .)
        cdf_z = J(n_z, 1, .)
        running = 0
        for (i = 1; i <= n_y; i++) {
            running = running + counts_y[i]
            cdf_y_counts[i] = running
        }

        running = 0
        for (i = 1; i <= n_z; i++) {
            running = running + counts_z[i]
            cdf_z[i] = running / n_z
        }

        running = 0
        for (i = 1; i <= n_x; i++) {
            draw_index = x_indices[b, i]
            if (draw_index < 1 | draw_index > n_x | draw_index >= .) {
                _error(3300)
            }

            pos = cic_upper_bound_sorted(sorted_z, x[draw_index]) - 1
            if (pos == 0) {
                u = 0
            }
            else {
                u = cdf_z[pos]
            }

            rank_y = ceil(u * n_y)
            if (rank_y < 1) {
                rank_y = 1
            }
            if (rank_y > n_y) {
                rank_y = n_y
            }
            y_position = cic_lower_bound_sorted(cdf_y_counts, rank_y)
            running = running + sorted_y[y_position]
        }

        theta_boot[b] = running / n_x
    }

    return(theta_boot)
}

/*
Inputs:
  theta_hat: original-sample point estimate.
  boot_values: B x 1 vector of bootstrap theta values.
  level: confidence level on R scale, e.g. 0.95.
Outputs:
  Row vector:
    1 se_boot
    2 bse_lower
    3 bse_upper
    4 bse_length
    5 bpc_lower
    6 bpc_upper
    7 bpc_length
    8 B
Assumptions:
  Uses sample standard deviation and R quantile type 7.
*/
real rowvector cic_bootstrap_summary_from_values(real scalar theta_hat,
                                                 real colvector boot_values,
                                                 real scalar level)
{
    real scalar B
    real scalar alpha
    real scalar z_alpha
    real scalar se_boot
    real scalar bse_lower
    real scalar bse_upper
    real scalar bse_length
    real scalar bpc_lower
    real scalar bpc_upper
    real scalar bpc_length

    B = rows(boot_values)
    if (B < 2 | level <= 0 | level >= 1) {
        _error(3498)
    }

    alpha = (1 - level) / 2
    z_alpha = invnormal(1 - alpha)
    se_boot = sqrt(sum((boot_values :- mean(boot_values)):^2) / (B - 1))
    bse_lower = theta_hat - z_alpha * se_boot
    bse_upper = theta_hat + z_alpha * se_boot
    bse_length = 2 * z_alpha * se_boot
    bpc_lower = cic_quantile_type7(boot_values, alpha)
    bpc_upper = cic_quantile_type7(boot_values, 1 - alpha)
    bpc_length = bpc_upper - bpc_lower

    return((se_boot, bse_lower, bse_upper, bse_length,
            bpc_lower, bpc_upper, bpc_length, B))
}

/*
Inputs:
  y, x, z: original samples.
  y_indices, x_indices, z_indices: B x n common bootstrap index matrices.
  level: confidence level on R scale, e.g. 0.95.
Outputs:
  Row vector:
    1 theta_hat
    2 se_boot
    3 bse_lower
    4 bse_upper
    5 bse_length
    6 bpc_lower
    7 bpc_upper
    8 bpc_length
    9 B
    10 n_y
    11 n_x
    12 n_z
Assumptions:
  This validation helper is deterministic conditional on supplied indices.
*/
real rowvector cic_bootstrap_fit_with_indices(real colvector y,
                                              real colvector x,
                                              real colvector z,
                                              real matrix y_indices,
                                              real matrix x_indices,
                                              real matrix z_indices,
                                              real scalar level)
{
    real scalar theta_hat
    real colvector boot_values
    real rowvector summary

    theta_hat = cic_rank_quantile_mean(y, x, z)
    boot_values = cic_bootstrap_values_with_indices(
        cic_sort_ascending(y),
        x,
        cic_sort_ascending(z),
        y_indices,
        x_indices,
        z_indices)
    summary = cic_bootstrap_summary_from_values(theta_hat, boot_values, level)

    return((theta_hat, summary, rows(y), rows(x), rows(z)))
}

/*
Inputs:
  sample_size: number of observations in a source sample.
  B: number of bootstrap replications.
Outputs:
  B x sample_size matrix of 1-based bootstrap draw indices.
Assumptions:
  Uses Mata's current RNG state; callers can set Stata's seed beforehand.
*/
real matrix cic_bootstrap_index_matrix(real scalar sample_size,
                                       real scalar B)
{
    if (sample_size < 1 | B < 1) {
        _error(3498)
    }

    return(floor(sample_size * runiform(B, sample_size)) :+ 1)
}

/*
Inputs:
  y, x, z: original samples.
  B: number of bootstrap replications.
  level: confidence level on R scale, e.g. 0.95.
Outputs:
  Same 12-column row vector as cic_bootstrap_fit_with_indices().
Assumptions:
  This is the live Stata bootstrap path. Exact equality with R is not expected
  from matching integer seeds because R and Stata RNGs differ.
*/
real rowvector cic_bootstrap_fit(real colvector y,
                                 real colvector x,
                                 real colvector z,
                                 real scalar B,
                                 real scalar level)
{
    return(cic_bootstrap_fit_with_indices(
        y,
        x,
        z,
        cic_bootstrap_index_matrix(rows(y), B),
        cic_bootstrap_index_matrix(rows(x), B),
        cic_bootstrap_index_matrix(rows(z), B),
        level))
}

end
