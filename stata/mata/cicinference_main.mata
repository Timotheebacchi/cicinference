version 16.0

mata:
mata set matastrict on

/*
Inputs:
  first_value, second_value, third_value: scalar sample sizes or weights.
Outputs:
  Minimum of the three scalar inputs.
Assumptions:
  Inputs are finite real scalars.
*/
real scalar cic_min3(real scalar first_value,
                     real scalar second_value,
                     real scalar third_value)
{
    real scalar current_min

    current_min = first_value
    if (second_value < current_min) {
        current_min = second_value
    }
    if (third_value < current_min) {
        current_min = third_value
    }

    return(current_min)
}

/*
Inputs:
  values: real column vector.
  first_index: first row index.
  last_index: last row index.
Outputs:
  values[first_index..last_index] as a column vector.
Assumptions:
  Indices are valid 1-based Mata row indices.
*/
real colvector cic_slice_rows(real colvector values,
                              real scalar first_index,
                              real scalar last_index)
{
    return(values[(first_index::last_index)])
}

/*
Inputs:
  y: outcome sample.
  x: treatment/endogenous-variable sample.
  z: instrument/exogenous-variable sample.
  level: confidence level on the R scale, e.g. 0.95.
Outputs:
  Row vector:
    1 theta_hat
    2 se
    3 ci_lower
    4 ci_upper
    5 ci_length
    6 eps_hat
    7 eta
    8 sigma_sq
    9 epsilon_n
    10 n_y
    11 n_x
    12 n_z
Assumptions:
  This implements only the cross-sectional no-split R path. It does not
  implement split, KDE, bootstrap, panel data, or the final ado interface.
*/
real rowvector cic_nosplit_fit(real colvector y,
                               real colvector x,
                               real colvector z,
                               real scalar level)
{
    real scalar n_y
    real scalar n_x
    real scalar n_z
    real scalar effective_n
    real scalar epsilon_n
    real scalar theta_hat
    real scalar eps_hat
    real scalar eta_nosplit
    real scalar sigma_sq
    real scalar se
    real scalar alpha
    real scalar z_alpha
    real scalar ci_lower
    real scalar ci_upper
    real scalar ci_length
    real scalar lbda1_3
    real scalar lbda2
    real matrix transform
    real colvector uhat
    real colvector qcdf_transform
    real colvector ysortdiff
    real colvector fyhat
    real colvector sorted_uhat
    real colvector fuhat

    n_y = rows(y)
    n_x = rows(x)
    n_z = rows(z)
    if (n_y < 4 | n_x < 4 | n_z < 4) {
        _error(2001)
    }
    if (level <= 0 | level >= 1) {
        _error(198)
    }

    effective_n = cic_min3(n_y, n_x, n_z)
    epsilon_n = 1 / log(n_x)

    transform = cic_rank_quantile_transform(y, x, z)
    uhat = transform[., 1]
    qcdf_transform = transform[., 2]
    theta_hat = mean(qcdf_transform)
    eps_hat = mean((qcdf_transform :- theta_hat):^2)

    ysortdiff = cic_adjacent_differences(y)
    fyhat = cic_empirical_grid(n_y)
    sorted_uhat = cic_sort_ascending(uhat)
    fuhat = cic_rect_density_sorted(sorted_uhat, fyhat, epsilon_n)
    eta_nosplit = cic_fast_eta_self(ysortdiff, fuhat, fyhat)

    lbda1_3 = effective_n * (n_y + n_z) / (n_y * n_z)
    lbda2 = effective_n / n_x
    sigma_sq = lbda1_3 * eta_nosplit + lbda2 * eps_hat
    se = sqrt(sigma_sq / effective_n)

    alpha = (1 - level) / 2
    z_alpha = invnormal(1 - alpha)
    ci_lower = theta_hat - z_alpha * se
    ci_upper = theta_hat + z_alpha * se
    ci_length = 2 * z_alpha * se

    return((theta_hat, se, ci_lower, ci_upper, ci_length,
            eps_hat, eta_nosplit, sigma_sq, epsilon_n,
            n_y, n_x, n_z))
}

/*
Inputs:
  y: outcome sample.
  x: treatment/endogenous-variable sample.
  z: instrument/exogenous-variable sample.
  level: confidence level on the R scale, e.g. 0.95.
Outputs:
  Row vector:
    1 theta_hat
    2 se
    3 ci_lower
    4 ci_upper
    5 ci_length
    6 eps_hat
    7 eta
    8 sigma_sq
    9 epsilon_n
    10 n_y
    11 n_x
    12 n_z
    13 n_half
Assumptions:
  Implements only the cross-sectional split R path. The split is deterministic:
  first n_half observations and last n_half observations are used, where
  n_half = min(floor(n_y/2), floor(n_x/2), floor(n_z/2)).
*/
real rowvector cic_split_fit(real colvector y,
                             real colvector x,
                             real colvector z,
                             real scalar level)
{
    real scalar n_y
    real scalar n_x
    real scalar n_z
    real scalar effective_n
    real scalar n_half
    real scalar epsilon_n
    real scalar theta_hat
    real scalar eps_hat
    real scalar eta_split
    real scalar sigma_sq
    real scalar sigma_sq_for_se
    real scalar se
    real scalar alpha
    real scalar z_alpha
    real scalar ci_lower
    real scalar ci_upper
    real scalar ci_length
    real scalar lbda1_3
    real scalar lbda2
    real matrix transform
    real colvector qcdf_transform
    real colvector y_first
    real colvector y_last
    real colvector x_first
    real colvector x_last
    real colvector z_first
    real colvector z_last
    real colvector uhat1
    real colvector uhat2
    real colvector ysort1diff
    real colvector ysort2diff
    real colvector fyhat_split
    real colvector fuhat1
    real colvector fuhat2

    n_y = rows(y)
    n_x = rows(x)
    n_z = rows(z)
    if (n_y < 4 | n_x < 4 | n_z < 4) {
        _error(2001)
    }
    if (level <= 0 | level >= 1) {
        _error(198)
    }

    effective_n = cic_min3(n_y, n_x, n_z)
    n_half = cic_min3(floor(n_y / 2), floor(n_x / 2), floor(n_z / 2))
    if (n_half < 2) {
        _error(2001)
    }
    epsilon_n = 1 / log(n_x)

    transform = cic_rank_quantile_transform(y, x, z)
    qcdf_transform = transform[., 2]
    theta_hat = mean(qcdf_transform)
    eps_hat = mean((qcdf_transform :- theta_hat):^2)

    y_first = cic_slice_rows(y, 1, n_half)
    y_last = cic_slice_rows(y, n_y - n_half + 1, n_y)
    x_first = cic_slice_rows(x, 1, n_half)
    x_last = cic_slice_rows(x, n_x - n_half + 1, n_x)
    z_first = cic_slice_rows(z, 1, n_half)
    z_last = cic_slice_rows(z, n_z - n_half + 1, n_z)

    uhat1 = cic_ecdf_at(z_first, x_first)
    uhat2 = cic_ecdf_at(z_last, x_last)
    ysort1diff = cic_adjacent_differences(y_first)
    ysort2diff = cic_adjacent_differences(y_last)
    fyhat_split = cic_empirical_grid(n_half)
    fuhat1 = cic_rect_density_sorted(cic_sort_ascending(uhat1), fyhat_split, epsilon_n)
    fuhat2 = cic_rect_density_sorted(cic_sort_ascending(uhat2), fyhat_split, epsilon_n)
    eta_split = cic_fast_eta(ysort1diff, fuhat1, ysort2diff, fuhat2, fyhat_split)

    lbda1_3 = effective_n * (n_y + n_z) / (n_y * n_z)
    lbda2 = effective_n / n_x
    sigma_sq = lbda1_3 * eta_split + lbda2 * eps_hat
    sigma_sq_for_se = sigma_sq
    if (sigma_sq_for_se < 0) {
        sigma_sq_for_se = 0
    }
    se = sqrt(sigma_sq_for_se / effective_n)

    alpha = (1 - level) / 2
    z_alpha = invnormal(1 - alpha)
    ci_lower = theta_hat - z_alpha * se
    ci_upper = theta_hat + z_alpha * se
    ci_length = 2 * z_alpha * se

    return((theta_hat, se, ci_lower, ci_upper, ci_length,
            eps_hat, eta_split, sigma_sq, epsilon_n,
            n_y, n_x, n_z, n_half))
}

/*
Inputs:
  y: outcome sample.
  x: treatment/endogenous-variable sample.
  z: instrument/exogenous-variable sample.
  level: confidence level on the R scale, e.g. 0.95.
Outputs:
  Row vector:
    1 theta_hat
    2 se
    3 ci_lower
    4 ci_upper
    5 ci_length
    6 eps_hat
    7 eta
    8 sigma_sq
    9 epsilon_n
    10 n_y
    11 n_x
    12 n_z
Assumptions:
  Implements only the cross-sectional KDE R path. It does not implement split,
  bootstrap, or panel data.
*/
real rowvector cic_kde_fit(real colvector y,
                           real colvector x,
                           real colvector z,
                           real scalar level)
{
    real scalar n_y
    real scalar n_x
    real scalar n_z
    real scalar effective_n
    real scalar epsilon_n
    real scalar theta_hat
    real scalar eps_hat
    real scalar eta_kde
    real scalar sigma_sq
    real scalar sigma_sq_for_se
    real scalar se
    real scalar alpha
    real scalar z_alpha
    real scalar ci_lower
    real scalar ci_upper
    real scalar ci_length
    real matrix transform
    real colvector uhat
    real colvector qcdf_transform
    real colvector bandwidths
    real colvector density_values

    n_y = rows(y)
    n_x = rows(x)
    n_z = rows(z)
    if (n_y < 4 | n_x < 4 | n_z < 4) {
        _error(2001)
    }
    if (level <= 0 | level >= 1) {
        _error(198)
    }

    effective_n = cic_min3(n_y, n_x, n_z)
    epsilon_n = 1 / log(n_x)

    transform = cic_rank_quantile_transform(y, x, z)
    uhat = transform[., 1]
    qcdf_transform = transform[., 2]
    theta_hat = mean(qcdf_transform)
    eps_hat = mean((qcdf_transform :- theta_hat):^2)

    bandwidths = epsilon_n * (uhat :* (J(rows(uhat), 1, 1) :- uhat))
    density_values = cic_epanechnikov_y_density(y, qcdf_transform, bandwidths)
    eta_kde = cic_kde_eta_from_density(density_values, uhat)
    sigma_sq = 2 * eta_kde + eps_hat
    sigma_sq_for_se = sigma_sq
    if (sigma_sq_for_se < 0) {
        sigma_sq_for_se = 0
    }
    se = sqrt(sigma_sq_for_se / effective_n)

    alpha = (1 - level) / 2
    z_alpha = invnormal(1 - alpha)
    ci_lower = theta_hat - z_alpha * se
    ci_upper = theta_hat + z_alpha * se
    ci_length = 2 * z_alpha * se

    return((theta_hat, se, ci_lower, ci_upper, ci_length,
            eps_hat, eta_kde, sigma_sq, epsilon_n,
            n_y, n_x, n_z))
}

end
