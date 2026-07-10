version 16.0

mata:
mata set matastrict on

/*
Inputs:
  sorted_values: ascending real column vector.
  value: scalar evaluation point.
Outputs:
  Insertion position before the first element >= value, matching C++ lower_bound.
Assumptions:
  sorted_values is already ascending and nonempty.
*/
real scalar cic_lower_bound_sorted(real colvector sorted_values, real scalar value)
{
    real scalar left
    real scalar right
    real scalar mid

    left = 1
    right = rows(sorted_values) + 1

    while (left < right) {
        mid = floor((left + right) / 2)
        if (sorted_values[mid] < value) {
            left = mid + 1
        }
        else {
            right = mid
        }
    }

    return(left)
}

/*
Inputs:
  grid: empirical CDF grid, typically (1:(n-1))/n.
  epsilon: scalar bandwidth multiplier.
Outputs:
  h_vals = epsilon * grid * (1 - grid), the no-split density bandwidths.
Assumptions:
  This is the pointwise = 1 bandwidth shape used by the no-split R path.
*/
real colvector cic_density_bandwidth(real colvector grid, real scalar epsilon)
{
    return(epsilon * (grid :* (J(rows(grid), 1, 1) :- grid)))
}

/*
Inputs:
  sorted_values: ascending observations used as density support.
  eval_points: density evaluation grid.
  bandwidths: interval half-widths, one per evaluation point.
Outputs:
  Counts of sorted_values in [eval_point - h, eval_point + h].
Assumptions:
  Matches rect_counts_rcpp(): lower_bound on the left endpoint and
  upper_bound on the right endpoint, so both exact endpoints are included.
*/
real colvector cic_rect_counts_sorted(real colvector sorted_values,
                                      real colvector eval_points,
                                      real colvector bandwidths)
{
    real scalar i
    real scalar lo
    real scalar hi
    real scalar n_eval
    real colvector counts

    n_eval = rows(eval_points)
    if (rows(bandwidths) != n_eval) {
        _error(3200)
    }

    counts = J(n_eval, 1, .)
    for (i = 1; i <= n_eval; i++) {
        if (eval_points[i] >= . | bandwidths[i] >= .) {
            counts[i] = .
        }
        else {
            lo = cic_lower_bound_sorted(sorted_values, eval_points[i] - bandwidths[i])
            hi = cic_upper_bound_sorted(sorted_values, eval_points[i] + bandwidths[i])
            counts[i] = hi - lo
        }
    }

    return(counts)
}

/*
Inputs:
  counts: interval counts.
  bandwidths: interval half-widths.
  sample_size: size of the density support sample.
Outputs:
  Rectangle density estimates counts / (2 * sample_size * bandwidth).
Assumptions:
  Mirrors counts_to_density(); caller supplies positive bandwidths for active grid.
*/
real colvector cic_counts_to_density(real colvector counts,
                                     real colvector bandwidths,
                                     real scalar sample_size)
{
    if (rows(counts) != rows(bandwidths)) {
        _error(3200)
    }
    if (sample_size <= 0) {
        _error(3498)
    }

    return(counts :/ (2 * sample_size * bandwidths))
}

/*
Inputs:
  sorted_values: ascending observations used as density support.
  eval_points: density evaluation grid.
  bandwidths: interval half-widths.
Outputs:
  Rectangle density estimates from interval counts and bandwidths.
Assumptions:
  This composes cic_rect_counts_sorted() and cic_counts_to_density().
*/
real colvector cic_rect_density_from_bandwidths(real colvector sorted_values,
                                                real colvector eval_points,
                                                real colvector bandwidths)
{
    real colvector counts

    counts = cic_rect_counts_sorted(sorted_values, eval_points, bandwidths)
    return(cic_counts_to_density(counts, bandwidths, rows(sorted_values)))
}

/*
Inputs:
  sorted_values: ascending observations used as density support.
  eval_points: empirical CDF grid.
  epsilon: scalar bandwidth multiplier.
Outputs:
  No-split rectangle density estimate evaluated on eval_points.
Assumptions:
  Uses the R no-split bandwidth h = epsilon * F * (1 - F).
*/
real colvector cic_rect_density_sorted(real colvector sorted_values,
                                       real colvector eval_points,
                                       real scalar epsilon)
{
    real colvector bandwidths

    bandwidths = cic_density_bandwidth(eval_points, epsilon)
    return(cic_rect_density_from_bandwidths(sorted_values, eval_points, bandwidths))
}

/*
Inputs:
  y: outcome sample.
  eval_points: points where the density is evaluated.
  bandwidths: bandwidths, one per evaluation point.
Outputs:
  Epanechnikov density estimates matching f_y_hat_epnechikov().
Assumptions:
  The implementation centers y before sorting and uses prefix sums for the
  first two centered moments, following the current Rcpp code.
*/
real colvector cic_epanechnikov_y_density(real colvector y,
                                          real colvector eval_points,
                                          real colvector bandwidths)
{
    real scalar n
    real scalar m
    real scalar i
    real scalar lo
    real scalar hi
    real scalar center
    real scalar sqrt5
    real scalar coef
    real scalar bandwidth
    real scalar window
    real scalar inv_h
    real scalar inv_5h2
    real scalar eval_centered
    real scalar s0
    real scalar s1
    real scalar s2
    real scalar quadratic_sum
    real scalar kernel_sum
    real colvector y_centered_sorted
    real colvector cs1
    real colvector cs2
    real colvector density

    n = rows(y)
    m = rows(eval_points)
    if (rows(bandwidths) != m) {
        _error(3200)
    }

    center = mean(y)
    y_centered_sorted = cic_sort_ascending(y :- center)
    cs1 = J(n + 1, 1, 0)
    cs2 = J(n + 1, 1, 0)
    for (i = 1; i <= n; i++) {
        cs1[i + 1] = cs1[i] + y_centered_sorted[i]
        cs2[i + 1] = cs2[i] + y_centered_sorted[i]^2
    }

    sqrt5 = sqrt(5)
    coef = 3 / (4 * sqrt5)
    density = J(m, 1, .)

    for (i = 1; i <= m; i++) {
        bandwidth = bandwidths[i]
        if (bandwidth <= 0) {
            density[i] = 0
        }
        else {
            window = bandwidth * sqrt5
            inv_h = 1 / bandwidth
            inv_5h2 = 1 / (5 * bandwidth^2)
            eval_centered = eval_points[i] - center
            lo = cic_lower_bound_sorted(y_centered_sorted, eval_centered - window)
            hi = cic_upper_bound_sorted(y_centered_sorted, eval_centered + window)
            s0 = hi - lo
            s1 = cs1[hi] - cs1[lo]
            s2 = cs2[hi] - cs2[lo]
            quadratic_sum = s2 - 2 * eval_centered * s1 + eval_centered^2 * s0
            kernel_sum = s0 - quadratic_sum * inv_5h2
            density[i] = kernel_sum * coef * inv_h / n
        }
    }

    return(density)
}

end
