version 16.0

mata:
mata set matastrict on

/*
Inputs:
  values: real column vector.
Outputs:
  Reverse cumulative sum: output[i] = sum(values[i..n]).
Assumptions:
  values is nonempty.
*/
real colvector cic_reverse_cumsum(real colvector values)
{
    real scalar i
    real scalar running_sum
    real colvector cumulative_values

    cumulative_values = J(rows(values), 1, .)
    running_sum = 0

    for (i = rows(values); i >= 1; i--) {
        running_sum = running_sum + values[i]
        cumulative_values[i] = running_sum
    }

    return(cumulative_values)
}

/*
Inputs:
  grid: increasing empirical CDF grid.
Outputs:
  Grid increments c(grid[1], diff(grid)).
Assumptions:
  grid is nonempty.
*/
real colvector cic_grid_increments(real colvector grid)
{
    real scalar i
    real colvector increments

    increments = J(rows(grid), 1, .)
    increments[1] = grid[1]
    for (i = 2; i <= rows(grid); i++) {
        increments[i] = grid[i] - grid[i - 1]
    }

    return(increments)
}

/*
Inputs:
  y_diff1, y_diff2: adjacent sorted-Y differences.
  density1, density2: density estimates on the same grid.
  grid: empirical CDF grid, usually (1:(n-1))/n.
Outputs:
  Fast eta variance component used by the R .fast_eta() helper.
Assumptions:
  Inputs have conformable lengths and correspond to the same empirical grid.
*/
real scalar cic_fast_eta(real colvector y_diff1,
                         real colvector density1,
                         real colvector y_diff2,
                         real colvector density2,
                         real colvector grid)
{
    real colvector u
    real colvector v
    real colvector reverse_u
    real colvector reverse_v
    real colvector delta_v
    real scalar term2

    if (rows(y_diff1) != rows(density1) |
        rows(y_diff2) != rows(density2) |
        rows(y_diff1) != rows(grid) |
        rows(y_diff2) != rows(grid)) {
        _error(3200)
    }

    u = y_diff1 :* density1
    v = y_diff2 :* density2
    term2 = sum(u :* grid) * sum(v :* grid)
    reverse_u = cic_reverse_cumsum(u)
    reverse_v = cic_reverse_cumsum(v)
    delta_v = cic_grid_increments(grid)

    return(sum(delta_v :* reverse_u :* reverse_v) - term2)
}

/*
Inputs:
  y_diff: adjacent sorted-Y differences.
  density: density estimates on the same grid.
  grid: empirical CDF grid.
Outputs:
  Same-sample no-split eta component, .fast_eta(y_diff, f, y_diff, f, grid).
Assumptions:
  This is a convenience wrapper only; it does not compute density estimates.
*/
real scalar cic_fast_eta_self(real colvector y_diff,
                              real colvector density,
                              real colvector grid)
{
    return(cic_fast_eta(y_diff, density, y_diff, density, grid))
}

/*
Inputs:
  sorted_values: ascending support vector.
  eval_points: points at which findInterval(eval_point, sorted_values)+1 is needed.
Outputs:
  k indices matching R findInterval(eval_points, sorted_values) + 1.
Assumptions:
  sorted_values is ascending; default R findInterval behavior is count(x_i <= t).
*/
real colvector cic_find_interval_plus_one(real colvector sorted_values,
                                          real colvector eval_points)
{
    real scalar i
    real colvector indices

    indices = J(rows(eval_points), 1, .)
    for (i = 1; i <= rows(eval_points); i++) {
        indices[i] = cic_upper_bound_sorted(sorted_values, eval_points[i])
    }

    return(indices)
}

/*
Inputs:
  density_values: density values f_i.
  uhat: unsorted rank transform values.
Outputs:
  KDE eta component matching .compute_eta_from_f().
Assumptions:
  Implements the current R code path with grid = (1:n)/n and
  findInterval(grid - 1e-12, sort(uhat)) + 1.
*/
real scalar cic_kde_eta_from_density(real colvector density_values,
                                     real colvector uhat)
{
    real scalar n
    real scalar i
    real scalar c2
    real colvector order_index
    real colvector sorted_uhat
    real colvector grid
    real colvector k
    real colvector inv_density_sorted
    real colvector reverse_inv_density
    real colvector t1

    n = rows(uhat)
    if (rows(density_values) != n) {
        _error(3200)
    }

    order_index = cic_order_ascending(uhat)
    sorted_uhat = uhat[order_index]
    grid = (1::n) / n
    k = cic_find_interval_plus_one(sorted_uhat, grid :- 1e-12)
    c2 = mean(uhat :/ density_values)
    inv_density_sorted = 1 :/ density_values[order_index]
    reverse_inv_density = cic_reverse_cumsum(inv_density_sorted) / n
    t1 = J(n, 1, 0)

    for (i = 1; i <= n; i++) {
        if (k[i] <= n) {
            t1[i] = reverse_inv_density[k[i]]
        }
    }

    return(mean((t1 :- c2):^2))
}

end
