#!/usr/bin/env Rscript

# Phase 2 reference-output generator for the Stata migration.
# It loads the local R package, creates deterministic datasets, runs the R
# reference estimators, and writes CSV files that Stata can import later.

options(stringsAsFactors = FALSE, warn = 1)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

find_repo_root <- function(start_dir) {
  current <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        dir.exists(file.path(current, "R"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find repository root containing DESCRIPTION and R/.", call. = FALSE)
    }
    current <- parent
  }
}

script_path <- get_script_path()
script_dir <- if (dir.exists(script_path)) script_path else dirname(script_path)
repo_root <- find_repo_root(script_dir)
reference_dir <- file.path(repo_root, "stata", "tests", "reference_outputs")
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "", quote = TRUE)
}

load_local_cicinference <- function(root) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    loaded <- tryCatch({
      suppressPackageStartupMessages(
        pkgload::load_all(
          path = root,
          reset = TRUE,
          recompile = TRUE,
          export_all = FALSE,
          helpers = FALSE,
          attach_testthat = FALSE,
          quiet = TRUE
        )
      )
      TRUE
    }, error = function(e) {
      message("pkgload::load_all() failed; falling back to sourceCpp/sys.source: ",
              conditionMessage(e))
      FALSE
    })
    if (loaded) {
      return("pkgload::load_all")
    }
  }

  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Rcpp is required to compile local C++ reference functions.", call. = FALSE)
  }

  Rcpp::sourceCpp(
    file.path(root, "src", "Density_cic_boostrap_functions.cpp"),
    rebuild = TRUE,
    verbose = FALSE
  )

  r_sources <- file.path(
    root,
    "R",
    c("utils.R", "data-generating.R", "cic.R", "print_summary.R")
  )
  for (source_file in r_sources) {
    sys.source(source_file, envir = .GlobalEnv)
  }
  "Rcpp::sourceCpp + sys.source"
}

load_method <- load_local_cicinference(repo_root)

get_reference_function <- function(name) {
  if (exists(name, envir = .GlobalEnv, mode = "function", inherits = FALSE)) {
    return(get(name, envir = .GlobalEnv, mode = "function"))
  }
  namespace_fun <- tryCatch(
    getFromNamespace(name, ns = "cicinference"),
    error = function(e) NULL
  )
  if (is.function(namespace_fun)) {
    return(namespace_fun)
  }
  stop(sprintf("Could not find required reference function `%s`.", name), call. = FALSE)
}

rect_counts_reference <- get_reference_function("rect_counts_rcpp")
counts_to_density_reference <- get_reference_function("counts_to_density")
fast_eta_reference <- get_reference_function(".fast_eta")
compute_eta_from_f_reference <- get_reference_function(".compute_eta_from_f")
f_y_hat_epanechnikov_reference <- get_reference_function("f_y_hat_epnechikov")

set_reference_seed <- function(seed) {
  set.seed(
    seed,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
}

clamp_to_open_range <- function(x, reference) {
  finite_reference <- reference[is.finite(reference)]
  if (length(finite_reference) < 2L) {
    stop("Need at least two finite reference values for clamping.", call. = FALSE)
  }
  reference_range <- range(finite_reference)
  span <- diff(reference_range)
  tolerance <- if (span > 0) span * 1e-10 else 1e-10
  pmin(pmax(x, reference_range[[1L]] + tolerance), reference_range[[2L]] - tolerance)
}

make_sim_dataset <- function(n, seed, b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05) {
  set_reference_seed(seed)
  data <- sim_dgp(n = n, b1 = b1, b2 = b2, d1 = d1, d2 = d2, seed = seed)
  data$X <- clamp_to_open_range(data$X, data$Z)
  data
}

make_ties_dataset <- function(n, seed) {
  data <- make_sim_dataset(n = n, seed = seed, b1 = 0.05, b2 = 0.08, d1 = 0.02, d2 = 0.05)
  data$Y <- round(data$Y, 1)
  data$Z <- round(data$Z, 1)
  data$X <- clamp_to_open_range(round(data$X, 1), data$Z)
  data
}

make_extreme_dataset <- function(n, seed) {
  set_reference_seed(seed)
  inner_probs <- seq(0.02, 0.98, length.out = n - 6L)
  y_probs <- c(1e-6, 1e-4, 1e-3, inner_probs, 0.999, 0.9999, 1 - 1e-6)
  y_values <- qY_dgp(y_probs, d1 = 0.12, d2 = 0.12)

  z_probs <- seq(0.01, 0.99, length.out = n)
  z_values <- stats::qnorm(z_probs) + stats::rnorm(n, mean = 0, sd = 0.01)

  x_probs <- seq(0.02, 0.98, length.out = n)
  x_values <- stats::qnorm(x_probs)

  list(
    Y = y_values[sample.int(n)],
    X = clamp_to_open_range(x_values[sample.int(n)], z_values),
    Z = z_values[sample.int(n)]
  )
}

make_panel_dataset <- function(n_yz, n_x, seed) {
  set_reference_seed(seed)
  panel_data <- sim_dgp(n = n_yz, seed = seed, panel_data = TRUE)
  x_data <- sim_dgp(n = n_x, seed = seed + 1000L)
  list(
    Y = panel_data$Y,
    X = clamp_to_open_range(x_data$X, panel_data$Z),
    Z = panel_data$Z
  )
}

make_missing_probe_dataset <- function(seed) {
  data <- make_sim_dataset(n = 64L, seed = seed)
  data$Y[[4L]] <- NA_real_
  data$X[[7L]] <- NA_real_
  data$Z[[9L]] <- NA_real_
  data
}

datasets <- list(
  list(
    dataset_id = "xs_small",
    kind = "cross_section",
    seed = 202601L,
    description = "Small simulated cross-sectional dataset with X clamped inside Z range.",
    data = make_sim_dataset(n = 64L, seed = 202601L)
  ),
  list(
    dataset_id = "xs_medium",
    kind = "cross_section",
    seed = 202602L,
    description = "Medium simulated cross-sectional dataset with non-default DGP parameters.",
    data = make_sim_dataset(n = 320L, seed = 202602L, b1 = 0.06, b2 = 0.08, d1 = 0.02, d2 = 0.04)
  ),
  list(
    dataset_id = "xs_edge_ties",
    kind = "cross_section_edge",
    seed = 202603L,
    description = "Edge dataset with many ties induced by rounding Y, X, and Z.",
    data = make_ties_dataset(n = 96L, seed = 202603L)
  ),
  list(
    dataset_id = "xs_edge_extreme",
    kind = "cross_section_edge",
    seed = 202604L,
    description = "Edge dataset with extreme Y quantile values and interior X ranks.",
    data = make_extreme_dataset(n = 96L, seed = 202604L)
  ),
  list(
    dataset_id = "panel_small",
    kind = "panel",
    seed = 202605L,
    description = "Small paired Y/Z panel dataset with same X sample size.",
    data = make_panel_dataset(n_yz = 64L, n_x = 64L, seed = 202605L)
  ),
  list(
    dataset_id = "panel_medium_unbalanced",
    kind = "panel",
    seed = 202606L,
    description = "Panel dataset with paired Y/Z length different from X length.",
    data = make_panel_dataset(n_yz = 160L, n_x = 120L, seed = 202606L)
  ),
  list(
    dataset_id = "xs_missing_probe",
    kind = "edge_probe",
    seed = 202607L,
    description = "Missing-value probe; behavior is recorded but not a required fit.",
    data = make_missing_probe_dataset(seed = 202607L)
  )
)

dataset_by_id <- stats::setNames(lapply(datasets, `[[`, "data"),
                                vapply(datasets, `[[`, character(1L), "dataset_id"))

input_rows <- do.call(rbind, lapply(datasets, function(dataset) {
  values <- dataset$data
  do.call(rbind, lapply(names(values), function(vector_name) {
    data.frame(
      dataset_id = dataset$dataset_id,
      kind = dataset$kind,
      vector = vector_name,
      index = seq_along(values[[vector_name]]),
      value = as.numeric(values[[vector_name]])
    )
  }))
}))

dataset_manifest <- do.call(rbind, lapply(datasets, function(dataset) {
  data.frame(
    dataset_id = dataset$dataset_id,
    kind = dataset$kind,
    seed = dataset$seed,
    n_y = length(dataset$data$Y),
    n_x = length(dataset$data$X),
    n_z = length(dataset$data$Z),
    has_missing = any(is.na(unlist(dataset$data, use.names = FALSE))),
    description = dataset$description
  )
}))

make_utility_component <- function(dataset, component, values) {
  data.frame(
    dataset_id = dataset$dataset_id,
    kind = dataset$kind,
    component = component,
    index = seq_along(values),
    value = as.numeric(values)
  )
}

make_nosplit_utility_references <- function(dataset) {
  data <- dataset$data
  if (anyNA(data$Y) || anyNA(data$X) || anyNA(data$Z)) {
    return(NULL)
  }

  sorted_y <- sort(data$Y, method = "radix")
  sorted_x <- sort(data$X, method = "radix")
  sorted_z <- sort(data$Z, method = "radix")
  uhat <- stats::ecdf(data$Z)(data$X)
  left_quantile_index <- pmax(ceiling(uhat * length(sorted_y)), 1L)
  qcdf_transform <- sorted_y[left_quantile_index]
  theta_from_transform <- mean(qcdf_transform)
  eps_hat <- mean((theta_from_transform - qcdf_transform)^2)
  fyhat <- seq_len(length(sorted_y) - 1L) / length(sorted_y)
  ysortdiff <- diff(sorted_y)
  epsilon_n <- 1 / log(length(data$X))
  sorted_uhat <- sort(uhat, method = "radix")
  h_vals <- epsilon_n * (fyhat * (1 - fyhat))
  rect_counts <- rect_counts_reference(sorted_uhat, fyhat, h_vals)
  fuhat <- counts_to_density_reference(rect_counts, h_vals, length(sorted_uhat))
  eta_u <- as.numeric(ysortdiff * fuhat)
  eta_reverse_cumsum <- rev(cumsum(rev(eta_u)))
  eta_delta_v <- c(fyhat[1], diff(fyhat))
  eta_term2 <- sum(eta_u * fyhat) * sum(eta_u * fyhat)
  eta_nosplit <- fast_eta_reference(ysortdiff, fuhat, ysortdiff, fuhat, fyhat)
  N <- min(length(data$Y), length(data$X), length(data$Z))
  lbda1_3 <- N * (length(data$Y) + length(data$Z)) / (length(data$Y) * length(data$Z))
  lbda2 <- N / length(data$X)
  sigma_sq_nosplit <- lbda1_3 * eta_nosplit + lbda2 * eps_hat
  se_nosplit <- sqrt(sigma_sq_nosplit / N)

  do.call(rbind, list(
    make_utility_component(dataset, "sorted_y", sorted_y),
    make_utility_component(dataset, "sorted_x", sorted_x),
    make_utility_component(dataset, "sorted_z", sorted_z),
    make_utility_component(dataset, "sorted_uhat", sorted_uhat),
    make_utility_component(dataset, "uhat_ecdf_z_at_x", uhat),
    make_utility_component(dataset, "left_quantile_index", left_quantile_index),
    make_utility_component(dataset, "qcdf_transform", qcdf_transform),
    make_utility_component(dataset, "theta_hat_from_transform", theta_from_transform),
    make_utility_component(dataset, "eps_hat", eps_hat),
    make_utility_component(dataset, "fyhat_grid", fyhat),
    make_utility_component(dataset, "ysortdiff", ysortdiff),
    make_utility_component(dataset, "epsilon_n_default", epsilon_n),
    make_utility_component(dataset, "density_bandwidth", h_vals),
    make_utility_component(dataset, "rect_counts", rect_counts),
    make_utility_component(dataset, "fuhat_density", fuhat),
    make_utility_component(dataset, "fast_eta_u", eta_u),
    make_utility_component(dataset, "fast_eta_reverse_cumsum", eta_reverse_cumsum),
    make_utility_component(dataset, "fast_eta_delta_v", eta_delta_v),
    make_utility_component(dataset, "fast_eta_term2", eta_term2),
    make_utility_component(dataset, "fast_eta_self", eta_nosplit),
    make_utility_component(dataset, "sigma_sq_nosplit", sigma_sq_nosplit),
    make_utility_component(dataset, "se_nosplit", se_nosplit)
  ))
}

make_split_utility_references <- function(dataset) {
  data <- dataset$data
  if (anyNA(data$Y) || anyNA(data$X) || anyNA(data$Z)) {
    return(NULL)
  }

  n1 <- length(data$Y)
  n2 <- length(data$X)
  n3 <- length(data$Z)
  n_half <- min(floor(n1 / 2), floor(n2 / 2), floor(n3 / 2))
  if (n_half < 2L) {
    return(NULL)
  }

  epsilon_n <- 1 / log(n2)
  N <- min(n1, n2, n3)
  lbda1_3 <- N * (n1 + n3) / (n1 * n3)
  lbda2 <- N / n2

  uhat_full <- stats::ecdf(data$Z)(data$X)
  qcdf_transform <- sort(data$Y, method = "radix")[pmax(ceiling(uhat_full * n1), 1L)]
  theta_hat <- mean(qcdf_transform)
  eps_hat <- mean((theta_hat - qcdf_transform)^2)

  z_first <- data$Z[seq_len(n_half)]
  z_last <- data$Z[seq.int(n3 - n_half + 1L, n3)]
  x_first <- data$X[seq_len(n_half)]
  x_last <- data$X[seq.int(n2 - n_half + 1L, n2)]
  y_first <- data$Y[seq_len(n_half)]
  y_last <- data$Y[seq.int(n1 - n_half + 1L, n1)]

  uhat1 <- stats::ecdf(z_first)(x_first)
  uhat2 <- stats::ecdf(z_last)(x_last)
  sorted_uhat1 <- sort(uhat1, method = "radix")
  sorted_uhat2 <- sort(uhat2, method = "radix")
  ysort1diff <- diff(sort(y_first, method = "radix"))
  ysort2diff <- diff(sort(y_last, method = "radix"))
  fyhat_split <- seq_len(n_half - 1L) / n_half
  h_vals <- epsilon_n * (fyhat_split * (1 - fyhat_split))
  counts1 <- rect_counts_reference(sorted_uhat1, fyhat_split, h_vals)
  counts2 <- rect_counts_reference(sorted_uhat2, fyhat_split, h_vals)
  fuhat1 <- counts_to_density_reference(counts1, h_vals, length(sorted_uhat1))
  fuhat2 <- counts_to_density_reference(counts2, h_vals, length(sorted_uhat2))
  eta_split <- fast_eta_reference(ysort1diff, fuhat1, ysort2diff, fuhat2, fyhat_split)
  sigma_sq_split <- lbda1_3 * eta_split + lbda2 * eps_hat
  se_split <- sqrt(max(sigma_sq_split, 0) / N)

  do.call(rbind, list(
    make_utility_component(dataset, "split_n_half", n_half),
    make_utility_component(dataset, "split_x_first", x_first),
    make_utility_component(dataset, "split_x_last", x_last),
    make_utility_component(dataset, "split_z_first", z_first),
    make_utility_component(dataset, "split_z_last", z_last),
    make_utility_component(dataset, "split_y_first", y_first),
    make_utility_component(dataset, "split_y_last", y_last),
    make_utility_component(dataset, "split_uhat1", uhat1),
    make_utility_component(dataset, "split_uhat2", uhat2),
    make_utility_component(dataset, "split_sorted_uhat1", sorted_uhat1),
    make_utility_component(dataset, "split_sorted_uhat2", sorted_uhat2),
    make_utility_component(dataset, "split_ysort1diff", ysort1diff),
    make_utility_component(dataset, "split_ysort2diff", ysort2diff),
    make_utility_component(dataset, "split_fyhat_grid", fyhat_split),
    make_utility_component(dataset, "split_density_bandwidth", h_vals),
    make_utility_component(dataset, "split_rect_counts1", counts1),
    make_utility_component(dataset, "split_rect_counts2", counts2),
    make_utility_component(dataset, "split_fuhat1", fuhat1),
    make_utility_component(dataset, "split_fuhat2", fuhat2),
    make_utility_component(dataset, "split_eta", eta_split),
    make_utility_component(dataset, "split_sigma_sq", sigma_sq_split),
    make_utility_component(dataset, "split_se", se_split)
  ))
}

make_kde_utility_references <- function(dataset) {
  data <- dataset$data
  if (anyNA(data$Y) || anyNA(data$X) || anyNA(data$Z)) {
    return(NULL)
  }

  n1 <- length(data$Y)
  n2 <- length(data$X)
  n3 <- length(data$Z)
  N <- min(n1, n2, n3)
  epsilon_n <- 1 / log(n2)

  uhat <- stats::ecdf(data$Z)(data$X)
  sorted_y <- sort(data$Y, method = "radix")
  qcdf_transform <- sorted_y[pmax(ceiling(uhat * n1), 1L)]
  theta_hat <- mean(qcdf_transform)
  eps_hat <- mean((theta_hat - qcdf_transform)^2)
  h <- epsilon_n * uhat * (1 - uhat)

  idx_sort <- order(uhat)
  u_sort <- uhat[idx_sort]
  grid <- seq_len(length(uhat)) / length(uhat)
  k <- findInterval(grid - 1e-12, u_sort) + 1L
  ok <- as.integer(k <= length(uhat))
  f_one <- f_y_hat_epanechnikov_reference(data$Y, qcdf_transform, h)
  eta_kde <- compute_eta_from_f_reference(f_one, uhat, idx_sort, k, ok == 1L, length(uhat))
  sigma_sq_kde <- 2 * eta_kde + eps_hat
  se_kde <- sqrt(max(sigma_sq_kde, 0) / N)

  do.call(rbind, list(
    make_utility_component(dataset, "kde_h", h),
    make_utility_component(dataset, "kde_idx_sort", idx_sort),
    make_utility_component(dataset, "kde_u_sort", u_sort),
    make_utility_component(dataset, "kde_grid", grid),
    make_utility_component(dataset, "kde_k", k),
    make_utility_component(dataset, "kde_ok", ok),
    make_utility_component(dataset, "kde_f_one", f_one),
    make_utility_component(dataset, "kde_eta", eta_kde),
    make_utility_component(dataset, "kde_sigma_sq", sigma_sq_kde),
    make_utility_component(dataset, "kde_se", se_kde)
  ))
}

utility_reference_datasets <- Filter(
  function(dataset) dataset$kind %in% c("cross_section", "cross_section_edge"),
  datasets
)
utility_nosplit_references <- do.call(
  rbind,
  c(
    lapply(utility_reference_datasets, make_nosplit_utility_references),
    lapply(utility_reference_datasets, make_split_utility_references),
    lapply(utility_reference_datasets, make_kde_utility_references)
  )
)

make_bootstrap_index_rows <- function(dataset, role, index_matrix) {
  do.call(rbind, lapply(seq_len(nrow(index_matrix)), function(replicate_id) {
    data.frame(
      dataset_id = dataset$dataset_id,
      kind = dataset$kind,
      role = role,
      replicate = replicate_id,
      draw = seq_len(ncol(index_matrix)),
      index = as.integer(index_matrix[replicate_id, ])
    )
  }))
}

bootstrap_values_from_indices <- function(data, y_index, x_index, z_index) {
  sorted_y <- sort(data$Y, method = "radix")
  sorted_z <- sort(data$Z, method = "radix")
  x_values <- data$X
  n1 <- length(sorted_y)
  n2 <- length(x_values)
  n3 <- length(sorted_z)

  vapply(seq_len(nrow(y_index)), function(replicate_id) {
    counts_y <- tabulate(y_index[replicate_id, ], nbins = n1)
    counts_z <- tabulate(z_index[replicate_id, ], nbins = n3)
    cdf_y_counts <- cumsum(counts_y)
    cdf_z <- cumsum(counts_z) / n3

    transformed <- vapply(x_index[replicate_id, ], function(x_draw) {
      xb <- x_values[[x_draw]]
      pos <- findInterval(xb, sorted_z)
      u <- if (pos == 0L) 0 else cdf_z[[pos]]
      rank_y <- ceiling(u * n1)
      rank_y <- min(max(rank_y, 1L), n1)
      sorted_y[[which(cdf_y_counts >= rank_y)[[1L]]]]
    }, numeric(1L))

    mean(transformed)
  }, numeric(1L))
}

make_bootstrap_common_references <- function(dataset, B, seed) {
  data <- dataset$data
  if (anyNA(data$Y) || anyNA(data$X) || anyNA(data$Z)) {
    return(NULL)
  }

  n1 <- length(data$Y)
  n2 <- length(data$X)
  n3 <- length(data$Z)
  set_reference_seed(seed)

  y_index <- t(replicate(B, sample.int(n1, n1, replace = TRUE)))
  x_index <- t(replicate(B, sample.int(n2, n2, replace = TRUE)))
  z_index <- t(replicate(B, sample.int(n3, n3, replace = TRUE)))
  boot_theta <- bootstrap_values_from_indices(data, y_index, x_index, z_index)

  uhat <- stats::ecdf(data$Z)(data$X)
  qcdf_transform <- sort(data$Y, method = "radix")[pmax(ceiling(uhat * n1), 1L)]
  theta_hat <- mean(qcdf_transform)
  alpha <- (0.05 / 2)
  z_value <- stats::qnorm(1 - alpha)
  se_boot <- stats::sd(boot_theta)
  bpc_lower <- stats::quantile(boot_theta, probs = alpha, names = FALSE, type = 7)
  bpc_upper <- stats::quantile(boot_theta, probs = 1 - alpha, names = FALSE, type = 7)

  list(
    index_rows = do.call(rbind, list(
      make_bootstrap_index_rows(dataset, "Y", y_index),
      make_bootstrap_index_rows(dataset, "X", x_index),
      make_bootstrap_index_rows(dataset, "Z", z_index)
    )),
    draw_rows = data.frame(
      dataset_id = dataset$dataset_id,
      kind = dataset$kind,
      replicate = seq_len(B),
      theta_boot = as.numeric(boot_theta)
    ),
    summary_rows = data.frame(
      dataset_id = dataset$dataset_id,
      kind = dataset$kind,
      B = B,
      seed = seed,
      level = 0.95,
      theta_hat = theta_hat,
      se_boot = se_boot,
      bse_lower = theta_hat - z_value * se_boot,
      bse_upper = theta_hat + z_value * se_boot,
      bse_length = 2 * z_value * se_boot,
      bpc_lower = bpc_lower,
      bpc_upper = bpc_upper,
      bpc_length = bpc_upper - bpc_lower,
      n_y = n1,
      n_x = n2,
      n_z = n3
    )
  )
}

bootstrap_reference_specs <- list(
  list(dataset_id = "xs_small", B = 40L, seed = 932601L),
  list(dataset_id = "xs_medium", B = 40L, seed = 932602L),
  list(dataset_id = "xs_edge_ties", B = 40L, seed = 932603L),
  list(dataset_id = "xs_edge_extreme", B = 40L, seed = 932604L)
)

bootstrap_reference_results <- lapply(bootstrap_reference_specs, function(spec) {
  dataset <- datasets[[match(spec$dataset_id, vapply(datasets, `[[`, character(1L), "dataset_id"))]]
  make_bootstrap_common_references(dataset, B = spec$B, seed = spec$seed)
})
bootstrap_index_references <- do.call(rbind, lapply(bootstrap_reference_results, `[[`, "index_rows"))
bootstrap_draw_references <- do.call(rbind, lapply(bootstrap_reference_results, `[[`, "draw_rows"))
bootstrap_summary_references <- do.call(rbind, lapply(bootstrap_reference_results, `[[`, "summary_rows"))

cross_section_methods <- c("no-split", "split", "kde", "bse", "bpc")
panel_methods <- c("no-split", "split")

fit_specs <- list(
  list(
    call_id = "xs_small_all_methods",
    dataset_id = "xs_small",
    panel_data = FALSE,
    methods = cross_section_methods,
    B = 200L,
    level = 0.95,
    fit_seed = 902601L,
    required = TRUE
  ),
  list(
    call_id = "xs_medium_all_methods",
    dataset_id = "xs_medium",
    panel_data = FALSE,
    methods = cross_section_methods,
    B = 200L,
    level = 0.95,
    fit_seed = 902602L,
    required = TRUE
  ),
  list(
    call_id = "xs_edge_ties_all_methods",
    dataset_id = "xs_edge_ties",
    panel_data = FALSE,
    methods = cross_section_methods,
    B = 200L,
    level = 0.95,
    fit_seed = 902603L,
    required = TRUE
  ),
  list(
    call_id = "xs_edge_extreme_all_methods",
    dataset_id = "xs_edge_extreme",
    panel_data = FALSE,
    methods = cross_section_methods,
    B = 200L,
    level = 0.95,
    fit_seed = 902604L,
    required = TRUE
  ),
  list(
    call_id = "panel_small_methods",
    dataset_id = "panel_small",
    panel_data = TRUE,
    methods = panel_methods,
    B = NA_integer_,
    level = 0.95,
    fit_seed = 902605L,
    required = TRUE
  ),
  list(
    call_id = "panel_medium_unbalanced_methods",
    dataset_id = "panel_medium_unbalanced",
    panel_data = TRUE,
    methods = panel_methods,
    B = NA_integer_,
    level = 0.95,
    fit_seed = 902606L,
    required = TRUE
  ),
  list(
    call_id = "xs_small_bootstrap_B50_probe",
    dataset_id = "xs_small",
    panel_data = FALSE,
    methods = c("bse", "bpc"),
    B = 50L,
    level = 0.95,
    fit_seed = 902607L,
    required = FALSE
  ),
  list(
    call_id = "xs_missing_probe_nosplit",
    dataset_id = "xs_missing_probe",
    panel_data = FALSE,
    methods = c("no-split"),
    B = NA_integer_,
    level = 0.95,
    fit_seed = 902608L,
    required = FALSE
  )
)

fit_spec_manifest <- do.call(rbind, lapply(fit_specs, function(spec) {
  output <- data.frame(
    call_id = spec$call_id,
    dataset_id = spec$dataset_id,
    panel_data = spec$panel_data,
    requested_methods = paste(spec$methods, collapse = ";"),
    requested_B = if (is.na(spec$B)) NA_integer_ else spec$B,
    effective_B = if (any(spec$methods %in% c("bse", "bpc"))) max(spec$B, 200L) else NA_integer_,
    level = spec$level,
    fit_seed = spec$fit_seed,
    required = spec$required
  )
}))

make_output_rows <- function(fit, spec, warnings) {
  ci <- fit$ci
  alpha <- (1 - fit$level) / 2
  z_value <- stats::qnorm(1 - alpha)
  se <- ifelse(ci$method == "bpc", NA_real_, ci$length / (2 * z_value))
  t_value <- ifelse(is.na(se) | se <= 0, NA_real_, fit$theta_hat / se)
  p_value <- ifelse(is.na(t_value), NA_real_,
                    2 * stats::pnorm(abs(t_value), lower.tail = FALSE))
  confint_matrix <- confint(fit)

  output <- data.frame(
    call_id = spec$call_id,
    dataset_id = spec$dataset_id,
    panel_data = spec$panel_data,
    requested_methods = paste(spec$methods, collapse = ";"),
    requested_B = if (is.na(spec$B)) NA_integer_ else spec$B,
    effective_B = if (any(spec$methods %in% c("bse", "bpc"))) max(spec$B, 200L) else NA_integer_,
    level = fit$level,
    fit_seed = spec$fit_seed,
    method = ci$method,
    theta_hat = fit$theta_hat,
    ci_lower = ci$lower,
    ci_upper = ci$upper,
    ci_length = ci$length,
    se_implied = se,
    t_value = t_value,
    p_value = p_value,
    confint_lower = confint_matrix[ci$method, "lower"],
    confint_upper = confint_matrix[ci$method, "upper"],
    coef_theta = as.numeric(coef(fit)),
    epsilon_n = fit$epsilon_n,
    n_y = unname(fit$n[["n1"]]),
    n_x = unname(fit$n[["n2"]]),
    n_z = unname(fit$n[["n3"]]),
    warning_messages = paste(unique(warnings), collapse = " | "),
    required = spec$required
  )
  output$finite_ci <- is.finite(output$theta_hat) &
    is.finite(output$ci_lower) &
    is.finite(output$ci_upper) &
    is.finite(output$ci_length)
  output
}

run_reference_fit <- function(spec) {
  data <- dataset_by_id[[spec$dataset_id]]
  warnings <- character()
  set_reference_seed(spec$fit_seed)

  fit <- tryCatch(
    withCallingHandlers(
      cic_inference(
        Y = data$Y,
        X = data$X,
        Z = data$Z,
        method = spec$methods,
        B = if (is.na(spec$B)) 1000L else spec$B,
        level = spec$level,
        panel_data = spec$panel_data,
        timings = FALSE
      ),
      warning = function(warning_condition) {
        warnings <<- c(warnings, conditionMessage(warning_condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error_condition) error_condition
  )

  if (inherits(fit, "error")) {
    status <- data.frame(
      call_id = spec$call_id,
      dataset_id = spec$dataset_id,
      status = "error",
      required = spec$required,
      warning_messages = paste(unique(warnings), collapse = " | "),
      error_message = conditionMessage(fit)
    )
    if (isTRUE(spec$required)) {
      stop(sprintf("Required reference fit failed for %s: %s",
                   spec$call_id, conditionMessage(fit)), call. = FALSE)
    }
    return(list(output = NULL, summary = NULL, status = status))
  }

  output <- make_output_rows(fit, spec, warnings)
  if (isTRUE(spec$required) && any(!output$finite_ci)) {
    stop(sprintf("Required reference fit produced non-finite CI output for %s.",
                 spec$call_id), call. = FALSE)
  }

  summary_lines <- utils::capture.output(summary(fit))
  summary <- data.frame(
    call_id = spec$call_id,
    line_no = seq_along(summary_lines),
    text = summary_lines
  )

  status <- data.frame(
    call_id = spec$call_id,
    dataset_id = spec$dataset_id,
    status = "success",
    required = spec$required,
    warning_messages = paste(unique(warnings), collapse = " | "),
    error_message = ""
  )

  list(output = output, summary = summary, status = status)
}

fit_results <- lapply(fit_specs, run_reference_fit)
reference_outputs <- do.call(rbind, lapply(fit_results, `[[`, "output"))
reference_summaries <- do.call(rbind, lapply(fit_results, `[[`, "summary"))
reference_status <- do.call(rbind, lapply(fit_results, `[[`, "status"))

dgp_function_references <- rbind(
  do.call(rbind, lapply(
    list(
      list(case_id = "qY_default_grid", d1 = 0, d2 = 0.05,
           t = c(0.01, 0.1, 0.5, 0.9, 0.99)),
      list(case_id = "qY_two_tail_grid", d1 = 0.08, d2 = 0.06,
           t = c(0.001, 0.05, 0.5, 0.95, 0.999))
    ),
    function(case) {
      data.frame(
        function_name = "qY_dgp",
        case_id = case$case_id,
        t = case$t,
        b1 = NA_real_,
        b2 = NA_real_,
        d1 = case$d1,
        d2 = case$d2,
        value = qY_dgp(case$t, d1 = case$d1, d2 = case$d2)
      )
    }
  )),
  do.call(rbind, lapply(
    list(
      list(case_id = "theta_default", b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05),
      list(case_id = "theta_two_tail", b1 = 0.06, b2 = 0.08, d1 = 0.02, d2 = 0.04),
      list(case_id = "theta_edge_heavy", b1 = 0.05, b2 = 0.05, d1 = 0.12, d2 = 0.12)
    ),
    function(case) {
      data.frame(
        function_name = "theta_true",
        case_id = case$case_id,
        t = NA_real_,
        b1 = case$b1,
        b2 = case$b2,
        d1 = case$d1,
        d2 = case$d2,
        value = theta_true(b1 = case$b1, b2 = case$b2, d1 = case$d1, d2 = case$d2)
      )
    }
  ))
)

generation_manifest <- data.frame(
  item = c(
    "repo_root",
    "load_method",
    "r_version",
    "rng_kind",
    "cross_section_methods",
    "panel_methods",
    "bootstrap_reference_included",
    "bootstrap_common_index_B",
    "bootstrap_common_index_datasets"
  ),
  value = c(
    repo_root,
    load_method,
    paste(R.version$major, R.version$minor, sep = "."),
    paste(RNGkind(), collapse = ";"),
    paste(cross_section_methods, collapse = ";"),
    paste(panel_methods, collapse = ";"),
    "TRUE",
    paste(unique(vapply(bootstrap_reference_specs, `[[`, integer(1L), "B")), collapse = ";"),
    paste(vapply(bootstrap_reference_specs, `[[`, character(1L), "dataset_id"), collapse = ";")
  )
)

write_csv(input_rows, file.path(reference_dir, "reference_inputs.csv"))
write_csv(dataset_manifest, file.path(reference_dir, "reference_dataset_manifest.csv"))
write_csv(fit_spec_manifest, file.path(reference_dir, "reference_fit_manifest.csv"))
write_csv(reference_outputs, file.path(reference_dir, "reference_outputs.csv"))
write_csv(reference_summaries, file.path(reference_dir, "reference_summaries.csv"))
write_csv(reference_status, file.path(reference_dir, "reference_status.csv"))
write_csv(dgp_function_references, file.path(reference_dir, "dgp_function_references.csv"))
write_csv(generation_manifest, file.path(reference_dir, "reference_generation_manifest.csv"))
write_csv(utility_nosplit_references, file.path(reference_dir, "utility_nosplit_references.csv"))
write_csv(bootstrap_index_references, file.path(reference_dir, "bootstrap_index_references.csv"))
write_csv(bootstrap_draw_references, file.path(reference_dir, "bootstrap_draw_references.csv"))
write_csv(bootstrap_summary_references, file.path(reference_dir, "bootstrap_summary_references.csv"))

cat("Reference generation complete.\n")
cat("Output directory:", reference_dir, "\n")
cat("Datasets:", nrow(dataset_manifest), "\n")
cat("Fit calls:", nrow(fit_spec_manifest), "\n")
cat("Reference output rows:", nrow(reference_outputs), "\n")
cat("Utility reference rows:", nrow(utility_nosplit_references), "\n")
cat("Bootstrap index rows:", nrow(bootstrap_index_references), "\n")
cat("Bootstrap draw rows:", nrow(bootstrap_draw_references), "\n")
cat("Statuses:\n")
print(table(reference_status$status, useNA = "ifany"))
