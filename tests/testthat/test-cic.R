library(testthat)
library(quantcdf.inference)

test_that("fit() exposes the new public API", {
  expect_equal(names(formals(fit))[1:4], c("sample1", "sample2", "sample3", "method"))

  d <- sim_dgp(200, seed = 1)
  out <- fit(d$Y, d$X, d$Z, method = c("no-split", "split", "ai"))

  expect_s3_class(out, "quantcdf_fit")
  expect_named(
    out,
    c("theta_hat", "ci", "level", "n", "method", "panel_data", "epsilon_n", "h_ai")
  )
  expect_equal(out$ci$method, c("no-split", "split", "ai"))
  expect_true(all(is.finite(out$ci$lower)))
  expect_true(all(is.finite(out$ci$upper)))
  expect_true(all(out$ci$lower < out$ci$upper))
})

test_that("fit() preserves the left-continuous quantile convention at zero", {
  sample1 <- c(10, 20, 30, 40)
  sample2 <- c(-1, 0, 1, 2)
  sample3 <- c(0, 1, 2, 3)

  out <- fit(sample1, sample2, sample3, method = "no-split")

  expect_equal(out$theta_hat, mean(c(10, 10, 20, 30)))
})

test_that("obsolete public names and method labels are rejected", {
  d <- sim_dgp(200, seed = 2)
  old_function <- paste0("cic", "_inference")
  old_method <- paste0("k", "de")

  expect_false(exists(old_function, mode = "function"))
  expect_error(fit(d$Y, d$X, d$Z, method = old_method))
  expect_error(fit(d$Y, d$X, d$Z, method = "no-split", h = 0.5),
               regexp = "Unused argument")
})

test_that("bandwidths follow the split/no-split and AI formulas", {
  d <- sim_dgp(300, seed = 3)

  out_no_split <- fit(d$Y, d$X, d$Z, method = "no-split")
  expect_equal(out_no_split$epsilon_n, 1 / log(length(d$X)))
  expect_true(is.na(out_no_split$h_ai))

  out_ai <- fit(d$Y, d$X, d$Z, method = "ai")
  expect_true(is.na(out_ai$epsilon_n))
  expect_equal(out_ai$h_ai, 1.06 * length(d$Y)^(-1 / 5) / stats::sd(d$Y))
})

test_that("AI method is robust for targeted generated-data seeds", {
  seeds <- c(1, 2, 42, 123, 2026, 9999)

  for (seed in seeds) {
    d <- sim_dgp(500, seed = seed)
    out <- fit(d$Y, d$X, d$Z, method = "ai")
    expect_true(is.finite(out$theta_hat), info = paste("seed", seed))
    expect_true(all(is.finite(out$ci$length)), info = paste("seed", seed))
    expect_true(all(out$ci$lower < out$ci$upper), info = paste("seed", seed))
  }
})

test_that("AI method is deterministic in the ambient RNG seed for fixed data", {
  d <- sim_dgp(500, seed = 2026)
  ambient_seeds <- c(1, 2, 42, 123, 2026, 9999)

  lengths <- vapply(ambient_seeds, function(seed) {
    set.seed(seed)
    fit(d$Y, d$X, d$Z, method = "ai")$ci$length
  }, numeric(1))

  expect_true(all(is.finite(lengths)))
  expect_equal(lengths, rep(lengths[[1]], length(lengths)))
})

test_that("AI method is finite on a larger deterministic seed collection", {
  seeds <- 100:125

  lengths <- vapply(seeds, function(seed) {
    d <- sim_dgp(250, seed = seed)
    fit(d$Y, d$X, d$Z, method = "ai")$ci$length
  }, numeric(1))

  expect_true(all(is.finite(lengths)))
  expect_true(all(lengths > 0))
})

test_that("coef(), confint(), and summary() use the quantcdf_fit class", {
  d <- sim_dgp(200, seed = 4)
  out <- fit(d$Y, d$X, d$Z, method = c("no-split", "bse"), B = 200)

  expect_equal(coef(out), out$theta_hat)
  expect_equal(rownames(confint(out)), c("no-split", "bse"))

  printed <- capture.output(summary(out))
  expect_true(any(grepl("Std. Error", printed, fixed = TRUE)))
  expect_true(any(grepl("Empirical Quantile Inference", printed, fixed = TRUE)))
})

test_that("panel_data remains limited to split and no-split", {
  d <- sim_dgp(200, seed = 5)

  expect_s3_class(
    fit(d$Y, d$X, d$Z, method = c("no-split", "split"), panel_data = TRUE),
    "quantcdf_fit"
  )
  expect_error(
    fit(d$Y, d$X, d$Z, method = "ai", panel_data = TRUE),
    regexp = "panel_data = TRUE"
  )
})

test_that("cic_fit() delegates the counterfactual point estimate to fit()", {
  d <- sim_dgp(300, seed = 6)
  Y11 <- d$Y + 1

  theta_fit <- fit(sample1 = d$Y, sample2 = d$X, sample3 = d$Z, method = "no-split")
  out <- cic_fit(Y11, d$X, d$Y, d$Z, method = "no-split")

  expect_s3_class(out, "quantcdf_cic_fit")
  expect_equal(out$theta_hat, mean(Y11) - theta_fit$theta_hat)
  expect_equal(out$components$counterfactual_theta_hat, theta_fit$theta_hat)
})

test_that("cic_fit() adds the Y11 contribution on the sqrt(N) variance scale", {
  d <- sim_dgp(300, seed = 7)
  Y11 <- d$Y + rnorm(length(d$Y), mean = 1, sd = 0.2)

  theta_fit <- fit(d$Y, d$X, d$Z, method = "no-split")
  out <- cic_fit(Y11, d$X, d$Y, d$Z, method = "no-split")

  N <- min(length(Y11), length(d$X), length(d$Y), length(d$Z))
  V11 <- N / length(Y11) * (mean(Y11^2) - mean(Y11)^2)
  expected_sigma <- theta_fit$ci$sigma_sq[theta_fit$ci$method == "no-split"] + V11

  expect_equal(out$ci$sigma_sq[out$ci$method == "no-split"], expected_sigma)
  expect_equal(out$ci$se[out$ci$method == "no-split"], sqrt(expected_sigma / N))
})

test_that("cic_fit() bootstraps the complete estimator directly", {
  d <- sim_dgp(250, seed = 8)
  Y11_constant <- rep(1, length(d$Y))
  Y11_variable <- rep(c(0, 2), length.out = length(d$Y))

  set.seed(10)
  out_constant <- cic_fit(Y11_constant, d$X, d$Y, d$Z, method = "bse", B = 250)
  set.seed(10)
  out_variable <- cic_fit(Y11_variable, d$X, d$Y, d$Z, method = "bse", B = 250)

  expect_true(is.finite(out_constant$ci$length))
  expect_true(is.finite(out_variable$ci$length))
  expect_gt(out_variable$ci$length, out_constant$ci$length)
})

test_that("cic_fit() percentile bootstrap returns estimator-scale intervals", {
  d <- sim_dgp(250, seed = 9)
  Y11 <- d$Y + 0.5

  set.seed(11)
  out <- cic_fit(Y11, d$X, d$Y, d$Z, method = c("bse", "bpc"), B = 250)

  expect_equal(out$ci$method, c("bse", "bpc"))
  expect_true(all(is.finite(out$ci$lower)))
  expect_true(all(is.finite(out$ci$upper)))
  expect_true(all(out$ci$lower < out$ci$upper))
  expect_true(is.na(out$ci$se[out$ci$method == "bpc"]))
})
