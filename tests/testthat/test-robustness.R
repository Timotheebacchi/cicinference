library(testthat)
library(quantcdf.inference)

test_that("fit() preserves requested method order", {
  d <- sim_dgp(180, seed = 101)
  out <- fit(d$Y, d$X, d$Z, method = c("ai", "no-split", "bpc"), B = 200)

  expect_equal(out$method, c("ai", "no-split", "bpc"))
  expect_equal(out$ci$method, c("ai", "no-split", "bpc"))
})

test_that("fit() accepts unequal sample sizes", {
  set.seed(102)
  sample1 <- rnorm(240)
  sample2 <- rnorm(180)
  sample3 <- rnorm(320)

  out <- fit(sample1, sample2, sample3, method = c("no-split", "split", "ai"))

  expect_equal(unname(out$n), c(240, 180, 320))
  expect_true(all(is.finite(out$ci$se)))
  expect_true(all(out$ci$length > 0))
})

test_that("epsilon_n can be supplied as a scalar or function", {
  d <- sim_dgp(220, seed = 103)

  scalar <- fit(d$Y, d$X, d$Z, method = "no-split", epsilon_n = 0.2)
  fun <- fit(d$Y, d$X, d$Z, method = "no-split", epsilon_n = function(n) 2 / log(n))

  expect_equal(scalar$epsilon_n, 0.2)
  expect_equal(fun$epsilon_n, 2 / log(length(d$X)))
  expect_true(is.finite(scalar$ci$se))
  expect_true(is.finite(fun$ci$se))
})

test_that("invalid epsilon_n values are rejected", {
  d <- sim_dgp(120, seed = 104)

  expect_error(fit(d$Y, d$X, d$Z, epsilon_n = 0), regexp = "epsilon_n")
  expect_error(fit(d$Y, d$X, d$Z, epsilon_n = -0.1), regexp = "epsilon_n")
  expect_error(fit(d$Y, d$X, d$Z, epsilon_n = Inf), regexp = "epsilon_n")
  expect_error(fit(d$Y, d$X, d$Z, epsilon_n = function(n) NA_real_), regexp = "epsilon_n")
})

test_that("bootstrap B validation matches public contract", {
  d <- sim_dgp(160, seed = 105)

  expect_error(fit(d$Y, d$X, d$Z, method = "bse", B = NA), regexp = "B must not")
  expect_error(fit(d$Y, d$X, d$Z, method = "bse", B = NULL), regexp = "B must not")
  expect_warning(
    out <- fit(d$Y, d$X, d$Z, method = "bse", B = 50),
    regexp = "B < 200"
  )
  expect_s3_class(out, "quantcdf_fit")
  expect_true(is.finite(out$ci$se))
})

test_that("bootstrap intervals are reproducible under a fixed RNG seed", {
  d <- sim_dgp(180, seed = 106)

  set.seed(1)
  first <- fit(d$Y, d$X, d$Z, method = c("bse", "bpc"), B = 220)
  set.seed(1)
  second <- fit(d$Y, d$X, d$Z, method = c("bse", "bpc"), B = 220)

  expect_equal(first$ci, second$ci)
})

test_that("bpc stores percentile intervals without an implied standard error", {
  d <- sim_dgp(180, seed = 107)
  out <- fit(d$Y, d$X, d$Z, method = "bpc", B = 220)

  expect_true(is.na(out$ci$se))
  expect_true(is.na(out$ci$sigma_sq))
  expect_true(out$ci$lower < out$ci$upper)
})

test_that("confidence level changes interval length monotonically", {
  d <- sim_dgp(220, seed = 108)

  low <- fit(d$Y, d$X, d$Z, method = "ai", level = 0.90)
  high <- fit(d$Y, d$X, d$Z, method = "ai", level = 0.99)

  expect_lt(low$ci$length, high$ci$length)
  expect_error(fit(d$Y, d$X, d$Z, level = 1), regexp = "level")
  expect_error(fit(d$Y, d$X, d$Z, level = 0), regexp = "level")
})

test_that("AI method rejects zero Sample 1 dispersion", {
  set.seed(109)
  sample1 <- rep(1, 80)
  sample2 <- rnorm(80)
  sample3 <- rnorm(80)

  expect_error(
    suppressWarnings(fit(sample1, sample2, sample3, method = "ai")),
    regexp = "positive finite empirical standard deviation"
  )
})

test_that("input coercion warnings are specific and usable", {
  d <- sim_dgp(100, seed = 110)

  expect_warning(
    out <- fit(matrix(d$Y, ncol = 1), d$X, d$Z, method = "no-split"),
    regexp = "Sample 1 should be a numeric vector"
  )
  expect_s3_class(out, "quantcdf_fit")
  expect_error(fit(as.character(d$Y), d$X, d$Z), regexp = "Sample 1 must be")
})

test_that("many ties in Sample 1 warn but still return an object for non-AI methods", {
  set.seed(111)
  sample1 <- rep(c(0, 1), each = 60)
  sample2 <- rnorm(120)
  sample3 <- rnorm(120)

  expect_warning(
    out <- fit(sample1, sample2, sample3, method = "no-split"),
    regexp = "many tied values"
  )
  expect_s3_class(out, "quantcdf_fit")
  expect_true(is.finite(out$theta_hat))
})

test_that("confint() warns when asked for a different stored level", {
  d <- sim_dgp(140, seed = 112)
  out <- fit(d$Y, d$X, d$Z, method = "no-split", level = 0.9)

  expect_warning(confint(out, level = 0.95), regexp = "differs")
})

test_that("cic_fit() supports all analytic methods together", {
  d <- sim_dgp(220, seed = 113)
  Y11 <- d$Y + 0.25

  out <- cic_fit(Y11, d$X, d$Y, d$Z, method = c("no-split", "split", "ai"))

  expect_equal(out$ci$method, c("no-split", "split", "ai"))
  expect_true(all(is.finite(out$ci$se)))
  expect_true(all(out$ci$sigma_sq >= 0))
})

test_that("cic_fit() handles unequal repeated-cross-section sizes", {
  set.seed(114)
  Y11 <- rnorm(260, 1)
  Y10 <- rnorm(180)
  Y01 <- rnorm(220)
  Y00 <- rnorm(300)

  out <- cic_fit(Y11, Y10, Y01, Y00, method = c("no-split", "ai"))

  expect_equal(unname(out$n), c(260, 180, 220, 300))
  expect_equal(out$components$N, 180)
  expect_true(all(is.finite(out$ci$length)))
})

test_that("cic_fit() bootstrap is reproducible under a fixed RNG seed", {
  d <- sim_dgp(180, seed = 115)
  Y11 <- d$Y + 0.3

  set.seed(2)
  first <- cic_fit(Y11, d$X, d$Y, d$Z, method = c("bse", "bpc"), B = 220)
  set.seed(2)
  second <- cic_fit(Y11, d$X, d$Y, d$Z, method = c("bse", "bpc"), B = 220)

  expect_equal(first$theta_hat, second$theta_hat)
  expect_equal(first$ci, second$ci)
})

test_that("sim_dgp() panel mode pairs Sample 1 and Sample 3 through the same rank", {
  d <- sim_dgp(120, seed = 116, panel_data = TRUE)

  expect_equal(length(d$Y), length(d$Z))
  expect_equal(stats::pnorm(d$Z), 1 - (1 + d$Y)^(-20), tolerance = 1e-12)
})

test_that("qY_dgp() and theta_true() return stable benchmark values", {
  expect_equal(qY_dgp(0.5, d1 = 0, d2 = 0), 0)
  expect_equal(theta_true(b1 = 0, b2 = 0, d1 = 0, d2 = 0), 0)
  expect_true(is.finite(theta_true(b1 = 0.05, b2 = 0.05, d1 = 0.02, d2 = 0.02)))
})
