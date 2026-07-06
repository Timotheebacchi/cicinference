library(testthat)
library(cic.newassumptions.newvarianceestimator)

# ── DGP helpers ───────────────────────────────────────────────────────────────
# Quantile function from the Monte Carlo section of the paper:
#   F_Y^{-1}(t) = -t^{-d1} + (1-t)^{-d2}
# with the convention that d1 = 0 => t^{-d1} = 1, (1-t)^{-d2} with d2=0 => 1
qY_dgp <- function(t, d1 = 0, d2 = 0.05) {
  term1 <- if (d1 == 0) 1 else t^(-d1)
  term2 <- if (d2 == 0) 1 else (1 - t)^(-d2)
  -term1 + term2
}

# True theta_0 = [B(1-b1, 1-b2-d2) - B(1-b1-d1, 1-b2)] / B(1-b1, 1-b2)
theta_true <- function(b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05) {
  (beta(1 - b1, 1 - b2 - d2) - beta(1 - b1 - d1, 1 - b2)) /
    beta(1 - b1, 1 - b2)
}

# Simulate one dataset from the DGP
sim_dgp <- function(n, b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  W <- runif(n)
  Y <- qY_dgp(W, d1, d2)
  Z <- rnorm(n)
  V <- rbeta(n, 1 - b1, 1 - b2)
  X <- qnorm(V)
  list(Y = Y, X = X, Z = Z)
}

# ── Tests : structure de l'output ─────────────────────────────────────────────
test_that("cic() retourne un objet de classe 'cic'", {
  d <- sim_dgp(200, seed = 1)
  fit <- cic(d$Y, d$X, d$Z, method = "no-split", B = 200)
  expect_s3_class(fit, "cic")
})

test_that("cic() retourne les bons champs", {
  d <- sim_dgp(200, seed = 2)
  fit <- cic(d$Y, d$X, d$Z, method = c("no-split", "bse", "bpc"), B = 200)
  expect_named(fit, c("theta_hat", "ci", "level", "n", "method", "panel_data", "h"))
  expect_equal(nrow(fit$ci), 3)
  expect_equal(fit$ci$method, c("no-split", "bse", "bpc"))
})

test_that("les tailles d'échantillon sont bien enregistrées", {
  d <- sim_dgp(200, seed = 3)
  fit <- cic(d$Y, d$X, d$Z, method = "no-split")
  expect_equal(unname(fit$n), c(200L, 200L, 200L))
})

test_that("le niveau de confiance est respecté", {
  d <- sim_dgp(200, seed = 4)
  fit90 <- cic(d$Y, d$X, d$Z, method = "no-split", level = 0.90)
  fit95 <- cic(d$Y, d$X, d$Z, method = "no-split", level = 0.95)
  # IC à 90% doit être plus court qu'à 95%
  expect_lt(fit90$ci$length[1], fit95$ci$length[1])
})

# ── Tests : validité des intervalles ──────────────────────────────────────────
test_that("lower < upper pour toutes les méthodes", {
  d <- sim_dgp(200, seed = 5)
  fit <- cic(d$Y, d$X, d$Z, method = c("no-split", "bse", "bpc"), B = 200)
  expect_true(all(fit$ci$lower < fit$ci$upper))
})

test_that("length == upper - lower", {
  d <- sim_dgp(200, seed = 6)
  fit <- cic(d$Y, d$X, d$Z, method = c("no-split", "bse", "bpc"), B = 200)
  expect_equal(fit$ci$length, fit$ci$upper - fit$ci$lower, tolerance = 1e-10)
})

test_that("theta_hat est dans l'IC (no-split) pour un grand n", {
  d <- sim_dgp(2000, seed = 7)
  fit <- cic(d$Y, d$X, d$Z, method = "no-split")
  expect_gte(fit$theta_hat, fit$ci$lower[1])
  expect_lte(fit$theta_hat, fit$ci$upper[1])
})

# ── Tests : convergence vers theta_0 ──────────────────────────────────────────
test_that("theta_hat converge vers theta_0 pour n grand", {
  # DGP de base : b1=0, b2=0.05, d1=0, d2=0.05
  theta_0 <- theta_true(b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05)
  d   <- sim_dgp(5000, b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05, seed = 42)
  fit <- cic(d$Y, d$X, d$Z, method = "no-split")
  # Pour n=5000, on tolère 3% d'écart
  expect_equal(fit$theta_hat, theta_0, tolerance = 0.03)
})

test_that("theta_0 est couvert par l'IC no-split pour n grand", {
  theta_0 <- theta_true(b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05)
  d   <- sim_dgp(5000, b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05, seed = 99)
  fit <- cic(d$Y, d$X, d$Z, method = "no-split")
  expect_gte(theta_0, fit$ci$lower[1])
  expect_lte(theta_0, fit$ci$upper[1])
})

# ── Tests : cas limites et erreurs attendues ───────────────────────────────────
test_that("B < 200 déclenche un warning et est corrigé à 200", {
  d <- sim_dgp(200, seed = 8)
  expect_warning(
    cic(d$Y, d$X, d$Z, method = "bse", B = 50),
    regexp = "B < 200"
  )
})

test_that("vecteurs non numériques déclenchent une erreur", {
  d <- sim_dgp(200, seed = 9)
  expect_error(cic(as.character(d$Y), d$X, d$Z))
})

test_that("taille impaire est acceptée pour les méthodes bootstrap", {
  d <- sim_dgp(201, seed = 17)
  expect_s3_class(cic(d$Y, d$X, d$Z, method = "bse", B = 200), "cic")
  expect_s3_class(cic(d$Y, d$X, d$Z, method = "bpc", B = 200), "cic")
})

test_that("taille impaire est acceptée pour la méthode no-split", {
  d <- sim_dgp(201, seed = 18)
  expect_s3_class(cic(d$Y, d$X, d$Z, method = "no-split"), "cic")
})

test_that("method invalide déclenche une erreur", {
  d <- sim_dgp(200, seed = 10)
  expect_error(cic(d$Y, d$X, d$Z, method = "mauvaise_methode"))
})

test_that("bandwidth manuel est bien utilisé", {
  d   <- sim_dgp(200, seed = 11)
  fit <- cic(d$Y, d$X, d$Z, method = "no-split", h = 0.5)
  expect_equal(fit$h, 0.5)
})

test_that("h influe sur l'estimation no-split", {
  d <- sim_dgp(200, seed = 20)
  fit_small <- cic(d$Y, d$X, d$Z, method = "no-split", h = 0.05)
  fit_large <- cic(d$Y, d$X, d$Z, method = "no-split", h = 0.5)
  expect_false(identical(fit_small$ci$length, fit_large$ci$length))
})

test_that("h nul ou négatif déclenche une erreur", {
  d <- sim_dgp(200, seed = 19)
  expect_error(cic(d$Y, d$X, d$Z, h = 0), regexp = "h")
  expect_error(cic(d$Y, d$X, d$Z, h = -0.5), regexp = "h")
})

# ── Tests : bootstrap ─────────────────────────────────────────────────────────
test_that("bse et bpc donnent des IC différents", {
  d   <- sim_dgp(200, seed = 12)
  fit <- cic(d$Y, d$X, d$Z, method = c("bse", "bpc"), B = 200)
  # Pas identiques (sauf coïncidence impossible)
  expect_false(fit$ci$lower[1] == fit$ci$lower[2])
})

test_that("split et kde sont disponibles", {
  d <- sim_dgp(200, seed = 23)
  fit <- cic(d$Y, d$X, d$Z, method = c("split", "kde"))
  expect_equal(fit$ci$method, c("split", "kde"))
  expect_true(all(is.finite(fit$ci$length)))
  expect_true(all(fit$ci$lower < fit$ci$upper))
})

test_that("plus de réplications bootstrap donne des IC plus stables", {
  d    <- sim_dgp(500, seed = 13)
  set.seed(1); fit_small <- cic(d$Y, d$X, d$Z, method = "bse", B = 200)
  set.seed(1); fit_large <- cic(d$Y, d$X, d$Z, method = "bse", B = 999)
  # Les deux doivent être du même ordre de grandeur (pas un test de précision,
  # juste qu'on n'explose pas)
  expect_lt(abs(fit_small$ci$length[1] - fit_large$ci$length[1]), 0.2)
})

test_that("coef.cic retourne theta_hat", {
  d <- sim_dgp(200, seed = 14)
  fit <- cic(d$Y, d$X, d$Z, method = "no-split")
  expect_equal(coef(fit), fit$theta_hat)
})

test_that("confint.cic retourne les intervalles de confiance", {
  d <- sim_dgp(200, seed = 15)
  fit <- cic(d$Y, d$X, d$Z, method = c("no-split", "bse"), B = 200)
  ci <- confint(fit)
  expected <- as.matrix(fit$ci[, c("lower", "upper")])
  rownames(expected) <- fit$ci$method
  expect_equal(ci, expected)
  expect_equal(rownames(ci), c("no-split", "bse"))
})

test_that("summary.cic montre une sortie familière aux économistes", {
  d <- sim_dgp(200, seed = 29)
  fit <- cic(d$Y, d$X, d$Z, method = "no-split")
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Std. Error", out, fixed = TRUE)))
  expect_true(any(grepl("t value", out, fixed = TRUE)))
  expect_true(any(grepl("Pr(>|t|)", out, fixed = TRUE)))
})

test_that("timings = TRUE affiche des jalons de calcul", {
  d <- sim_dgp(200, seed = 30)
  msgs <- character()
  fit <- withCallingHandlers(
    cic(d$Y, d$X, d$Z, method = c("no-split", "bse"), B = 200, timings = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_s3_class(fit, "cic")
  expect_true(any(grepl("cic timing \\[bootstrap\\]", msgs)))
})


# ── Tests : Input Sanitization (Code Audit) ──────────────────────────────────
test_that("B = NA déclenche une erreur", {
  d <- sim_dgp(200, seed = 21)
  expect_error(cic(d$Y, d$X, d$Z, method = "bse", B = NA),
               regexp = "B must not be NULL or NA")
})

test_that("B = NULL déclenche une erreur", {
  d <- sim_dgp(200, seed = 22)
  expect_error(cic(d$Y, d$X, d$Z, method = "bse", B = NULL),
               regexp = "B must not be NULL or NA")
})

# Basic input validation tests (replaced extensive diagnostics tests)
test_that("non-vector inputs trigger a warning", {
  d <- sim_dgp(100, seed = 40)
  # pass a matrix instead of a vector
  Ym <- matrix(d$Y, ncol = 1)
  expect_warning(cic(Ym, d$X, d$Z), regexp = "should be a numeric vector|coercing to numeric vector")
})

test_that("inputs with many ties produce a warning", {
  n <- 200
  Y <- rep(1, n)
  X <- rnorm(n)
  Z <- rnorm(n)
  expect_warning(cic(Y, X, Z), regexp = "tied values|many ties|uniqueness")
})

Here is the complete test block translated into English, ready to be added to your testthat file to validate the new panel_data behavior.
R

# ── Tests: Panel Data (panel_data) ─────────────────────────────────────
test_that("panel_data = TRUE works for 'no-split' and 'split'", {
  d <- sim_dgp(200, seed = 101)
  
  # Test the 'no-split' method in panel mode
  fit_panel_nosplit <- cic(d$Y, d$X, d$Z, method = "no-split", panel_data = TRUE)
  expect_s3_class(fit_panel_nosplit, "cic")
  expect_true(fit_panel_nosplit$panel_data)
  expect_true(all(is.finite(fit_panel_nosplit$ci$length)))
  
  # Test the 'split' method in panel mode (new integration)
  fit_panel_split <- cic(d$Y, d$X, d$Z, method = "split", panel_data = TRUE)
  expect_s3_class(fit_panel_split, "cic")
  expect_true(fit_panel_split$panel_data)
  expect_true(all(is.finite(fit_panel_split$ci$length)))
})

test_that("panel_data = TRUE accepts simultaneous calls for 'no-split' and 'split'", {
  d <- sim_dgp(200, seed = 102)
  
  # Verification of the logical bug fix: calling both methods simultaneously should pass smoothly
  expect_silent(
    fit_both <- cic(d$Y, d$X, d$Z, method = c("no-split", "split"), panel_data = TRUE)
  )
  expect_equal(nrow(fit_both$ci), 2)
  expect_equal(fit_both$ci$method, c("no-split", "split"))
})

test_that("panel_data = TRUE throws an error with unsupported methods", {
  d <- sim_dgp(200, seed = 103)
  
  # A single unsupported method
  expect_error(
    cic(d$Y, d$X, d$Z, method = "kde", panel_data = TRUE),
    regexp = "panel_data = TRUE is currently supported only for method"
  )
  
  # A mix containing an unsupported method
  expect_error(
    cic(d$Y, d$X, d$Z, method = c("no-split", "bse"), panel_data = TRUE),
    regexp = "panel_data = TRUE is currently supported only for method"
  )
})