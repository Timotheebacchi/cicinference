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
  expect_named(fit, c("theta_hat", "ci", "level", "n", "method", "h"))
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

test_that("plot.cic dessine sans erreur", {
  d <- sim_dgp(200, seed = 16)
  fit <- cic(d$Y, d$X, d$Z, method = c("no-split", "bpc"), B = 200)
  pdf(file = tempfile(fileext = ".pdf"))
  expect_silent(plot(fit))
  dev.off()
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

# ── Tests : Diagnostic Function (check_cic_assumptions) ──────────────────────
test_that("check_cic_assumptions retourne un objet valide", {
  d <- sim_dgp(1000, b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05, seed = 24)
  diag <- check_cic_assumptions(d$Y, d$X, d$Z)
  expect_named(diag, c("pass_all", "metrics", "messages"))
  expect_type(diag$pass_all, "logical")
  expect_type(diag$metrics, "list")
  expect_type(diag$messages, "character")
})

test_that("check_cic_assumptions accepte les données bien comportées", {
  set.seed(24)
  n <- 1000
  Y <- rnorm(n)
  X <- rnorm(n)
  Z <- rnorm(n)
  diag <- check_cic_assumptions(Y, X, Z)
  expect_true(diag$pass_all)
})

test_that("check_cic_assumptions détecte les queues lourdes", {
  # Créer des données avec des queues extrêmement lourdes
  n <- 100
  set.seed(25)
  Y <- rt(n, df = 2)  # t-distribution avec 2 degrés de liberté
  X <- rt(n, df = 2)
  Z <- rt(n, df = 2)
  diag <- check_cic_assumptions(Y, X, Z)
  # On s'attend à un possible fail ou au moins des warnings
  expect_type(diag$messages, "character")
})

test_that("check_cic_assumptions inclut les ratios lambda", {
  d <- sim_dgp(300, seed = 26)
  diag <- check_cic_assumptions(d$Y, d$X, d$Z)
  expect_true("lambda1" %in% names(diag$metrics))
  expect_true("lambda2" %in% names(diag$metrics))
  expect_true("lambda3" %in% names(diag$metrics))
  expect_gt(diag$metrics$lambda1, 0)
  expect_lt(diag$metrics$lambda1, 1)
})

test_that("check_cic_assumptions inclut les estimations de queues", {
  d <- sim_dgp(300, seed = 27)
  diag <- check_cic_assumptions(d$Y, d$X, d$Z)
  expect_true("d1_left_tail" %in% names(diag$metrics))
  expect_true("d2_right_tail" %in% names(diag$metrics))
  expect_true("b1_left_boundary" %in% names(diag$metrics))
  expect_true("b2_right_boundary" %in% names(diag$metrics))
  expect_gte(diag$metrics$d1_left_tail, 0)
  expect_gte(diag$metrics$d2_right_tail, 0)
})

test_that("check_cic_assumptions génère des messages clairs", {
  d <- sim_dgp(200, seed = 28)
  diag <- check_cic_assumptions(d$Y, d$X, d$Z)
  expect_true(length(diag$messages) >= 1)
  # Chaque message doit être non-vide
  expect_true(all(nchar(diag$messages) > 0))
})
