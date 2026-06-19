# ── print.cic ──────────────────────────────────────────────────────────────────

#' @export
print.cic <- function(x, digits = 4, ...) {
  cat("Changes-in-Changes Estimator\n")
  cat(rep("-", 44), "\n", sep = "")
  cat(sprintf("Estimate (theta_hat) : %.4f\n", x$theta_hat))
  cat(sprintf("Confidence level     : %.0f%%\n", x$level * 100))
  cat(sprintf("Sample sizes         : n1 = %d, n2 = %d, n3 = %d\n",
              x$n["n1"], x$n["n2"], x$n["n3"]))
  if (!is.na(x$h))
    cat(sprintf("Bandwidth (h)        : %.4f\n", x$h))
  cat(rep("-", 44), "\n", sep = "")
  cat("Confidence intervals:\n\n")

  # Format the CI table
  ci <- x$ci
  ci$lower  <- round(ci$lower,  digits)
  ci$upper  <- round(ci$upper,  digits)
  ci$length <- round(ci$length, digits)

  label_map <- c(split = "Split-sample", bpc = "Bootstrap (pct.)",
                 bse   = "Bootstrap (SE)")
  ci$method <- label_map[ci$method]

  print(ci, row.names = FALSE, right = FALSE)
  invisible(x)
}

# ── summary.cic ────────────────────────────────────────────────────────────────

#' @export
summary.cic <- function(object, digits = 4, ...) {
  cat("Changes-in-Changes — Summary\n")
  cat(rep("=", 44), "\n", sep = "")
  cat(sprintf("Point estimate : %.4f\n\n", object$theta_hat))

  cat("Sample sizes:\n")
  cat(sprintf("  Y sample (n1) : %d\n", object$n["n1"]))
  cat(sprintf("  X sample (n2) : %d\n", object$n["n2"]))
  cat(sprintf("  Z sample (n3) : %d\n", object$n["n3"]))

  cat(sprintf("\nConfidence level : %.0f%%\n", object$level * 100))

  if (!is.na(object$h))
    cat(sprintf("Bandwidth (h)    : %.4f  [Silverman rule-of-thumb]\n",
                object$h))

  cat(rep("-", 44), "\n", sep = "")
  cat("Confidence intervals:\n\n")

  ci <- object$ci
  label_map <- c(split = "Split-sample    ",
                 bpc   = "Bootstrap (pct.)",
                 bse   = "Bootstrap (SE)  ")
  for (i in seq_len(nrow(ci))) {
    cat(sprintf("  %s  [%.*f, %.*f]  (length %.4f)\n",
                label_map[ci$method[i]],
                digits, ci$lower[i],
                digits, ci$upper[i],
                ci$length[i]))
  }
  cat(rep("=", 44), "\n", sep = "")
  invisible(object)
}
