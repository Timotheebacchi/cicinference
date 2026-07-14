# Internal helpers for S3 output.

.quantcdf_table <- function(object) {
  ci <- object$ci
  se <- ci$se
  stat <- ifelse(is.na(se) | se <= 0, NA_real_, object$theta_hat / se)
  p_value <- ifelse(is.na(stat), NA_real_, 2 * stats::pnorm(abs(stat), lower.tail = FALSE))

  out <- data.frame(
    Method = ci$method,
    Std.Error = se,
    t.value = stat,
    p.value = p_value,
    Lower.CI = ci$lower,
    Upper.CI = ci$upper,
    Length = ci$length,
    check.names = FALSE
  )

  names(out)[2:7] <- c("Std. Error", "t value", "Pr(>|t|)", "Lower CI", "Upper CI", "Length")
  out
}

#' @export
summary.quantcdf_fit <- function(object, digits = 4, ...) {
  title <- if (inherits(object, "quantcdf_cic_fit")) {
    "Changes-in-Changes - Summary"
  } else {
    "Empirical Quantile Inference - Summary"
  }
  cat(title, "\n", sep = "")
  cat(rep("=", 44), "\n", sep = "")

  cat("Sample sizes:\n")
  if (inherits(object, "quantcdf_cic_fit")) {
    cat(sprintf("  Y11 sample (n11) : %d\n", object$n["n11"]))
    cat(sprintf("  Y10 sample (n10) : %d\n", object$n["n10"]))
    cat(sprintf("  Y01 sample (n01) : %d\n", object$n["n01"]))
    cat(sprintf("  Y00 sample (n00) : %d\n", object$n["n00"]))
  } else {
    cat(sprintf("  Sample 1 (n1) : %d\n", object$n["n1"]))
    cat(sprintf("  Sample 2 (n2) : %d\n", object$n["n2"]))
    cat(sprintf("  Sample 3 (n3) : %d\n", object$n["n3"]))
  }

  cat(sprintf("\nConfidence level : %.0f%%\n", object$level * 100))

  if (!is.na(object$epsilon_n)) {
    cat(sprintf("Bandwidth epsilon_n : %.4f  [split/no-split default is 1/log(n2)]\n",
                object$epsilon_n))
  }
  if (!is.na(object$h_ai)) {
    cat(sprintf("AI bandwidth h_AI   : %.4f\n", object$h_ai))
  }

  cat(rep("-", 44), "\n", sep = "")
  cat(sprintf("Point estimate : %.4f\n\n", object$theta_hat))
  cat("Inference table:\n\n")

  tab <- .quantcdf_table(object)
  numeric_cols <- vapply(tab, is.numeric, logical(1))
  tab[numeric_cols] <- lapply(tab[numeric_cols], round, digits = digits)

  label_map <- c(
    `no-split` = "No-split        ",
    split = "Split           ",
    ai = "Athey-Imbens   ",
    bpc = "Bootstrap (pct.)",
    bse = "Bootstrap (SE)  "
  )
  tab$Method <- unname(label_map[tab$Method])
  print(tab, row.names = FALSE, right = FALSE)
  cat(rep("=", 44), "\n", sep = "")
  invisible(object)
}

#' @export
coef.quantcdf_fit <- function(object, ...) {
  object$theta_hat
}

#' @export
confint.quantcdf_fit <- function(object, parm = NULL, level = object$level, ...) {
  if (!is.null(level) && !identical(level, object$level)) {
    warning(
      "`level` differs from the confidence level used to fit the object. ",
      "Returning stored intervals from `object$level`."
    )
  }
  ci <- object$ci[, c("lower", "upper")]
  rownames(ci) <- object$ci$method
  as.matrix(ci)
}
