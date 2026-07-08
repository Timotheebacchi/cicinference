# ── Internal helpers ───────────────────────────────────────────────────────────

.cic_inference_table <- function(object) {
  ci <- object$ci
  alpha <- (1 - object$level) / 2
  z_a <- stats::qnorm(1 - alpha)

  se <- ifelse(ci$method == "bpc", NA_real_, ci$length / (2 * z_a))
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


# ── summary.cic ────────────────────────────────────────────────────────────────

#' @export
summary.cic_inference <- function(object, digits = 4, ...) {
  cat("Changes-in-Changes - Summary\n")
  cat(rep("=", 44), "\n", sep = "")
  

  cat("Sample sizes:\n")
  cat(sprintf("  Y sample (n1) : %d\n", object$n["n1"]))
  cat(sprintf("  X sample (n2) : %d\n", object$n["n2"]))
  cat(sprintf("  Z sample (n3) : %d\n", object$n["n3"]))

  cat(sprintf("\nConfidence level : %.0f%%\n", object$level * 100))

  if (!is.na(object$epsilon_n))
    cat(sprintf("Bandwidth (epsilon_n)    : %.4f  [1/log(n2) plug-in default]\n",
                object$epsilon_n))
  else
    cat("Bandwidth (epsilon_n)    : NA  [not applicable for bootstrap methods]\n")

  cat(rep("-", 44), "\n", sep = "")
  cat(sprintf("Point estimate : %.4f\n\n", object$theta_hat))
  cat("Inference table:\n\n")

  tab <- .cic_inference_table(object)
  numeric_cols <- vapply(tab, is.numeric, logical(1))
  tab[numeric_cols] <- lapply(tab[numeric_cols], round, digits = digits)

  label_map <- c(`no-split` = "No-split        ",
                 split = "Split           ",
                 kde   = "KDE             ",
                 bpc   = "Bootstrap (pct.)",
                 bse   = "Bootstrap (SE)  ")
  tab$Method <- label_map[tab$Method]
  print(tab, row.names = FALSE, right = FALSE)
  cat(rep("=", 44), "\n", sep = "")
  invisible(object)
}

#' @export
coef.cic_inference <- function(object, ...) {
  object$theta_hat
}

#' @export
confint.cic_inference <- function(object, parm = NULL, level = object$level, ...) {
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

