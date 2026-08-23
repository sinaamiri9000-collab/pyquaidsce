#' @importFrom stats coef vcov sigma fitted residuals printCoefmat pnorm
NULL

#' Extract Model Coefficients
#' @param object A \code{quaidsce} object.
#' @param ... Additional arguments (ignored).
#' @return A named numeric vector of estimated coefficients.
#' @export
#' @method coef quaidsce
coef.quaidsce <- function(object, ...) {
  object$coefficients
}

#' Extract Variance-Covariance Matrix
#' @param object A \code{quaidsce} object.
#' @param ... Additional arguments (ignored).
#' @return A numeric matrix of parameter variances and covariances.
#' @export
#' @method vcov quaidsce
vcov.quaidsce <- function(object, ...) {
  object$vcov
}

#' Extract Residual Covariance Matrix
#' @param object A \code{quaidsce} object.
#' @param ... Additional arguments (ignored).
#' @return The residual covariance matrix Sigma.
#' @export
#' @method sigma quaidsce
sigma.quaidsce <- function(object, ...) {
  object$sigma
}

#' Extract Fitted Values
#' @param object A \code{quaidsce} object.
#' @param ... Additional arguments (ignored).
#' @return Matrix of fitted budget shares.
#' @export
#' @method fitted quaidsce
fitted.quaidsce <- function(object, ...) {
  object$fitted.values
}

#' Extract Model Residuals
#' @param object A \code{quaidsce} object.
#' @param ... Additional arguments (ignored).
#' @return Matrix of budget share residuals.
#' @export
#' @method residuals quaidsce
residuals.quaidsce <- function(object, ...) {
  object$residuals
}

#' Extract Demand Elasticities
#'
#' @param object A \code{quaidsce} object.
#' @param type Character vector indicating which elasticities to return:
#'   \code{"all"} (default), \code{"income"} (or expenditure), \code{"uncompensated"} (Marshallian),
#'   or \code{"compensated"} (Hicksian).
#' @param ... Additional arguments.
#' @return A list or matrix of elasticity estimates at sample means.
#' @export
elasticities <- function(object, ...) {
  UseMethod("elasticities")
}

#' @rdname elasticities
#' @export
#' @method elasticities quaidsce
elasticities.quaidsce <- function(object, type = c("all", "income", "uncompensated", "compensated"), ...) {
  type <- match.arg(type)
  if (type == "income") {
    return(object$elasticities$income)
  } else if (type == "uncompensated") {
    return(object$elasticities$uncompensated)
  } else if (type == "compensated") {
    return(object$elasticities$compensated)
  } else {
    return(object$elasticities)
  }
}

#' Extract Elasticity Standard Errors
#'
#' @param object A \code{quaidsce} object.
#' @param ... Additional arguments.
#' @return Information on elasticity standard errors (available when bootstrap replications were estimated).
#' @export
elasticity_se <- function(object, ...) {
  UseMethod("elasticity_se")
}

#' @rdname elasticity_se
#' @export
#' @method elasticity_se quaidsce
elasticity_se.quaidsce <- function(object, ...) {
  if (is.null(object$bootstrap)) {
    message("Note: Elasticity standard errors require bootstrap estimation (set `reps > 0` in `quaidsce()`).")
    return(invisible(NULL))
  }
  return(object$bootstrap$se)
}

#' Print Method for quaidsce
#' @param x A \code{quaidsce} object.
#' @param ... Additional arguments (ignored).
#' @export
#' @method print quaidsce
print.quaidsce <- function(x, ...) {
  cat("\nCensored QUAIDS Demand System Estimation (rquaidsce)\n")
  cat(paste(rep("-", 65), collapse = ""), "\n")
  cat(sprintf("Number of observations : %d\n", x$nobs))
  cat(sprintf("Number of equations    : %d (%s)\n", length(x$shares), paste(x$shares, collapse = ", ")))
  cat(sprintf("Demographics included  : %d (%s)\n", x$ndemo, ifelse(x$ndemo > 0, paste(x$demographics, collapse = ", "), "none")))
  cat(sprintf("Estimation method      : %s\n", toupper(x$method)))
  cat(sprintf("Specification          : %s\n", ifelse(x$quadratic, "QUAIDS (Quadratic)", "AIDS (Linear)")))
  cat(sprintf("Censoring correction   : %s\n", ifelse(x$censor, sprintf("Shonkwiler-Yen (first_stage_predict = '%s')", x$first_stage_predict), "None (Uncensored)")))
  cat(sprintf("Log-likelihood         : %.4f\n", x$llf))
  cat(sprintf("Alpha_0 (anot)         : %.4f\n", x$anot))
  cat(sprintf("Convergence status     : %s (outer: %d, inner: %d)\n", ifelse(x$converged, "Converged", "NOT converged"), x$n_outer, x$n_gn))

  if (!is.null(x$bootstrap)) {
    cat(sprintf("Bootstrap replications : %d / %d successful\n", x$bootstrap$reps_ok, x$bootstrap$reps_requested))
  }
  cat("\nCoefficients (first 8 displayed):\n")
  disp_k <- min(8, length(x$coefficients))
  print(round(x$coefficients[seq_len(disp_k)], 5))
  if (length(x$coefficients) > disp_k) {
    cat(sprintf("... and %d more coefficients (use summary() or coef() to view all)\n", length(x$coefficients) - disp_k))
  }
  invisible(x)
}

#' Summary Method for quaidsce
#' @param object A \code{quaidsce} object.
#' @param ... Additional arguments (ignored).
#' @return An object of class \code{summary.quaidsce}.
#' @export
#' @method summary quaidsce
summary.quaidsce <- function(object, ...) {
  b <- object$coefficients
  V <- object$vcov

  se <- sqrt(diag(V))
  z_stat <- b / se
  p_val <- 2 * (1 - stats::pnorm(abs(z_stat)))

  coef_table <- cbind(
    Estimate = b,
    `Std. Error` = se,
    `z value` = z_stat,
    `Pr(>|z|)` = p_val
  )

  res <- list(
    call = object$call,
    coefficients = coef_table,
    elasticities = object$elasticities,
    nobs = object$nobs,
    llf = object$llf,
    anot = object$anot,
    ndemo = object$ndemo,
    converged = object$converged,
    shares = object$shares,
    method = object$method,
    quadratic = object$quadratic,
    censor = object$censor,
    first_stage_predict = object$first_stage_predict,
    bootstrap = object$bootstrap
  )
  class(res) <- "summary.quaidsce"
  return(res)
}

#' Print Method for summary.quaidsce
#' @param x A \code{summary.quaidsce} object.
#' @param digits Number of digits to print.
#' @param ... Additional arguments.
#' @export
#' @method print summary.quaidsce
print.summary.quaidsce <- function(x, digits = 4, ...) {
  cat("\n================================================================================\n")
  cat("             CENSORED QUAIDS DEMAND SYSTEM ESTIMATION RESULTS                   \n")
  cat("================================================================================\n")
  cat(sprintf("Number of obs      = %10d      Log-likelihood     = %12.4f\n", x$nobs, x$llf))
  cat(sprintf("Demographics       = %10d      Alpha_0            = %12.4f\n", x$ndemo, x$anot))
  cat(sprintf("Model Type         = %10s      Censoring          = %12s\n",
              ifelse(x$quadratic, "QUAIDS", "AIDS"),
              ifelse(x$censor, sprintf("SY (%s)", x$first_stage_predict), "None")))
  cat(sprintf("Estimation Method  = %10s      Convergence        = %12s\n",
              toupper(x$method),
              ifelse(x$converged, "Converged", "FAILED")))

  if (!is.null(x$bootstrap)) {
    cat(sprintf("Bootstrap Replications = %d / %d (standard errors from bootstrap)\n", x$bootstrap$reps_ok, x$bootstrap$reps_requested))
  } else {
    cat("Standard Errors        = Analytical conditional VCE\n")
  }

  cat("\n--------------------------------------------------------------------------------\n")
  cat("Parameter Estimates:\n")
  cat("--------------------------------------------------------------------------------\n")
  stats::printCoefmat(x$coefficients, digits = digits, P.values = TRUE, has.Pvalue = TRUE)

  cat("\n--------------------------------------------------------------------------------\n")
  cat("Expenditure (Income) Elasticities (at sample means):\n")
  cat("--------------------------------------------------------------------------------\n")
  print(round(x$elasticities$income, digits))

  cat("\n--------------------------------------------------------------------------------\n")
  cat("Uncompensated (Marshallian) Price Elasticities [Row: Good, Column: Price]:\n")
  cat("--------------------------------------------------------------------------------\n")
  print(round(x$elasticities$uncompensated, digits))

  cat("\n--------------------------------------------------------------------------------\n")
  cat("Compensated (Hicksian) Price Elasticities [Row: Good, Column: Price]:\n")
  cat("--------------------------------------------------------------------------------\n")
  print(round(x$elasticities$compensated, digits))
  cat("================================================================================\n")
  invisible(x)
}
