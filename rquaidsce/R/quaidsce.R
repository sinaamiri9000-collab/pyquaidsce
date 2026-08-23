#' Estimate Censored QUAIDS Demand System
#'
#' Estimates a Censored Quadratic Almost Ideal Demand System (QUAIDS) with Ray (1983)
#' demographic scaling and Shonkwiler & Yen (1999) two-step censoring correction,
#' powered by the \code{pyquaidsce} Python backend.
#'
#' @param data A data.frame containing all required variables.
#' @param shares Character vector of expenditure share variable names (at least 3).
#' @param prices Character vector of price variable names.
#' @param lnprices Character vector of log-price variable names (alternative to prices).
#' @param expenditure Character name of total expenditure variable.
#' @param lnexpenditure Character name of log-expenditure variable.
#' @param demographics Optional character vector of demographic variable names.
#' @param anot Numeric scalar normalization for alpha_0 (default = 10.0).
#' @param quadratic Logical, whether to include quadratic expenditure term (default = TRUE).
#' @param censor Logical, whether to apply Shonkwiler-Yen censoring correction (default = TRUE).
#' @param method Estimation method: 'ifgnls' (default), 'fgnls', or 'nls'.
#' @param algorithm Numerical optimizer: 'gn' (Gauss-Newton, default) or 'lm' (Levenberg-Marquardt).
#' @param start Starting values: 'zero' (default) or 'linear'.
#' @param initial Optional numeric vector of custom starting free parameters (or a previous \code{quaidsce} object).
#' @param sigma_initial Optional numeric matrix of custom starting residual covariance Sigma (dimension (n_goods - 1) x (n_goods - 1)).
#' @param stop_rule Convergence rule: 'standard' (default) or 'tight'.
#' @param bootstrap_start Bootstrap starting scheme: 'zero' (default) or 'warm'.
#' @param first_stage_predict First-stage probit predictor: 'xb' (default, theoretical linear index) or 'pr' (legacy Stata-compatible prediction).
#' @param strict_stata Logical, whether to match Stata numerical quirks (default = FALSE).
#' @param vce_sigma Residual covariance for analytical VCE: 'objective' (default) or 'final'.
#' @param control_function Optional character name of endogeneity control function residual.
#' @param selection_control_function Optional character name of selection-stage control function.
#' @param selection_prices Optional character vector specifying subset of prices in selection probits.
#' @param selection_covariates Optional character vector specifying subset of covariates in selection probits.
#' @param selection_expenditure Logical, whether expenditure enters selection probits (default = TRUE).
#' @param reps Integer, number of bootstrap replications (default = 0).
#' @param seed Optional integer, random number seed for bootstrap.
#' @param n_jobs Integer, number of parallel bootstrap workers (default = 1).
#' @param mp_context Multiprocessing context method (default = "spawn").
#' @param rep_timeout Optional numeric timeout per bootstrap replication in seconds.
#' @param tol Optional outer convergence tolerance (default = 1e-13). If NULL, default tolerance is used.
#' @param max_outer Maximum outer iterations (default = 200).
#' @param max_iter Maximum inner iterations (default = 300).
#' @param chunk Chunk size for memory management (default = 2000).
#' @param verbose Logical, whether to print progress (default = FALSE).
#' @return An object of class \code{quaidsce}.
#' @export
#' @examples
#' \dontrun{
#' # Assuming df is a data.frame with budget shares w1..w4 and prices p1..p4
#' fit <- quaidsce(
#'   data = df,
#'   shares = c("w1", "w2", "w3", "w4"),
#'   prices = c("p1", "p2", "p3", "p4"),
#'   expenditure = "exp_total",
#'   demographics = c("hh_size", "urban"),
#'   anot = 10.0
#' )
#' summary(fit)
#' elasticities(fit)
#' }
quaidsce <- function(data,
                     shares,
                     prices = NULL,
                     lnprices = NULL,
                     expenditure = NULL,
                     lnexpenditure = NULL,
                     demographics = NULL,
                     anot = 10.0,
                     quadratic = TRUE,
                     censor = TRUE,
                     method = "ifgnls",
                     algorithm = "gn",
                     start = "zero",
                     initial = NULL,
                     sigma_initial = NULL,
                     stop_rule = "standard",
                     bootstrap_start = "zero",
                     first_stage_predict = "xb",
                     strict_stata = FALSE,
                     vce_sigma = "objective",
                     control_function = NULL,
                     selection_control_function = NULL,
                     selection_prices = NULL,
                     selection_covariates = NULL,
                     selection_expenditure = TRUE,
                     reps = 0L,
                     seed = NULL,
                     n_jobs = 1L,
                     mp_context = "spawn",
                     rep_timeout = NULL,
                     tol = NULL,
                     max_outer = 200L,
                     max_iter = 300L,
                     chunk = 2000L,
                     verbose = FALSE) {

  # 1. Validate data.frame and variable inputs
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (length(shares) < 3) {
    stop("At least 3 expenditure share variables must be specified.", call. = FALSE)
  }

  if (is.null(prices) && is.null(lnprices)) {
    stop("Must specify either 'prices' or 'lnprices'.", call. = FALSE)
  }
  if (!is.null(prices) && !is.null(lnprices)) {
    stop("Cannot specify both 'prices' and 'lnprices'.", call. = FALSE)
  }

  p_vars <- if (!is.null(prices)) prices else lnprices
  if (length(p_vars) != length(shares)) {
    stop(sprintf("Number of prices (%d) must equal number of shares (%d).", length(p_vars), length(shares)), call. = FALSE)
  }

  if (is.null(expenditure) && is.null(lnexpenditure)) {
    stop("Must specify either 'expenditure' or 'lnexpenditure'.", call. = FALSE)
  }
  if (!is.null(expenditure) && !is.null(lnexpenditure)) {
    stop("Cannot specify both 'expenditure' and 'lnexpenditure'.", call. = FALSE)
  }
  exp_var <- if (!is.null(expenditure)) expenditure else lnexpenditure

  all_vars <- unique(c(shares, p_vars, exp_var, demographics, control_function, selection_control_function, selection_prices, selection_covariates))
  missing_cols <- setdiff(all_vars, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("The following specified variables were not found in data: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  # Coerce factors/characters to numeric if necessary in the required subset
  clean_data <- data[, all_vars, drop = FALSE]
  for (col in names(clean_data)) {
    if (is.factor(clean_data[[col]]) || is.character(clean_data[[col]])) {
      clean_data[[col]] <- as.numeric(clean_data[[col]])
    }
  }

  # Handle initial argument if passed as a quaidsce object, list, or numeric vector
  init_vec <- NULL
  if (!is.null(initial)) {
    if (inherits(initial, "quaidsce") || is.list(initial)) {
      if (!is.null(initial$theta)) {
        init_vec <- as.numeric(initial$theta)
      } else if (!is.null(initial$py_object) && !is.null(initial$py_object$theta)) {
        init_vec <- as.numeric(initial$py_object$theta)
      } else {
        stop("The provided object does not contain structural parameters '$theta'.", call. = FALSE)
      }
    } else if (is.numeric(initial)) {
      init_vec <- as.numeric(initial)
    } else {
      stop("'initial' must be a numeric vector of free parameters or a fitted quaidsce object.", call. = FALSE)
    }
  }

  # 2. Get pyquaidsce backend
  py_pkg <- .get_pyquaidsce()

  # 3. Prepare arguments for Python estimator
  args <- list(
    data = clean_data,
    shares = as.list(shares),
    prices = if (!is.null(prices)) as.list(prices) else NULL,
    lnprices = if (!is.null(lnprices)) as.list(lnprices) else NULL,
    expenditure = expenditure,
    lnexpenditure = lnexpenditure,
    demographics = if (!is.null(demographics)) as.list(demographics) else NULL,
    anot = as.numeric(anot),
    quadratic = as.logical(quadratic),
    censor = as.logical(censor),
    method = as.character(method),
    algorithm = as.character(algorithm),
    start = as.character(start),
    initial = init_vec,
    sigma_initial = if (!is.null(sigma_initial)) as.matrix(sigma_initial) else NULL,
    stop_rule = as.character(stop_rule),
    bootstrap_start = as.character(bootstrap_start),
    first_stage_predict = as.character(first_stage_predict),
    strict_stata = as.logical(strict_stata),
    vce_sigma = as.character(vce_sigma),
    control_function = control_function,
    selection_control_function = selection_control_function,
    selection_prices = if (!is.null(selection_prices)) as.list(selection_prices) else NULL,
    selection_covariates = if (!is.null(selection_covariates)) as.list(selection_covariates) else NULL,
    selection_expenditure = as.logical(selection_expenditure),
    reps = as.integer(reps),
    seed = if (!is.null(seed)) as.integer(seed) else NULL,
    n_jobs = as.integer(n_jobs),
    mp_context = if (!is.null(mp_context)) as.character(mp_context) else NULL,
    rep_timeout = if (!is.null(rep_timeout)) as.numeric(rep_timeout) else NULL,
    tol = if (!is.null(tol)) as.numeric(tol) else 1e-13,
    max_outer = as.integer(max_outer),
    max_iter = as.integer(max_iter),
    chunk = as.integer(chunk),
    verbose = as.logical(verbose)
  )

  # 4. Call Python estimator
  res_py <- do.call(py_pkg$quaidsce, args)

  # 5. Extract and convert results to native R objects
  b_vec <- as.numeric(res_py$b)
  names(b_vec) <- as.character(res_py$names)

  V_mat <- as.matrix(res_py$V)
  rownames(V_mat) <- names(b_vec)
  colnames(V_mat) <- names(b_vec)

  share_names <- as.character(res_py$share_names)

  # Elasticities
  elas_inc <- as.numeric(res_py$elas$income)
  names(elas_inc) <- share_names

  elas_uncomp <- as.matrix(res_py$elas$uncompensated)
  rownames(elas_uncomp) <- share_names
  colnames(elas_uncomp) <- share_names

  elas_comp <- as.matrix(res_py$elas$compensated)
  rownames(elas_comp) <- share_names
  colnames(elas_comp) <- share_names

  # Extract fitted shares & residuals if available
  fitted_mat <- tryCatch(as.matrix(res_py$fitted_shares()), error = function(e) NULL)
  resid_mat <- NULL
  if (!is.null(fitted_mat)) {
    rownames(fitted_mat) <- seq_len(nrow(fitted_mat))
    colnames(fitted_mat) <- share_names
    w_actual <- as.matrix(clean_data[, share_names, drop = FALSE])
    resid_mat <- w_actual - fitted_mat
  }

  # Bootstrap info
  boot_info <- NULL
  if (!is.null(res_py$boot)) {
    boot_info <- list(
      reps_requested = as.integer(res_py$boot$reps_requested),
      reps_ok = as.integer(res_py$boot$reps_ok),
      se = as.numeric(res_py$boot$se),
      failures = as.character(res_py$boot$failures),
      b_star = tryCatch(as.matrix(res_py$boot$b_star), error = function(e) NULL),
      V = tryCatch(as.matrix(res_py$boot$V), error = function(e) NULL)
    )
    names(boot_info$se) <- names(b_vec)
  }

  # Build return object
  out <- list(
    coefficients = b_vec,
    theta = as.numeric(res_py$theta),
    vcov = V_mat,
    sigma = as.matrix(res_py$sigma),
    elasticities = list(
      income = elas_inc,
      uncompensated = elas_uncomp,
      compensated = elas_comp
    ),
    fitted.values = fitted_mat,
    residuals = resid_mat,
    nobs = as.integer(res_py$nobs),
    llf = as.numeric(res_py$llf),
    anot = as.numeric(res_py$anot),
    ndemo = as.integer(res_py$spec$ndemo),
    converged = as.logical(res_py$converged),
    n_outer = as.integer(res_py$n_outer),
    n_gn = as.integer(res_py$n_gn),
    shares = share_names,
    prices = as.character(res_py$price_names),
    demographics = as.character(res_py$demo_names),
    method = as.character(method),
    quadratic = as.logical(quadratic),
    censor = as.logical(censor),
    first_stage_predict = as.character(first_stage_predict),
    bootstrap = boot_info,
    call = match.call(),
    py_object = res_py
  )

  class(out) <- "quaidsce"
  return(out)
}
