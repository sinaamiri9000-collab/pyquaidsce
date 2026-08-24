# ==============================================================================
# RQUAIDSCE 1.4.0 MASTER COMPREHENSIVE 24-SCENARIO TEST SUITE
# ==============================================================================
# Dataset: bd_uruguay.csv (14 food groups, 6,848 observations)
# ==============================================================================
# Instructions:
# Open in RStudio and press Ctrl+A (Select All) then Ctrl+Enter (Run).
# ==============================================================================

suppressPackageStartupMessages({
  library(rquaidsce)
})

cat("\n================================================================================\n")
cat("RQUAIDSCE 1.4.0 MASTER COMPREHENSIVE 24-SCENARIO TEST SUITE\n")
cat("================================================================================\n\n")

# ---- STEP 0: Load dataset ----------------------------------------------------
data_candidates <- c(
  "C:/Users/sina/Downloads/bd_uruguay.csv",
  "bd_uruguay.csv",
  "../benchmarks/bd_uruguay.csv"
)

df <- NULL
for (p in data_candidates) {
  if (file.exists(p)) {
    df <- read.csv(p)
    cat(sprintf("Loaded dataset from: %s (N = %d observations)\n", p, nrow(df)))
    break
  }
}

if (is.null(df)) {
  stop("Could not find bd_uruguay.csv in candidate paths.")
}

# Clean shares & build log variables
for (i in 1:14) {
  col <- paste0("w", i, "_red")
  if (col %in% names(df)) {
    df[[col]] <- pmax(df[[col]], 0.0)
  }
  p_col <- paste0("P_med", i)
  if (p_col %in% names(df)) {
    df[[paste0("ln_", p_col)]] <- log(df[[p_col]])
  }
}

if (!("ln_gasto_total" %in% names(df))) {
  df$ln_gasto_total <- log(df$gasto_total)
}

# Generate endogeneity residual proxy for test 18
ols_mod <- lm(ln_gasto_total ~ npersonas + edad + Sex + log_ing, data = df)
df$vhat <- residuals(ols_mod)

shares_14 <- paste0("w", 1:14, "_red")
prices_14 <- paste0("P_med", 1:14)
lnprices_14 <- paste0("ln_P_med", 1:14)
demographics <- c("npersonas", "edad", "Sex", "prim_comp", "sec_comp", "sup_comp", "log_ing")
anot_val <- 2.9766226

results <- list()

record_test <- function(test_id, title, expr_func) {
  cat(sprintf("\n[%02d/24] RUNNING: %s ...\n", test_id, title))
  t0 <- proc.time()
  res <- tryCatch({
    fit <- expr_func()
    el <- (proc.time() - t0)["elapsed"]
    llf <- if (!is.null(fit$llf)) fit$llf else fit$ll
    conv <- if (!is.null(fit$converged)) fit$converged else TRUE
    st <- if (conv) "PASS" else "WARN"
    msg <- sprintf("LLF = %.4f, Time = %.2fs", llf, el)
    cat(sprintf("       --> [%s] %s\n", st, msg))
    results[[length(results) + 1]] <<- list(id = test_id, title = title, status = st, time = el, details = msg)
    fit
  }, error = function(e) {
    el <- (proc.time() - t0)["elapsed"]
    msg <- conditionMessage(e)
    cat(sprintf("       --> [FAIL] %s\n", msg))
    results[[length(results) + 1]] <<- list(id = test_id, title = title, status = "FAIL", time = el, details = msg)
    NULL
  })
  return(res)
}

# ------------------------------------------------------------------------------
# Tier 1: Model Specification & Variable Formats
# ------------------------------------------------------------------------------
t1_fit <- record_test(1, "Baseline 14-Good IFGNLS (first_stage_predict='xb')", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, method = "ifgnls", first_stage_predict = "xb")
})

record_test(2, "Direct Log-Prices & Log-Expenditure (lnprices & lnexpenditure)", function() {
  quaidsce(data = df, shares = shares_14, lnprices = lnprices_14, lnexpenditure = "ln_gasto_total",
           demographics = demographics, anot = anot_val, method = "ifgnls")
})

record_test(3, "Linear AIDS Model (quadratic=FALSE)", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, quadratic = FALSE, method = "ifgnls")
})

record_test(4, "Uncensored QUAIDS (censor=FALSE, Poi 2012 specification)", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, censor = FALSE, method = "ifgnls")
})

record_test(5, "Translog Constant Variation (anot = 10.0)", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = 10.0, method = "ifgnls")
})

# 3-good subsystem
df_3g <- df
sub_shares <- c("w1_red", "w2_red", "w3_red")
tot_sub <- rowSums(df_3g[, sub_shares])
for (s in sub_shares) {
  df_3g[[paste0(s, "_sub")]] <- df_3g[[s]] / tot_sub
}
df_3g$sub_exp <- pmax(df_3g$gasto1 + df_3g$gasto2 + df_3g$gasto3, 1.0)

record_test(6, "3-Good Subsystem Estimation", function() {
  quaidsce(data = df_3g, shares = paste0(sub_shares, "_sub"),
           prices = c("P_med1", "P_med2", "P_med3"), expenditure = "sub_exp",
           demographics = demographics, anot = anot_val, method = "ifgnls")
})

# ------------------------------------------------------------------------------
# Tier 2: Solvers, Algorithms & Numerical Tolerances
# ------------------------------------------------------------------------------
t7_fit <- record_test(7, "Nonlinear Least Squares (method='nls')", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, method = "nls")
})

t8_fit <- record_test(8, "Feasible Generalized NLS (method='fgnls')", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, method = "fgnls")
})

record_test(9, "Iterated FGNLS (method='ifgnls')", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, method = "ifgnls")
})

record_test(10, "Levenberg-Marquardt Optimizer (algorithm='lm')", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, algorithm = "lm", method = "ifgnls")
})

record_test(11, "Strict Gradient Stopping Rule (stop_rule='tight')", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, stop_rule = "tight", method = "ifgnls")
})

record_test(12, "Alternative VCE Sigma Formula (vce_sigma='final')", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, vce_sigma = "final", method = "ifgnls")
})

record_test(13, "Chunk Size Variation (chunk=500)", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, chunk = 500L, method = "ifgnls")
})

# ------------------------------------------------------------------------------
# Tier 3: Censoring, Probit Specifications & Control Functions
# ------------------------------------------------------------------------------
record_test(14, "Legacy Stata Probit CDF Predictor (first_stage_predict='pr')", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, first_stage_predict = "pr", method = "ifgnls")
})

record_test(15, "Selection Price Subset (selection_prices=P_med1..3)", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val,
           selection_prices = c("P_med1", "P_med2", "P_med3"), method = "ifgnls")
})

record_test(16, "Selection Independent Covariates (npersonas, edad)", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val,
           selection_covariates = c("npersonas", "edad"), method = "ifgnls")
})

record_test(17, "Selection Omit Log Expenditure (selection_expenditure=FALSE)", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val,
           selection_expenditure = FALSE, method = "ifgnls")
})

record_test(18, "Endogeneity Control Function (control_function='vhat')", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val,
           control_function = "vhat", first_stage_predict = "xb", method = "ifgnls")
})

# ------------------------------------------------------------------------------
# Tier 4: Warm-Starting & Matrix Transfers
# ------------------------------------------------------------------------------
record_test(19, "Linearized AIDS Starting Values (start='linear')", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, start = "linear", method = "ifgnls")
})

theta_init <- if (!is.null(t1_fit)) t1_fit$theta else NULL
sigma_init <- if (!is.null(t1_fit)) t1_fit$sigma else NULL

record_test(20, "Custom Structural Free Parameters Vector (initial=theta_init)", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, initial = theta_init, method = "ifgnls")
})

record_test(21, "Custom Residual Covariance Matrix (sigma_initial=sigma_init)", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, initial = theta_init,
           sigma_initial = sigma_init, method = "ifgnls")
})

record_test(22, "Chained Multi-Stage Warm-Start (NLS theta/sigma -> IFGNLS)", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val,
           initial = if (!is.null(t7_fit)) t7_fit$theta else NULL,
           sigma_initial = if (!is.null(t7_fit)) t7_fit$sigma else NULL,
           method = "ifgnls")
})

# ------------------------------------------------------------------------------
# Tier 5: Bootstrap, Multiprocessing & Postestimation
# ------------------------------------------------------------------------------
record_test(23, "Parallel Bootstrap Standard Errors (reps=10, n_jobs=4, seed=123456)", function() {
  quaidsce(data = df, shares = shares_14, prices = prices_14, expenditure = "gasto_total",
           demographics = demographics, anot = anot_val, reps = 10L, n_jobs = 4L,
           seed = 123456L, mp_context = "spawn", method = "ifgnls")
})

record_test(24, "Complete S3 Methods (summary, coef, vcov, residuals, fitted, elasticities)", function() {
  if (is.null(t1_fit)) stop("Base fit t1_fit is missing.")
  sm <- summary(t1_fit)
  cf <- coef(t1_fit)
  vc <- vcov(t1_fit)
  rs <- residuals(t1_fit)
  ft <- fitted(t1_fit)
  sg <- sigma(t1_fit)
  el_inc <- t1_fit$elasticities$income
  el_unc <- t1_fit$elasticities$uncompensated
  stopifnot(length(cf) > 0, nrow(vc) == length(cf), nrow(rs) == nrow(df), nrow(ft) == nrow(df))
  stopifnot(length(el_inc) == 14, nrow(el_unc) == 14, ncol(el_unc) == 14)
  t1_fit
})

# ------------------------------------------------------------------------------
# Final Scorecard
# ------------------------------------------------------------------------------
cat("\n================================================================================\n")
cat("RQUAIDSCE 1.4.0 MASTER TEST SUITE SCORECARD\n")
cat("================================================================================\n")
cat(sprintf("%-3s %-52s %-8s %-10s %s\n", "#", "Test Scenario Title", "Status", "Time (s)", "Details"))
cat(paste(rep("-", 80), collapse = ""), "\n")
n_pass <- 0
total_time <- 0
for (item in results) {
  total_time <- total_time + item$time
  if (item$status == "PASS") n_pass <- n_pass + 1
  cat(sprintf("%-3d %-52s %-8s %6.2fs    %s\n", item$id, substr(item$title, 1, 50), item$status, item$time, substr(item$details, 1, 30)))
}
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("Summary: %d/24 tests PASSED successfully in %.2f seconds.\n", n_pass, total_time))
cat("================================================================================\n\n")
