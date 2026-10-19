uncensored_iv_data <- function() {
  set.seed(918)
  n <- 240
  lp <- matrix(rnorm(n * 3, 0, 0.15), ncol = 3)
  z <- rnorm(n, 0, 0.15)
  instrument <- rnorm(n, 0, 0.4)
  residual <- rnorm(n, 0, 0.15)
  lnexp <- 2 + drop(lp %*% c(0.12, -0.08, 0.06)) +
    0.1 * z + 0.5 * instrument + residual
  alpha <- c(0.32, 0.30, 0.38)
  beta <- c(0.015, -0.01, -0.005)
  control <- c(0.08, -0.05, -0.03)
  real_exp <- lnexp - 1.6 - drop(lp %*% alpha)
  shares <- matrix(alpha, n, 3, byrow = TRUE) +
    outer(real_exp, beta) + outer(residual, control)
  error <- matrix(rnorm(n * 2, 0, 0.001), ncol = 2)
  shares <- shares + cbind(error, -rowSums(error))
  data.frame(w1 = shares[, 1], w2 = shares[, 2], w3 = shares[, 3],
             p1 = exp(lp[, 1]), p2 = exp(lp[, 2]), p3 = exp(lp[, 3]),
             expenditure = exp(lnexp), z = z, instrument = instrument)
}

uncensored_iv_args <- function(data, quadratic = TRUE, reps = 0L) {
  list(data = data, shares = c("w1", "w2", "w3"),
       prices = c("p1", "p2", "p3"), expenditure = "expenditure",
       demographics = "z", ivexp = "instrument", censor = FALSE,
       quadratic = quadratic, anot = 1.6, start = "linear", method = "ifgnls",
       reps = reps, seed = 573L, n_jobs = 1L, bootstrap_start = "warm",
       max_outer = 200L, verbose = FALSE)
}

uncensored_iv_python_fit <- function(args) {
  py <- reticulate::import("pyquaidsce")
  for (name in c("shares", "prices", "demographics", "ivexp")) {
    args[[name]] <- as.list(args[[name]])
  }
  do.call(py$quaidsce, args)
}

test_that("uncensored ivexp is passed to Python for AIDS and QUAIDS", {
  data <- uncensored_iv_data()
  for (quadratic in c(FALSE, TRUE)) {
    args <- uncensored_iv_args(data, quadratic)
    fit <- do.call(quaidsce, args)
    direct <- uncensored_iv_python_fit(args)
    expect_true(fit$converged)
    expect_false(fit$censor)
    expect_equal(fit$ivexp, "instrument")
    expect_equal(as.numeric(coef(fit)), as.numeric(direct$b), tolerance = 1e-10)
    expect_equal(unname(vcov(fit)), unname(as.matrix(direct$V)), tolerance = 1e-10)
    expect_equal(names(fit$standard.errors), names(coef(fit)))
    direct_se <- as.numeric(direct$se)
    direct_se[!is.finite(direct_se)] <- NA_real_
    expect_equal(as.numeric(fit$standard.errors), direct_se, tolerance = 1e-10)
    expect_equal(fit$reduced_form$residuals, as.numeric(direct$reduced_form$residuals), tolerance = 1e-10)
    expect_true(is.finite(fit$reduced_form$excluded.f))
    controls <- grep("^cfcoef:", names(coef(fit)))
    expect_length(controls, 3)
    expect_equal(sum(coef(fit)[controls]), 0, tolerance = 1e-12)
    expect_false(any(grepl("^(tau|delta):", names(coef(fit)))))
    tab <- summary(fit)$coefficients
    elasticities <- grepl("^ELAS_", rownames(tab))
    expect_equal(sum(elasticities), 21)
    expect_true(all(is.na(tab[elasticities, "Std. Error"])))
    expect_true(all(is.na(tab[elasticities, "z value"])))
    expect_true(all(is.na(tab[elasticities, "Pr(>|z|)"])))
    expect_true(all(is.finite(tab[controls, "Std. Error"])))
    expect_equal(unname(fit$elasticities$uncompensated),
                 unname(as.matrix(direct$elas$uncompensated)), tolerance = 1e-10)
  }
})

test_that("uncensored ivexp bootstrap SE and inference match Python", {
  args <- uncensored_iv_args(uncensored_iv_data(), reps = 4L)
  fit <- do.call(quaidsce, args)
  direct <- uncensored_iv_python_fit(args)
  expect_true(fit$converged)
  expect_equal(fit$bootstrap$reps_ok, 4L)
  expect_true(all(is.finite(fit$standard.errors)))
  expect_equal(as.numeric(fit$standard.errors), as.numeric(direct$se), tolerance = 1e-10)
  expect_equal(unname(vcov(fit)), unname(as.matrix(direct$V)), tolerance = 1e-10)
  tab <- summary(fit)$coefficients
  expect_equal(unname(tab[, "Std. Error"]), as.numeric(direct$se), tolerance = 1e-10)
  expect_equal(unname(tab[, "z value"]), as.numeric(direct$b / direct$se), tolerance = 1e-10)
  expect_equal(unname(tab[, "Pr(>|z|)"]),
               unname(2 * pnorm(abs(tab[, "z value"]), lower.tail = FALSE)), tolerance = 1e-12)
  old <- fit
  old$standard.errors <- NULL
  expect_equal(summary(old)$coefficients, tab)
})
