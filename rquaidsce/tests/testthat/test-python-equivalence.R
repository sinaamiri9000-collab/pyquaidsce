test_that("R wrapper produces exact numerical match with direct Python estimator", {
  set.seed(99)
  N <- 300
  p1 <- exp(rnorm(N, 0, 0.1))
  p2 <- exp(rnorm(N, 0, 0.1))
  p3 <- exp(rnorm(N, 0, 0.1))
  exp_tot <- exp(rnorm(N, 4, 0.2))

  w1 <- pmax(0, 0.5 + 0.01 * log(p1) + rnorm(N, 0, 0.02))
  w2 <- pmax(0, 0.3 - 0.01 * log(p2) + rnorm(N, 0, 0.02))
  w3 <- pmax(0, 0.2 + rnorm(N, 0, 0.02))
  w1[sample(N, 20)] <- 0
  w2[sample(N, 20)] <- 0
  w3[sample(N, 20)] <- 0
  tot <- pmax(w1 + w2 + w3, 0.01)
  w1 <- w1 / tot
  w2 <- w2 / tot
  w3 <- pmax(0, 1 - (w1 + w2))

  df <- data.frame(w1 = w1, w2 = w2, w3 = w3, p1 = p1, p2 = p2, p3 = p3, exp_tot = exp_tot, z1 = rnorm(N))

  # Estimate via R package
  fit_r <- quaidsce(
    data = df,
    shares = c("w1", "w2", "w3"),
    prices = c("p1", "p2", "p3"),
    expenditure = "exp_tot",
    demographics = "z1",
    anot = 3.5,
    method = "fgnls",
    first_stage_predict = "xb",
    verbose = FALSE
  )

  # Direct Python call via reticulate
  py <- reticulate::import("pyquaidsce")
  res_py <- py$quaidsce(
    data = df,
    shares = as.list(c("w1", "w2", "w3")),
    prices = as.list(c("p1", "p2", "p3")),
    expenditure = "exp_tot",
    demographics = as.list("z1"),
    anot = 3.5,
    method = "fgnls",
    first_stage_predict = "xb",
    verbose = FALSE
  )

  # Cross-platform robust numerical equivalence assertions
  expect_equal(as.numeric(fit_r$coefficients), as.numeric(res_py$b), tolerance = 1e-10)
  expect_equal(as.numeric(fit_r$vcov), as.numeric(res_py$V), tolerance = 1e-10)
  expect_equal(as.numeric(fit_r$elasticities$income), as.numeric(res_py$elas$income), tolerance = 1e-10)
  expect_equal(as.numeric(fit_r$elasticities$uncompensated), as.numeric(res_py$elas$uncompensated), tolerance = 1e-10)
  expect_equal(as.numeric(fit_r$elasticities$compensated), as.numeric(res_py$elas$compensated), tolerance = 1e-10)
  expect_equal(as.numeric(fit_r$llf), as.numeric(res_py$llf), tolerance = 1e-8)
})
