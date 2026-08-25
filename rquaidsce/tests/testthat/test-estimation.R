test_that("quaidsce basic estimation works with standard arguments", {
  set.seed(42)
  N <- 500
  p1 <- exp(rnorm(N, 0, 0.2))
  p2 <- exp(rnorm(N, 0, 0.2))
  p3 <- exp(rnorm(N, 0, 0.2))
  p4 <- exp(rnorm(N, 0, 0.2))
  exp_tot <- exp(rnorm(N, 5, 0.5))
  hh_size <- sample(1:6, N, replace = TRUE)
  urban <- sample(0:1, N, replace = TRUE)

  # Synthetic shares with valid proportions
  w1 <- pmax(0, 0.3 + 0.05 * log(p1) - 0.02 * log(p2) + rnorm(N, 0, 0.03))
  w2 <- pmax(0, 0.25 - 0.02 * log(p1) + 0.04 * log(p2) + rnorm(N, 0, 0.03))
  w3 <- pmax(0, 0.25 + 0.01 * log(p3) + rnorm(N, 0, 0.03))
  w4 <- pmax(0, 0.20 + 0.01 * log(p4) + rnorm(N, 0, 0.03))
  
  # Induce censoring in all shares
  w1[sample(N, 30)] <- 0
  w2[sample(N, 30)] <- 0
  w3[sample(N, 30)] <- 0
  w4[sample(N, 30)] <- 0
  
  tot_w <- pmax(w1 + w2 + w3 + w4, 0.05)
  w1 <- w1 / tot_w
  w2 <- w2 / tot_w
  w3 <- w3 / tot_w
  w4 <- pmax(0, 1 - (w1 + w2 + w3))

  df <- data.frame(
    w1 = w1, w2 = w2, w3 = w3, w4 = w4,
    p1 = p1, p2 = p2, p3 = p3, p4 = p4,
    exp_tot = exp_tot,
    hh_size = hh_size,
    urban = urban
  )

  fit <- quaidsce(
    data = df,
    shares = c("w1", "w2", "w3", "w4"),
    prices = c("p1", "p2", "p3", "p4"),
    expenditure = "exp_tot",
    demographics = c("hh_size", "urban"),
    anot = 4.0,
    first_stage_predict = "xb",
    method = "fgnls",
    verbose = FALSE
  )

  expect_s3_class(fit, "quaidsce")
  expect_true(fit$converged)
  expect_equal(fit$nobs, N)
  expect_equal(length(fit$shares), 4)
  expect_true(is.numeric(fit$coefficients))
  expect_true(is.matrix(fit$vcov))
  expect_true(is.matrix(fit$sigma))
  expect_equal(length(fit$elasticities$income), 4)
  expect_equal(dim(fit$elasticities$uncompensated), c(4, 4))
  expect_equal(dim(fit$elasticities$compensated), c(4, 4))
})

test_that("input validation catches invalid specifications", {
  df <- data.frame(w1 = c(0.5, 0.5), w2 = c(0.5, 0.5), p1 = c(1, 1), p2 = c(1, 1), exp = c(10, 10))

  # Less than 3 shares
  expect_error(quaidsce(df, shares = c("w1", "w2"), prices = c("p1", "p2"), expenditure = "exp"),
               "At least 3 expenditure share variables")

  # Mismatched prices and shares
  expect_error(quaidsce(df, shares = c("w1", "w2", "w3"), prices = c("p1", "p2"), expenditure = "exp"),
               "Number of prices")

  # Missing column in data
  expect_error(quaidsce(df, shares = c("w1", "w2", "w3"), prices = c("p1", "p2", "p3"), expenditure = "nonexistent"),
               "were not found in data")
})

test_that("ivexp is forwarded and reduced-form diagnostics are returned", {
  set.seed(2026)
  N <- 350
  z <- rnorm(N)
  instrument <- rnorm(N)
  p1 <- exp(rnorm(N, 0, 0.15))
  p2 <- exp(rnorm(N, 0, 0.15))
  p3 <- exp(rnorm(N, 0, 0.15))
  exp_tot <- exp(4 + 0.7 * instrument + 0.1 * z + rnorm(N, 0, 0.25))
  raw <- cbind(
    pmax(0, 0.45 + rnorm(N, 0, 0.08)),
    pmax(0, 0.35 + rnorm(N, 0, 0.08)),
    pmax(0, 0.20 + rnorm(N, 0, 0.08))
  )
  raw[sample(N, 25), 1] <- 0
  raw[sample(N, 25), 2] <- 0
  raw[sample(N, 25), 3] <- 0
  shares <- raw / rowSums(raw)
  df <- data.frame(
    w1 = shares[, 1], w2 = shares[, 2], w3 = shares[, 3],
    p1 = p1, p2 = p2, p3 = p3, exp_tot = exp_tot,
    z = z, instrument = instrument
  )

  fit <- quaidsce(
    df,
    shares = c("w1", "w2", "w3"),
    prices = c("p1", "p2", "p3"),
    expenditure = "exp_tot",
    demographics = "z",
    ivexp = "instrument",
    anot = 4,
    method = "nls",
    first_stage_predict = "xb",
    verbose = FALSE
  )
  expect_equal(fit$ivexp, "instrument")
  expect_true(is.list(fit$reduced_form))
  expect_true(is.numeric(fit$reduced_form$coefficients))
  expect_true(is.matrix(fit$reduced_form$vcov))
  expect_gt(fit$reduced_form$excluded.f, 10)
  expect_equal(length(fit$reduced_form$residuals), N)
})
