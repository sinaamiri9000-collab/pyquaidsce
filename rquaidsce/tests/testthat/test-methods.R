test_that("S3 methods work as expected", {
  set.seed(123)
  N <- 400
  p1 <- exp(rnorm(N, 0, 0.1))
  p2 <- exp(rnorm(N, 0, 0.1))
  p3 <- exp(rnorm(N, 0, 0.1))
  exp_tot <- exp(rnorm(N, 5, 0.3))
  demo1 <- rnorm(N)

  w1 <- pmax(0, 0.4 + 0.02 * log(p1) + rnorm(N, 0, 0.05))
  w2 <- pmax(0, 0.35 - 0.02 * log(p2) + rnorm(N, 0, 0.05))
  w3 <- pmax(0, 0.25 + rnorm(N, 0, 0.05))
  w1[sample(N, 30)] <- 0
  w2[sample(N, 30)] <- 0
  w3[sample(N, 30)] <- 0
  tot <- pmax(w1 + w2 + w3, 0.01)
  w1 <- w1 / tot
  w2 <- w2 / tot
  w3 <- 1 - (w1 + w2)

  df <- data.frame(w1 = w1, w2 = w2, w3 = w3, p1 = p1, p2 = p2, p3 = p3, exp_tot = exp_tot, demo1 = demo1)

  fit <- quaidsce(df, shares = c("w1", "w2", "w3"), prices = c("p1", "p2", "p3"), expenditure = "exp_tot", demographics = "demo1", anot = 4.0, method = "fgnls", verbose = FALSE)

  # coef
  cf <- coef(fit)
  expect_true(is.numeric(cf))
  expect_equal(cf, fit$coefficients)

  # vcov
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(length(cf), length(cf)))

  # sigma
  S <- sigma(fit)
  expect_true(is.matrix(S))

  # elasticities extractor
  el_all <- elasticities(fit, type = "all")
  expect_type(el_all, "list")
  expect_named(el_all, c("income", "uncompensated", "compensated"))

  el_inc <- elasticities(fit, type = "income")
  expect_equal(length(el_inc), 3)

  el_uncomp <- elasticities(fit, type = "uncompensated")
  expect_equal(dim(el_uncomp), c(3, 3))

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.quaidsce")
  expect_true(is.matrix(s$coefficients))
  expect_equal(colnames(s$coefficients), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

  # print methods should not fail
  expect_output(print(fit))
  expect_output(print(s))
})
