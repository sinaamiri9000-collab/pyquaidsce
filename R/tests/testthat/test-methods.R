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
  elasticity_rows <- grepl("^ELAS_", rownames(s$coefficients))
  expect_true(all(is.na(s$coefficients[elasticity_rows, "Std. Error"])))
  expect_true(all(is.na(s$coefficients[elasticity_rows, "Pr(>|z|)"])))

  # Independently connect named coefficients to [good, price] matrix cells.
  for (i in seq_len(3)) {
    for (j in seq_len(3)) {
      for (block in c("ELAS_UNCOMP", "ELAS_COMP")) {
        type <- if (block == "ELAS_UNCOMP") "uncompensated" else "compensated"
        key <- sprintf("%s:e_%d_%d", block, i, j)
        expected <- elasticities(fit, type = type)[i, j]
        expect_equal(unname(cf[key]), unname(expected), tolerance = 1e-12)
        expect_equal(unname(s$coefficients[key, "Estimate"]),
                     unname(expected), tolerance = 1e-12)
      }
    }
  }

  # print methods should not fail
  expect_output(print(fit))
  expect_output(print(s))
})

summary_fixture <- function() {
  b <- c("beta:beta_1" = 0.2, "ELAS_INC:e_1" = 1.2)
  V <- diag(c(0.01, 0))
  dimnames(V) <- list(names(b), names(b))
  structure(list(coefficients = b, vcov = V,
                 standard.errors = setNames(c(0.1, NA_real_), names(b))),
            class = "quaidsce")
}

test_that("summary keeps undefined elasticity inference missing", {
  fit <- summary_fixture()
  before <- vcov(fit)
  tab <- summary(fit)$coefficients
  expect_equal(unname(tab[1, "Std. Error"]), 0.1)
  expect_equal(unname(tab[1, "z value"]), 2)
  expect_equal(unname(tab[1, "Pr(>|z|)"]), 0.0455002638963584)
  expect_true(all(is.na(tab[2, c("Std. Error", "z value", "Pr(>|z|)")])))
  expect_equal(vcov(fit), before)
  expect_equal(coef(fit), c("beta:beta_1" = 0.2, "ELAS_INC:e_1" = 1.2))
})

test_that("summary handles objects saved before standard.errors was added", {
  fit <- summary_fixture()
  fit$standard.errors <- NULL
  tab <- summary(fit)$coefficients
  expect_equal(unname(tab[1, "Std. Error"]), 0.1)
  expect_true(all(is.na(tab[2, c("Std. Error", "z value", "Pr(>|z|)")])))
})

test_that("zero or nonfinite SE do not produce spurious z and p values", {
  fit <- summary_fixture()
  for (se in c(0, NaN, Inf)) {
    fit$standard.errors[2] <- se
    tab <- summary(fit)$coefficients
    expect_true(is.na(tab[2, "z value"]))
    expect_true(is.na(tab[2, "Pr(>|z|)"]))
  }
})
