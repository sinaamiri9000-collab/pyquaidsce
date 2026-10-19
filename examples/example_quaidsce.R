# Make sure both the R interface and Python package are installed.
# See R/README.md for installation instructions.
# Run from the repository root: source("examples/example_quaidsce.R")

library(pyquaidsce)

# Read the dataset from your checkout, or download it from the repository.
data_file <- "examples/data/household_demand.csv"
if (!file.exists(data_file)) {
  data_file <- "https://raw.githubusercontent.com/sinaamiri9000-collab/pyquaidsce/publication/examples/data/household_demand.csv"
}
df <- read.csv(data_file)
shares <- c("w1", "w2", "w3", "w4")
prices <- c("p1", "p2", "p3", "p4")
demographics <- c("hh_size", "urban")

# 1. Basic run: no bootstrap and no expenditure endogeneity correction.
basic <- quaidsce(
  df, shares = shares, prices = prices, expenditure = "expenditure",
  demographics = demographics, anot = 8, reps = 0L
)
print(summary(basic))
print(elasticities(basic))

# 2. Use income as an instrument to account for endogenous expenditure.
# The package estimates the expenditure equation and uses its residual.
iv <- quaidsce(
  df, shares = shares, prices = prices, expenditure = "expenditure",
  demographics = demographics, anot = 8, reps = 0L, ivexp = c("income")
)
print(summary(iv))

# 3. Run bootstrap with the expenditure correction.
# reps = 0 skips bootstrap. reps = 10 attempts 10 replications.
# reps = 1 attempts one, but at least two must succeed for standard errors.
# We use 10 here for practice; increase reps for your research results.
boot <- quaidsce(
  df, shares = shares, prices = prices, expenditure = "expenditure",
  demographics = demographics, anot = 8, ivexp = c("income"),
  reps = 10L, seed = 12345L
)
print(summary(boot))
print(elasticities(boot))

# 4. Run bootstrap on two workers at the same time.
# warm starts each replication from the full-sample estimates.
# Each replication may run for at most 900 seconds. A longer replication
# counts as unsuccessful, and the program moves on to the next one.
parallel <- quaidsce(
  df, shares = shares, prices = prices, expenditure = "expenditure",
  demographics = demographics, anot = 8, ivexp = c("income"),
  reps = 10L, n_jobs = 2L, seed = 12345L, bootstrap_start = "warm",
  rep_timeout = 900
)
print(elasticities(parallel))

# 5. Reuse starting values from the same model fitted in example 2.
# Use theta, not coefficients: theta contains only the free demand parameters.
restarted <- quaidsce(
  df, shares = shares, prices = prices, expenditure = "expenditure",
  demographics = demographics, anot = 8, ivexp = c("income"), reps = 0L,
  initial = iv$theta, sigma_initial = iv$sigma
)
print(elasticities(restarted))

# 6. Let the package calculate starting values from a simpler AIDS fit.
linear_start <- quaidsce(
  df, shares = shares, prices = prices, expenditure = "expenditure",
  demographics = demographics, anot = 8, start = "linear", reps = 0L
)
print(elasticities(linear_start))

# 7. Fit AIDS by removing the quadratic expenditure term.
aids <- quaidsce(
  df, shares = shares, prices = prices, expenditure = "expenditure",
  demographics = demographics, anot = 8, quadratic = FALSE, reps = 0L
)
print(elasticities(aids))
