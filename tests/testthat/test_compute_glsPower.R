context("compute_glsPower function")

## Basic functionality #####

test_that("compute_glsPower works with simple SWD design", {
  desmat <- construct_DesMat(Cl = c(2, 2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1)
  expect_is(result, "list")
  expect_true("power" %in% names(result))
  expect_true("Params" %in% names(result))
  expect_true("ProjMatrix" %in% names(result))
  expect_is(result$power, "numeric")
  expect_true(result$power > 0)
  expect_true(result$power <= 1)
})

test_that("compute_glsPower works with verbose=0", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, verbose = 0)
  expect_is(result, "numeric")
  expect_length(result, 1)
  expect_true(is.finite(result))
})

test_that("compute_glsPower works with parallel design", {
  desmat <- construct_DesMat(Cl = c(3, 3), dsntype = "parallel", timepoints = 3)
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, tau = 0.2)
  expect_is(result, "list")
  expect_true(result$power > 0)
})

## Effect size tests #####

test_that("compute_glsPower returns higher power for larger effect sizes", {
  desmat <- construct_DesMat(Cl = c(3, 3, 3), dsntype = "SWD")
  power_small <- compute_glsPower(DesMat = desmat, EffSize = 0.2, sigma = 1, verbose = 0)
  power_large <- compute_glsPower(DesMat = desmat, EffSize = 0.8, sigma = 1, verbose = 0)
  expect_true(power_large > power_small)
})

test_that("compute_glsPower returns power near 0 for effect size 0", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0, sigma = 1, verbose = 0)
  expect_true(result < 0.1)
})

## Covariance structure tests #####

test_that("compute_glsPower works with tau parameter", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, tau = 0.5)
  expect_is(result, "list")
  expect_true(result$power > 0)
})

test_that("compute_glsPower works with eta parameter", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, tau = 0.3, eta = 0.2)
  expect_is(result, "list")
  expect_true(result$power > 0)
})

test_that("compute_glsPower works with AR correlation", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, tau = 0.3, AR = c(0.5))
  expect_is(result, "list")
  expect_true(result$power > 0)
})

test_that("compute_glsPower works with rho parameter", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, tau = 0.3, eta = 0.2, rho = 0.5)
  expect_is(result, "list")
  expect_true(result$power > 0)
})

test_that("compute_glsPower works with gamma parameter", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, gamma = 0.3)
  expect_is(result, "list")
  expect_true(result$power > 0)
})

test_that("compute_glsPower works with psi parameter", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD", N = c(10, 10, 10, 10), INDIV_LVL = TRUE)
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, psi = 0.3, INDIV_LVL = TRUE)
  expect_is(result, "list")
  expect_true(result$power > 0)
})

test_that("compute_glsPower works with pre-computed CovMat", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  covmat <- construct_CovMat(sumCl = 4, timepoints = 3, sigma = 1, tau = 0.2)
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, CovMat = covmat)
  expect_is(result, "list")
  expect_true(result$power > 0)
})

## dfAdjust tests #####

test_that("compute_glsPower works with dfAdjust='none'", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, dfAdjust = "none")
  expect_is(result, "list")
  expect_equal(result$Params$dfAdjust, "none")
  expect_true(is.infinite(result$Params$denomDF))
})

test_that("compute_glsPower works with dfAdjust='between-within'", {
  desmat <- construct_DesMat(Cl = c(5, 5, 5), dsntype = "SWD")
  result <- suppressWarnings(compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, dfAdjust = "between-within"))
  expect_is(result, "list")
  expect_equal(result$Params$dfAdjust, "between-within")
  expect_true(is.finite(result$Params$denomDF))
})

test_that("compute_glsPower works with dfAdjust='containment'", {
  desmat <- construct_DesMat(Cl = c(2, 2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, dfAdjust = "containment")
  expect_is(result, "list")
  expect_equal(result$Params$dfAdjust, "containment")
  expect_true(is.finite(result$Params$denomDF))
})

test_that("compute_glsPower works with dfAdjust='residual'", {
  desmat <- construct_DesMat(Cl = c(2, 2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, dfAdjust = "residual")
  expect_is(result, "list")
  expect_equal(result$Params$dfAdjust, "residual")
  expect_true(is.finite(result$Params$denomDF))
})

## Significance level tests #####

test_that("compute_glsPower works with custom significance level", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, sig.level = 0.10)
  expect_is(result, "list")
  expect_equal(result$Params$sig.level, 0.10)
  expect_true(result$power > 0)
})

test_that("compute_glsPower returns higher power for higher significance level", {
  desmat <- construct_DesMat(Cl = c(3, 3), dsntype = "SWD")
  power_01 <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, sig.level = 0.01, verbose = 0)
  power_05 <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, sig.level = 0.05, verbose = 0)
  expect_true(power_05 > power_01)
})

## Individual level tests #####

test_that("compute_glsPower works with INDIV_LVL=TRUE", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD", N = c(10, 10, 10, 10), INDIV_LVL = TRUE)
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, psi = 0.2, INDIV_LVL = TRUE)
  expect_is(result, "list")
  expect_true(result$power > 0)
})

## Information content tests #####

test_that("compute_glsPower works with INFO_CONTENT=TRUE", {
  desmat <- construct_DesMat(Cl = c(3, 3, 3), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, tau = 0.3, INFO_CONTENT = TRUE)
  expect_is(result, "list")
  expect_true("InformationContent" %in% names(result))
  expect_is(result$InformationContent, "list")
})

test_that("compute_glsPower handles INFO_CONTENT failure gracefully", {
  desmat <- construct_DesMat(Cl = c(1, 1), dsntype = "parallel", timepoints = 2)
  expect_warning(
    result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, INFO_CONTENT = TRUE),
    "Information content calculation failed"
  )
  expect_is(result, "list")
  expect_true("power" %in% names(result))
  expect_false("InformationContent" %in% names(result))
})

## Verbose mode tests #####

test_that("compute_glsPower works with verbose=2", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, verbose = 2)
  expect_is(result, "list")
  expect_true("DesignMatrix" %in% names(result))
  expect_true("CovarianceMatrix" %in% names(result))
  expect_true("VarianceMatrix" %in% names(result))
})

## Different design types #####

test_that("compute_glsPower works with parallel_baseline design", {
  desmat <- construct_DesMat(Cl = c(3, 3), dsntype = "parallel_baseline", timepoints = 3)
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, tau = 0.2)
  expect_is(result, "list")
  expect_true(result$power > 0)
})

test_that("compute_glsPower works with crossover design", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "crossover", timepoints = c(2, 2))
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, tau = 0.2)
  expect_is(result, "list")
  expect_true(result$power > 0)
})

## Error handling tests #####

test_that("compute_glsPower throws error with missing DesMat", {
  expect_error(
    compute_glsPower(DesMat = NULL, EffSize = 0.5, sigma = 1)
  )
})

test_that("compute_glsPower throws error with missing EffSize", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  expect_error(
    compute_glsPower(DesMat = desmat, sigma = 1)
  )
})

test_that("compute_glsPower throws error with missing sigma", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  expect_error(
    compute_glsPower(DesMat = desmat, EffSize = 0.5)
  )
})

## Edge cases #####

test_that("compute_glsPower works with single cluster per group", {
  desmat <- construct_DesMat(Cl = c(1, 1), dsntype = "parallel", timepoints = 2)
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, verbose = 0)
  expect_is(result, "numeric")
  expect_true(is.finite(result))
})

test_that("compute_glsPower works with many timepoints", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD", timepoints = 10)
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, verbose = 0)
  expect_is(result, "numeric")
  expect_true(is.finite(result))
})

test_that("compute_glsPower works with tau=0", {
  desmat <- construct_DesMat(Cl = c(2, 2), dsntype = "SWD")
  result <- compute_glsPower(DesMat = desmat, EffSize = 0.5, sigma = 1, tau = 0)
  expect_is(result, "list")
  expect_true(result$power > 0)
})
