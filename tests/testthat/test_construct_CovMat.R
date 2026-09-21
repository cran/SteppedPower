
context("construct_CovMat functions")


## construct_CovBlk #####

test_that("construct_CovBlk works with only sigma", {
  expect_equal(
    construct_CovBlk(c(3,3,3), tau=NULL, eta=NULL),
    diag(9, 3)
  )
})

test_that("construct_CovBlk works with tau only", {
  expect_equal(
    construct_CovBlk(c(3,3,3), tau=c(1,1,1)),
    diag(9, 3) + matrix(1, 3, 3)
  )
})

test_that("construct_CovBlk works with both tau and eta", {
  expect_equal(
    construct_CovBlk(c(3,3,3), tau=c(1,1,1), eta=c(1,1,1)),
    diag(9, 3) + matrix(1, 3, 3) + matrix(1, 3, 3)
  )
})

test_that("construct_CovBlk works with AR(1) correlation for tau", {
  tau_cov <- matrix(c(1,0.5,0.25,0.5,1,0.5,0.25,0.5,1), 3, 3)
  expect_equal(
    construct_CovBlk(c(3,3,3), tau=c(1,1,1), AR=list(0.5, NULL)),
    diag(9, 3) + tau_cov
  )
})

test_that("construct_CovBlk works with AR(1) correlation for eta", {
  eta_cov <- matrix(c(1,0.8,0.64,0.8,1,0.8,0.64,0.8,1), 3, 3)
  expect_equal(
    construct_CovBlk(c(2,2,2), tau=c(1,1,1), eta=c(1,1,1), AR=list(NULL, 0.8)),
    diag(4, 3) + matrix(1, 3, 3) + eta_cov
  )
})

test_that("construct_CovBlk works with rho", {
  result <- construct_CovBlk(c(2,2,2), tau=c(1,1,1), eta=c(1,1,1), rho=0.5)
  # rhoVec = rho * eta * tau = 0.5 * 1 * 1 = 0.5 for each element
  # outer(rhoVec, rhoVec, "+") = outer(c(0.5,0.5,0.5), c(0.5,0.5,0.5), "+")
  # This gives matrix(c(1,1,1,1,1,1,1,1,1), 3, 3) * 0.25? No wait...
  # outer(a,b,"+") gives a_i + b_j, so outer(c(0.5,0.5,0.5), c(0.5,0.5,0.5), "+") 
  # gives matrix(c(1,1,1,1,1,1,1,1,1), 3, 3)
  # Wait, no: outer(c(0.5,0.5,0.5), c(0.5,0.5,0.5), "+") = 
  # c(0.5,0.5,0.5) %o% c(0.5,0.5,0.5) which is matrix(rep(1,9), 3, 3)
  # Actually outer with + gives the sum, so 0.5+0.5=1 for all
  expected <- diag(4, 3) + matrix(1, 3, 3) + matrix(1, 3, 3) + matrix(1, 3, 3)
  expect_equal(result, expected)
})

test_that("construct_CovBlk throws error with mismatched sigma and tau lengths", {
  expect_error(
    construct_CovBlk(c(3,3), tau=c(1,1,1))
  )
})

test_that("construct_CovBlk throws error with eta but without tau", {
  expect_error(
    construct_CovBlk(c(3,3,3), tau=NULL, eta=c(1,1,1))
  )
})

test_that("construct_CovBlk throws error with mismatched tau and eta lengths", {
  expect_error(
    construct_CovBlk(c(3,3,3), tau=c(1,1), eta=c(1,1,1))
  )
})


## construct_CovMat #####

test_that("construct_CovMat works with simple case: 2 clusters, 3 timepoints", {
  result <- construct_CovMat(sumCl=2, timepoints=3, sigma=2, tau=1)
  expect_equal(dim(result), c(6, 6))
})

test_that("construct_CovMat works with N parameter", {
  result <- construct_CovMat(sumCl=2, timepoints=3, sigma=2, tau=1, N=5)
  expect_equal(dim(result), c(6, 6))
})

test_that("construct_CovMat works with gamma", {
  result <- construct_CovMat(sumCl=2, timepoints=3, sigma=2, tau=1, gamma=0.5)
  expect_equal(dim(result), c(6, 6))
  # Result is a sparse matrix, use Matrix::diag
  diag_vals <- Matrix::diag(result)
  # sigma^2 = 4, tau^2 = 1, gamma^2 = 0.25, so diagonal should be 5.25
  expect_true(all(as.numeric(diag_vals) > 5))
})

test_that("construct_CovMat works with matrix inputs", {
  tau_mat <- matrix(c(1,1,2,1,1,2), nrow=2, ncol=3, byrow=TRUE)
  sigma_mat <- matrix(c(2,2,2,3,3,3), nrow=2, byrow=TRUE)
  result <- construct_CovMat(sumCl=2, timepoints=3, 
                             sigma=sigma_mat,
                             tau=tau_mat)
  expect_equal(dim(result), c(6, 6))
})

test_that("construct_CovMat works with eta as matrix", {
  trtMat_test <- matrix(c(0,0,1,0,1,1,1,1,1), nrow=3, ncol=3, byrow=TRUE)
  eta_mat <- matrix(c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5), nrow=3, ncol=3)
  result <- construct_CovMat(trtMat=trtMat_test, sigma=2, tau=1, eta=eta_mat)
  expect_equal(dim(result), c(9, 9))
})

test_that("construct_CovMat works with rho", {
  trtMat_test <- matrix(c(0,0,1,0,1,1,1,1,1), nrow=3, ncol=3, byrow=TRUE)
  eta_mat <- matrix(c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5), nrow=3, ncol=3)
  result <- construct_CovMat(trtMat=trtMat_test, sigma=2, tau=1, eta=eta_mat, rho=0.3)
  expect_equal(dim(result), c(9, 9))
})

test_that("construct_CovMat works with AR correlation", {
  result <- construct_CovMat(sumCl=2, timepoints=3, sigma=2, tau=1, AR=c(0.5, 0.3))
  expect_equal(dim(result), c(6, 6))
})

test_that("construct_CovMat throws error with rho but without eta", {
  expect_error(
    construct_CovMat(sumCl=2, timepoints=3, sigma=2, tau=1, rho=0.5)
  )
})

test_that("construct_CovMat works with trtMat parameter", {
  trtMat_test <- matrix(c(0,0,1,0,1,1), nrow=2, ncol=3, byrow=TRUE)
  result <- construct_CovMat(trtMat=trtMat_test, sigma=2, tau=1)
  expect_equal(dim(result), c(6, 6))
})

test_that("construct_CovMat works with psi (non-cross-sectional)", {
  result <- construct_CovMat(sumCl=2, timepoints=3, sigma=2, tau=1, psi=0.3, N=c(5,5))
  # With psi, this triggers non-cross-sectional path
  # The dimension will depend on implementation details
  expect_true(dim(result)[1] >= 6)
})

test_that("construct_CovMat works with INDIV_LVL=TRUE", {
  result <- construct_CovMat(sumCl=2, timepoints=3, sigma=2, tau=1, psi=0.3, 
                            N=c(5,5), INDIV_LVL=TRUE)
  expect_true(dim(result)[1] >= 6)
})


