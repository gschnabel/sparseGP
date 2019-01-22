library(sparseGP)
library(mvtnorm)
context("test histogram covariance kernel functionality")

hyperpars <- list(
  binloc = seq(0,10,by=1),
  binvar = rep(1,10)
)

test_that("matrix block and element mode agree", {
  gridx <- seq(0,9.9,by=0.1)
  covmat <- getMatrixHisto(gridx, gridx, hyperpars, FALSE)
  diagcov <- getMatrixHisto(gridx, gridx, hyperpars, TRUE)
  expect_equal(diagcov,diag(covmat))
})


test_that("sampling from histogram kernel really yields a histogram", {
  gridx <- seq(0,9.9,by=0.1)
  covmat <- getMatrixHisto(gridx, gridx, hyperpars, FALSE)
  smpl <- rmvnorm(1,rep(0,length(gridx)), covmat)
  fact <- cut(gridx, hyperpars$binloc, right=FALSE)
  expect_true(all(tapply(smpl, fact, sd) < 1e-06))
})


test_that("variance/covariance outside boundaries equals zero",{

  res1 <- getMatrixHisto(cbind(-1), cbind(5), hyperpars, TRUE)
  res2 <- getMatrixHisto(cbind(5), cbind(11), hyperpars, TRUE)
  res3 <- getMatrixHisto(cbind(-1), cbind(-1), hyperpars, TRUE)
  res4 <- getMatrixHisto(cbind(11), cbind(11), hyperpars, TRUE)
  expect_equal(res1,0)
  expect_equal(res2,0)
  expect_equal(res3,0)
  expect_equal(res4,0)
})




