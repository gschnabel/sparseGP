library(sparseGP)
library(mvtnorm)
context("test polygon chain covariance kernel functionality")

hyperpars <- list(
  loc = seq(0,10,by=1),
  var = rep(1,11)
)

test_that("matrix block and element mode agree", {
  gridx <- seq(-5,15,by=0.1)
  covmat <- getMatrixPolyChain(gridx, gridx, hyperpars, FALSE)
  diagcov <- getMatrixPolyChain(gridx, gridx, hyperpars, TRUE)
  expect_equal(diagcov,diag(covmat))
})


test_that("sampling from polychain kernel really yields a polygon chain", {
  gridx <- seq(0,10,by=0.1)
  covmat <- getMatrixPolyChain(gridx, gridx, hyperpars, FALSE)
  smpl <- rmvnorm(1,rep(0,length(gridx)), covmat)
  dim(smpl) <- NULL
  idx <- match(hyperpars$loc, gridx)
  linApprox <- approx(gridx[idx], smpl[idx], xout=gridx)$y
  expect_equal(smpl, linApprox, tolerance=sqrt(.Machine$double.eps)*100)
})


test_that("variance/covariance outside boundaries equals zero",{

  res1 <- getMatrixPolyChain(cbind(-1), cbind(5), hyperpars, TRUE)
  res2 <- getMatrixPolyChain(cbind(5), cbind(11), hyperpars, TRUE)
  res3 <- getMatrixPolyChain(cbind(-1), cbind(-1), hyperpars, TRUE)
  res4 <- getMatrixPolyChain(cbind(11), cbind(11), hyperpars, TRUE)
  expect_equal(res1,0)
  expect_equal(res2,0)
  expect_equal(res3,0)
  expect_equal(res4,0)
})


test_that("non-rectangular covariance does not yield an error", {

  gridx <- c(1,3,5)
  gridy <- c(2,8)
  expect_error(getMatrixPolyChain(gridx,gridy, hyperpars, FALSE), NA)
})









