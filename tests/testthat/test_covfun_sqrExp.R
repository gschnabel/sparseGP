library(sparseGP)
library(numDeriv)
context("test squared exp covfun functionality")

gp <- list(fun = getMatrix,
           dfunInd = deriveMatrixInd,
           dfunHyp = deriveMatrixHyp,
           x = matrix(runif(9,0,5),3,3),
           y = matrix(runif(9,0,5),3,3),
           obs = seq_len(3),
           B = rep(2.5,3),
           hyperpars = list(sigma=2,len=rep(2,3),nugget=1.2))

covfun <- function(x,y,hyperpars,elementMode) {

  dist <- 0
  for (i in seq_len(ncol(x))) {
    z <- outer(x[,i],y[,i],`-`)/hyperpars$len[i]
    curdist <- -0.5 * z*z
    dist <- dist + curdist
  }
  nuggetInc <- ifelse(dist==0,hyperpars$nugget^2,0)

  res <- hyperpars$sigma^2 * exp(dist) + nuggetInc
  if (!isTRUE(elementMode)) res else diag(res)
}


adjustedCovfun <- function(x,y,z,hyperpars,elementMode,what) {
  switch(what,
         "sigma" = { hyperpars$sigma <- hyperpars$sigma + x },
         "len1" = { hyperpars$len[1] <- hyperpars$len[1] + x },
         "len2" = { hyperpars$len[2] <- hyperpars$len[2] + x })
  covfun(y,z,hyperpars,elementMode)
}

numDeriveCovfun <- function(x,y,hyperpars,elementMode,what) {
  res <- jacobian(adjustedCovfun,0,y=x,z=y,hyperpars=hyperpars,
                  elementMode=elementMode,what=what)
  if (!elementMode)
    dim(res) <- c(nrow(x),nrow(y))
  else
    res <- as.vector(res)
  res
}


test_that("covariance function is computed correctly",{

  realCovmat <- with(gp,covfun(x,y,hyperpars,FALSE))
  covmat <- with(gp,getMatrix(x,y,hyperpars,FALSE))
  expect_equal(covmat,realCovmat)
  # with nugget parameter
  realCovmat <- with(gp,covfun(x,x,hyperpars,FALSE))
  covmat <- with(gp,getMatrix(x,x,hyperpars,FALSE))
  expect_equal(covmat,realCovmat)
})


test_that("derivative of covfun wrt induction points is correctly computed",{

  refDeriveInd <- with(gp,jacobian(covfun,x,y=y,hyperpars=hyperpars,elementMode=FALSE))
  for (i in seq_len(nrow(gp$x))) {
    dimIndex <- floor((i-1) / nrow(gp$x))
    rowIndex <- (i-1) %% nrow(gp$x) + 1
    deriveInd <- with(gp,deriveMatrixInd(x,y,hyperpars,FALSE,dimIndex))
    tmp <- refDeriveInd[,i]
    dim(tmp) <- c(nrow(gp$x),ncol(gp$x))
    curRef <-tmp[rowIndex,]
    curRes <- deriveInd[rowIndex,]
    expect_equal(curRes,curRef)
  }
  # with nugget parameter
  refDeriveInd <- with(gp,jacobian(covfun,x,y=x,hyperpars=hyperpars,elementMode=FALSE))
  for (i in seq_len(nrow(gp$x))) {
    dimIndex <- floor((i-1) / nrow(gp$x))
    rowIndex <- (i-1) %% nrow(gp$x) + 1
    deriveInd <- with(gp,deriveMatrixInd(x,x,hyperpars,FALSE,dimIndex))
    tmp <- refDeriveInd[,i]
    dim(tmp) <- c(nrow(gp$x),ncol(gp$x))
    curRef <-tmp[rowIndex,]
    curRes <- deriveInd[rowIndex,]
    expect_equal(curRes,curRef)
  }
})



test_that("derivative of covfun wrt hyperperameters is correctly computed",{

  refMatSigma <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"sigma"))
  matSigma <- with(gp,deriveMatrixHyp(x,y,hyperpars,list(sigma=TRUE),FALSE,TRUE))
  refMatLen1 <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"len1"))
  matLen1 <- with(gp,deriveMatrixHyp(x,y,hyperpars,list(len=c(TRUE,FALSE,FALSE)),FALSE,TRUE))
  refMatLen2 <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"len2"))
  matLen2 <- with(gp,deriveMatrixHyp(x,y,hyperpars,list(len=c(FALSE,TRUE,FALSE)),FALSE,TRUE))

  expect_equal(matSigma,refMatSigma)
  expect_equal(matLen1,refMatLen1)
  expect_equal(matLen2,refMatLen2)

  # diagonal elements
  refDiagSigma <- with(gp,numDeriveCovfun(x,x,hyperpars,TRUE,"sigma"))
  diagSigma <- as.vector(with(gp,deriveMatrixHyp(x,x,hyperpars,list(sigma=TRUE),TRUE,TRUE)))
  refDiagLen1 <- with(gp,numDeriveCovfun(x,x,hyperpars,TRUE,"len1"))
  diagLen1 <- as.vector(with(gp,deriveMatrixHyp(x,x,hyperpars,list(len=c(TRUE,FALSE,FALSE)),TRUE,TRUE)))
  refDiagLen2 <- with(gp,numDeriveCovfun(x,x,hyperpars,TRUE,"len2"))
  diagLen2 <- as.vector(with(gp,deriveMatrixHyp(x,x,hyperpars,list(len=c(FALSE,TRUE,FALSE)),TRUE,TRUE)))

  expect_equal(diagSigma,refDiagSigma)
  expect_equal(diagLen1,refDiagLen1)
  expect_equal(diagLen2,refDiagLen2)
})








