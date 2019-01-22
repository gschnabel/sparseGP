library(sparseGP)
library(numDeriv)
context("test sparseGP functionality")


gp <- list(x = matrix(runif(9,0,5),3,3),
           y = seq_len(3),
           err = rep(2.5,3),
           hyperpars = list(
             gpFuns = list(
               covfun = getMatrix,
               deriveCovfunInd = deriveMatrixInd,
               deriveCovfunHyp = deriveMatrixHyp
             ),
             hyperpars = list(
               sigma=2,
               len=rep(2,3),
               nugget=1.5
             ),
             ind = matrix(runif(9,0,5),3,3)
           ))


gpPostCov <- function(gp, px, py, exact=FALSE) {
  Kuu <- with(gp$hyperpars, gpFuns$covfun(ind, ind, hyperpars,FALSE))
  Kupx <- with(gp$hyperpars, gpFuns$covfun(ind,px,hyperpars,FALSE))
  Kupy <- with(gp$hyperpars, gpFuns$covfun(ind,py,hyperpars,FALSE))

  Kpxpy <- with(gp$hyperpars, gpFuns$covfun(px,py,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  Qpxpy <- t(Kupx)%*%invKuu%*%Kupy
  Kuf <- with(gp$hyperpars, gpFuns$covfun(ind,gp$x,hyperpars,FALSE))
  Kff <- with(gp$hyperpars, gpFuns$covfun(gp$x,gp$x,hyperpars,FALSE))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + as.vector(gp$err)^2)
  invDelta <- chol2inv(chol(delta))
  covmat <- t(Kupx) %*% solve(Kuu + Kuf %*% invDelta %*% t(Kuf)) %*% Kupy
  if (exact)
  {
    covmat <- Kpxpy - Qpxpy + covmat
  }
  else
  {
    equalMask <- outer(seq(nrow(px)),seq(nrow(py)),function(i,j) {
      apply(px[i,]-py[j,]==0, 1, all)
    })
    covmat[equalMask] <- covmat[equalMask] + Kpxpy[equalMask] - Qpxpy[equalMask]
  }
  covmat
}



gpMarlike <- function(gp) {
  Kuu <- with(gp$hyperpars, gpFuns$covfun(ind, ind, hyperpars,FALSE))
  Kuf <- with(gp$hyperpars, gpFuns$covfun(ind,gp$x,hyperpars,FALSE))
  Kff <- with(gp$hyperpars, gpFuns$covfun(gp$x,gp$x,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$err^2)
  invDelta <- chol2inv(chol(delta))
  Z <- (t(Kuf)%*%invKuu%*%Kuf + delta)
  res <- (-0.5) * (2*sum(log(diag(chol(Z)))) + t(gp$y) %*% solve(Z) %*% gp$y + length(gp$err)*log(2*pi))
  as.vector(res)
}

adjustedGpMarlike <- function(x,gp,what) {
  switch(what,
         "sigma" = { gp$hyperpars$hyperpars$sigma <- gp$hyperpars$hyperpars$sigma + x },
         "len1" = { gp$hyperpars$hyperpars$len[1] <- gp$hyperpars$hyperpars$len[1] + x },
         "len2" = { gp$hyperpars$hyperpars$len[2] <- gp$hyperpars$hyperpars$len[2] + x },
         "ind11" = { gp$hyperpars$ind[1,1] <- gp$hyperpars$ind[1,1] + x },
         "ind13" = { gp$hyperpars$ind[1,3] <- gp$hyperpars$ind[1,3] + x },
         "ind22" = { gp$hyperpars$ind[2,2] <- gp$hyperpars$ind[2,2] + x })
  gpMarlike(gp)
}

numDeriveGpMarlike <- function(gp,what) {
  res <- jacobian(adjustedGpMarlike,0,gp=gp,what=what)
  as.vector(res)
}



test_that("marginal likelihood is correctly computed",{

  refMarlike <- gpMarlike(gp)
  marlike <- with(gp, fitSparseGP(gp$x, gp$y, gp$err, hyperpars))$marlike
  expect_equal(marlike,refMarlike)
})


test_that("derivative of marginal likelihood wrt induction points is correctly computed",{

  refDerive11 <- numDeriveGpMarlike(gp,"ind11")
  refDerive13 <- numDeriveGpMarlike(gp,"ind13")
  refDerive22 <- numDeriveGpMarlike(gp,"ind22")

  myfit <- with(gp,fitSparseGP(x,y,err,hyperpars))
  derive <- with(gp,deriveSparseGPMarlikeInd(myfit))
  expect_equal(derive[1,1],refDerive11)
  expect_equal(derive[1,3],refDerive13)
  expect_equal(derive[2,2],refDerive22)
})


test_that("derivative of marginal likelihood wrt hyperparameters is correctly computed",{

  refDeriveSigma <- numDeriveGpMarlike(gp,"sigma")
  refDeriveLen1 <- numDeriveGpMarlike(gp,"len1")
  refDeriveLen2 <- numDeriveGpMarlike(gp,"len2")

  myfit <- with(gp,fitSparseGP(x,y,err,hyperpars))
  deriveSigma <- with(gp,deriveSparseGPMarlikeHyp(myfit,list(hyperpars=list(sigma=TRUE))))
  deriveLen1 <- with(gp,deriveSparseGPMarlikeHyp(myfit,list(hyperpars=list(len=c(TRUE,FALSE,FALSE)))))
  deriveLen2 <- with(gp,deriveSparseGPMarlikeHyp(myfit,list(hyperpars=list(len=c(FALSE,TRUE,FALSE)))))

  expect_equal(deriveSigma,refDeriveSigma)
  expect_equal(deriveLen1,refDeriveLen1)
  expect_equal(deriveLen2,refDeriveLen2)
})


test_that("partial derivative of marginal likelihood wrt induction points is correctly computed",{

  myfit <- with(gp,fitSparseGP(x,y,err,hyperpars))
  derive <- deriveSparseGPMarlikeInd(myfit)
  derive12 <- deriveSparseGPMarlikeIndSel(myfit,c(TRUE,TRUE,FALSE))
  derive23 <- deriveSparseGPMarlikeIndSel(myfit,c(FALSE,TRUE,TRUE))
  expect_equal(derive[1:2,],derive12)
  expect_equal(derive[2:3,],derive23)
})


test_that("exact posterior covariance computed correctly",{
  myfit <- with(gp, fitSparseGP(gp$x, gp$y, gp$err, hyperpars))
  refCovmat <- gpPostCov(myfit,gp$x,gp$x,exact=TRUE)
  actCovmat <- getPostCovFromSparseGP(gp$x,gp$x,myfit,FALSE,exact=TRUE)
  expect_equal(actCovmat, refCovmat)

  px <- matrix(runif(15,0,5),5,3)
  py <- matrix(runif(9,0,5),3,3)
  refCovmat <- gpPostCov(myfit,px,py,exact=TRUE)
  actCovmat <- getPostCovFromSparseGP(px, py, myfit,FALSE,exact=TRUE)
  expect_equal(actCovmat, refCovmat)
})


test_that("FITC posterior covariance computed correctly",{
  myfit <- with(gp, fitSparseGP(gp$x, gp$y, gp$err, hyperpars))
  refCovmat <- gpPostCov(myfit,gp$x,gp$x,exact=FALSE)
  actCovmat <- getPostCovFromSparseGP(gp$x,gp$x,myfit,FALSE,exact=FALSE)
  expect_equal(actCovmat, refCovmat)

  px <- matrix(runif(15,0,5),5,3)
  py <- matrix(runif(9,0,5),3,3)
  refCovmat <- gpPostCov(myfit,px,py,exact=FALSE)
  actCovmat <- getPostCovFromSparseGP(px, py, myfit,FALSE,exact=FALSE)
  expect_equal(actCovmat, refCovmat)
})


test_that("exact posterior covariance agree in element- and matrix-mode",{
  myfit <- with(gp, fitSparseGP(gp$x, gp$y, gp$err, hyperpars))

  px <- gp$x
  py <- gp$x
  grid <- expand.grid(i = seq(nrow(px)), j = seq(nrow(py)))
  refCovmat <- getPostCovFromSparseGP(px,py,myfit,FALSE,exact=TRUE)
  actCovmat <- getPostCovFromSparseGP(px[grid$i,],py[grid$j,],myfit,TRUE,exact=TRUE)
  dim(actCovmat) <- c(nrow(px),nrow(py))
  expect_equal(actCovmat, refCovmat)

  px <- matrix(runif(15,0,5),5,3)
  py <- matrix(runif(9,0,5),3,3)
  grid <- expand.grid(i = seq(nrow(px)), j = seq(nrow(py)))
  refCovmat <- getPostCovFromSparseGP(px,py,myfit,FALSE,exact=TRUE)
  actCovmat <- getPostCovFromSparseGP(px[grid$i,],py[grid$j,],myfit,TRUE,exact=TRUE)
  dim(actCovmat) <- c(nrow(px),nrow(py))
  expect_equal(actCovmat, refCovmat)
})


test_that("FITC posterior covariance agree in element- and matrix-mode",{
  myfit <- with(gp, fitSparseGP(gp$x, gp$y, gp$err, hyperpars))

  px <- gp$x
  py <- gp$x
  grid <- expand.grid(i = seq(nrow(px)), j = seq(nrow(py)))
  refCovmat <- getPostCovFromSparseGP(px,py,myfit,FALSE,exact=FALSE)
  actCovmat <- getPostCovFromSparseGP(px[grid$i,],py[grid$j,],myfit,TRUE,exact=FALSE)
  dim(actCovmat) <- c(nrow(px),nrow(py))
  expect_equal(actCovmat, refCovmat)

  px <- matrix(runif(15,0,5),5,3)
  py <- matrix(runif(9,0,5),3,3)
  grid <- expand.grid(i = seq(nrow(px)), j = seq(nrow(py)))
  refCovmat <- getPostCovFromSparseGP(px,py,myfit,FALSE,exact=FALSE)
  actCovmat <- getPostCovFromSparseGP(px[grid$i,],py[grid$j,],myfit,TRUE,exact=FALSE)
  dim(actCovmat) <- c(nrow(px),nrow(py))
  expect_equal(actCovmat, refCovmat)
})



