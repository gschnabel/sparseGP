library(sparseGP)
library(numDeriv)
context("test sparseGP_freeDelta with denseLinearMapping functionality")


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
             deltaFuns = createDenseLinearMapping(),
             deltapars = list(
               D = rep(0,3),
               S = matrix(c(3,1,1,1,5,1,1,1,7),3,3),
               P = data.frame(VAR=c(20,25,30))
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
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + as.vector(gp$err)^2) + gp$hyperpars$deltapars$delta
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
  deltaPart <- diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + as.vector(gp$err)^2
  delta <-  diag(deltaPart) + with(gp$hyperpars$deltapars, diag(D) + S %*% diag(P$VAR) %*% t(S))
  invDelta <- chol2inv(chol(delta))
  cholX <- chol(Kuu + Kuf%*%invDelta%*%t(Kuf))
  invX <- chol2inv(cholX)
  Z <- (t(Kuf)%*%invKuu%*%Kuf + delta)
  maha <- t(gp$y) %*% solve(Z) %*% gp$y
  logDet <-  2*sum(log(diag(chol(Z))))
  res <- (-0.5) * (logDet + maha + length(gp$hyperpars$deltapars$D)*log(2*pi))
  as.vector(res)

  #logDet #debug
  # list(res=as.vector(res),
  #      maha=maha,
  #      logDet=logDet,
  #      deltaPart=deltaPart, delta=delta,
  #      KufInvDeltaKfu=Kuf%*%invDelta%*%t(Kuf),
  #      yInvDeltay = t(gp$y) %*% invDelta %*% gp$y,
  #      invDeltay = invDelta %*% gp$y,
  #      y = gp$y) # debug
}

adjustedGpMarlike <- function(x,gp,what) {
  switch(what,
         "sigma" = { gp$hyperpars$hyperpars$sigma <- gp$hyperpars$hyperpars$sigma + x },
         "len1" = { gp$hyperpars$hyperpars$len[1] <- gp$hyperpars$hyperpars$len[1] + x },
         "len2" = { gp$hyperpars$hyperpars$len[2] <- gp$hyperpars$hyperpars$len[2] + x },
         "ind11" = { gp$hyperpars$ind[1,1] <- gp$hyperpars$ind[1,1] + x },
         "ind13" = { gp$hyperpars$ind[1,3] <- gp$hyperpars$ind[1,3] + x },
         "ind22" = { gp$hyperpars$ind[2,2] <- gp$hyperpars$ind[2,2] + x },
         # deltapars
         "pvar1" = { gp$hyperpars$deltapars$P$VAR[1] <- gp$hyperpars$deltapars$P$VAR[1] + x },
         "pvar2" = { gp$hyperpars$deltapars$P$VAR[2] <- gp$hyperpars$deltapars$P$VAR[2] + x },
         "pvar3" = { gp$hyperpars$deltapars$P$VAR[3] <- gp$hyperpars$deltapars$P$VAR[3] + x } )
  gpMarlike(gp)
}

numDeriveGpMarlike <- function(gp,what) {
  res <- jacobian(adjustedGpMarlike,0,gp=gp,what=what)
  as.vector(res)
}



test_that("marginal likelihood is correctly computed",{

  refMarlike <- gpMarlike(gp)
  marlike <- with(gp, fitSparseGPfD(gp$x, gp$y, gp$err, hyperpars))$marlike
  expect_equal(marlike,refMarlike)

  # refMarlike <- gpMarlike(gp)
  # marlike <- with(gp, fitSparseGPfD(gp$x, gp$y, gp$err, hyperpars))$marlike
  # expect_equal(marlike,refMarlike$res)
  #
  # Kuf <- with(gp$hyperpars, gpFuns$covfun(ind,gp$x,hyperpars,FALSE))
  # myfit <- with(gp, fitSparseGPfD(gp$x, gp$y, gp$err, hyperpars))
  # blub <- gp$hyperpars$deltaFuns$bothMultInv(t(Kuf), myfit$delta, gp$hyperpars$deltapars)
  # delta <- with(gp, fitSparseGPfD(gp$x, gp$y, gp$err, hyperpars))$delta
  # cholX <- with(gp, fitSparseGPfD(gp$x, gp$y, gp$err, hyperpars))$cholX
  # myfit$KfuInvDeltaKuf
  # blub
})


test_that("derivative of marginal likelihood wrt induction points is correctly computed",{

  refDerive11 <- numDeriveGpMarlike(gp,"ind11")
  refDerive13 <- numDeriveGpMarlike(gp,"ind13")
  refDerive22 <- numDeriveGpMarlike(gp,"ind22")

  myfit <- with(gp,fitSparseGPfD(x,y,err,hyperpars))
  derive <- with(gp,deriveSparseGPMarlikeIndfD(myfit))
  expect_equal(derive[1,1],refDerive11)
  expect_equal(derive[1,3],refDerive13)
  expect_equal(derive[2,2],refDerive22)
})



test_that("partial derivative of marginal likelihood wrt induction points is correctly computed",{

  myfit <- with(gp,fitSparseGPfD(x,y,err,hyperpars))
  derive <- deriveSparseGPMarlikeIndfD(myfit)
  derive12 <- deriveSparseGPMarlikeIndSelfD(myfit,c(TRUE,TRUE,FALSE))
  derive23 <- deriveSparseGPMarlikeIndSelfD(myfit,c(FALSE,TRUE,TRUE))
  expect_equal(derive[1:2,],derive12)
  expect_equal(derive[2:3,],derive23)
})



test_that("derivative of marginal likelihood wrt hyperparameters is correctly computed",{

  refDeriveSigma <- numDeriveGpMarlike(gp,"sigma")
  refDeriveLen1 <- numDeriveGpMarlike(gp,"len1")
  refDeriveLen2 <- numDeriveGpMarlike(gp,"len2")

  myfit <- with(gp,fitSparseGPfD(x,y,err,hyperpars))
  deriveSigma <- with(gp,deriveSparseGPMarlikeHypfD(myfit,list(hyperpars=list(sigma=TRUE))))
  deriveLen1 <- with(gp,deriveSparseGPMarlikeHypfD(myfit,list(hyperpars=list(len=c(TRUE,FALSE,FALSE)))))
  deriveLen2 <- with(gp,deriveSparseGPMarlikeHypfD(myfit,list(hyperpars=list(len=c(FALSE,TRUE,FALSE)))))

  expect_equal(deriveSigma,refDeriveSigma)
  expect_equal(deriveLen1,refDeriveLen1)
  expect_equal(deriveLen2,refDeriveLen2)
})



test_that("derivative of marginal likelihood wrt delta parameters is correctly computed",{

  refDeriveDelta1 <- numDeriveGpMarlike(gp,"pvar1")
  refDeriveDelta2 <- numDeriveGpMarlike(gp,"pvar2")
  refDeriveDelta3 <- numDeriveGpMarlike(gp,"pvar3")

  myfit <- with(gp,fitSparseGPfD(x,y,err,hyperpars))
  deriveDelta1 <- with(gp,deriveSparseGPMarlikeDeltafD(myfit,list(deltapars=list(P=list(VAR=c(TRUE,FALSE,FALSE))))))
  deriveDelta2 <- with(gp,deriveSparseGPMarlikeDeltafD(myfit,list(deltapars=list(P=list(VAR=c(FALSE,TRUE,FALSE))))))
  deriveDelta3 <- with(gp,deriveSparseGPMarlikeDeltafD(myfit,list(deltapars=list(P=list(VAR=c(FALSE,FALSE,TRUE))))))
  dim(deriveDelta1) <- dim(deriveDelta2) <- dim(deriveDelta3) <- NULL

  deriveDeltaAll <- with(gp,deriveSparseGPMarlikeDeltafD(myfit,list(deltapars=list(P=list(VAR=c(TRUE,TRUE,TRUE))))))
  dim(deriveDeltaAll) <- NULL

  expect_equal(deriveDelta1,refDeriveDelta1)
  expect_equal(deriveDelta2,refDeriveDelta2)
  expect_equal(deriveDelta3,refDeriveDelta3)
  expect_equal(deriveDeltaAll,c(refDeriveDelta1,refDeriveDelta2,refDeriveDelta3))
})




# these tests are taken from sparseGP with diagonal delta
# functionality needs still to be implemented in sparseGP_freeDelta
# note: maybe the functions can be taken to be one to one the same?
notrun <- function() {

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

}

