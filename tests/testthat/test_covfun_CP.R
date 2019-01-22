library(sparseGP)
library(numDeriv)
context("test changepoint covfun functionality")

gp <- list(fun = getMatrixCP,
           dfunInd = deriveMatrixIndCP,
           dfunHyp = deriveMatrixHypCP,
           x = matrix(runif(9,0,5),3,3),
           y = matrix(runif(9,0,5),3,3),
           obs = seq_len(3),
           B = rep(2.5,3),
           hyperpars = list(gp1Funs = list(covfun=getMatrix,
                                           deriveCovfunInd=deriveMatrixInd,
                                           deriveCovfunHyp=deriveMatrixHyp),
                            gp2Funs = list(covfun=getMatrix,
                                           deriveCovfunInd=deriveMatrixInd,
                                           deriveCovfunHyp=deriveMatrixHyp),
                            gp1=list(sigma=2,len=rep(2,3),nugget=0),
                            gp2=list(sigma=3,len=rep(4,3),nugget=0),
                            cp=list(kc=0.5, xc=2.5,
                                    ord = 1:2,
                                    ang = c(runif(1,0,pi),runif(1,0,2*pi))),
                            nugget=3))



invLogit <- function(x, k, x0) {
  1/(1+exp(-k*(x-x0)))
}

deriveInvLogit <- function(x, k, x0) {

  invLogit(x, k, x0)^2 * exp(-k*(x-x0)) * k
}

sphereCoordToCartR <- function(ang, ord) {

  ang <- ang[ord]
  res <- rep(1,length(ang)+1)
  for (i in seq_len(length(ang)+1)) {
    if (i>1)
      res[i] <- res[i] * prod(sin((ang[1:(i-1)])))
    if (i<(length(ang)+1))
      res[i] <- res[i] * cos(ang[i])
  }
  res
}


covfun <- function(x, y, hyperpars, elementMode) {

  cov1 <- hyperpars$gp1Funs$covfun(x,y,hyperpars$gp1,FALSE)
  cov2 <- hyperpars$gp2Funs$covfun(x,y,hyperpars$gp2,FALSE)
  #hyperpars$cp$dir <- hyperpars$cp$dir / sqrt(sum(hyperpars$cp$dir^2))
  dir <- sphereCoordToCartR(hyperpars$cp$ang[hyperpars$cp$ord])

  proj1 <- colSums(apply(x,1, function(x) x*dir))
  proj2 <- colSums(apply(y,1, function(x) x*dir))
  s1 <- invLogit(proj1, hyperpars$cp$kc, hyperpars$cp$xc)
  s2 <- invLogit(proj2, hyperpars$cp$kc, hyperpars$cp$xc)
  r1 <- 1 - s1
  r2 <- 1 - s2

  res <- matrix(0,nrow(x),nrow(y))
  outer(seq_len(nrow(x)),seq_len(nrow(y)),function(i,j) {
    nuggetInc <- ifelse(apply(x[i,,drop=FALSE]==y[j,,drop=FALSE],1,all),hyperpars$nugget^2,0)
    res[(j-1)*nrow(res)+i] <<- s1[i]*s2[j]*cov2[(j-1)*nrow(res)+i] + r1[i]*r2[j]*cov1[(j-1)*nrow(res)+i] + nuggetInc
  })
  if (!elementMode) res else diag(res)
}




adjustedCovfun <- function(x,y,z,hyperpars,elementMode,what) {
  switch(what,
         "sigma1" = { hyperpars$gp1$sigma <- hyperpars$gp1$sigma + x },
         "sigma2" = { hyperpars$gp2$sigma <- hyperpars$gp2$sigma + x },
         "len11" = { hyperpars$gp1$len[1] <- hyperpars$gp1$len[1] + x },
         "len12" = { hyperpars$gp1$len[2] <- hyperpars$gp1$len[2] + x },
         "len21" = { hyperpars$gp2$len[1] <- hyperpars$gp2$len[1] + x },
         "len22" = { hyperpars$gp2$len[2] <- hyperpars$gp2$len[2] + x },
         "kc" = { hyperpars$cp$kc <- hyperpars$cp$kc + x },
         "xc" = { hyperpars$cp$xc <- hyperpars$cp$xc + x },
         "ang1" = { hyperpars$cp$ang[1] <- hyperpars$cp$ang[1] + x },
         "ang2" = { hyperpars$cp$ang[2] <- hyperpars$cp$ang[2] + x })
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




test_that("Covariance function computes correctly", {

  realCovmat <- with(gp,covfun(x,y,hyperpars,FALSE))
  covmat <- with(gp,getMatrixCP(x,y,hyperpars,FALSE))
  expect_equal(covmat,realCovmat)
})



test_that("Covariance function with nugget parameter computes correctly", {

  realCovmat <- with(gp,covfun(x,x,hyperpars,FALSE))
  covmat <- with(gp,getMatrixCP(x,x,hyperpars,FALSE))
  expect_equal(covmat,realCovmat)
})



test_that("Element-wise and block-wise calculation of covariance agree", {

  covmat <- with(gp,getMatrixCP(x,y,hyperpars,FALSE))
  diagCovmat <- as.vector(with(gp,getMatrixCP(x,y,hyperpars,TRUE)))
  expect_equal(diag(covmat),diagCovmat)
  # test with nugget parameter
  covmat <- with(gp,getMatrixCP(x,x,hyperpars,FALSE))
  diagCovmat <- as.vector(with(gp,getMatrixCP(x,x,hyperpars,TRUE)))
  expect_equal(diag(covmat),diagCovmat)
})



test_that("derivative of covfun wrt induction points is correctly computed",{

  refDeriveInd <- with(gp,jacobian(covfun,x,y=y,hyperpars=hyperpars,elementMode=FALSE))
  for (i in seq_len(nrow(gp$x))) {
    dimIndex <- floor((i-1) / nrow(gp$x))
    rowIndex <- (i-1) %% nrow(gp$x) + 1
    deriveInd <- with(gp,deriveMatrixIndCP(x,y,hyperpars,FALSE,dimIndex))
    tmp <- refDeriveInd[,i]
    dim(tmp) <- c(nrow(gp$x),ncol(gp$x))
    curRef <-tmp[rowIndex,]
    curRes <- deriveInd[rowIndex,]
    expect_equal(curRes,curRef)
  }
})



test_that("derivative of covfun wrt induction points with nugget parameter is correctly computed",{

  refDeriveInd <- with(gp,jacobian(covfun,x,y=x,hyperpars=hyperpars,elementMode=FALSE))
  for (i in seq_len(nrow(gp$x))) {
    dimIndex <- floor((i-1) / nrow(gp$x))
    rowIndex <- (i-1) %% nrow(gp$x) + 1
    deriveInd <- with(gp,deriveMatrixIndCP(x,x,hyperpars,FALSE,dimIndex))
    tmp <- refDeriveInd[,i]
    dim(tmp) <- c(nrow(gp$x),ncol(gp$x))
    curRef <-tmp[rowIndex,]
    curRes <- deriveInd[rowIndex,]
    expect_equal(curRes,curRef)
  }
})



test_that("Element-wise and block-wise calculation of derivative wrt induction points agree", {

  deriveCovmat <- with(gp,deriveMatrixIndCP(x,y,hyperpars,FALSE,1))
  diagDeriveCovmat <- as.vector(with(gp,deriveMatrixIndCP(x,y,hyperpars,TRUE,1)))
  expect_equal(diag(deriveCovmat),diagDeriveCovmat)
})



test_that("Element-wise and block-wise calculation of derivative wrt induction points with nugget parameter agree", {

  deriveCovmat <- with(gp,deriveMatrixIndCP(x,x,hyperpars,FALSE,1))
  diagDeriveCovmat <- as.vector(with(gp,deriveMatrixIndCP(x,x,hyperpars,TRUE,1)))
  expect_equal(diag(deriveCovmat),diagDeriveCovmat)
})



test_that("derivative of covfun wrt gp1$sigma hyperperameters is correctly computed",{

  refMatSigma1 <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"sigma1"))
  matSigma1 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(sigma=TRUE)),FALSE,TRUE))
  matSigma1x <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(sigma=TRUE)),FALSE,FALSE))$gp1$sigma
  expect_equal(matSigma1, refMatSigma1)
  expect_equal(matSigma1x, refMatSigma1)
})



test_that("derivative of covfun wrt gp1$sigma hyperperameter with nugget is correctly computed",{

  refMatSigma1 <- with(gp,numDeriveCovfun(x,x,hyperpars,FALSE,"sigma1"))
  matSigma1 <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp1=list(sigma=TRUE)),FALSE,TRUE))
  matSigma1x <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp1=list(sigma=TRUE)),FALSE,FALSE))$gp1$sigma
  expect_equal(matSigma1, refMatSigma1)
  expect_equal(matSigma1x, refMatSigma1)
})



test_that("elementwise and blockwise derivative of covfun wrt gp1$sigma hyperperameters agree",{

  matSigma1Elem <- as.vector(with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(sigma=TRUE)),TRUE,TRUE)))
  matSigma1Elemx <- as.vector(with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(sigma=TRUE)),TRUE,FALSE))$gp1$sigma)
  matSigma1Block <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(sigma=TRUE)),FALSE,TRUE))
  matSigma1Blockx <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(sigma=TRUE)),FALSE,FALSE))$gp1$sigma

  expect_equal(matSigma1Elem, diag(matSigma1Block))
  expect_equal(matSigma1Elemx, diag(matSigma1Blockx))
  expect_equal(matSigma1Elem, matSigma1Elemx)
})



test_that("derivative of covfun wrt gp1$len1 hyperperameter is correctly computed",{

  refMatLen11 <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"len11"))
  matLen11 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(len=c(TRUE,FALSE,FALSE))),FALSE,TRUE))
  matLen11x <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(len=c(TRUE,FALSE,FALSE))),FALSE,FALSE))$gp1$len[[1]]
  expect_equal(matLen11, refMatLen11)
  expect_equal(matLen11x, refMatLen11)
})



test_that("derivative of covfun wrt gp1$len1 hyperperameter with nugget is correctly computed",{

  refMatLen11 <- with(gp,numDeriveCovfun(x,x,hyperpars,FALSE,"len11"))
  matLen11 <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp1=list(len=c(TRUE,FALSE,FALSE))),FALSE,TRUE))
  matLen11x <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp1=list(len=c(TRUE,FALSE,FALSE))),FALSE,FALSE))$gp1$len[[1]]
  expect_equal(matLen11, refMatLen11)
  expect_equal(matLen11x, refMatLen11)
})



test_that("derivative of covfun wrt gp1$len2 hyperperameter is correctly computed",{

  refMatLen12 <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"len12"))
  matLen12 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),FALSE,TRUE))
  matLen12x <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),FALSE,FALSE))$gp1$len[[1]]
  expect_equal(matLen12, refMatLen12)
  expect_equal(matLen12x, refMatLen12)
})



test_that("derivative of covfun wrt gp1$len2 hyperperameter with nugget is correctly computed",{

  refMatLen12 <- with(gp,numDeriveCovfun(x,x,hyperpars,FALSE,"len12"))
  matLen12 <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),FALSE,TRUE))
  matLen12x <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),FALSE,FALSE))$gp1$len[[1]]
  expect_equal(matLen12, refMatLen12)
  expect_equal(matLen12x, refMatLen12)
})


test_that("elementwise and blockwise derivative of covfun wrt gp1$len2 hyperperameters agree",{

  matLen12Elem <- as.vector(with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),TRUE,TRUE)))
  matLen12Elemx <- as.vector(with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),TRUE,FALSE))$gp1$len[[1]])
  matLen12Block <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),FALSE,TRUE))
  matLen12Blockx <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),FALSE,FALSE))$gp1$len[[1]]

  expect_equal(matLen12Elem, diag(matLen12Block))
  expect_equal(matLen12Elemx, diag(matLen12Blockx))
  expect_equal(matLen12Elem, matLen12Elemx)
})



test_that("elementwise and blockwise derivative of covfun wrt gp1$len2 hyperperameters with nugget agree",{

  matLen12Elem <- as.vector(with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),TRUE,TRUE)))
  matLen12Elemx <- as.vector(with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),TRUE,FALSE))$gp1$len[[1]])
  matLen12Block <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),FALSE,TRUE))
  matLen12Blockx <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),FALSE,FALSE))$gp1$len[[1]]

  expect_equal(matLen12Elem, diag(matLen12Block))
  expect_equal(matLen12Elemx, diag(matLen12Blockx))
  expect_equal(matLen12Elem, matLen12Elemx)
})



test_that("derivative of covfun wrt gp2$sigma hyperperameter is correctly computed",{

  refMatSigma2 <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"sigma2"))
  matSigma2 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp2=list(sigma=TRUE)),FALSE,TRUE))
  matSigma2x <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp2=list(sigma=TRUE)),FALSE,FALSE))$gp2$sigma
  expect_equal(matSigma2, refMatSigma2)
  expect_equal(matSigma2x, refMatSigma2)
})



test_that("derivative of covfun wrt gp2$sigma hyperperameter with nugget is correctly computed",{

  refMatSigma2 <- with(gp,numDeriveCovfun(x,x,hyperpars,FALSE,"sigma2"))
  matSigma2 <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp2=list(sigma=TRUE)),FALSE,TRUE))
  matSigma2x <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp2=list(sigma=TRUE)),FALSE,FALSE))$gp2$sigma
  expect_equal(matSigma2, refMatSigma2)
  expect_equal(matSigma2x, refMatSigma2)
})



test_that("derivative of covfun wrt gp2$len1 hyperperameter is correctly computed",{

  refMatLen21 <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"len21"))
  matLen21 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp2=list(len=c(TRUE,FALSE,FALSE))),FALSE,TRUE))
  matLen21x <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp2=list(len=c(TRUE,FALSE,FALSE))),FALSE,FALSE))$gp2$len[[1]]
  expect_equal(matLen21, refMatLen21)
  expect_equal(matLen21x, refMatLen21)
})



test_that("derivative of covfun wrt gp2$len1 hyperperameter with nugget is correctly computed",{

  refMatLen21 <- with(gp,numDeriveCovfun(x,x,hyperpars,FALSE,"len21"))
  matLen21 <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp2=list(len=c(TRUE,FALSE,FALSE))),FALSE,TRUE))
  matLen21x <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp2=list(len=c(TRUE,FALSE,FALSE))),FALSE,FALSE))$gp2$len[[1]]
  expect_equal(matLen21, refMatLen21)
  expect_equal(matLen21x, refMatLen21)
})



test_that("derivative of covfun wrt gp2$len2 hyperperameter is correctly computed",{

  refMatLen22 <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"len22"))
  matLen22 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp2=list(len=c(FALSE,TRUE,FALSE))),FALSE,TRUE))
  matLen22x <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp2=list(len=c(FALSE,TRUE,FALSE))),FALSE,FALSE))$gp$len[[1]]
  expect_equal(matLen22, refMatLen22)
  expect_equal(matLen22x, refMatLen22)
})



test_that("derivative of covfun wrt gp2$len2 hyperperameter with nugget is correctly computed",{

  refMatLen22 <- with(gp,numDeriveCovfun(x,x,hyperpars,FALSE,"len22"))
  matLen22 <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp2=list(len=c(FALSE,TRUE,FALSE))),FALSE,TRUE))
  matLen22x <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(gp2=list(len=c(FALSE,TRUE,FALSE))),FALSE,FALSE))$gp$len[[1]]
  expect_equal(matLen22, refMatLen22)
  expect_equal(matLen22x, refMatLen22)
})



test_that("derivative of covfun wrt cp$kc hyperperameter is correctly computed",{

  refMatKc <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"kc"))
  matKc <-  with(gp,deriveMatrixHypCP(x,y,hyperpars,list(cp=list(kc=TRUE)),FALSE,TRUE))
  matKcx <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(cp=list(kc=TRUE)),FALSE,FALSE))$cp$kc
  expect_equal(matKc, refMatKc)
  expect_equal(matKcx, refMatKc)
})



test_that("derivative of covfun wrt cp$kc hyperperameter with nugget is correctly computed",{

  refMatKc <- with(gp,numDeriveCovfun(x,x,hyperpars,FALSE,"kc"))
  matKc <-  with(gp,deriveMatrixHypCP(x,x,hyperpars,list(cp=list(kc=TRUE)),FALSE,TRUE))
  matKcx <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(cp=list(kc=TRUE)),FALSE,FALSE))$cp$kc
  expect_equal(matKc, refMatKc)
  expect_equal(matKcx, refMatKc)
})



test_that("derivative of covfun wrt cp$xc hyperperameter is correctly computed",{

  refMatXc <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"xc"))
  matXc <-  with(gp,deriveMatrixHypCP(x,y,hyperpars,list(cp=list(xc=TRUE)),FALSE,TRUE))
  matXcx <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(cp=list(xc=TRUE)),FALSE,FALSE))$cp$xc
  expect_equal(matXc, refMatXc)
  expect_equal(matXcx, refMatXc)
})



test_that("derivative of covfun wrt cp$xc hyperperameter with nugget is correctly computed",{

  refMatXc <- with(gp,numDeriveCovfun(x,x,hyperpars,FALSE,"xc"))
  matXc <-  with(gp,deriveMatrixHypCP(x,x,hyperpars,list(cp=list(xc=TRUE)),FALSE,TRUE))
  matXcx <- with(gp,deriveMatrixHypCP(x,x,hyperpars,list(cp=list(xc=TRUE)),FALSE,FALSE))$cp$xc
  expect_equal(matXc, refMatXc)
  expect_equal(matXcx, refMatXc)
})



test_that("derivative of covfun wrt cp$ang1 hyperperameter is correctly computed",{

  refMatAng1 <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"ang1"))
  matAng1 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(cp=list(ang=c(TRUE,FALSE))),FALSE,TRUE))
  expect_equal(matAng1,refMatAng1)
})



test_that("derivative of covfun wrt cp$ang2 hyperperameter is correctly computed",{

  refMatAng2 <- with(gp,numDeriveCovfun(x,y,hyperpars,FALSE,"ang2"))
  matAng2 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(cp=list(ang=c(FALSE,TRUE))),FALSE,TRUE))
  expect_equal(matAng2,refMatAng2)
})


test_that("multiple calculation and single calculation of covfun derivative wrt hyperparameters agree", {

  matXc <-  with(gp,deriveMatrixHypCP(x,y,hyperpars,list(cp=list(xc=TRUE)),FALSE,TRUE))
  matKc <-  with(gp,deriveMatrixHypCP(x,y,hyperpars,list(cp=list(kc=TRUE)),FALSE,TRUE))
  matLen12 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp1=list(len=c(FALSE,TRUE,FALSE))),FALSE,TRUE))
  matLen21 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(gp2=list(len=c(TRUE,FALSE,FALSE))),FALSE,TRUE))
  matAng1 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(cp=list(ang=c(TRUE,FALSE))),FALSE,TRUE))
  matAng2 <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(cp=list(ang=c(FALSE,TRUE))),FALSE,TRUE))

  compRes <- with(gp,deriveMatrixHypCP(x,y,hyperpars,list(cp=list(xc=TRUE,kc=TRUE, ang=c(TRUE,TRUE)),
                                                          gp1=list(len=c(FALSE,TRUE,FALSE)),
                                                          gp2=list(len=c(TRUE,FALSE,FALSE))),FALSE,FALSE))
  expect_equal(matLen12, compRes$gp1$len[[1]])
  expect_equal(matLen21, compRes$gp2$len[[1]])
  expect_equal(matXc, compRes$cp$xc)
  expect_equal(matKc, compRes$cp$kc)
  expect_equal(matAng1, compRes$cp$ang[[1]])
  expect_equal(matAng2, compRes$cp$ang[[2]])
})




