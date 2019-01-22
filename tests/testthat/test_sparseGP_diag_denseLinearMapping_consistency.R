library(sparseGP)
library(numDeriv)
context("test coherence of sparseGP_freeDelta with sparseGP")


gp_diag <- list(x = matrix(runif(9,0,5),3,3),
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


gp_free <- gp_diag
gp_free$hyperpars$deltaFuns <- createDenseLinearMapping()
gp_free$hyperpars$deltapars <- list(
               D = rep(0,3),
               S = matrix(c(3,1,1,1,5,1,1,1,7),3,3),
               P = data.frame(VAR=c(20,25,30))
             )


test_that("fitSparseGP and fitSparseGPfD yield the same result if S==0", {

  gp_free_mod <- gp_free
  gp_free_mod$hyperpars$deltapars$S[,] <- 0
  myfit_diag <- with(gp_diag, fitSparseGP(x, y, err, hyperpars))
  myfit_free <- with(gp_free_mod, fitSparseGPfD(x, y, err, hyperpars))
  expect_equal(myfit_diag$marlike, myfit_free$marlike)
  expect_equal(myfit_diag$maha, myfit_free$maha)
  expect_equal(myfit_diag$logDet, myfit_free$logDet)
})



test_that("fitSparseGP and fitSparseGPfD yield different results if S!=0", {

  gp_free$hyperpars$deltaFuns$resetCache()
  myfit_diag <- with(gp_diag, fitSparseGP(x, y, err, hyperpars))
  myfit_free <- with(gp_free, fitSparseGPfD(x, y, err, hyperpars))
  expect_false(isTRUE(all.equal(myfit_diag$marlike, myfit_free$marlike)))
  expect_false(isTRUE(all.equal(myfit_diag$maha,    myfit_free$maha)))
  expect_false(isTRUE(all.equal(myfit_diag$logDet,  myfit_free$logDet)))
})


test_that("gradients wrt hyperparameters of sparseGP and sparseGP_freedelte coincide if S==0",{

  gp_free$hyperpars$deltaFuns$resetCache()
  gp_free_mod <- gp_free
  gp_free_mod$hyperpars$deltapars$S[,] <- 0
  myfit_diag <- with(gp_diag, fitSparseGP(x, y, err, hyperpars))
  myfit_free <- with(gp_free_mod, fitSparseGPfD(x, y, err, hyperpars))

  resGrad1 <- gradientSparseGPMarlike(myfit_diag, list(hyperpars=list(sigma=TRUE, len=TRUE)))
  resGrad2 <- gradientSparseGPMarlikefD(myfit_free, list(hyperpars=list(sigma=TRUE, len=TRUE)))
  expect_equal(resGrad1, resGrad2)
})



test_that("gradients wrt hyperpars of sparseGP and sparseGP_freedelte different if S!=0",{

  gp_free$hyperpars$deltaFuns$resetCache()
  gp_free_mod <- gp_free
  myfit_diag <- with(gp_diag, fitSparseGP(x, y, err, hyperpars))
  myfit_free <- with(gp_free_mod, fitSparseGPfD(x, y, err, hyperpars))

  resGrad1 <- gradientSparseGPMarlike(myfit_diag, list(hyperpars=list(sigma=TRUE, len=TRUE)))
  resGrad2 <- gradientSparseGPMarlikefD(myfit_free, list(hyperpars=list(sigma=TRUE, len=TRUE)))
  expect_false(isTRUE(all.equal(resGrad1, resGrad2)))
})


test_that("gradient wrt deltapars vanishes if S==0", {

  gp_free$hyperpars$deltaFuns$resetCache()
  gp_free_mod <- gp_free
  gp_free_mod$hyperpars$deltapars$S[,] <- 0
  myfit_free <- with(gp_free_mod, fitSparseGPfD(x, y, err, hyperpars))
  resGrad <- gradientSparseGPMarlikefD(myfit_free, list(deltapars=list(P=list(VAR=TRUE))))
  expect_true(all(unlist(resGrad)==0))
})


test_that("hyperpars structure associated with sparseGP_freeDelta is really a specialization of sparseGP",{

  gp_free$hyperpars$deltaFuns$resetCache()
  myfit_diag <- with(gp_diag, fitSparseGP(x, y, err, hyperpars))
  myfit_free <- with(gp_free, fitSparseGP(x, y, err, hyperpars))
  expect_equal(myfit_diag$marlike, myfit_free$marlike)

  resGrad1 <- gradientSparseGPMarlike(myfit_diag, list(hyperpars=list(sigma=TRUE, len=TRUE)))
  resGrad2 <- gradientSparseGPMarlike(myfit_free, list(hyperpars=list(sigma=TRUE, len=TRUE)))
  expect_equal(resGrad1, resGrad2)
})



