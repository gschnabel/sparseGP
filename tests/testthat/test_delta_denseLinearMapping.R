library(sparseGP)
library(numDeriv)
context("test delta_denseLinearMapping functionality")


deltapars <- list(
  D = 1:10,
  S = matrix(runif(10*3),10,3),
  P = data.frame(VAR=c(3,5,7))
)

mapFuns <- createDenseLinearMapping()



test_that("logDet is correctly calculated", {
  refRes <- with(deltapars,determinant(diag(D) + S%*%diag(P$VAR)%*%t(S)),logarithm=TRUE)$modulus
  myRes <- mapFuns$logDet(rep(0,10), deltapars)
  expect_equivalent(myRes, refRes)
})



test_that("diagonal of inverse delta matrix is correctly calculated", {
  myDiagInv <- mapFuns$diagInv(rep(0,10), deltapars)
  refDiagInv <- with(deltapars,diag(chol2inv(chol((diag(D) + S%*%diag(P$VAR)%*%t(S))))))
  expect_equal(myDiagInv, refDiagInv)
})


test_that("multiplication from left with inverse delta matrix correctly calculatd", {
  x <- matrix(runif(10*5),10,5)
  refRes <- with(deltapars,chol2inv(chol((diag(D) + S%*%diag(P$VAR)%*%t(S)))) %*% x)
  myRes <- mapFuns$leftMultInv(x, rep(0,10), deltapars)
  expect_equal(myRes, refRes)
})



test_that("multiplication from right with inverse delta matrix correctly calculatd", {
  x <- matrix(runif(10*5),5,10)
  refRes <- with(deltapars,x %*% chol2inv(chol((diag(D) + S%*%diag(P$VAR)%*%t(S)))))
  myRes <- mapFuns$rightMultInv(x, rep(0,10), deltapars)
  expect_equal(myRes, refRes)
})


test_that("multiplying inverse delta matrix from left and right with the same matrix is correctly calculated", {
  x <- matrix(runif(10*5),10,5)
  refRes <- with(deltapars,t(x) %*% chol2inv(chol((diag(D) + S%*%diag(P$VAR)%*%t(S)))) %*% x)
  myRes <- mapFuns$bothMultInv(x, rep(0,10), deltapars)
  expect_equal(myRes, refRes)
})



test_that("trace of inverse delta matrix with derivative of delta matrix correctly computed", {
  invDelta <- with(deltapars,chol2inv(chol((diag(D) + S%*%diag(P$VAR)%*%t(S)))))
  refRes <- numeric(nrow(deltapars$P))
  for (i in seq_along(refRes)) {
    deriveDelta <- matrix(0,nrow(deltapars$P),nrow(deltapars$P))
    diag(deriveDelta)[i] <- 1
    deriveDelta <- deltapars$S%*%deriveDelta%*%t(deltapars$S)
    refRes[i] <- sum(diag(invDelta %*% deriveDelta))
  }
  myRes <- mapFuns$trInvDerive(rep(0,10), deltapars, list(P=list(VAR=c(TRUE,TRUE,TRUE))))
  expect_equal(myRes, refRes)
})



test_that("trace of inverse delta matrix with derivative of delta matrix matches numerical derivative", {

  wrapfun <- function(x, sel, addD) {
    newDeltapars <- deltapars
    newDeltapars$P$VAR[sel] <- deltapars$P$VAR[sel] + x
    res <- with(newDeltapars,2*sum(log(diag(chol((diag(D+addD) + S%*%diag(P$VAR)%*%t(S)))))))
    res
  }
  refRes <- numeric(nrow(deltapars$P))
  for (i in seq_along(refRes)) {
    refRes[i] <- as.vector(jacobian(wrapfun,x=0,sel=i, addD=seq_along(deltapars$D)))
  }
  myRes <- mapFuns$trInvDerive(seq_along(deltapars$D), deltapars, list(P=list(VAR=c(TRUE,TRUE,TRUE))))
  expect_equal(myRes, refRes)
})


test_that("trace of M x Ut x invDeriveDelta x U correctly computed",{
  tmp <- matrix(runif(5*5),5,5)
  M <- t(tmp) %*% tmp
  U <- matrix(runif(10*5),10,5)
  invDelta <- with(deltapars,chol2inv(chol((diag(D) + S%*%diag(P$VAR)%*%t(S)))))
  cholM <- chol(M)
  refRes <- numeric(nrow(deltapars$P))
  for (i in seq_along(refRes)) {
    deriveDelta <- matrix(0,nrow(deltapars$P),nrow(deltapars$P))
    diag(deriveDelta)[i] <- 1
    deriveDelta <- deltapars$S%*%deriveDelta%*%t(deltapars$S)
    refRes[i] <- sum(diag(M %*% t(U) %*% invDelta %*% deriveDelta %*% invDelta %*% U))
  }
  myRes <- mapFuns$trMxtUinvDeriveInvU(cholM %*% t(U), rep(0,10),deltapars, list(P=list(VAR=c(TRUE,TRUE,TRUE))))
  expect_equal(myRes, refRes)
})



test_that("multiplication of invDelta x deriveDelta x invDelta with vectors t(x) and y from left and right correctly computed", {
  x <- runif(10)
  y <- runif(10)
  invDelta <- with(deltapars,chol2inv(chol((diag(D) + S%*%diag(P$VAR)%*%t(S)))))
  refRes <- numeric(nrow(deltapars$P))
  for (i in seq_along(refRes)) {
    deriveDelta <- matrix(0,nrow(deltapars$P),nrow(deltapars$P))
    diag(deriveDelta)[i] <- 1
    deriveDelta <- deltapars$S%*%deriveDelta%*%t(deltapars$S)
    refRes[i] <- t(x) %*% invDelta %*% deriveDelta %*% invDelta %*% y
  }
  myRes <- mapFuns$multInvDeriveInvByVecs(x, y, rep(0,10), deltapars, list(P=list(VAR=c(TRUE,TRUE,TRUE))))
  expect_equal(myRes, refRes)
})


test_that("updating the cache works correctly", {
  # fill cache
  refRes <- with(deltapars,determinant(diag(D) + S%*%diag(P$VAR)%*%t(S)),logarithm=TRUE)$modulus
  myRes <- mapFuns$logDet(rep(0,10), deltapars)
  expect_equivalent(myRes, refRes)
  # change of P to clear cache
  deltapars$P[1] <- deltapars$P[1] + 1
  refRes <- with(deltapars,determinant(diag(D) + S%*%diag(P$VAR)%*%t(S)),logarithm=TRUE)$modulus
  myRes <- mapFuns$logDet(rep(0,10), deltapars)
  expect_equivalent(myRes, refRes)
  # change of D to clear cache
  deltapars$D[1] <- deltapars$D[1] + 1
  refRes <- with(deltapars,determinant(diag(D) + S%*%diag(P$VAR)%*%t(S)),logarithm=TRUE)$modulus
  myRes <- mapFuns$logDet(rep(0,10), deltapars)
  expect_equivalent(myRes, refRes)
})

