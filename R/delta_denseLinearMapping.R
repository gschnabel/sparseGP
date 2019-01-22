

createDenseLinearMapping <- function() {

  cache <- list()

  updateCache <- function(D, deltapars) {
    effD <- as.vector(D) + as.vector(deltapars$D)
    if (is.null(cache$PVAR) ||
        !identical(cache$PVAR, deltapars$P$VAR) ||
        !identical(cache$effD, effD)) {
      X <- t(deltapars$S) %*% (1/effD*deltapars$S)
      diag(X) <- diag(X) + 1/deltapars$P[["VAR"]]
      cholX <- chol(X)
      cache <<- list(
        PVAR = deltapars$P$VAR,
        effD = effD,
        cholX = cholX,
        invLtStinvD = forwardsolve(cholX, t((1/effD)*deltapars$S), transpose=TRUE, upper.tri=TRUE)
      )
    }
  }
  resetCache <- function() { cache <<- NULL }

  logDet <- function(D, deltapars) {
    updateCache(D, deltapars)
    sum(log(cache$PVAR)) + sum(log(cache$effD)) + 2*sum(log(diag(cache$cholX)))
  }

  diagInv <- function(D, deltapars) {
    updateCache(D, deltapars)
    1/cache$effD - colSums(cache$invLtStinvD*cache$invLtStinvD)
  }

  leftMultInv <- function(x, D, deltapars) {
    updateCache(D, deltapars)
    (1/cache$effD) * x - crossprod(cache$invLtStinvD, cache$invLtStinvD %*% x)
  }

  rightMultInv <- function(x, D, deltapars) {
    updateCache(D, deltapars)
    t((1/cache$effD)*t(x)) - tcrossprod(x, cache$invLtStinvD) %*% cache$invLtStinvD
  }

  bothMultInv <- function(x, D, deltapars) {
    updateCache(D, deltapars)
    crossprod(x, (1/cache$effD) * x) - crossprod(cache$invLtStinvD %*% x)
  }

  # derivatives
  trInvDerive <- function(D, deltapars, selection) {
    updateCache(D, deltapars)
    U <- cache$invLtStinvD %*% deltapars$S
    diagInvDelta1 <- colSums(deltapars$S * ((1/cache$effD) * deltapars$S))
    diagInvDelta2 <- colSums(U*U)
    (diagInvDelta1 - diagInvDelta2)[selection$P$VAR]
  }

  trMxtUinvDeriveInvU <- function(cholMUt, D, deltapars, selection) {
    updateCache(D, deltapars)
    U <- rightMultInv(cholMUt, D, deltapars) %*% deltapars$S
    diagVec <- colSums(U*U)
    diagVec[selection$P$VAR]
  }

  multInvDeriveInvByVecs = function(x, y, D, deltapars, selection) {
    updateCache(D, deltapars)
    u1 <- t(deltapars$S) %*% leftMultInv(x, D, deltapars)
    u2 <- t(deltapars$S) %*% leftMultInv(y, D, deltapars)
    (u1 * u2)[selection$P$VAR]
  }

  list(logDet=logDet,
       diagInv=diagInv,
       leftMultInv=leftMultInv,
       rightMultInv=rightMultInv,
       bothMultInv=bothMultInv,
       # derivatives
       trInvDerive=trInvDerive,
       trMxtUinvDeriveInvU=trMxtUinvDeriveInvU,
       multInvDeriveInvByVecs=multInvDeriveInvByVecs,
       # cache
       resetCache=resetCache)
}



