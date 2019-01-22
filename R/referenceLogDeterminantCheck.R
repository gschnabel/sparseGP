

# quick and dirty implementation of the log det component
# of the gradient of the sparse GP marginal likelihood

logDet <- function(gp) {
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(gp,fun(y,y,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  Z <- (t(Kuf)%*%invKuu%*%Kuf + delta)
  2*sum(log(diag(chol(Z))))
}

adjustedLogDet <- function(x,gp,ptIndex,dimIndex) {
  gp$x[ptIndex,dimIndex] <- gp$x[ptIndex,dimIndex] + x
  logDet(gp)
}

numDeriveLogDet <- function(gp,ptIndex,dimIndex) {
  res <- jacobian(adjustedLogDet,0,gp=gp,ptIndex=ptIndex,dimIndex=dimIndex)
  as.vector(res)
}

deriveLogDet1 <- function(gp,ptIndex,dimIndex) {
  
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(gp,fun(y,y,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  Z <- t(Kuf)%*%invKuu%*%Kuf + delta
  invZ <- chol2inv(chol(Z))
  
  deriveKuu <- with(gp,dfun(x,x,hyperpars,FALSE,dimIndex-1))
  deriveKuf <- with(gp,dfun(x,y,hyperpars,FALSE,dimIndex-1))
  
  curDeriveKuu <- filterMat(deriveKuu,ptIndex)
  curDeriveKuf <- filterMat(deriveKuf,ptIndex)
  curDeriveDelta <- diag(diag(-2*t(curDeriveKuf)%*%invKuu%*%Kuf + 2*t(Kuf)%*%invKuu%*%curDeriveKuu%*%invKuu%*%Kuf))
  
  deriveZ <- (2*(t(curDeriveKuf)%*%invKuu%*%Kuf - t(Kuf)%*%invKuu%*%curDeriveKuu%*%invKuu%*%Kuf) + curDeriveDelta)
  sum(diag(invZ %*% deriveZ))
}


deriveLogDet2 <- function(gp,ptIndex,dimIndex) {
  
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(gp,fun(y,y,hyperpars,FALSE))

  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  invDelta <- chol2inv(chol(delta))
  Z <- t(Kuf)%*%invKuu%*%Kuf + delta
  invZ <- chol2inv(chol(Z))
  X <- Kuu + Kuf%*%invDelta%*%t(Kuf)
  invX <- chol2inv(chol(X))
  
  deriveKuu <- with(gp,dfun(x,x,hyperpars,FALSE,dimIndex-1))
  deriveKuf <- with(gp,dfun(x,y,hyperpars,FALSE,dimIndex-1))
  
  curDeriveKuu <- filterMat(deriveKuu,ptIndex)
  curDeriveKuf <- filterMat(deriveKuf,ptIndex)
  curDeriveDelta <- diag(diag(-2*t(curDeriveKuf)%*%invKuu%*%Kuf + 2*t(Kuf)%*%invKuu%*%curDeriveKuu%*%invKuu%*%Kuf))
  
  dinvDelta <- -invDelta%*%curDeriveDelta%*%invDelta
  deriveX1 <- 2*curDeriveKuu + 2* curDeriveKuf%*%invDelta%*%t(Kuf)
  deriveX <- 2*curDeriveKuu + 2* curDeriveKuf%*%invDelta%*%t(Kuf) + Kuf%*%dinvDelta%*%t(Kuf)
  
  deriveLogDetKuu <- 2*sum(diag(invKuu %*% curDeriveKuu))
  deriveLogDetDelta <- sum(diag(invDelta %*% curDeriveDelta))
  deriveLogDetX <- sum(diag(invX %*% deriveX))
  
  # for debugging
  deriveLogDetX1 <- sum(diag(invX %*% deriveX1))
  diagKfuInvXKuf <- diag(t(Kuf)%*%invX%*%Kuf)
  
  print(list(deriveLogDetKuu=deriveLogDetKuu,
             deriveLogDetDelta=deriveLogDetDelta,
             deriveLogDetX=deriveLogDetX))
  
  -deriveLogDetKuu + deriveLogDetDelta + deriveLogDetX
}






