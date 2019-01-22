
# quick and dirty implementation of the mahalanobis component
# of the gradient of the sparse GP marginal likelihood

filterMat <- function(mat,i=NULL,j=NULL) {
  fMat <- mat
  if (!is.null(i)) fMat[-i,] <- 0
  if (!is.null(j)) fMat[,-j] <- 0
  fMat
}

# check derivative of invKuu

invKuu <- function(gp) {
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  chol2inv(chol(Kuu))
}

adjustedInvKuu <- function(x,gp,ptIndex,dimIndex) {
  gp$x[ptIndex,dimIndex] <- gp$x[ptIndex,dimIndex] + x
  invKuu(gp)
}

numDeriveInvKuu <- function(gp,ptIndex,dimIndex) {
  res <- jacobian(adjustedInvKuu,0,gp=gp,ptIndex=ptIndex,dimIndex=dimIndex)
  dim(res) <- c(nrow(gp$x),nrow(gp$x))
  res
}

deriveInvKuu <- function(gp,ptIndex,dimIndex) {
  deriveKuu <- with(gp,dfun(x,x,hyperpars,FALSE,dimIndex-1))
  deriveKuf <- with(gp,dfun(x,y,hyperpars,FALSE,dimIndex-1))
  curDeriveKuu <- filterMat(deriveKuu,ptIndex)
  curDeriveKuf <- filterMat(deriveKuf,ptIndex)
  -invKuu(gp) %*% (curDeriveKuu+t(curDeriveKuu)) %*% invKuu(gp)
}
### works ###

# check derivative of one part of delta
deltaPart1 <- function(gp,orgp) {
  Kuu <- with(orgp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(orgp,fun(y,y,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  delta
}

adjustedDeltaPart1 <- function(x,gp,ptIndex,dimIndex) {
  orgp <- gp
  gp$x[ptIndex,dimIndex] <- gp$x[ptIndex,dimIndex] + x
  deltaPart1(gp,orgp)
}

numDeriveDeltaPart1 <- function(gp,ptIndex,dimIndex) {
  res <- jacobian(adjustedDeltaPart1,0,gp=gp,ptIndex=ptIndex,dimIndex=dimIndex)
  dim(res) <- c(nrow(gp$y),nrow(gp$y))
  res
}

deriveDeltaPart1 <- function(gp,ptIndex,dimIndex) {
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(gp,fun(y,y,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  deriveKuu <- with(gp,dfun(x,x,hyperpars,FALSE,dimIndex-1))
  deriveKuf <- with(gp,dfun(x,y,hyperpars,FALSE,dimIndex-1))
  curDeriveKuu <- filterMat(deriveKuu,ptIndex)
  curDeriveKuf <- filterMat(deriveKuf,ptIndex)
  print(list(curDeriveKuf=curDeriveKuf,
             invKuu=invKuu, Kuf=Kuf)) #debug
  
  diag(diag(-2*t(curDeriveKuf)%*%invKuu%*%Kuf))
}
### works ###


# check derivative of delta
delta <- function(gp) {
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(gp,fun(y,y,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  delta
}

adjustedDelta <- function(x,gp,ptIndex,dimIndex) {
  gp$x[ptIndex,dimIndex] <- gp$x[ptIndex,dimIndex] + x
  delta(gp)
}

numDeriveDelta <- function(gp,ptIndex,dimIndex) {
  res <- jacobian(adjustedDelta,0,gp=gp,ptIndex=ptIndex,dimIndex=dimIndex)
  dim(res) <- c(nrow(gp$y),nrow(gp$y))
  res
}

deriveDelta <- function(gp,ptIndex,dimIndex) {
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(gp,fun(y,y,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  deriveKuu <- with(gp,dfun(x,x,hyperpars,FALSE,dimIndex-1))
  deriveKuf <- with(gp,dfun(x,y,hyperpars,FALSE,dimIndex-1))
  curDeriveKuu <- filterMat(deriveKuu,ptIndex)
  curDeriveKuf <- filterMat(deriveKuf,ptIndex)
  
  diag(diag(-t(curDeriveKuf)%*%invKuu%*%Kuf - t(Kuf)%*%invKuu%*%curDeriveKuf + 
              t(Kuf)%*%invKuu%*%(curDeriveKuu + t(curDeriveKuu))%*%invKuu%*%Kuf))
}
### works ###


# check derivative of invmaha

invmaha <- function(gp) {
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(gp,fun(y,y,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  Z <- t(Kuf)%*%invKuu%*%Kuf + delta
  as.vector(t(gp$obs) %*% Z %*% gp$obs)
}

adjustedInvmaha <- function(x,gp,ptIndex,dimIndex) {
  gp$x[ptIndex,dimIndex] <- gp$x[ptIndex,dimIndex] + x
  invmaha(gp)
}

numDeriveInvmaha <- function(gp,ptIndex,dimIndex) {
  jacobian(adjustedInvmaha,0,gp=gp,ptIndex=ptIndex,dimIndex=dimIndex)
}

deriveInvmaha <- function(gp,ptIndex,dimIndex) {
  
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(gp,fun(y,y,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  
  deriveKuu <- with(gp,dfun(x,x,hyperpars,FALSE,dimIndex-1))
  deriveKuf <- with(gp,dfun(x,y,hyperpars,FALSE,dimIndex-1))
  
  curDeriveKuu <- filterMat(deriveKuu,ptIndex)
  curDeriveKuf <- filterMat(deriveKuf,ptIndex)
  curDeriveDelta <- diag(diag(-t(curDeriveKuf)%*%invKuu%*%Kuf - t(Kuf)%*%invKuu%*%curDeriveKuf + 
                                t(Kuf)%*%invKuu%*%(curDeriveKuu + t(curDeriveKuu))%*%invKuu%*%Kuf))
  
  t(gp$obs) %*%
    (2*t(curDeriveKuf)%*%invKuu%*%Kuf + 
       (-2)*t(Kuf)%*%invKuu%*%curDeriveKuu%*%invKuu%*%Kuf + 
       curDeriveDelta) %*% gp$obs
}
### works ###

# check that maha works

maha <- function(gp) {
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(gp,fun(y,y,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  Z <- chol2inv(chol(t(Kuf)%*%invKuu%*%Kuf + delta))
  as.vector(t(gp$obs) %*% Z %*% gp$obs)
}

adjustedMaha <- function(x,gp,ptIndex,dimIndex) {
  gp$x[ptIndex,dimIndex] <- gp$x[ptIndex,dimIndex] + x
  maha(gp)
}

numDeriveMaha <- function(gp,ptIndex,dimIndex) {
  jacobian(adjustedMaha,0,gp=gp,ptIndex=ptIndex,dimIndex=dimIndex)
}



deriveMaha1 <- function(gp,ptIndex,dimIndex) {
  
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(gp,fun(y,y,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  Z <- chol2inv(chol(t(Kuf)%*%invKuu%*%Kuf + delta))
  
  deriveKuu <- with(gp,dfun(x,x,hyperpars,FALSE,dimIndex-1))
  deriveKuf <- with(gp,dfun(x,y,hyperpars,FALSE,dimIndex-1))
  
  curDeriveKuu <- filterMat(deriveKuu,ptIndex)
  curDeriveKuf <- filterMat(deriveKuf,ptIndex)
  curDeriveDelta <- diag(diag(-2*t(curDeriveKuf)%*%invKuu%*%Kuf + 2*t(Kuf)%*%invKuu%*%curDeriveKuu%*%invKuu%*%Kuf))
  
  -t(gp$obs) %*% Z %*%
    (2*(t(curDeriveKuf)%*%invKuu%*%Kuf - t(Kuf)%*%invKuu%*%curDeriveKuu%*%invKuu%*%Kuf) + curDeriveDelta) %*%
    Z %*% gp$obs
}


deriveMaha2 <- function(gp,ptIndex,dimIndex) {
  
  Kuu <- with(gp,fun(x,x,hyperpars,FALSE))
  Kuf <- with(gp,fun(x,y,hyperpars,FALSE))
  Kff <- with(gp,fun(y,y,hyperpars,FALSE))
  invKuu <- chol2inv(chol(Kuu))
  
  delta <- diag(diag(Kff-t(Kuf)%*%invKuu%*%Kuf) + gp$B)
  invDelta <- chol2inv(chol(delta))
  invX <- chol2inv(chol(Kuu + Kuf%*%invDelta%*%t(Kuf)))
  
  deriveKuu <- with(gp,dfun(x,x,hyperpars,FALSE,dimIndex-1))
  deriveKuf <- with(gp,dfun(x,y,hyperpars,FALSE,dimIndex-1))
  
  curDeriveKuu <- filterMat(deriveKuu,ptIndex)
  curDeriveKuf <- filterMat(deriveKuf,ptIndex)
  curDeriveDelta <- diag(diag(-2*t(curDeriveKuf)%*%invKuu%*%Kuf + 2*t(Kuf)%*%invKuu%*%curDeriveKuu%*%invKuu%*%Kuf))

  dinvDelta <- -invDelta%*%curDeriveDelta%*%invDelta
  deriveX <- 2*curDeriveKuu + 2* curDeriveKuf%*%invDelta%*%t(Kuf) + Kuf%*%dinvDelta%*%t(Kuf)
  
  invX_Kuf_invDelta <- invX%*%Kuf%*%invDelta
  
  term1 <- t(gp$obs)%*%dinvDelta%*%gp$obs
  term2 <- t(gp$obs)%*% (-2*dinvDelta %*% t(Kuf) %*% invX_Kuf_invDelta) %*% gp$obs
  term3 <- t(gp$obs)%*% (-2*invDelta %*% t(curDeriveKuf) %*% invX_Kuf_invDelta) %*% gp$obs
  term4 <- t(gp$obs)%*% (t(invX_Kuf_invDelta) %*% deriveX %*% invX_Kuf_invDelta) %*% gp$obs
  
  print(list(mahaTerm2=term2))
  
  term1 + term2 + term3 + term4
}
### works ###




