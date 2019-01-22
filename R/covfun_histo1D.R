
getMatrixHisto <- function(x, y, hyperpars, elementMode) {

  stopifnot(length(hyperpars$binloc)-1==length(hyperpars$binvar))
  xidx <- findInterval(x, hyperpars$binloc) + 1
  yidx <- findInterval(y, hyperpars$binloc) + 1
  effvar <- c(0,hyperpars$binvar,0)

  if (elementMode) {
    stopifnot(length(x)==length(y))
    res <- numeric(length(x))
    res[xidx==yidx] <- effvar[xidx[xidx==yidx]]
  }
  else
  {
    res <- outer(effvar[xidx], effvar[yidx])
    res[outer(xidx,yidx,`!=`)] <- 0
  }
  res
}
