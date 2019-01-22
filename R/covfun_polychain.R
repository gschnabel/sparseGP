
getMatrixPolyChain <- function(x, y, hyperpars, elementMode) {

  stopifnot(length(hyperpars$loc)==length(hyperpars$var))
  extloc <- hyperpars$loc
  extvar <- hyperpars$var

  xidx <- findInterval(x, extloc, rightmost.closed=TRUE)
  yidx <- findInterval(y, extloc, rightmost.closed=TRUE)
  xidx[xidx==0 | xidx==length(extloc)] <- NA_integer_
  yidx[yidx==0 | yidx==length(extloc)] <- NA_integer_

  xdiffvec <- extloc[xidx+1]-extloc[xidx]
  xloFact <- (extloc[xidx+1] - x) / xdiffvec
  xhiFact <- (x - extloc[xidx]) / xdiffvec

  xloVar <- extvar[xidx]
  xhiVar <- extvar[xidx+1]

  ydiffvec <- extloc[yidx+1]-extloc[yidx]
  yloFact <- (extloc[yidx+1] - y) / ydiffvec
  yhiFact <- (y - extloc[yidx]) / ydiffvec
  yloFact[!is.finite(yloFact)] <- 0 # to treat case Inf/Inf = NaN
  yhiFact[!is.finite(yhiFact)] <- 0

  yloVar <- extvar[yidx]
  yhiVar <- extvar[yidx+1]

  if (elementMode) {
    stopifnot(length(x)==length(y))
    res <- numeric(length(x))
    sel <- which(xidx==yidx)
    res[sel] <- (xloFact[sel]*yloFact[sel]*xloVar[sel] +
                   xhiFact[sel]*yhiFact[sel]*xhiVar[sel])

    sel <- which(xidx==yidx+1)
    sel <- sel[!is.na(sel)]
    res[sel] <- xloFact[sel]*yhiFact[sel]*xloVar[sel]

    sel <- which(xidx+1==yidx)
    sel <- sel[!is.na(sel)]
    res[sel] <- xhiFact[sel]*yloFact[sel]*yloVar[sel]
  }
  else
  {
    res <- matrix(0,length(x),length(y))
    f <- outer(xloFact*xloVar, yloFact) + outer(xhiFact*xhiVar, yhiFact)
    sel <- which(outer(xidx,yidx,`==`))
    res[sel] <- f[sel]

    f <- outer(xloFact*xloVar, yhiFact)
    sel <- which(outer(xidx,yidx+1,`==`))
    res[sel] <- f[sel]

    f <- outer(xhiFact, yloFact*yloVar)
    sel <- which(outer(xidx+1,yidx,`==`))
    res[sel] <- f[sel]
  }
  res
}

