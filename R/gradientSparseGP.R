

gradientSparseGPMarlike <- function(model, selection) {

  selnameChain <- character(0)
  redcursel <- list()
  depth <- 0
  recurseSelection <- function(curhyp, cursel) {
    stopifnot(all(names(cursel)%in%names(curhyp)))
    depth <<- depth + 1
    selnames <- names(cursel)
    partgrad <- cursel
    for (sn in selnames) {
      selnameChain[depth] <<- sn
      if (is.list(cursel[[sn]])) {
        redcursel[[selnameChain]] <<- list()
        partgrad[[sn]] <- recurseSelection(curhyp[[sn]], cursel[[sn]])
        redcursel[[selnameChain]] <<- NULL
      } else {
        idxvec <- if (sn!="ind") seq_along(curhyp[[sn]]) else seq_len(nrow(curhyp[[sn]]))
        idxvec <- idxvec[cursel[[sn]]]
        if (sn!="ind") {
          partgrad[[sn]] <- numeric(length(idxvec))
          logvec <- logical(length(curhyp[[sn]]))
          lastidx <- NULL
          for (i in seq_along(idxvec)) {
            curidx <- idxvec[i]
            logvec[lastidx] <- FALSE
            logvec[curidx] <- TRUE
            lastidx <- curidx
            redcursel[[selnameChain]] <<- logvec
            partgrad[[sn]][i] <- deriveSparseGPMarlikeHyp(model,redcursel)
          }
        } else { # induction points
          logvec <- logical(nrow(curhyp[[sn]]))
          logvec[idxvec] <- TRUE
          partgrad[[sn]] <- deriveSparseGPMarlikeIndSel(model,logvec)
        }
      }
      redcursel[[selnameChain]] <<- NULL
    }
    selnameChain <<- selnameChain[-depth]
    depth <<- depth - 1
    partgrad
  }

  hyperpars <- model$hyperpars
  recurseSelection(hyperpars,selection)
}





gradientSparseGPMarlikefD <- function(model, selection) {

  selnameChain <- character(0)
  redcursel <- list()
  depth <- 0
  recurseSelection <- function(curhyp, cursel) {
    stopifnot(all(names(cursel)%in%names(curhyp)))
    depth <<- depth + 1
    selnames <- names(cursel)
    partgrad <- cursel
    for (sn in selnames) {
      selnameChain[depth] <<- sn
      if (is.list(cursel[[sn]])) {
        redcursel[[selnameChain]] <<- list()
        partgrad[[sn]] <- recurseSelection(curhyp[[sn]], cursel[[sn]])
        redcursel[[selnameChain]] <<- NULL
      } else {
        idxvec <- if (sn!="ind") seq_along(curhyp[[sn]]) else seq_len(nrow(curhyp[[sn]]))
        idxvec <- idxvec[cursel[[sn]]]
        if (sn!="ind") {
          if (selnameChain[1]=="deltapars")
          {
            redcursel[[selnameChain]] <- idxvec
            partgrad[[sn]] <- as.vector(deriveSparseGPMarlikeDeltafD(model,redcursel))
          }
          else
          {
            partgrad[[sn]] <- numeric(length(idxvec))
            logvec <- logical(length(curhyp[[sn]]))
            lastidx <- NULL
            for (i in seq_along(idxvec)) {
              curidx <- idxvec[i]
              logvec[lastidx] <- FALSE
              logvec[curidx] <- TRUE
              lastidx <- curidx
              redcursel[[selnameChain]] <<- logvec
              partgrad[[sn]][i] <- deriveSparseGPMarlikeHypfD(model,redcursel)
            }
          }
        } else { # induction points
          logvec <- logical(nrow(curhyp[[sn]]))
          logvec[idxvec] <- TRUE
          partgrad[[sn]] <- deriveSparseGPMarlikeIndSel(model,logvec)
        }
      }
      redcursel[[selnameChain]] <<- NULL
    }
    selnameChain <<- selnameChain[-depth]
    depth <<- depth - 1
    partgrad
  }

  hyperpars <- model$hyperpars
  recurseSelection(hyperpars,selection)
}
