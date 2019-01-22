


#' Normalize selection list
#'
#' Converts every element of a selection list for hyperparameters into a boolean vector
#'
#' @param hyperpars List that contains hyperparameter specifications
#' @param selection List that contains selection information for hyperparameters
#'
#' @return
#' Normalized list, i.e. integer indexes are replaced by boolean vectors
#'
#' @example
#' hyperpars <- list(sigma=10, len=c(15,20,25,30))
#' selection <- list(sigma=TRUE, len=2:3)
#' normSelection <- normalizeSelection(selection)
#'
#' @export
normalizeSelection <- function (hyperpars, selection)
{
  stopifnot(is.list(hyperpars))
  hypnames <- names(hyperpars)
  selnames <- names(selection)
  selnames <- selnames[nzchar(selnames)]
  for (v in hypnames) {
    hyperpars[[v]] <- if (is.list(hyperpars[[v]])) {
      normalizeSelection(hyperpars[[v]],selection[[v]])
    }
    else {
      if (is.null(selection[[v]])) {
        rep(FALSE, if (v!="ind") length(hyperpars[[v]]) else nrow(hyperpars[[v]]))
      }
      else {
        logvec <- logical(if (v!="ind") length(hyperpars[[v]]) else nrow(hyperpars[[v]]))
        if (is.logical(selection[[v]])) {
          logvec | selection[[v]]
        } else {
          logvec[selection[[v]]] <- TRUE
          logvec
        }
      }
    }
  }
  hyperpars
}

filterHyperpars <- function(hyperpars, selection) {

  stopifnot(is.list(hyperpars) ||
              is.environment(hyperpars))
  hypnames <- names(hyperpars)
  selnames <- names(selection)
  stopifnot(all(nzchar(hypnames)))

  for (v in hypnames) {
    hyperpars[[v]] <- if (is.null(selection[[v]])) {
      NULL
    }
    else if (is.list(hyperpars[[v]])) {
      filterHyperpars(hyperpars[[v]],selection[[v]])
    }
    else {
      if (is.null(dim(hyperpars[[v]])))
        hyperpars[[v]][selection[[v]]]
      else # in case we have a matrix
        hyperpars[[v]][selection[[v]],]
    }
  }
  if (length(hyperpars)==0)
    hyperpars <- NULL

  hyperpars
}


replaceHyperpars <- function(hyperpars, newHyperpars, selection=NULL) {

  stopifnot(is.list(hyperpars))
  stopifnot(all(nzchar(hyperpars)))
  hypnames <- names(hyperpars)
  selnames <- names(selection)

  for (v in hypnames) {
    if (is.list(hyperpars[[v]])) {
      hyperpars[[v]] <- replaceHyperpars(hyperpars[[v]], newHyperpars[[v]], selection[[v]])
    }
    else if (!is.null(newHyperpars[[v]])) {
      stopifnot(!is.null(hyperpars[[v]]),!is.null(newHyperpars[[v]]))
      if (is.null(selection[[v]])) {
        stopifnot(length(hyperpars[[v]])==length(newHyperpars[[v]]))
        hyperpars[[v]][] <- newHyperpars[[v]]
      }
      else {
        if (v!="ind")
          hyperpars[[v]][selection[[v]]] <- newHyperpars[[v]]
        else
          hyperpars[[v]][selection[[v]],] <- newHyperpars[[v]]
      }
    }
  }
  hyperpars
}
