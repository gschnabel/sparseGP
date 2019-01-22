# under construction

trafoGP <- function(x, y, hyperpars, elementMode) {

  hyperpars$gpFuns$covfun(hyperpars$trafoFun(x),
                          hyperpars$trafoFun(y),
                          hyperpars$hyperpars,
                          elementMode)
}


deriveTrafoGPHyp <- function(x, y, hyperpars, elementMode) {

  hyperpars$gpFuns$deriveCovfunHyp(hyperpars$trafoFun(x),
                                   hyperpars$trafoFun(y),
                                   hyperpars$hyperpars,
                                   elementMode)
}

deriveTrafoGPInd <- function(x, y, hyperpars, elementMode) {

  hyperpars$gpFuns$deriveCovfunInd(hyperpars$trafoFun(x),
                                   hyperpars$trafoFun(y),
                                   hyperpars$hyperpars,
                                   elementMode)
  #TODO apply chain rule and add to the induction point
}
