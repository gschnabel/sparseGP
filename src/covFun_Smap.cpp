// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]


// [[Rcpp::export]]
arma::mat getMatrixSmap(Rcpp::DataFrame & x, Rcpp::DataFrame & y, Rcpp::List & hyperpars, bool & elementMode) {

  Rcpp::List chanGpFuns = hyperpars["gpFuns"];
  Rcpp::List chanHyperpars = hyperpars["hyperpars"];
  arma::mat sensMat = hyperpars["Smat"];

  Rcpp::CharacterVector sensRowNames = hyperpars["allChannels"];
  Rcpp::CharacterVector sensColNames = hyperpars["baseChannels"];
  Rcpp::CharacterVector xChans = x[0];
  Rcpp::CharacterVector yChans = y[0];
  arma::uvec sensRowIndices = Rcpp::as<arma::uvec>(Rcpp::match(xChans, sensRowNames)) - 1;
  arma::uvec sensColIndices = Rcpp::as<arma::uvec>(Rcpp::match(yChans, sensRowNames)) - 1;

  Rcpp::NumericMatrix redx(x.nrows(), x.size()-1);
  Rcpp::NumericMatrix redy(y.nrows(), y.size()-1);

  // skip first column that contains the
  for (int i=0; i<(x.size()-1);i++) {
    redx( Rcpp::_,i)=Rcpp::NumericVector(x[i+1]);
  }
  for (int i=0; i<(y.size()-1);i++) {
    redy( Rcpp::_,i)=Rcpp::NumericVector(y[i+1]);
  }

  arma::mat resultCov;
  arma::mat tmpCov;
  arma::uvec colIdx;
  Rcpp::List workHyperpars;
  Rcpp::List workGpFuns;

  if (!elementMode) resultCov.zeros(x.nrows(),y.nrows());
  else resultCov.zeros(1, x.nrows());

  for (int k=0; k<chanGpFuns.size(); k++) {
    colIdx = k;
    workGpFuns = chanGpFuns[k];
    Rcpp::Function covfun = workGpFuns["covfun"];

    workHyperpars = chanHyperpars[k];
    tmpCov = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(redx,redy,workHyperpars,elementMode)));

    if (!elementMode)
    {
      tmpCov = tmpCov.each_col() % sensMat.submat(sensRowIndices, colIdx);
      tmpCov = tmpCov.each_row() % arma::trans(sensMat.submat(sensColIndices, colIdx));
      resultCov += tmpCov;
    }
    else // element mode
    {
      tmpCov = tmpCov % arma::trans(sensMat.submat(sensRowIndices, colIdx));
      tmpCov = tmpCov % arma::trans(sensMat.submat(sensColIndices, colIdx));
      resultCov += tmpCov;
    }
  }
  return resultCov;
}


// // [[Rcpp::export]]
// SEXP deriveMatrixHypSmap(arma::mat & x, arma::mat & y, Rcpp::List & hyperpars, Rcpp::List & selection,
//                        bool elementMode, bool singleDerive=false) {
//
// }



