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
arma::mat getMatrix(arma::mat & x, arma::mat & y, Rcpp::List & hyperpars, bool & elementMode) {

  double sigma = hyperpars["sigma"];
  double sigmaSqr = sigma*sigma;
  arma::rowvec len = hyperpars["len"];
  double nugget = hyperpars["nugget"];
  double nuggetSqr = nugget*nugget;

  int n = x.n_rows;
  int m = y.n_rows;
  arma::rowvec z(len.size());
  double maha;
  arma::mat res;

  if (elementMode) {
    res.set_size(1,n);
    for (int i=0; i < n; i++) {
      z = (x.row(i) - y.row(i)) / len;
      maha = dot(z,z);
      res(0,i) = sigmaSqr * exp(-0.5*maha);
      if (maha==0) res(0,i) += nuggetSqr;
    }
  }
  else
  {
    res.set_size(n,m);
    #pragma omp parallel for private(z,maha)
    for (int j=0; j < m; j++) {
      for (int i=0; i < n; i++) {
        z = (x.row(i)-y.row(j))/len;
        maha = dot(z,z);
        res(i,j) = sigmaSqr*exp(-0.5*maha);
        if (maha==0) res(i,j) += nuggetSqr;
      }
    }
  }
  return res;
}

// [[Rcpp::export]]
arma::mat deriveMatrixInd(arma::mat & x, arma::mat & y, Rcpp::List & hyperpars, bool & elementMode, int & dimIndex) {

  double sigma = hyperpars["sigma"];
  double sigmaSqr = sigma*sigma;
  arma::rowvec len = hyperpars["len"];
  double nugget = hyperpars["nugget"];
  double nuggetSqr = nugget*nugget;

  int n = x.n_rows;
  int m = y.n_rows;
  int cDimIndex = dimIndex;
  arma::rowvec z(len.size());
  double maha;
  double innerDeriv;
  arma::mat res;

  if (dimIndex<0 || len.size()<=dimIndex) {
    throw std::range_error("Input dimension lower than dimIndex");
    return arma::mat(0,0);
  }

  if (elementMode) {
    res.set_size(1,n);
    for (int i=0; i < n; i++) {
      z = (x.row(i) - y.row(i)) / len;
      maha = dot(z,z);
      innerDeriv = (-z(cDimIndex)/len(cDimIndex));
      res(0,i) = sigmaSqr * exp(-0.5*maha) * innerDeriv;
    }
  }
  else
  {
    res.set_size(n,m);
    #pragma omp parallel for private(z,maha,innerDeriv)
    for (int j=0; j < m; j++) {
      for (int i=0; i < n; i++) {
        z = (x.row(i)-y.row(j))/len;
        maha = dot(z,z);
        innerDeriv = (-z(cDimIndex)/len(cDimIndex));
        res(i,j) = sigmaSqr * exp(-0.5*maha) * innerDeriv;
      }
    }
  }
  return res;
}



// [[Rcpp::export]]
SEXP deriveMatrixHyp(arma::mat & x, arma::mat & y, Rcpp::List & hyperpars, Rcpp::List & selection,
                     bool elementMode, bool singleDerive=false) {

  Rcpp::List res;

  // extract all hyperparameters
  double sigma = hyperpars["sigma"];
  double sigmaSqr = sigma*sigma;
  arma::rowvec len = hyperpars["len"];

  int n = x.n_rows;
  int m = y.n_rows;
  arma::rowvec z(len.size());
  double maha;
  double innerDeriv;
  arma::mat workMat;

  // this building block occurs in both derivative wrt sigma and length scale
  Rcpp::List modifiedHyperpars = Rcpp::clone(hyperpars);
  modifiedHyperpars["nugget"] = 0;
  arma::mat Kxy = getMatrix(x,y,modifiedHyperpars,elementMode);

  if (selection.size()==0) {
    if (singleDerive)
    {
      workMat.reset();
      return Rcpp::wrap(workMat);
    }
    else
    {
      return res;
    }
  }

  Rcpp::CharacterVector hyperNames = selection.names();
  for (int idx=0; idx<hyperNames.size(); idx++)
  {
    if (hyperNames.at(idx) == ("sigma"))
    {
      bool sigmaSelection = selection["sigma"];
      if (sigmaSelection)
      {
        if (elementMode) {
          workMat.set_size(1,n);
          for (int i=0; i < n; i++) {
            innerDeriv = 2/sigma;
            workMat(0,i) = Kxy(0,i) * innerDeriv;
          }
          if (singleDerive)
            return Rcpp::wrap(workMat);
          else
            res["sigma"] = Rcpp::wrap(workMat);
        }
        else
        {
          workMat.set_size(n,m);
#pragma omp parallel for private(innerDeriv)
          for (int j=0; j < m; j++) {
            for (int i=0; i < n; i++) {
              innerDeriv = 2/sigma;
              workMat(i,j) = Kxy(i,j) * innerDeriv;
            }
          }
          if (singleDerive)
            return Rcpp::wrap(workMat);
          else
            res["sigma"] = Rcpp::wrap(workMat);
        }
      }
    }

    if (hyperNames.at(idx) == ("len"))
    {
      Rcpp::LogicalVector lenSelection = selection["len"];
      std::vector<Rcpp::NumericMatrix> lenDeriveMats;
      for (int k=0; k<len.size(); k++) {

        if (lenSelection(k))
        {
          if (elementMode)
          {
            workMat.set_size(1,n);
            for (int i=0; i < n; i++)
            {
              innerDeriv = ((x(i,k) - y(i,k))/len(k));
              innerDeriv *= innerDeriv / len(k);
              workMat(0,i) = Kxy(0,i) * innerDeriv;
            }
            if (singleDerive)
              return Rcpp::wrap(workMat);
          }
          else
          {
            workMat.set_size(n,m);
#pragma omp parallel for private(innerDeriv)
            for (int j=0; j < m; j++) {
              for (int i=0; i < n; i++) {
                innerDeriv = ((x(i,k) - y(j,k))/len(k));
                innerDeriv *= innerDeriv / len(k);
                workMat(i,j) = Kxy(i,j) * innerDeriv;
              }
            }
            if (singleDerive)
              return Rcpp::wrap(workMat);
          }
          lenDeriveMats.push_back(Rcpp::wrap(workMat));
        }
      }
      res["len"] = Rcpp::List(Rcpp::wrap(lenDeriveMats));
    }
  }

  if (singleDerive) {
    workMat.reset();
    return(Rcpp::wrap(workMat));
  }
  else
    return res;
}


