// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]


// helper functions

inline double invLogit(const double & x, const double & kc, const double & xc)
{
  return(1/(1+exp(-kc*(x-xc))));
}


inline double deriveInvLogit(const double & x, const double & kc, const double & xc)
{
  double expval = exp(-kc*(x-xc));
  if (!std::isinf(expval))
    return((expval/(1+expval))/(1+expval));
  else
    return 0.0;
}


// remark: the length of the vector order must always equals the number of dimensions minus 1
// assumption: radius from origin is fixed to 1
// x1 = cos(theta1)
// x2 = sin(theta1)*cos(theta2)
// x3 = sin(theta1)*sin(theta2)*cos(theta3)
// ...
// x(n-1) = sin(theta1)*sin(theta2)*...*sin(theta(n-2))*cos(theta(n-1))
// xn = sin(theta1)*sin(theta(n-1))
arma::vec sphereCoordToCart(const arma::vec & ang, const arma::uvec & order) {

  arma::vec angOrd = ang(order-1);
  arma::vec res(ang.n_elem+1,arma::fill::ones);

  for (int i=0; i<=ang.n_elem; i++)
  {
    if (i>0)
      res(i)  *= arma::prod(arma::sin(angOrd(arma::span(0,i-1))));
    if (i < ang.n_elem)
      res(i) *= cos(angOrd(i));
  }
  return(res);
}



arma::mat jacobianSphereCoordToCart(const arma::vec & ang, const arma::uvec & order) {

  arma::vec angOrd = ang(order-1);
  arma::vec fx = sphereCoordToCart(ang, order);
  arma::mat res(ang.n_elem+1, ang.n_elem, arma::fill::zeros);
  res = res.each_col() + fx;

  for (int i=0; i<=ang.n_elem; i++)
    for (int j=0; j<ang.n_elem; j++)
  {
      if (j > i)
      {
        res(i,j) = 0;
      }
      else if (j < i)
      {
        res(i,j) *= cos(angOrd(j)) / sin(angOrd(j));
      }
      else // i==j
      {
        res(i,j) *= -sin(angOrd(j)) / cos(angOrd(j));
      }
  }
  return(res);
}



arma::rowvec prepare_svec(const arma::mat & x,
                       const double & kc, const double & xc, const arma::rowvec & dir)
{
  arma::rowvec result(x.n_rows);
  double projdist;
  for (int i=0; i < x.n_rows; i++)
  {
    projdist = dot(arma::trans(x.row(i)), dir);
    result(i) = invLogit(projdist, kc, xc);
  }
  return result;
}



arma::rowvec prepare_dsvec(const arma::mat & x,
                          const double & kc, const double & xc, const arma::rowvec & dir)
{
  arma::rowvec result(x.n_rows);
  double projdist;
  for (int i=0; i < x.n_rows; i++)
  {
    projdist = dot(arma::trans(x.row(i)), dir);
    result(i) = deriveInvLogit(projdist, kc, xc);
  }
  return result;
}



void deriveWrtCP(const arma::rowvec & s1vec, const arma::rowvec & s2vec,
                      const arma::rowvec & ds1vecMod, const arma::rowvec & ds2vecMod,
                      const arma::mat & Kxy1, const arma::mat & Kxy2, bool & elementMode,
                      arma::mat & workMat)
{
  arma::rowvec r1vec = 1 - s1vec;
  arma::rowvec r2vec = 1 - s2vec;
  if (elementMode)
  {
    workMat  = -(ds1vecMod % r2vec + r1vec % ds2vecMod) % Kxy1;
    workMat +=  (ds1vecMod % s2vec + s1vec % ds2vecMod) % Kxy2;
  }
  else
  {
    workMat  = -(arma::trans(ds1vecMod) * r2vec + arma::trans(r1vec) * ds2vecMod) % Kxy1;
    workMat +=  (arma::trans(ds1vecMod) * s2vec + arma::trans(s1vec) * ds2vecMod) % Kxy2;
  }
}


// [[Rcpp::export]]
arma::mat getMatrixCP(arma::mat & x, arma::mat & y, Rcpp::List & hyperpars, bool & elementMode) {

  Rcpp::List hyperparsGP1 = hyperpars["gp1"];
  Rcpp::List hyperparsGP2 = hyperpars["gp2"];
  Rcpp::List hyperparsCP = hyperpars["cp"];
  Rcpp::List funsGP1 = hyperpars["gp1Funs"];
  Rcpp::List funsGP2 = hyperpars["gp2Funs"];
  Rcpp::Function covfunGP1 = funsGP1["covfun"];
  Rcpp::Function covfunGP2 = funsGP2["covfun"];

  double kc = hyperparsCP["kc"];
  double xc = hyperparsCP["xc"];

  arma::rowvec dir;
  if (x.n_cols > 1)
  {
    arma::uvec angOrd = hyperparsCP["ord"];
    arma::vec ang = hyperparsCP["ang"];
    dir = arma::trans(sphereCoordToCart(ang,angOrd));
  }
  else
  {
    dir.ones(1);
  }

  double nugget = hyperpars["nugget"];
  double nuggetSqr = nugget*nugget;

  arma::mat res;
  arma::mat Kxy1 = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfunGP1(x,y,hyperparsGP1,elementMode)));
  arma::mat Kxy2 = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfunGP2(x,y,hyperparsGP2,elementMode)));

  arma::rowvec s1vec = prepare_svec(x,kc,xc,dir);
  arma::rowvec s2vec = prepare_svec(y,kc,xc,dir);
  arma::rowvec r1vec = 1 - s1vec;
  arma::rowvec r2vec = 1 - s2vec;

  if (elementMode)
  {
    res = s1vec % s2vec % Kxy2 + r1vec % r2vec % Kxy1;
    for (int i=0; i < x.n_rows; i++)
      if (!any(x.row(i)-y.row(i))) res(0,i) += nuggetSqr;
  }
  else // elementMode==FALSE
  {
    res = (arma::trans(r1vec) * r2vec) % Kxy1 + (arma::trans(s1vec) * s2vec) % Kxy2;
    #pragma omp parallel for
    for (int j=0; j < y.n_rows; j++)
      for (int i=0; i < x.n_rows; i++)
        if (!any(x.row(i)-y.row(j))) res(i,j) += nuggetSqr;
  }
  return res;
}


// [[Rcpp::export]]
arma::mat deriveMatrixIndCP(arma::mat & x, arma::mat & y, Rcpp::List & hyperpars, bool & elementMode, int & dimIndex) {

  Rcpp::List hyperparsGP1 = hyperpars["gp1"];
  Rcpp::List hyperparsGP2 = hyperpars["gp2"];
  Rcpp::List hyperparsCP = hyperpars["cp"];

  Rcpp::List funsGP1 = hyperpars["gp1Funs"];
  Rcpp::List funsGP2 = hyperpars["gp2Funs"];

  Rcpp::Function covfunGP1 = funsGP1["covfun"];
  Rcpp::Function covfunGP2 = funsGP2["covfun"];
  Rcpp::Function deriveCovfunIndGP1 = funsGP1["deriveCovfunInd"];
  Rcpp::Function deriveCovfunIndGP2 = funsGP2["deriveCovfunInd"];

  double kc = hyperparsCP["kc"];
  double xc = hyperparsCP["xc"];

  arma::rowvec dir;
  if (x.n_cols > 1)
  {
    arma::uvec angOrd = hyperparsCP["ord"];
    arma::vec ang = hyperparsCP["ang"];
    dir = arma::trans(sphereCoordToCart(ang,angOrd));
  }
  else
  {
    dir.ones(1);
  }

  double nugget = hyperpars["nugget"];
  double nuggetSqr = nugget*nugget;

  int n = x.n_rows;
  int m = y.n_rows;

  double s1, s2, r1, r2, projdist1, projdist2;
  double derive_s1, derive_r1;
  double expval1, expval2;

  arma::mat res;
  arma::mat Kxy1 = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfunGP1(x,y,hyperparsGP1,elementMode)));
  arma::mat Kxy2 = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfunGP2(x,y,hyperparsGP2,elementMode)));
  arma::mat dKxy1 = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfunIndGP1(x,y,hyperparsGP1,elementMode,dimIndex)));
  arma::mat dKxy2 = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfunIndGP2(x,y,hyperparsGP2,elementMode,dimIndex)));

  arma::rowvec s1vec = prepare_svec(x,kc,xc,dir);
  arma::rowvec s2vec = prepare_svec(y,kc,xc,dir);
  arma::rowvec r1vec = 1 - s1vec;
  arma::rowvec r2vec = 1 - s2vec;
  arma::rowvec dsvec = prepare_dsvec(x,kc,xc,dir) * kc * dir(dimIndex);

  if (elementMode)
  {
    res  =  r1vec % r2vec % dKxy1;
    res -= dsvec % r2vec % Kxy1;
    res += s1vec % s2vec % dKxy2;
    res += dsvec % s2vec % Kxy2;
  }
  else // elementMode==FALSE
  {
    res =  (arma::trans(r1vec)*r2vec) % dKxy1;
    res -= (arma::trans(dsvec)*r2vec) % Kxy1;
    res += (arma::trans(s1vec)*s2vec) % dKxy2;
    res += (arma::trans(dsvec)*s2vec) % Kxy2;
  }
  return res;
}


void correctResultGP(const arma::mat & x, const arma::mat & y,
                     const double & kc, const double & xc, const arma::rowvec & dir,
                     bool & elementMode,
                     const int comp, Rcpp::List resultGP,
                     arma::rowvec & s1vec, arma::rowvec & s2vec) {

  for( Rcpp::List::iterator it = resultGP.begin(); it != resultGP.end(); ++it )
  {
    switch( TYPEOF(*it) ) {
    case VECSXP: {
      correctResultGP(x, y, kc, xc, dir, elementMode, comp,
                      *it, s1vec, s2vec);
      break;
    }
    case REALSXP: {
      if (s1vec.n_elem==0)
      {
        s1vec = prepare_svec(x,kc,xc,dir);
        s2vec = prepare_svec(y,kc,xc,dir);
        if (comp==1) {
          s1vec = 1 - s1vec;
          s2vec = 1 - s2vec;
        }
      }
      arma::mat cov = Rcpp::as<arma::mat>(*it);
      if (elementMode)
      {
        cov %= s1vec % s2vec;
      }
      else
      {
        cov %= (arma::trans(s1vec)*s2vec);
      }
      (*it) = Rcpp::wrap(cov);
      break;
    }
    default: {
      Rcpp::stop("incompatible SEXP encountered; only accepts lists with REALSXPs");
    }
    }
  }
}


// [[Rcpp::export]]
SEXP deriveMatrixHypCP(arma::mat & x, arma::mat & y, Rcpp::List & hyperpars, Rcpp::List & selection,
                     bool elementMode, bool singleDerive=false) {

  Rcpp::List res;

  Rcpp::List hyperparsGP1 = hyperpars["gp1"];
  Rcpp::List hyperparsGP2 = hyperpars["gp2"];
  Rcpp::List hyperparsCP = hyperpars["cp"];

  Rcpp::List funsGP1 = hyperpars["gp1Funs"];
  Rcpp::List funsGP2 = hyperpars["gp2Funs"];

  Rcpp::Function covfunGP1 = funsGP1["covfun"];
  Rcpp::Function covfunGP2 = funsGP2["covfun"];

  Rcpp::Function deriveCovfunHypGP1 = funsGP1["deriveCovfunHyp"];
  Rcpp::Function deriveCovfunHypGP2 = funsGP2["deriveCovfunHyp"];

  double kc = hyperparsCP["kc"];
  double xc = hyperparsCP["xc"];
  arma::uvec angOrd;
  arma::vec ang;
  arma::rowvec dir;
  if (x.n_cols > 1)
  {
    angOrd = Rcpp::as<arma::uvec>(hyperparsCP["ord"]);
    ang = Rcpp::as<arma::vec>(hyperparsCP["ang"]);
    dir = arma::trans(sphereCoordToCart(ang,angOrd));
  }
  else
  {
    dir.ones(1);
  }

  double nugget = hyperpars["nugget"];
  double nuggetSqr = nugget*nugget;

  Rcpp::List resultGP1;
  Rcpp::List resultGP2;

  arma::mat workMat;
  arma::rowvec s1vec;
  arma::rowvec s2vec;
  arma::rowvec r1vec;
  arma::rowvec r2vec;

  if (selection.size()==0)
  {
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
    if (hyperNames.at(idx) == ("gp1"))
    {
      Rcpp::List selectionGP1 = selection["gp1"];
      SEXP resultGP1 = deriveCovfunHypGP1(x,y,hyperparsGP1,selectionGP1,elementMode,singleDerive);
      if (singleDerive)
      {
        arma::mat matrixResultGP1 = Rcpp::as<arma::mat>(resultGP1);
        if (matrixResultGP1.n_elem > 0)
        {
          r1vec = 1 - prepare_svec(x,kc,xc,dir);
          r2vec = 1 - prepare_svec(y,kc,xc,dir);
          if (elementMode)
          {
            matrixResultGP1 %= r1vec % r2vec;
          }
          else
          {
            matrixResultGP1 %= (arma::trans(r1vec) * r2vec);
          }

          return(Rcpp::wrap(matrixResultGP1));
        }
      }
      else
      {
        res["gp1"] = resultGP1;
        correctResultGP(x, y, kc, xc, dir, elementMode, 1, resultGP1, r1vec, r2vec);
      }
    }

    if (hyperNames.at(idx) == ("gp2"))
    {
      Rcpp::List selectionGP2 = selection["gp2"];
      SEXP resultGP2 = deriveCovfunHypGP2(x,y,hyperparsGP2,selectionGP2,elementMode,singleDerive);
      if (singleDerive)
      {
        arma::mat matrixResultGP2 = Rcpp::as<arma::mat>(resultGP2);
        if (matrixResultGP2.n_elem > 0)
        {
          s1vec = prepare_svec(x,kc,xc,dir);
          s2vec = prepare_svec(y,kc,xc,dir);
          if (elementMode)
          {
            matrixResultGP2 %= s1vec % s2vec;
          }
          else
          {
            matrixResultGP2 %= (arma::trans(s1vec) * s2vec);
          }
          return(Rcpp::wrap(matrixResultGP2));
        }
      }
      else
      {
        res["gp2"] = resultGP2;
        correctResultGP(x, y, kc, xc, dir, elementMode, 2, resultGP2, s1vec, s2vec);
      }
    }


    if (hyperNames.at(idx)=="cp")
    {
      Rcpp::List resultCP;
      Rcpp::List selectionCP = selection["cp"];

      if (s1vec.n_elem==0) s1vec = prepare_svec(x,kc,xc,dir);
      if (s2vec.n_elem==0) s2vec = prepare_svec(y,kc,xc,dir);
      if (r1vec.n_elem==0) r1vec = 1-s1vec;
      if (r2vec.n_elem==0) r2vec = 1-s2vec;
      arma::mat ds1vec = prepare_dsvec(x,kc,xc,dir);
      arma::mat ds2vec = prepare_dsvec(y,kc,xc,dir);
      arma::mat Kxy1 = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfunGP1(x,y,hyperparsGP1,elementMode)));
      arma::mat Kxy2 = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfunGP2(x,y,hyperparsGP2,elementMode)));
      arma::mat ds1vecMod, ds2vecMod;

      Rcpp::CharacterVector hyperNamesCP = selectionCP.names();
      for (int idx=0; idx<hyperNamesCP.size(); idx++)
      {
        if (hyperNamesCP.at(idx) == "kc" || hyperNamesCP.at(idx) == "xc")
        {
          if (hyperNamesCP.at(idx) == "kc")
          {
            ds1vecMod = ds1vec % (arma::trans(sum(x.each_row() % dir,1)) - xc);
            ds2vecMod = ds2vec % (arma::trans(sum(y.each_row() % dir,1)) - xc);
          }
          else if (hyperNamesCP.at(idx) == "xc")
          {
            ds1vecMod = -ds1vec * kc;
            ds2vecMod = -ds2vec * kc;
          }

          deriveWrtCP(s1vec, s2vec, ds1vecMod, ds2vecMod, Kxy1, Kxy2, elementMode, workMat);

          if (singleDerive)
          {
            return Rcpp::wrap(workMat);
          }
          else
          {
            resultCP[Rcpp::as<std::string>(hyperNamesCP.at(idx))] = Rcpp::wrap(workMat);
          }
        }

        // special treatment for derivative of projection direction
        if (hyperNamesCP.at(idx) == "ang")
        {
          std::vector<Rcpp::NumericMatrix> deriveAngList;
          Rcpp::LogicalVector angSel = selectionCP["ang"];
          arma::mat jacmat = jacobianSphereCoordToCart(ang, angOrd);
          for (int i=0; i < angSel.length(); i++)
          {
            if (angSel(i))
            {
              ds1vecMod = (kc * ds1vec) % arma::trans(x * jacmat.col(i));
              ds2vecMod = (kc * ds2vec) % arma::trans(y * jacmat.col(i));

              deriveWrtCP(s1vec, s2vec, ds1vecMod, ds2vecMod, Kxy1, Kxy2, elementMode, workMat);

              if (singleDerive)
              {
                return Rcpp::wrap(workMat);
              }
              else
              {
                deriveAngList.push_back(Rcpp::wrap(workMat));
              }
            }
          }
          resultCP["ang"] = Rcpp::List(Rcpp::wrap(deriveAngList));
        }

      }
      res["cp"] = resultCP;
    }
  }

  if (singleDerive)
  {
    workMat.reset();
    return(Rcpp::wrap(workMat));
  }

  return res;
}



