// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
Rcpp::List fitSparseGP(arma::mat & x, arma::colvec & y, arma::colvec & err, Rcpp::List & hyperpars) {

  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::List hyperparsGP = hyperpars["hyperpars"];
  arma::mat ind = hyperpars["ind"];

  arma::mat Kuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,x,hyperparsGP,bool(FALSE))));
  arma::mat Kuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,ind,hyperparsGP,bool(FALSE))));
  arma::vec diagKff = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::NumericVector>(covfun(x,x,hyperparsGP,bool(TRUE))));

  // prepare invKuu
  arma::mat cholKuu = arma::chol(Kuu);

  // prepare delta
  arma::mat z = solve(trimatl(arma::trans(cholKuu)),Kuf);
  arma::colvec delta = diagKff - arma::trans(sum(z%z)) + err%err;
  arma::colvec invDelta = 1/delta;

  // prepare X
  //arma::mat cholX = chol(Kuu + Kuf * arma::diagmat(invDelta) * Kuf.t());
  arma::mat cholX = chol(Kuu + (Kuf.each_row() % arma::trans(invDelta)) * Kuf.t());
  arma::colvec invDelta_y = invDelta % y;
  arma::colvec Kuf_invDelta_y = Kuf * invDelta_y;
  arma::colvec invX_Kuf_invDelta_y = solve(arma::trimatu(cholX),
                                           solve(arma::trimatl(arma::trans(cholX)),Kuf_invDelta_y));

  // because it is cheap, calculate the log-marginal likelihood here
  // calculate the mahalanobis distance
  double yInvDeltay = dot(y,invDelta_y);
  arma::colvec sqrtMaha = solve(arma::trimatl(arma::trans(cholX)), Kuf_invDelta_y);
  double partMaha = dot(sqrtMaha,sqrtMaha);
  double maha = yInvDeltay - partMaha;

  // calculate the log determinant
  double logDetKuu = 2*arma::sum(log(cholKuu.diag()));
  double logDetDelta = arma::sum(log(delta));
  double logDetX = 2*arma::sum(log(cholX.diag()));
  double logDet =  logDetDelta + logDetX - logDetKuu;

  // all together
  double marlike = (-0.5) * (maha + logDet + x.n_rows*log2pi);

  return Rcpp::List::create(Rcpp::Named("invX_Kuf_invDelta_y") = invX_Kuf_invDelta_y,
                            Rcpp::Named("cholKuu") = cholKuu,
                            Rcpp::Named("cholX") = cholX,
                            Rcpp::Named("delta") = delta,
                            Rcpp::Named("covfun") = covfun,
                            Rcpp::Named("hyperpars") = hyperpars,
                            Rcpp::Named("x") = x,
                            Rcpp::Named("y") = y,
                            Rcpp::Named("err") = err,
                            Rcpp::Named("logDet") = logDet,
                            Rcpp::Named("maha") = maha,
                            Rcpp::Named("marlike") = marlike);
}


// [[Rcpp::export]]
Rcpp::DataFrame predictSparseGP(Rcpp::List & model, arma::mat & px, bool reterr = true) {

  Rcpp::List hyperpars = model["hyperpars"];
  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::List hyperparsGP = hyperpars["hyperpars"];
  arma::mat ind = hyperpars["ind"];

  arma::colvec invX_Kuf_invDelta_y = model["invX_Kuf_invDelta_y"];

  // prediction
  arma::mat Kup = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,px,hyperparsGP,bool(FALSE))));
  arma::vec py = arma::trans(Kup) * invX_Kuf_invDelta_y;

  arma::colvec predVar(px.n_rows,arma::fill::zeros);
  arma::mat covmat;
  if (reterr) {
    arma::mat cholKuu = model["cholKuu"];
    arma::mat cholX = model["cholX"];

    arma::vec diagKpp = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::NumericVector>(covfun(px,px,hyperparsGP,bool(TRUE))));
    arma::mat z = arma::solve(arma::trans(cholKuu),Kup);
    arma::mat diagQpp = sum(arma::trans(z%z),1);
    z = arma::solve(arma::trans(cholX),Kup);
    predVar = diagKpp - diagQpp + sum(arma::trans(z%z),1);
  }

  return Rcpp::DataFrame::create(Rcpp::Named("x")=px,
                                                  Rcpp::Named("y")=py,
                                                  Rcpp::Named("err")=arma::sqrt(predVar));
}


// [[Rcpp::export]]
arma::mat getPostCovFromSparseGP(arma::mat & px, arma::mat & py, Rcpp::List & model,
                                 bool elementMode, bool exact = false) {

  Rcpp::List hyperpars = model["hyperpars"];
  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::List hyperparsGP = hyperpars["hyperpars"];
  arma::mat ind = hyperpars["ind"];

  arma::mat cholKuu = model["cholKuu"];
  arma::mat cholX = model["cholX"];

  arma::mat Kupx = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,px,hyperparsGP,bool(FALSE))));
  arma::mat Kupy = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,py,hyperparsGP,bool(FALSE))));

  arma::mat ux = arma::solve(arma::trans(cholKuu),Kupx);
  arma::mat uy = arma::solve(arma::trans(cholKuu),Kupy);

  arma::mat zx = arma::solve(arma::trans(cholX),Kupx);
  arma::mat zy = arma::solve(arma::trans(cholX),Kupy);

  arma::mat covmat;

  if (!elementMode)
  {
    if (!exact)
    {
      covmat = arma::trans(zx) * zy;
      // WARNING: activation of omp causes crash (called from R); possibly interaction with covfun call?
      //#pragma omp parallel
      {
        arma::uvec isDiag(px.n_rows);
        arma::rowvec diagQpp;
        arma::rowvec diagKff;
        arma::uvec diagIdx;
        arma::mat diagPx;
        //#pragma omp for
        for (int j=0; j < py.n_rows; j++)
        {
          for (int i=0; i < px.n_rows; i++)
          {
            isDiag(i) = !any(px.row(i)-py.row(j));
          }
          diagIdx = arma::find(isDiag);
          if (diagIdx.n_elem > 0)
          {
            diagQpp = arma::trans(uy.col(j)) * ux.cols(diagIdx);
            diagPx = px.rows(diagIdx);
            diagKff = Rcpp::as<arma::rowvec>(Rcpp::as<Rcpp::NumericMatrix>(covfun(diagPx,diagPx,hyperparsGP,bool(TRUE))));
            covmat(j*px.n_rows + diagIdx) += diagKff - diagQpp;
          }
        }
      }
    }
    else
    {
      arma::mat Kpxpy = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(px,py,hyperparsGP,bool(FALSE))));
      covmat = Kpxpy - ux.t()*uy + zx.t()*zy;
    }
  }
  else // elementMode
  {
    arma::uvec diagIdx;
    arma::mat diagPx;

    if (!exact)
    {
      arma::vec isDiag(px.n_rows, arma::fill::zeros);
      covmat = sum(zx % zy,0);
      for (int i=0; i < px.n_rows; i++)
        isDiag(i) = !any(px.row(i)-py.row(i));
      diagIdx = arma::find(isDiag);
      diagPx = px.rows(diagIdx);
      covmat(diagIdx) += Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(
        covfun(diagPx, diagPx, hyperparsGP,bool(TRUE))));
      covmat.cols(diagIdx) -= sum(ux.cols(diagIdx) % uy.cols(diagIdx),0);
    }
    else // exact
    {
      arma::rowvec Kpxpy = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(px,py,hyperparsGP,bool(TRUE))));
      covmat = Kpxpy - sum(ux % uy,0) + sum(zx % zy,0);
    }
  }

  return covmat;
}

// [[Rcpp::export]]
Rcpp::List preparePostSamplingFromSparseGP(Rcpp::List & model, arma::mat & px, bool exact = false) {

  Rcpp::List hyperpars = model["hyperpars"];
  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::List hyperparsGP = hyperpars["hyperpars"];
  arma::mat ind = hyperpars["ind"];

  arma::mat cholKuu = model["cholKuu"];
  arma::mat cholX = model["cholX"];

  arma::colvec invX_Kuf_invDelta_y = model["invX_Kuf_invDelta_y"];

  // prediction
  arma::mat Kup = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,px,hyperparsGP,bool(FALSE))));
  arma::vec py = arma::trans(Kup) * invX_Kuf_invDelta_y;

  arma::vec diagKpp = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::NumericVector>(covfun(px,px,hyperparsGP,bool(TRUE))));
  arma::mat u = arma::solve(arma::trans(cholKuu),Kup);
  arma::mat z = arma::solve(arma::trans(cholX),Kup);
  arma::vec diagQpp = sum(arma::trans(u%u),1);

  arma::mat cholPostKpp;
  if (exact) {
    arma::mat Kpp = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(px,px,hyperparsGP,bool(FALSE))));
    cholPostKpp = arma::chol(Kpp - u.t()*u + z.t()*z);
  }

  arma::vec diagKppMinusQpp = diagKpp - diagQpp;
  arma::vec sqrtDiagKppMinusQpp(diagKpp.size(),arma::fill::zeros);
  arma::uvec notNegSel = find(diagKppMinusQpp>=0);
  sqrtDiagKppMinusQpp.elem(notNegSel) = arma::sqrt(diagKppMinusQpp.elem(notNegSel));

  return Rcpp::List::create(Rcpp::Named("diagKpp") = arma::trans(diagKpp),
                            Rcpp::Named("diagQpp") = arma::trans(diagQpp),
                            Rcpp::Named("predMean") = arma::trans(py),
                            Rcpp::Named("sqrtDiagKppMinusQpp") = arma::trans(sqrtDiagKppMinusQpp),
                            Rcpp::Named("z") = z,
                            Rcpp::Named("cholPostKpp") = cholPostKpp);
}

// [[Rcpp::export]]
arma::mat getSamplesFromSparseGP(Rcpp::List sampleSpec, int num, bool exact=false) {

  arma::rowvec sqrtDiagKppMinusQpp = sampleSpec["sqrtDiagKppMinusQpp"];
  arma::mat cholPostKpp = sampleSpec["cholPostKpp"];
  arma::mat z = sampleSpec["z"];
  arma::rowvec predMean = sampleSpec["predMean"];
  arma::mat obsmat;
  if (!exact) {
    obsmat = arma::randn(num, sqrtDiagKppMinusQpp.n_elem);
    obsmat = obsmat.each_row() % sqrtDiagKppMinusQpp;
    obsmat += arma::randn(num, z.n_rows) * z;
  }
  else
  {
    obsmat = arma::randn(num, cholPostKpp.n_rows) * cholPostKpp;
  }

  return (obsmat.each_row() + predMean);
}



// [[Rcpp::export]]
arma::mat deriveSparseGPMarlikeInd(Rcpp::List & model) {

  Rcpp::List hyperpars = model["hyperpars"];
  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::Function deriveCovfun = gpFuns["deriveCovfunInd"];
  Rcpp::List hyperparsGP = hyperpars["hyperpars"];
  arma::mat ind = hyperpars["ind"];

  arma::mat x = model["x"];
  arma::colvec y = model["y"];
  arma::colvec delta = model["delta"];

  arma::mat cholX = model["cholX"];
  arma::mat cholKuu = model["cholKuu"];
  arma::colvec invX_Kuf_invDelta_y = model["invX_Kuf_invDelta_y"];

  arma::mat Kuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,ind,hyperparsGP,bool(FALSE))));
  arma::vec diagKff = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::NumericVector>(covfun(x,x,hyperparsGP,bool(TRUE))));
  arma::colvec invDelta = 1/delta;

  // prepare auxiliary quantities to maintain O(d*m^2*n) for complete gradient
  // at expense of additional memory requirement
  arma::mat Kuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,x,hyperparsGP,bool(FALSE))));
  arma::mat invKuuKuf = solve(arma::trimatu(cholKuu),solve(arma::trimatl(arma::trans(cholKuu)),Kuf));

  arma::colvec Kfu_invX_Kuf_invDelta_y = arma::trans(Kuf) * invX_Kuf_invDelta_y;
  arma::mat diagKfuInvXKuf = solve(arma::trimatl(arma::trans(cholX)),Kuf);
  diagKfuInvXKuf = arma::trans(sum(diagKfuInvXKuf % diagKfuInvXKuf));
  arma::mat invKuu = arma::inv(arma::trimatu(cholKuu));
  invKuu = arma::trimatu(invKuu) * arma::trimatl(arma::trans(invKuu));
  arma::mat invX = arma::inv(arma::trimatu(cholX));
  invX = arma::trimatu(invX) * arma::trimatl(arma::trans(invX));

  // variables that change according to current derivative
  arma::mat gradMarginalLikelihood(ind.n_rows,ind.n_cols,arma::fill::zeros);

  arma::colvec deriveLogDetKuu(ind.n_rows);
  arma::colvec deriveLogDetXMat(ind.n_rows);
  arma::colvec deriveLogDetDeltaMat(ind.n_rows);
  arma::colvec deriveLogDetMat(ind.n_rows);
  arma::colvec deriveMahaMat(ind.n_rows);

  // memory intense
  arma::mat deriveDeltaMat; //(ind.n_rows,x.n_rows);
  arma::mat deriveInvDeltaMat; //(ind.n_rows, x.n_rows);
  arma::mat deriveX1Mat; //(ind.n_rows,ind.n_rows);

  for (int j=0; j<ind.n_cols; j++) {

    arma::mat deriveKuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(ind,x,hyperparsGP,bool(FALSE),j)));
    arma::mat deriveKuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(ind,ind,hyperparsGP,bool(FALSE),j)));

    deriveLogDetKuu = 2*sum(invKuu % deriveKuu,1);
    deriveDeltaMat = 2*(-deriveKuf + (deriveKuu * invKuuKuf)) % invKuuKuf;
    deriveInvDeltaMat = deriveDeltaMat.each_row() % arma::trans(-invDelta % invDelta);

    deriveX1Mat = 2*(deriveKuu + (deriveKuf.each_row() % arma::trans(invDelta)) * arma::trans(Kuf));

    deriveLogDetXMat = sum(invX % deriveX1Mat,1);
    deriveLogDetXMat += sum(deriveInvDeltaMat.each_row() % arma::trans(diagKfuInvXKuf),1);
    deriveLogDetDeltaMat = sum(deriveDeltaMat.each_row() % arma::trans(invDelta),1);

    deriveLogDetMat = -deriveLogDetKuu + deriveLogDetDeltaMat + deriveLogDetXMat;

    deriveMahaMat = sum(deriveInvDeltaMat.each_row() % arma::trans(y % y),1);
    deriveMahaMat += (deriveInvDeltaMat.each_row() % arma::trans(y)) * arma::trans(Kuf) * (-2)*invX_Kuf_invDelta_y;
    deriveMahaMat += (deriveKuf * (y % invDelta)) % ((-2)*invX_Kuf_invDelta_y);
    deriveMahaMat += (deriveX1Mat * invX_Kuf_invDelta_y) % invX_Kuf_invDelta_y;
    deriveMahaMat += (deriveInvDeltaMat.each_row() % arma::trans(Kfu_invX_Kuf_invDelta_y)) * Kfu_invX_Kuf_invDelta_y;

    gradMarginalLikelihood.col(j) = -0.5 * (deriveLogDetMat + deriveMahaMat);
  }
  return gradMarginalLikelihood;
}




// [[Rcpp::export]]
double deriveSparseGPMarlikeHyp(Rcpp::List & model, Rcpp::List & selection) {

  Rcpp::List hyperpars = model["hyperpars"];
  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::Function deriveCovfun = gpFuns["deriveCovfunHyp"];
  Rcpp::List selectionGP = selection["hyperpars"];
  Rcpp::List hyperparsGP = hyperpars["hyperpars"];
  arma::mat ind = hyperpars["ind"];

  arma::mat x = model["x"];
  arma::colvec y = model["y"];
  arma::colvec delta = model["delta"];

  arma::mat cholX = model["cholX"];
  arma::mat cholKuu = model["cholKuu"];
  arma::colvec invX_Kuf_invDelta_y = model["invX_Kuf_invDelta_y"];

  arma::mat Kuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,ind,hyperparsGP,bool(FALSE))));
  arma::vec diagKff = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::NumericVector>(covfun(x,x,hyperparsGP,bool(TRUE))));
  arma::colvec invDelta = 1/delta;

  arma::mat deriveKuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(ind,ind,hyperparsGP,selectionGP,false,true)));
  arma::mat deriveKuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(ind,x,hyperparsGP,selectionGP,false,true)));
  arma::mat deriveKff = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(x,x,hyperparsGP,selectionGP,true,true)));

  // prepare auxiliary quantities to maintain O(d*m^2*n) for complete gradient
  // at expense of additional memory requirement
  arma::mat Kuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,x,hyperparsGP,bool(FALSE))));
  arma::mat invKuuKuf = solve(arma::trimatu(cholKuu),solve(arma::trimatl(arma::trans(cholKuu)),Kuf));
  arma::colvec Kfu_invX_Kuf_invDelta_y = arma::trans(Kuf) * invX_Kuf_invDelta_y;
  arma::mat diagKfuInvXKuf = solve(arma::trimatl(arma::trans(cholX)),Kuf);
  diagKfuInvXKuf = arma::trans(sum(diagKfuInvXKuf % diagKfuInvXKuf));
  arma::mat invKuu = arma::inv(arma::trimatu(cholKuu));
  invKuu = arma::trimatu(invKuu) * arma::trimatl(arma::trans(invKuu));
  arma::mat invX = arma::inv(arma::trimatu(cholX));
  invX = arma::trimatu(invX) * arma::trimatl(arma::trans(invX));
  // variables that change according to current derivative
  double gradMarginalLikelihood;

  arma::colvec deriveLogDetKuu(ind.n_rows);
  arma::colvec deriveLogDetXMat(ind.n_rows);
  arma::colvec deriveLogDetDeltaMat(ind.n_rows);
  arma::colvec deriveMahaMat(ind.n_rows);
  double deriveLogDet(ind.n_rows);

  // memory intense
  arma::mat deriveDeltaMat; //(ind.n_rows,x.n_rows);
  arma::mat deriveInvDeltaMat; //(ind.n_rows, x.n_rows);
  arma::mat deriveX1Mat; //(ind.n_rows,ind.n_rows);

  deriveLogDetKuu = arma::accu(invKuu % deriveKuu);
  deriveDeltaMat = deriveKff - 2*sum(deriveKuf%invKuuKuf,0) + sum(invKuuKuf % (deriveKuu * invKuuKuf),0);
  deriveInvDeltaMat = deriveDeltaMat % arma::trans(-invDelta % invDelta);

  deriveX1Mat = deriveKuu + deriveKuf * (arma::trans(Kuf.each_row() % (2*invDelta.t())));
  deriveLogDetXMat = accu(invX % deriveX1Mat);
  deriveLogDetXMat += sum(deriveInvDeltaMat % arma::trans(diagKfuInvXKuf),1);
  deriveLogDetDeltaMat = sum(deriveDeltaMat % arma::trans(invDelta),1);

  deriveLogDet = arma::as_scalar(-deriveLogDetKuu + deriveLogDetDeltaMat + deriveLogDetXMat);

  deriveMahaMat = sum(deriveInvDeltaMat.each_row() % arma::trans(y % y),1);
  deriveMahaMat += (deriveInvDeltaMat.each_row() % arma::trans(y)) * arma::trans(Kuf) * (-2)*invX_Kuf_invDelta_y;
  deriveMahaMat += (-2)*dot(deriveKuf * (y % invDelta), invX_Kuf_invDelta_y);
  deriveMahaMat += dot(deriveX1Mat * invX_Kuf_invDelta_y, invX_Kuf_invDelta_y);
  deriveMahaMat += (deriveInvDeltaMat.each_row() % arma::trans(Kfu_invX_Kuf_invDelta_y)) * Kfu_invX_Kuf_invDelta_y;

  gradMarginalLikelihood = -0.5 * arma::as_scalar(deriveLogDet + deriveMahaMat);

  return gradMarginalLikelihood;
}





// [[Rcpp::export]]
arma::mat deriveSparseGPMarlikeIndSel(Rcpp::List & model, Rcpp::LogicalVector selection) {

  int indCnt = 0;
  for (int i=0; i<selection.size(); i++)
    if (selection[i]==TRUE) indCnt++;

  arma::uvec indIdx(indCnt,arma::fill::zeros);
  indCnt = 0;
  for (int i=0; i<selection.size(); i++)
    if (selection[i]==TRUE) indIdx(indCnt++) = i;

  Rcpp::List hyperpars = model["hyperpars"];
  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::Function deriveCovfun = gpFuns["deriveCovfunInd"];
  Rcpp::List hyperparsGP = hyperpars["hyperpars"];
  arma::mat ind = hyperpars["ind"];

  arma::mat x = model["x"];
  arma::colvec y = model["y"];
  arma::colvec delta = model["delta"];

  arma::mat cholX = model["cholX"];
  arma::mat cholKuu = model["cholKuu"];
  arma::colvec invX_Kuf_invDelta_y = model["invX_Kuf_invDelta_y"];

  arma::mat Kuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,ind,hyperparsGP,bool(FALSE))));
  arma::vec diagKff = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::NumericVector>(covfun(x,x,hyperparsGP,bool(TRUE))));
  arma::colvec invDelta = 1/delta;

  // prepare auxiliary quantities to maintain O(d*m^2*n) for complete gradient
  // at expense of additional memory requirement
  arma::mat Kuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,x,hyperparsGP,bool(FALSE))));
  arma::mat invKuuKuf = solve(arma::trimatu(cholKuu),solve(arma::trimatl(arma::trans(cholKuu)),Kuf));

  arma::colvec Kfu_invX_Kuf_invDelta_y = arma::trans(Kuf) * invX_Kuf_invDelta_y;
  arma::mat diagKfuInvXKuf = solve(arma::trimatl(arma::trans(cholX)),Kuf);
  diagKfuInvXKuf = arma::trans(sum(diagKfuInvXKuf % diagKfuInvXKuf));
  arma::mat invKuu = arma::inv(arma::trimatu(cholKuu));
  invKuu = arma::trimatu(invKuu) * arma::trimatl(arma::trans(invKuu));
  arma::mat invX = arma::inv(arma::trimatu(cholX));
  invX = arma::trimatu(invX) * arma::trimatl(arma::trans(invX));

  // variables that change according to current derivative
  arma::mat gradMarginalLikelihood(indCnt,ind.n_cols,arma::fill::zeros);

  arma::colvec deriveLogDetKuu(ind.n_rows);
  arma::colvec deriveLogDetXMat(ind.n_rows);
  arma::colvec deriveLogDetDeltaMat(ind.n_rows);
  arma::colvec deriveLogDetMat(ind.n_rows);
  arma::colvec deriveMahaMat(ind.n_rows);

  // memory intense
  arma::mat deriveDeltaMat; //(ind.n_rows,x.n_rows);
  arma::mat deriveInvDeltaMat; //(ind.n_rows, x.n_rows);
  arma::mat deriveX1Mat; //(ind.n_rows,ind.n_rows);

  arma::mat redInd = ind.rows(indIdx);  //TODO

  for (int j=0; j<ind.n_cols; j++) {

    arma::mat deriveKuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(redInd,x,hyperparsGP,bool(FALSE),j)));
    arma::mat deriveKuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(redInd,ind,hyperparsGP,bool(FALSE),j)));

    deriveLogDetKuu = 2*sum(invKuu.rows(indIdx) % deriveKuu,1);
    deriveDeltaMat = 2*(-deriveKuf + (deriveKuu * invKuuKuf)) % invKuuKuf.rows(indIdx);
    deriveInvDeltaMat = deriveDeltaMat.each_row() % arma::trans(-invDelta % invDelta);

    deriveX1Mat = 2*(deriveKuu + (deriveKuf.each_row() % arma::trans(invDelta)) * arma::trans(Kuf));

    deriveLogDetXMat = sum(invX.rows(indIdx) % deriveX1Mat,1);
    deriveLogDetXMat += sum(deriveInvDeltaMat.each_row() % arma::trans(diagKfuInvXKuf),1);
    deriveLogDetDeltaMat = sum(deriveDeltaMat.each_row() % arma::trans(invDelta),1);

    deriveLogDetMat = -deriveLogDetKuu + deriveLogDetDeltaMat + deriveLogDetXMat;

    deriveMahaMat = sum(deriveInvDeltaMat.each_row() % arma::trans(y % y),1);
    deriveMahaMat += (deriveInvDeltaMat.each_row() % arma::trans(y)) * arma::trans(Kuf) * (-2)*invX_Kuf_invDelta_y;

    deriveMahaMat += (deriveKuf * (y % invDelta)) % ((-2)*invX_Kuf_invDelta_y.rows(indIdx));
    deriveMahaMat += (deriveX1Mat * invX_Kuf_invDelta_y) % invX_Kuf_invDelta_y.rows(indIdx);
    deriveMahaMat += (deriveInvDeltaMat.each_row() % arma::trans(Kfu_invX_Kuf_invDelta_y)) * Kfu_invX_Kuf_invDelta_y;

    gradMarginalLikelihood.col(j) = -0.5 * (deriveLogDetMat + deriveMahaMat);
  }
  return gradMarginalLikelihood;
}




