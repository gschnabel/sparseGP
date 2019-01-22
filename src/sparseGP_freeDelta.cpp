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
Rcpp::List fitSparseGPfD(arma::mat & x, arma::colvec & y, arma::colvec & err, Rcpp::List & hyperpars) {

  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::List deltaFuns = hyperpars["deltaFuns"];

  Rcpp::List hyperparsGP = hyperpars["hyperpars"];
  Rcpp::List deltapars = hyperpars["deltapars"];

  arma::mat ind = hyperpars["ind"];

  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::Function compLogDetDelta = deltaFuns["logDet"];
  Rcpp::Function leftMultInvDelta = deltaFuns["leftMultInv"];
  Rcpp::Function bothMultInvDelta = deltaFuns["bothMultInv"];

  arma::mat Kuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,x,hyperparsGP,bool(FALSE))));
  arma::mat Kuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,ind,hyperparsGP,bool(FALSE))));
  arma::vec diagKff = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::NumericVector>(covfun(x,x,hyperparsGP,bool(TRUE))));

  // prepare invKuu
  arma::mat cholKuu = arma::chol(Kuu);

  // prepare delta
  arma::mat z = solve(trimatl(arma::trans(cholKuu)),Kuf);
  arma::colvec deltaPart = diagKff - arma::trans(sum(z%z)) + err%err;

  // prepare X
  //debug
  //arma::mat KfuInvDeltaKuf = Rcpp::as<arma::mat>(bothMultInvDelta(Kuf.t(), deltaPart, deltapars));

  arma::mat cholX = chol(Kuu + Rcpp::as<arma::mat>(bothMultInvDelta(Kuf.t(), deltaPart, deltapars)));
  arma::colvec invDelta_y = Rcpp::as<arma::colvec>(leftMultInvDelta(y, deltaPart, deltapars));
  //arma::colvec invDelta_y = invDelta % y;
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
  double logDetDelta = Rcpp::as<double>(compLogDetDelta(deltaPart, deltapars));
  //double logDetDelta = arma::sum(log(deltaPart));
  double logDetX = 2*arma::sum(log(cholX.diag()));
  double logDet =  logDetDelta + logDetX - logDetKuu;

  // all together
  double marlike = (-0.5) * (maha + logDet + x.n_rows*log2pi);

  return Rcpp::List::create(Rcpp::Named("invX_Kuf_invDelta_y") = invX_Kuf_invDelta_y,
                            Rcpp::Named("cholKuu") = cholKuu,
                            Rcpp::Named("cholX") = cholX,
                            Rcpp::Named("delta") = deltaPart,
                            Rcpp::Named("hyperpars") = hyperpars,
                            Rcpp::Named("x") = x,
                            Rcpp::Named("y") = y,
                            Rcpp::Named("err") = err,
                            Rcpp::Named("logDet") = logDet,
                            Rcpp::Named("maha") = maha,
                            Rcpp::Named("marlike") = marlike);
                            //debug
                            // Rcpp::Named("KfuInvDeltaKuf") = KfuInvDeltaKuf,
                            // Rcpp::Named("yInvDeltay") = yInvDeltay,
                            // Rcpp::Named("partMaha") = partMaha,
                            // Rcpp::Named("invDelta_y") = invDelta_y);
}


// [[Rcpp::export]]
Rcpp::DataFrame predictSparseGPfD(Rcpp::List & model, arma::mat & px, bool reterr = true) {

  Rcpp::Function covfun = model["covfun"];
  Rcpp::List hyperpars = model["hyperpars"];
  arma::mat ind = model["indpoints"];

  arma::colvec invX_Kuf_invDelta_y = model["invX_Kuf_invDelta_y"];

  // prediction
  arma::mat Kup = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,px,hyperpars,bool(FALSE))));
  arma::vec py = arma::trans(Kup) * invX_Kuf_invDelta_y;

  arma::colvec predVar(px.n_rows,arma::fill::zeros);
  if (reterr) {
    arma::mat cholKuu = model["cholKuu"];
    arma::mat cholX = model["cholX"];

    arma::vec diagKpp = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::NumericVector>(covfun(px,px,hyperpars,bool(TRUE))));
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
arma::mat deriveSparseGPMarlikeIndfD(Rcpp::List & model) {

  Rcpp::List hyperpars = model["hyperpars"];
  Rcpp::List deltaFuns = hyperpars["deltaFuns"];
  Rcpp::List deltapars = hyperpars["deltapars"];
  Rcpp::List hyperparsGP = hyperpars["hyperpars"];

  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::Function deriveCovfun = gpFuns["deriveCovfunInd"];
  Rcpp::Function compLogDetDelta = deltaFuns["logDet"];
  Rcpp::Function diagInvDelta = deltaFuns["diagInv"];
  Rcpp::Function leftMultInvDelta = deltaFuns["leftMultInv"];
  Rcpp::Function rightMultInvDelta = deltaFuns["rightMultInv"];
  Rcpp::Function bothMultInvDelta = deltaFuns["bothMultInv"];

  arma::colvec deltaPart = model["delta"];
  arma::mat ind = hyperpars["ind"];
  arma::mat x = model["x"];
  arma::colvec y = model["y"];

  arma::mat cholX = model["cholX"];
  arma::mat cholKuu = model["cholKuu"];
  arma::colvec invX_Kuf_invDelta_y = model["invX_Kuf_invDelta_y"];

  // prepare auxiliary quantities to maintain O(d*m^2*n) for complete gradient
  // at expense of additional memory requirement
  arma::mat Kuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,x,hyperparsGP,bool(FALSE))));
  arma::mat invKuuKuf = solve(arma::trimatu(cholKuu),solve(arma::trimatl(arma::trans(cholKuu)),Kuf));

  arma::colvec Kfu_invX_Kuf_invDelta_y = arma::trans(Kuf) * invX_Kuf_invDelta_y;
  arma::mat invKuu = arma::inv(arma::trimatu(cholKuu));
  invKuu = arma::trimatu(invKuu) * arma::trimatl(arma::trans(invKuu));
  arma::mat invX = arma::inv(arma::trimatu(cholX));
  invX = arma::trimatu(invX) * arma::trimatl(arma::trans(invX));
  arma::mat cholInvX_Kuf_invDelta = solve(arma::trimatl(arma::trans(cholX)),
                                      Rcpp::as<arma::mat>(rightMultInvDelta(Kuf, deltaPart, deltapars)));
  arma::rowvec diag_invDelta_Kfu_invX_Kuf_invDelta = sum(cholInvX_Kuf_invDelta % cholInvX_Kuf_invDelta,0);
  arma::colvec invDelta_y = Rcpp::as<arma::mat>(leftMultInvDelta(y, deltaPart, deltapars));
  arma::colvec invDelta_Kfu_invX_Kuf_invDelta_y = Rcpp::as<arma::mat>(leftMultInvDelta(arma::trans(Kuf)*invX_Kuf_invDelta_y, deltaPart, deltapars));


  // variables that change according to current derivative
  arma::mat gradMarginalLikelihood(ind.n_rows,ind.n_cols,arma::fill::zeros);

  arma::colvec deriveLogDetKuu(ind.n_rows);
  arma::colvec deriveLogDetXMat(ind.n_rows);
  arma::colvec deriveLogDetDeltaMat(ind.n_rows);
  arma::colvec deriveLogDetMat(ind.n_rows);
  arma::colvec deriveMahaMat(ind.n_rows);

  // memory intense
  arma::mat deriveDeltaMat; //(ind.n_rows,x.n_rows);
  arma::mat deriveX1Mat; //(ind.n_rows,ind.n_rows);

  for (int j=0; j<ind.n_cols; j++) {

    arma::mat deriveKuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(ind,x,hyperparsGP,bool(FALSE),j)));
    arma::mat deriveKuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(ind,ind,hyperparsGP,bool(FALSE),j)));

    deriveLogDetKuu = 2*sum(invKuu % deriveKuu,1);
    deriveDeltaMat = 2*(-deriveKuf + (deriveKuu * invKuuKuf)) % invKuuKuf;
    deriveX1Mat = 2*(deriveKuu + (Rcpp::as<arma::mat>(rightMultInvDelta(deriveKuf, deltaPart, deltapars))) * arma::trans(Kuf));

    deriveLogDetXMat = sum(invX % deriveX1Mat,1);
    deriveLogDetXMat -= sum(deriveDeltaMat.each_row() % diag_invDelta_Kfu_invX_Kuf_invDelta, 1);
    deriveLogDetDeltaMat = sum(deriveDeltaMat.each_row() % Rcpp::as<arma::rowvec>(diagInvDelta(deltaPart, deltapars)),1);
    deriveLogDetMat = -deriveLogDetKuu + deriveLogDetDeltaMat + deriveLogDetXMat;

    deriveMahaMat = sum(deriveDeltaMat.each_row() % arma::trans(invDelta_y % (-invDelta_y)),1);
    deriveMahaMat += (deriveDeltaMat.each_row() % arma::trans(2*invDelta_Kfu_invX_Kuf_invDelta_y)) * invDelta_y;
    deriveMahaMat += (deriveKuf * invDelta_y) % ((-2)*invX_Kuf_invDelta_y);
    deriveMahaMat += (deriveX1Mat * invX_Kuf_invDelta_y) % invX_Kuf_invDelta_y;
    deriveMahaMat += (deriveDeltaMat.each_row() % arma::trans(-invDelta_Kfu_invX_Kuf_invDelta_y)) * invDelta_Kfu_invX_Kuf_invDelta_y;

    gradMarginalLikelihood.col(j) = -0.5 * (deriveLogDetMat + deriveMahaMat);
  }
  return gradMarginalLikelihood;
}




// [[Rcpp::export]]
double deriveSparseGPMarlikeHypfD(Rcpp::List & model, Rcpp::List & selection) {

  Rcpp::List hyperpars = model["hyperpars"];
  Rcpp::List deltapars = hyperpars["deltapars"];

  Rcpp::List deltaFuns = hyperpars["deltaFuns"];
  Rcpp::Function compLogDetDelta = deltaFuns["logDet"];
  Rcpp::Function diagInvDelta = deltaFuns["diagInv"];
  Rcpp::Function leftMultInvDelta = deltaFuns["leftMultInv"];
  Rcpp::Function rightMultInvDelta = deltaFuns["rightMultInv"];
  Rcpp::Function bothMultInvDelta = deltaFuns["bothMultInv"];

  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::Function deriveCovfun = gpFuns["deriveCovfunHyp"];

  Rcpp::List selectionGP = selection["hyperpars"];
  Rcpp::List hyperparsGP = hyperpars["hyperpars"];
  arma::mat ind = hyperpars["ind"];

  arma::mat x = model["x"];
  arma::colvec y = model["y"];
  arma::colvec deltaPart = model["delta"];

  arma::mat cholX = model["cholX"];
  arma::mat cholKuu = model["cholKuu"];
  arma::colvec invX_Kuf_invDelta_y = model["invX_Kuf_invDelta_y"];

  arma::mat Kuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,ind,hyperparsGP,bool(FALSE))));
  arma::vec diagKff = Rcpp::as<arma::vec>(Rcpp::as<Rcpp::NumericVector>(covfun(x,x,hyperparsGP,bool(TRUE))));

  arma::mat deriveKuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(ind,ind,hyperparsGP,selectionGP,false,true)));
  arma::mat deriveKuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(ind,x,hyperparsGP,selectionGP,false,true)));
  arma::mat deriveKff = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(x,x,hyperparsGP,selectionGP,true,true)));

  // prepare auxiliary quantities to maintain O(d*m^2*n) for complete gradient
  // at expense of additional memory requirement
  arma::mat Kuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,x,hyperparsGP,bool(FALSE))));
  arma::mat invKuuKuf = solve(arma::trimatu(cholKuu),solve(arma::trimatl(arma::trans(cholKuu)),Kuf));
  arma::mat invKuu = arma::inv(arma::trimatu(cholKuu));
  invKuu = arma::trimatu(invKuu) * arma::trimatl(arma::trans(invKuu));
  arma::mat invX = arma::inv(arma::trimatu(cholX));
  invX = arma::trimatu(invX) * arma::trimatl(arma::trans(invX));

  arma::mat cholInvX_Kuf_invDelta = solve(arma::trimatl(arma::trans(cholX)),
                                          Rcpp::as<arma::mat>(rightMultInvDelta(Kuf, deltaPart, deltapars)));
  arma::rowvec diag_invDelta_Kfu_invX_Kuf_invDelta = sum(cholInvX_Kuf_invDelta % cholInvX_Kuf_invDelta,0);
  cholInvX_Kuf_invDelta.reset();

  arma::colvec invDelta_y = Rcpp::as<arma::mat>(leftMultInvDelta(y, deltaPart, deltapars));
  arma::colvec invDelta_Kfu_invX_Kuf_invDelta_y = Rcpp::as<arma::mat>(leftMultInvDelta(arma::trans(Kuf)*invX_Kuf_invDelta_y, deltaPart, deltapars));


  // variables that change according to current derivative
  double gradMarginalLikelihood;

  arma::colvec deriveLogDetKuu(ind.n_rows);
  arma::colvec deriveLogDetXMat(ind.n_rows);
  arma::colvec deriveLogDetDeltaMat(ind.n_rows);
  arma::colvec deriveMahaMat(ind.n_rows);
  double deriveLogDet(ind.n_rows);

  // memory intense
  arma::mat deriveDeltaMat; //(ind.n_rows,x.n_rows);
  arma::mat deriveX1Mat; //(ind.n_rows,ind.n_rows);

  deriveLogDetKuu = arma::accu(invKuu % deriveKuu);
  deriveDeltaMat = deriveKff - 2*sum(deriveKuf%invKuuKuf,0) + sum(invKuuKuf % (deriveKuu * invKuuKuf),0);

  deriveX1Mat = deriveKuu + 2*((Rcpp::as<arma::mat>(rightMultInvDelta(deriveKuf, deltaPart, deltapars)) * arma::trans(Kuf)));
  deriveLogDetXMat = accu(invX % deriveX1Mat);
  deriveLogDetXMat -= sum(deriveDeltaMat.each_row() % diag_invDelta_Kfu_invX_Kuf_invDelta, 1);
  deriveLogDetDeltaMat = sum(deriveDeltaMat.each_row() % Rcpp::as<arma::rowvec>(diagInvDelta(deltaPart, deltapars)),1);
  deriveLogDet = arma::as_scalar(-deriveLogDetKuu + deriveLogDetDeltaMat + deriveLogDetXMat);

  deriveMahaMat = sum(deriveDeltaMat.each_row() % arma::trans(invDelta_y % (-invDelta_y)),1);
  deriveMahaMat += (deriveDeltaMat.each_row() % arma::trans(2*invDelta_Kfu_invX_Kuf_invDelta_y)) * invDelta_y;
  deriveMahaMat += (-2) * dot(deriveKuf * invDelta_y, invX_Kuf_invDelta_y);
  deriveMahaMat += dot(deriveX1Mat * invX_Kuf_invDelta_y, invX_Kuf_invDelta_y);
  deriveMahaMat += (deriveDeltaMat.each_row() % arma::trans(-invDelta_Kfu_invX_Kuf_invDelta_y)) * invDelta_Kfu_invX_Kuf_invDelta_y;

  gradMarginalLikelihood = -0.5 * arma::as_scalar(deriveLogDet + deriveMahaMat);

  return gradMarginalLikelihood;
}


// [[Rcpp::export]]
arma::colvec deriveSparseGPMarlikeDeltafD(Rcpp::List & model, Rcpp::List & selection) {

  Rcpp::List hyperpars = model["hyperpars"];
  Rcpp::List deltapars = hyperpars["deltapars"];

  Rcpp::List deltaFuns = hyperpars["deltaFuns"];
  Rcpp::Function trInvDerive = deltaFuns["trInvDerive"];
  Rcpp::Function trMxtUinvDeriveInvU = deltaFuns["trMxtUinvDeriveInvU"];
  Rcpp::Function multInvDeriveInvByVecs = deltaFuns["multInvDeriveInvByVecs"];

  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::List hyperparsGP = hyperpars["hyperpars"];
  Rcpp::List selectionDelta = selection["deltapars"];

  arma::mat ind = hyperpars["ind"];
  arma::mat x = model["x"];
  arma::colvec y = model["y"];
  arma::colvec deltaPart = model["delta"];

  arma::mat cholX = model["cholX"];
  arma::colvec invX_Kuf_invDelta_y = model["invX_Kuf_invDelta_y"];

  // prepare auxiliary quantities to maintain O(d*m^2*n) for complete gradient
  // at expense of additional memory requirement
  arma::mat Kuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,x,hyperparsGP,bool(FALSE))));
  arma::mat invX = arma::inv(arma::trimatu(cholX));
  invX = arma::trimatu(invX) * arma::trimatl(arma::trans(invX));
  arma::colvec Kfu_invX_Kuf_invDelta_y = arma::trans(Kuf)*invX_Kuf_invDelta_y;

  // variables that change according to current derivative
  arma::colvec gradMarginalLikelihood;

  arma::colvec deriveLogDetXMat;
  arma::colvec deriveLogDetDeltaMat;
  arma::colvec deriveMahaMat;
  arma::colvec deriveLogDet;

  // memory intense
  arma::mat tInvCholXKfu = arma::solve(arma::trimatl(arma::trans(cholX)), Kuf);
  deriveLogDetXMat = -(Rcpp::as<arma::colvec>(trMxtUinvDeriveInvU(tInvCholXKfu, deltaPart, deltapars, selectionDelta)));
  tInvCholXKfu.reset();
  deriveLogDetDeltaMat = (Rcpp::as<arma::colvec>(trInvDerive(deltaPart, deltapars, selectionDelta)));
  deriveLogDet = (deriveLogDetDeltaMat + deriveLogDetXMat);

  deriveMahaMat = (-(Rcpp::as<arma::colvec>(multInvDeriveInvByVecs(y, y, deltaPart, deltapars, selectionDelta))));
  deriveMahaMat += 2*(Rcpp::as<arma::colvec>(multInvDeriveInvByVecs(y, Kfu_invX_Kuf_invDelta_y, deltaPart, deltapars, selectionDelta)));
  deriveMahaMat -= (Rcpp::as<arma::colvec>(multInvDeriveInvByVecs(Kfu_invX_Kuf_invDelta_y, Kfu_invX_Kuf_invDelta_y, deltaPart, deltapars, selectionDelta)));

  gradMarginalLikelihood = -0.5 * (deriveLogDet + deriveMahaMat);
  return gradMarginalLikelihood;
}






// [[Rcpp::export]]
arma::mat deriveSparseGPMarlikeIndSelfD(Rcpp::List & model, Rcpp::LogicalVector selection) {

  int indCnt = 0;
  for (int i=0; i<selection.size(); i++)
    if (selection[i]==TRUE) indCnt++;

  arma::uvec indIdx(indCnt,arma::fill::zeros);
  indCnt = 0;
  for (int i=0; i<selection.size(); i++)
    if (selection[i]==TRUE) indIdx(indCnt++) = i;

  Rcpp::List hyperpars = model["hyperpars"];
  Rcpp::List deltaFuns = hyperpars["deltaFuns"];
  Rcpp::List deltapars = hyperpars["deltapars"];
  Rcpp::List hyperparsGP = hyperpars["hyperpars"];

  Rcpp::List gpFuns = hyperpars["gpFuns"];
  Rcpp::Function covfun = gpFuns["covfun"];
  Rcpp::Function deriveCovfun = gpFuns["deriveCovfunInd"];
  Rcpp::Function compLogDetDelta = deltaFuns["logDet"];
  Rcpp::Function diagInvDelta = deltaFuns["diagInv"];
  Rcpp::Function leftMultInvDelta = deltaFuns["leftMultInv"];
  Rcpp::Function rightMultInvDelta = deltaFuns["rightMultInv"];
  Rcpp::Function bothMultInvDelta = deltaFuns["bothMultInv"];

  arma::colvec deltaPart = model["delta"];
  arma::mat ind = hyperpars["ind"];
  arma::mat x = model["x"];
  arma::colvec y = model["y"];

  arma::mat cholX = model["cholX"];
  arma::mat cholKuu = model["cholKuu"];
  arma::colvec invX_Kuf_invDelta_y = model["invX_Kuf_invDelta_y"];

  // prepare auxiliary quantities to maintain O(d*m^2*n) for complete gradient
  // at expense of additional memory requirement
  arma::mat Kuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(covfun(ind,x,hyperparsGP,bool(FALSE))));
  arma::mat invKuuKuf = solve(arma::trimatu(cholKuu),solve(arma::trimatl(arma::trans(cholKuu)),Kuf));

  arma::colvec Kfu_invX_Kuf_invDelta_y = arma::trans(Kuf) * invX_Kuf_invDelta_y;
  arma::mat invKuu = arma::inv(arma::trimatu(cholKuu));
  invKuu = arma::trimatu(invKuu) * arma::trimatl(arma::trans(invKuu));
  arma::mat invX = arma::inv(arma::trimatu(cholX));
  invX = arma::trimatu(invX) * arma::trimatl(arma::trans(invX));
  arma::mat cholInvX_Kuf_invDelta = solve(arma::trimatl(arma::trans(cholX)),
                                          Rcpp::as<arma::mat>(rightMultInvDelta(Kuf, deltaPart, deltapars)));
  arma::rowvec diag_invDelta_Kfu_invX_Kuf_invDelta = sum(cholInvX_Kuf_invDelta % cholInvX_Kuf_invDelta,0);
  arma::colvec invDelta_y = Rcpp::as<arma::mat>(leftMultInvDelta(y, deltaPart, deltapars));
  arma::colvec invDelta_Kfu_invX_Kuf_invDelta_y = Rcpp::as<arma::mat>(leftMultInvDelta(arma::trans(Kuf)*invX_Kuf_invDelta_y, deltaPart, deltapars));


  // variables that change according to current derivative
  arma::mat gradMarginalLikelihood(indCnt,ind.n_cols,arma::fill::zeros);

  arma::colvec deriveLogDetKuu(indCnt);
  arma::colvec deriveLogDetXMat(indCnt);
  arma::colvec deriveLogDetDeltaMat(indCnt);
  arma::colvec deriveLogDetMat(indCnt);
  arma::colvec deriveMahaMat(indCnt);

  // memory intense
  arma::mat deriveDeltaMat;
  arma::mat deriveX1Mat;
  arma::mat redInd = ind.rows(indIdx);

  for (int j=0; j<ind.n_cols; j++) {

    arma::mat deriveKuf = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(redInd,x,hyperparsGP,bool(FALSE),j)));
    arma::mat deriveKuu = Rcpp::as<arma::mat>(Rcpp::as<Rcpp::NumericMatrix>(deriveCovfun(redInd,ind,hyperparsGP,bool(FALSE),j)));

    deriveLogDetKuu = 2*sum(invKuu.rows(indIdx) % deriveKuu,1);
    deriveDeltaMat = 2*(-deriveKuf + (deriveKuu * invKuuKuf)) % invKuuKuf.rows(indIdx);
    deriveX1Mat = 2*(deriveKuu + (Rcpp::as<arma::mat>(rightMultInvDelta(deriveKuf, deltaPart, deltapars))) * arma::trans(Kuf));

    deriveLogDetXMat = sum(invX.rows(indIdx) % deriveX1Mat,1);
    deriveLogDetXMat -= sum(deriveDeltaMat.each_row() % diag_invDelta_Kfu_invX_Kuf_invDelta, 1);
    deriveLogDetDeltaMat = sum(deriveDeltaMat.each_row() % Rcpp::as<arma::rowvec>(diagInvDelta(deltaPart, deltapars)),1);
    deriveLogDetMat = -deriveLogDetKuu + deriveLogDetDeltaMat + deriveLogDetXMat;

    deriveMahaMat = sum(deriveDeltaMat.each_row() % arma::trans(invDelta_y % (-invDelta_y)),1);
    deriveMahaMat += (deriveDeltaMat.each_row() % arma::trans(2*invDelta_Kfu_invX_Kuf_invDelta_y)) * invDelta_y;
    deriveMahaMat += (deriveKuf * invDelta_y) % ((-2)*invX_Kuf_invDelta_y.rows(indIdx));
    deriveMahaMat += (deriveX1Mat * invX_Kuf_invDelta_y) % invX_Kuf_invDelta_y.rows(indIdx);
    deriveMahaMat += (deriveDeltaMat.each_row() % arma::trans(-invDelta_Kfu_invX_Kuf_invDelta_y)) * invDelta_Kfu_invX_Kuf_invDelta_y;

    gradMarginalLikelihood.col(j) = -0.5 * (deriveLogDetMat + deriveMahaMat);
  }
  return gradMarginalLikelihood;
}




