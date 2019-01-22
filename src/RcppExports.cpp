// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// getMatrixCP
arma::mat getMatrixCP(arma::mat& x, arma::mat& y, Rcpp::List& hyperpars, bool& elementMode);
RcppExport SEXP _sparseGP_getMatrixCP(SEXP xSEXP, SEXP ySEXP, SEXP hyperparsSEXP, SEXP elementModeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type hyperpars(hyperparsSEXP);
    Rcpp::traits::input_parameter< bool& >::type elementMode(elementModeSEXP);
    rcpp_result_gen = Rcpp::wrap(getMatrixCP(x, y, hyperpars, elementMode));
    return rcpp_result_gen;
END_RCPP
}
// deriveMatrixIndCP
arma::mat deriveMatrixIndCP(arma::mat& x, arma::mat& y, Rcpp::List& hyperpars, bool& elementMode, int& dimIndex);
RcppExport SEXP _sparseGP_deriveMatrixIndCP(SEXP xSEXP, SEXP ySEXP, SEXP hyperparsSEXP, SEXP elementModeSEXP, SEXP dimIndexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type hyperpars(hyperparsSEXP);
    Rcpp::traits::input_parameter< bool& >::type elementMode(elementModeSEXP);
    Rcpp::traits::input_parameter< int& >::type dimIndex(dimIndexSEXP);
    rcpp_result_gen = Rcpp::wrap(deriveMatrixIndCP(x, y, hyperpars, elementMode, dimIndex));
    return rcpp_result_gen;
END_RCPP
}
// deriveMatrixHypCP
SEXP deriveMatrixHypCP(arma::mat& x, arma::mat& y, Rcpp::List& hyperpars, Rcpp::List& selection, bool elementMode, bool singleDerive);
RcppExport SEXP _sparseGP_deriveMatrixHypCP(SEXP xSEXP, SEXP ySEXP, SEXP hyperparsSEXP, SEXP selectionSEXP, SEXP elementModeSEXP, SEXP singleDeriveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type hyperpars(hyperparsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type selection(selectionSEXP);
    Rcpp::traits::input_parameter< bool >::type elementMode(elementModeSEXP);
    Rcpp::traits::input_parameter< bool >::type singleDerive(singleDeriveSEXP);
    rcpp_result_gen = Rcpp::wrap(deriveMatrixHypCP(x, y, hyperpars, selection, elementMode, singleDerive));
    return rcpp_result_gen;
END_RCPP
}
// getMatrixSmap
arma::mat getMatrixSmap(Rcpp::DataFrame& x, Rcpp::DataFrame& y, Rcpp::List& hyperpars, bool& elementMode);
RcppExport SEXP _sparseGP_getMatrixSmap(SEXP xSEXP, SEXP ySEXP, SEXP hyperparsSEXP, SEXP elementModeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type hyperpars(hyperparsSEXP);
    Rcpp::traits::input_parameter< bool& >::type elementMode(elementModeSEXP);
    rcpp_result_gen = Rcpp::wrap(getMatrixSmap(x, y, hyperpars, elementMode));
    return rcpp_result_gen;
END_RCPP
}
// getMatrix
arma::mat getMatrix(arma::mat& x, arma::mat& y, Rcpp::List& hyperpars, bool& elementMode);
RcppExport SEXP _sparseGP_getMatrix(SEXP xSEXP, SEXP ySEXP, SEXP hyperparsSEXP, SEXP elementModeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type hyperpars(hyperparsSEXP);
    Rcpp::traits::input_parameter< bool& >::type elementMode(elementModeSEXP);
    rcpp_result_gen = Rcpp::wrap(getMatrix(x, y, hyperpars, elementMode));
    return rcpp_result_gen;
END_RCPP
}
// deriveMatrixInd
arma::mat deriveMatrixInd(arma::mat& x, arma::mat& y, Rcpp::List& hyperpars, bool& elementMode, int& dimIndex);
RcppExport SEXP _sparseGP_deriveMatrixInd(SEXP xSEXP, SEXP ySEXP, SEXP hyperparsSEXP, SEXP elementModeSEXP, SEXP dimIndexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type hyperpars(hyperparsSEXP);
    Rcpp::traits::input_parameter< bool& >::type elementMode(elementModeSEXP);
    Rcpp::traits::input_parameter< int& >::type dimIndex(dimIndexSEXP);
    rcpp_result_gen = Rcpp::wrap(deriveMatrixInd(x, y, hyperpars, elementMode, dimIndex));
    return rcpp_result_gen;
END_RCPP
}
// deriveMatrixHyp
SEXP deriveMatrixHyp(arma::mat& x, arma::mat& y, Rcpp::List& hyperpars, Rcpp::List& selection, bool elementMode, bool singleDerive);
RcppExport SEXP _sparseGP_deriveMatrixHyp(SEXP xSEXP, SEXP ySEXP, SEXP hyperparsSEXP, SEXP selectionSEXP, SEXP elementModeSEXP, SEXP singleDeriveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type hyperpars(hyperparsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type selection(selectionSEXP);
    Rcpp::traits::input_parameter< bool >::type elementMode(elementModeSEXP);
    Rcpp::traits::input_parameter< bool >::type singleDerive(singleDeriveSEXP);
    rcpp_result_gen = Rcpp::wrap(deriveMatrixHyp(x, y, hyperpars, selection, elementMode, singleDerive));
    return rcpp_result_gen;
END_RCPP
}
// fitSparseGPfD
Rcpp::List fitSparseGPfD(arma::mat& x, arma::colvec& y, arma::colvec& err, Rcpp::List& hyperpars);
RcppExport SEXP _sparseGP_fitSparseGPfD(SEXP xSEXP, SEXP ySEXP, SEXP errSEXP, SEXP hyperparsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type err(errSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type hyperpars(hyperparsSEXP);
    rcpp_result_gen = Rcpp::wrap(fitSparseGPfD(x, y, err, hyperpars));
    return rcpp_result_gen;
END_RCPP
}
// predictSparseGPfD
Rcpp::DataFrame predictSparseGPfD(Rcpp::List& model, arma::mat& px, bool reterr);
RcppExport SEXP _sparseGP_predictSparseGPfD(SEXP modelSEXP, SEXP pxSEXP, SEXP reterrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< bool >::type reterr(reterrSEXP);
    rcpp_result_gen = Rcpp::wrap(predictSparseGPfD(model, px, reterr));
    return rcpp_result_gen;
END_RCPP
}
// deriveSparseGPMarlikeIndfD
arma::mat deriveSparseGPMarlikeIndfD(Rcpp::List& model);
RcppExport SEXP _sparseGP_deriveSparseGPMarlikeIndfD(SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(deriveSparseGPMarlikeIndfD(model));
    return rcpp_result_gen;
END_RCPP
}
// deriveSparseGPMarlikeHypfD
double deriveSparseGPMarlikeHypfD(Rcpp::List& model, Rcpp::List& selection);
RcppExport SEXP _sparseGP_deriveSparseGPMarlikeHypfD(SEXP modelSEXP, SEXP selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type selection(selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(deriveSparseGPMarlikeHypfD(model, selection));
    return rcpp_result_gen;
END_RCPP
}
// deriveSparseGPMarlikeDeltafD
arma::colvec deriveSparseGPMarlikeDeltafD(Rcpp::List& model, Rcpp::List& selection);
RcppExport SEXP _sparseGP_deriveSparseGPMarlikeDeltafD(SEXP modelSEXP, SEXP selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type selection(selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(deriveSparseGPMarlikeDeltafD(model, selection));
    return rcpp_result_gen;
END_RCPP
}
// deriveSparseGPMarlikeIndSelfD
arma::mat deriveSparseGPMarlikeIndSelfD(Rcpp::List& model, Rcpp::LogicalVector selection);
RcppExport SEXP _sparseGP_deriveSparseGPMarlikeIndSelfD(SEXP modelSEXP, SEXP selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type selection(selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(deriveSparseGPMarlikeIndSelfD(model, selection));
    return rcpp_result_gen;
END_RCPP
}
// fitSparseGP
Rcpp::List fitSparseGP(arma::mat& x, arma::colvec& y, arma::colvec& err, Rcpp::List& hyperpars);
RcppExport SEXP _sparseGP_fitSparseGP(SEXP xSEXP, SEXP ySEXP, SEXP errSEXP, SEXP hyperparsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type err(errSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type hyperpars(hyperparsSEXP);
    rcpp_result_gen = Rcpp::wrap(fitSparseGP(x, y, err, hyperpars));
    return rcpp_result_gen;
END_RCPP
}
// predictSparseGP
Rcpp::DataFrame predictSparseGP(Rcpp::List& model, arma::mat& px, bool reterr);
RcppExport SEXP _sparseGP_predictSparseGP(SEXP modelSEXP, SEXP pxSEXP, SEXP reterrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< bool >::type reterr(reterrSEXP);
    rcpp_result_gen = Rcpp::wrap(predictSparseGP(model, px, reterr));
    return rcpp_result_gen;
END_RCPP
}
// getPostCovFromSparseGP
arma::mat getPostCovFromSparseGP(arma::mat& px, arma::mat& py, Rcpp::List& model, bool elementMode, bool exact);
RcppExport SEXP _sparseGP_getPostCovFromSparseGP(SEXP pxSEXP, SEXP pySEXP, SEXP modelSEXP, SEXP elementModeSEXP, SEXP exactSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type py(pySEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< bool >::type elementMode(elementModeSEXP);
    Rcpp::traits::input_parameter< bool >::type exact(exactSEXP);
    rcpp_result_gen = Rcpp::wrap(getPostCovFromSparseGP(px, py, model, elementMode, exact));
    return rcpp_result_gen;
END_RCPP
}
// preparePostSamplingFromSparseGP
Rcpp::List preparePostSamplingFromSparseGP(Rcpp::List& model, arma::mat& px, bool exact);
RcppExport SEXP _sparseGP_preparePostSamplingFromSparseGP(SEXP modelSEXP, SEXP pxSEXP, SEXP exactSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type px(pxSEXP);
    Rcpp::traits::input_parameter< bool >::type exact(exactSEXP);
    rcpp_result_gen = Rcpp::wrap(preparePostSamplingFromSparseGP(model, px, exact));
    return rcpp_result_gen;
END_RCPP
}
// getSamplesFromSparseGP
arma::mat getSamplesFromSparseGP(Rcpp::List sampleSpec, int num, bool exact);
RcppExport SEXP _sparseGP_getSamplesFromSparseGP(SEXP sampleSpecSEXP, SEXP numSEXP, SEXP exactSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type sampleSpec(sampleSpecSEXP);
    Rcpp::traits::input_parameter< int >::type num(numSEXP);
    Rcpp::traits::input_parameter< bool >::type exact(exactSEXP);
    rcpp_result_gen = Rcpp::wrap(getSamplesFromSparseGP(sampleSpec, num, exact));
    return rcpp_result_gen;
END_RCPP
}
// deriveSparseGPMarlikeInd
arma::mat deriveSparseGPMarlikeInd(Rcpp::List& model);
RcppExport SEXP _sparseGP_deriveSparseGPMarlikeInd(SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(deriveSparseGPMarlikeInd(model));
    return rcpp_result_gen;
END_RCPP
}
// deriveSparseGPMarlikeHyp
double deriveSparseGPMarlikeHyp(Rcpp::List& model, Rcpp::List& selection);
RcppExport SEXP _sparseGP_deriveSparseGPMarlikeHyp(SEXP modelSEXP, SEXP selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type selection(selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(deriveSparseGPMarlikeHyp(model, selection));
    return rcpp_result_gen;
END_RCPP
}
// deriveSparseGPMarlikeIndSel
arma::mat deriveSparseGPMarlikeIndSel(Rcpp::List& model, Rcpp::LogicalVector selection);
RcppExport SEXP _sparseGP_deriveSparseGPMarlikeIndSel(SEXP modelSEXP, SEXP selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type selection(selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(deriveSparseGPMarlikeIndSel(model, selection));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sparseGP_getMatrixCP", (DL_FUNC) &_sparseGP_getMatrixCP, 4},
    {"_sparseGP_deriveMatrixIndCP", (DL_FUNC) &_sparseGP_deriveMatrixIndCP, 5},
    {"_sparseGP_deriveMatrixHypCP", (DL_FUNC) &_sparseGP_deriveMatrixHypCP, 6},
    {"_sparseGP_getMatrixSmap", (DL_FUNC) &_sparseGP_getMatrixSmap, 4},
    {"_sparseGP_getMatrix", (DL_FUNC) &_sparseGP_getMatrix, 4},
    {"_sparseGP_deriveMatrixInd", (DL_FUNC) &_sparseGP_deriveMatrixInd, 5},
    {"_sparseGP_deriveMatrixHyp", (DL_FUNC) &_sparseGP_deriveMatrixHyp, 6},
    {"_sparseGP_fitSparseGPfD", (DL_FUNC) &_sparseGP_fitSparseGPfD, 4},
    {"_sparseGP_predictSparseGPfD", (DL_FUNC) &_sparseGP_predictSparseGPfD, 3},
    {"_sparseGP_deriveSparseGPMarlikeIndfD", (DL_FUNC) &_sparseGP_deriveSparseGPMarlikeIndfD, 1},
    {"_sparseGP_deriveSparseGPMarlikeHypfD", (DL_FUNC) &_sparseGP_deriveSparseGPMarlikeHypfD, 2},
    {"_sparseGP_deriveSparseGPMarlikeDeltafD", (DL_FUNC) &_sparseGP_deriveSparseGPMarlikeDeltafD, 2},
    {"_sparseGP_deriveSparseGPMarlikeIndSelfD", (DL_FUNC) &_sparseGP_deriveSparseGPMarlikeIndSelfD, 2},
    {"_sparseGP_fitSparseGP", (DL_FUNC) &_sparseGP_fitSparseGP, 4},
    {"_sparseGP_predictSparseGP", (DL_FUNC) &_sparseGP_predictSparseGP, 3},
    {"_sparseGP_getPostCovFromSparseGP", (DL_FUNC) &_sparseGP_getPostCovFromSparseGP, 5},
    {"_sparseGP_preparePostSamplingFromSparseGP", (DL_FUNC) &_sparseGP_preparePostSamplingFromSparseGP, 3},
    {"_sparseGP_getSamplesFromSparseGP", (DL_FUNC) &_sparseGP_getSamplesFromSparseGP, 3},
    {"_sparseGP_deriveSparseGPMarlikeInd", (DL_FUNC) &_sparseGP_deriveSparseGPMarlikeInd, 1},
    {"_sparseGP_deriveSparseGPMarlikeHyp", (DL_FUNC) &_sparseGP_deriveSparseGPMarlikeHyp, 2},
    {"_sparseGP_deriveSparseGPMarlikeIndSel", (DL_FUNC) &_sparseGP_deriveSparseGPMarlikeIndSel, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_sparseGP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
