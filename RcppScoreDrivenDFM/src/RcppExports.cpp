#include "cppad/cppad.hpp"
#include "RcppEigen.h"

using namespace Eigen;
using namespace CppAD;

typedef AD<double> a_double;
typedef Matrix<a_double, Eigen::Dynamic, 1> a_vector;

template<typename T>
using Tvec =  Matrix<T, Eigen::Dynamic, 1>;

template<typename T>
using Tmat =  Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using Trowvec =  Matrix<T, 1, Eigen::Dynamic>;

// // //

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// filter_list
Rcpp::List filter_list(const Rcpp::List& model);
RcppExport SEXP _RcppScoreDrivenDFM_filter_list(SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(filter_list(model));
    return rcpp_result_gen;
END_RCPP
}
// loglik_rcpp
double loglik_rcpp(const VectorXd& pars, const MatrixXd& y, MatrixXd& T_loc, MatrixXd& T_scale, MatrixXd& T_shape, MatrixXd& Z_loc, MatrixXd& Z_scale, MatrixXd& Z_shape, const bool& stoch_loc, const bool& stoch_vol, const bool& stoch_shape, const MatrixXi& NAs, const VectorXd& guess_loc, const VectorXd& guess_scale);
RcppExport SEXP _RcppScoreDrivenDFM_loglik_rcpp(SEXP parsSEXP, SEXP ySEXP, SEXP T_locSEXP, SEXP T_scaleSEXP, SEXP T_shapeSEXP, SEXP Z_locSEXP, SEXP Z_scaleSEXP, SEXP Z_shapeSEXP, SEXP stoch_locSEXP, SEXP stoch_volSEXP, SEXP stoch_shapeSEXP, SEXP NAsSEXP, SEXP guess_locSEXP, SEXP guess_scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type T_loc(T_locSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type T_scale(T_scaleSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type T_shape(T_shapeSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type Z_loc(Z_locSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type Z_scale(Z_scaleSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type Z_shape(Z_shapeSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_loc(stoch_locSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_vol(stoch_volSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_shape(stoch_shapeSEXP);
    Rcpp::traits::input_parameter< const MatrixXi& >::type NAs(NAsSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type guess_loc(guess_locSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type guess_scale(guess_scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_rcpp(pars, y, T_loc, T_scale, T_shape, Z_loc, Z_scale, Z_shape, stoch_loc, stoch_vol, stoch_shape, NAs, guess_loc, guess_scale));
    return rcpp_result_gen;
END_RCPP
}
// loglik_deriv
VectorXd loglik_deriv(const VectorXd& pars, const MatrixXd& y, MatrixXd& T_loc, MatrixXd& T_scale, MatrixXd& T_shape, MatrixXd& Z_loc, MatrixXd& Z_scale, MatrixXd& Z_shape, const bool& stoch_loc, const bool& stoch_vol, const bool& stoch_shape, const MatrixXi& NAs, const VectorXd& guess_loc, const VectorXd& guess_scale);
RcppExport SEXP _RcppScoreDrivenDFM_loglik_deriv(SEXP parsSEXP, SEXP ySEXP, SEXP T_locSEXP, SEXP T_scaleSEXP, SEXP T_shapeSEXP, SEXP Z_locSEXP, SEXP Z_scaleSEXP, SEXP Z_shapeSEXP, SEXP stoch_locSEXP, SEXP stoch_volSEXP, SEXP stoch_shapeSEXP, SEXP NAsSEXP, SEXP guess_locSEXP, SEXP guess_scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type T_loc(T_locSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type T_scale(T_scaleSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type T_shape(T_shapeSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type Z_loc(Z_locSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type Z_scale(Z_scaleSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type Z_shape(Z_shapeSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_loc(stoch_locSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_vol(stoch_volSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_shape(stoch_shapeSEXP);
    Rcpp::traits::input_parameter< const MatrixXi& >::type NAs(NAsSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type guess_loc(guess_locSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type guess_scale(guess_scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_deriv(pars, y, T_loc, T_scale, T_shape, Z_loc, Z_scale, Z_shape, stoch_loc, stoch_vol, stoch_shape, NAs, guess_loc, guess_scale));
    return rcpp_result_gen;
END_RCPP
}
// loglik_hess
MatrixXd loglik_hess(const VectorXd& pars, const MatrixXd& y, MatrixXd& T_loc, MatrixXd& T_scale, MatrixXd& T_shape, MatrixXd& Z_loc, MatrixXd& Z_scale, MatrixXd& Z_shape, const bool& stoch_loc, const bool& stoch_vol, const bool& stoch_shape, const MatrixXi& NAs, const VectorXd& guess_loc, const VectorXd& guess_scale);
RcppExport SEXP _RcppScoreDrivenDFM_loglik_hess(SEXP parsSEXP, SEXP ySEXP, SEXP T_locSEXP, SEXP T_scaleSEXP, SEXP T_shapeSEXP, SEXP Z_locSEXP, SEXP Z_scaleSEXP, SEXP Z_shapeSEXP, SEXP stoch_locSEXP, SEXP stoch_volSEXP, SEXP stoch_shapeSEXP, SEXP NAsSEXP, SEXP guess_locSEXP, SEXP guess_scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type T_loc(T_locSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type T_scale(T_scaleSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type T_shape(T_shapeSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type Z_loc(Z_locSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type Z_scale(Z_scaleSEXP);
    Rcpp::traits::input_parameter< MatrixXd& >::type Z_shape(Z_shapeSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_loc(stoch_locSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_vol(stoch_volSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_shape(stoch_shapeSEXP);
    Rcpp::traits::input_parameter< const MatrixXi& >::type NAs(NAsSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type guess_loc(guess_locSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type guess_scale(guess_scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_hess(pars, y, T_loc, T_scale, T_shape, Z_loc, Z_scale, Z_shape, stoch_loc, stoch_vol, stoch_shape, NAs, guess_loc, guess_scale));
    return rcpp_result_gen;
END_RCPP
}
// loglik_list
Rcpp::List loglik_list(const Tvec<double>& pars, const Rcpp::List& model, const bool& give_grad, const bool& give_hess);
RcppExport SEXP _RcppScoreDrivenDFM_loglik_list(SEXP parsSEXP, SEXP modelSEXP, SEXP give_gradSEXP, SEXP give_hessSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Tvec<double>& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const bool& >::type give_grad(give_gradSEXP);
    Rcpp::traits::input_parameter< const bool& >::type give_hess(give_hessSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_list(pars, model, give_grad, give_hess));
    return rcpp_result_gen;
END_RCPP
}
// hstep_filter
Rcpp::List hstep_filter(const MatrixXd& y, const VectorXd& a_loc_0, const VectorXd& a_scale_0, const VectorXd& a_shape_0, const MatrixXd& Z_loc, const MatrixXd& Z_scale, const MatrixXd& Z_shape, const MatrixXd& T_loc, const MatrixXd& T_scale, const MatrixXd& T_shape, const MatrixXd& K_loc, const MatrixXd& K_scale, const MatrixXd& K_shape, const MatrixXd& parameters, const MatrixXi& NAs);
RcppExport SEXP _RcppScoreDrivenDFM_hstep_filter(SEXP ySEXP, SEXP a_loc_0SEXP, SEXP a_scale_0SEXP, SEXP a_shape_0SEXP, SEXP Z_locSEXP, SEXP Z_scaleSEXP, SEXP Z_shapeSEXP, SEXP T_locSEXP, SEXP T_scaleSEXP, SEXP T_shapeSEXP, SEXP K_locSEXP, SEXP K_scaleSEXP, SEXP K_shapeSEXP, SEXP parametersSEXP, SEXP NAsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type a_loc_0(a_loc_0SEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type a_scale_0(a_scale_0SEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type a_shape_0(a_shape_0SEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Z_loc(Z_locSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Z_scale(Z_scaleSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Z_shape(Z_shapeSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type T_loc(T_locSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type T_scale(T_scaleSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type T_shape(T_shapeSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type K_loc(K_locSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type K_scale(K_scaleSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type K_shape(K_shapeSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const MatrixXi& >::type NAs(NAsSEXP);
    rcpp_result_gen = Rcpp::wrap(hstep_filter(y, a_loc_0, a_scale_0, a_shape_0, Z_loc, Z_scale, Z_shape, T_loc, T_scale, T_shape, K_loc, K_scale, K_shape, parameters, NAs));
    return rcpp_result_gen;
END_RCPP
}
// filter_reg_list
Rcpp::List filter_reg_list(const Rcpp::List& model);
RcppExport SEXP _RcppScoreDrivenDFM_filter_reg_list(SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(filter_reg_list(model));
    return rcpp_result_gen;
END_RCPP
}
// loglik_reg_rcpp
double loglik_reg_rcpp(const Tvec<double>& pars, const Tvec<double>& y, const bool& stoch_loc, const bool& stoch_vol, const bool& stoch_shape, const Tvec<double>& X_loc, const Tvec<double>& X_scale, const Tvec<double>& X_shape, const Tvec<double>& NAs);
RcppExport SEXP _RcppScoreDrivenDFM_loglik_reg_rcpp(SEXP parsSEXP, SEXP ySEXP, SEXP stoch_locSEXP, SEXP stoch_volSEXP, SEXP stoch_shapeSEXP, SEXP X_locSEXP, SEXP X_scaleSEXP, SEXP X_shapeSEXP, SEXP NAsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Tvec<double>& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const Tvec<double>& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_loc(stoch_locSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_vol(stoch_volSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_shape(stoch_shapeSEXP);
    Rcpp::traits::input_parameter< const Tvec<double>& >::type X_loc(X_locSEXP);
    Rcpp::traits::input_parameter< const Tvec<double>& >::type X_scale(X_scaleSEXP);
    Rcpp::traits::input_parameter< const Tvec<double>& >::type X_shape(X_shapeSEXP);
    Rcpp::traits::input_parameter< const Tvec<double>& >::type NAs(NAsSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_reg_rcpp(pars, y, stoch_loc, stoch_vol, stoch_shape, X_loc, X_scale, X_shape, NAs));
    return rcpp_result_gen;
END_RCPP
}
// loglik_reg_list
Rcpp::List loglik_reg_list(const Tvec<double>& pars, const Rcpp::List& model);
RcppExport SEXP _RcppScoreDrivenDFM_loglik_reg_list(SEXP parsSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Tvec<double>& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_reg_list(pars, model));
    return rcpp_result_gen;
END_RCPP
}
// loglik_reg_deriv
VectorXd loglik_reg_deriv(const VectorXd& pars, const VectorXd& y, const bool& stoch_loc, const bool& stoch_vol, const bool& stoch_shape, const VectorXd& X_loc, const VectorXd& X_scale, const VectorXd& X_shape, const VectorXd& NAs);
RcppExport SEXP _RcppScoreDrivenDFM_loglik_reg_deriv(SEXP parsSEXP, SEXP ySEXP, SEXP stoch_locSEXP, SEXP stoch_volSEXP, SEXP stoch_shapeSEXP, SEXP X_locSEXP, SEXP X_scaleSEXP, SEXP X_shapeSEXP, SEXP NAsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_loc(stoch_locSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_vol(stoch_volSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_shape(stoch_shapeSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type X_loc(X_locSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type X_scale(X_scaleSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type X_shape(X_shapeSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type NAs(NAsSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_reg_deriv(pars, y, stoch_loc, stoch_vol, stoch_shape, X_loc, X_scale, X_shape, NAs));
    return rcpp_result_gen;
END_RCPP
}
// loglik_reg_hess
MatrixXd loglik_reg_hess(const VectorXd& pars, const VectorXd& y, const bool& stoch_loc, const bool& stoch_vol, const bool& stoch_shape, const VectorXd& X_loc, const VectorXd& X_scale, const VectorXd& X_shape, const VectorXd& NAs);
RcppExport SEXP _RcppScoreDrivenDFM_loglik_reg_hess(SEXP parsSEXP, SEXP ySEXP, SEXP stoch_locSEXP, SEXP stoch_volSEXP, SEXP stoch_shapeSEXP, SEXP X_locSEXP, SEXP X_scaleSEXP, SEXP X_shapeSEXP, SEXP NAsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_loc(stoch_locSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_vol(stoch_volSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stoch_shape(stoch_shapeSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type X_loc(X_locSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type X_scale(X_scaleSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type X_shape(X_shapeSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type NAs(NAsSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_reg_hess(pars, y, stoch_loc, stoch_vol, stoch_shape, X_loc, X_scale, X_shape, NAs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RcppScoreDrivenDFM_filter_list", (DL_FUNC) &_RcppScoreDrivenDFM_filter_list, 1},
    {"_RcppScoreDrivenDFM_loglik_rcpp", (DL_FUNC) &_RcppScoreDrivenDFM_loglik_rcpp, 14},
    {"_RcppScoreDrivenDFM_loglik_deriv", (DL_FUNC) &_RcppScoreDrivenDFM_loglik_deriv, 14},
    {"_RcppScoreDrivenDFM_loglik_hess", (DL_FUNC) &_RcppScoreDrivenDFM_loglik_hess, 14},
    {"_RcppScoreDrivenDFM_loglik_list", (DL_FUNC) &_RcppScoreDrivenDFM_loglik_list, 4},
    {"_RcppScoreDrivenDFM_hstep_filter", (DL_FUNC) &_RcppScoreDrivenDFM_hstep_filter, 15},
    {"_RcppScoreDrivenDFM_filter_reg_list", (DL_FUNC) &_RcppScoreDrivenDFM_filter_reg_list, 1},
    {"_RcppScoreDrivenDFM_loglik_reg_rcpp", (DL_FUNC) &_RcppScoreDrivenDFM_loglik_reg_rcpp, 9},
    {"_RcppScoreDrivenDFM_loglik_reg_list", (DL_FUNC) &_RcppScoreDrivenDFM_loglik_reg_list, 2},
    {"_RcppScoreDrivenDFM_loglik_reg_deriv", (DL_FUNC) &_RcppScoreDrivenDFM_loglik_reg_deriv, 9},
    {"_RcppScoreDrivenDFM_loglik_reg_hess", (DL_FUNC) &_RcppScoreDrivenDFM_loglik_reg_hess, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_RcppScoreDrivenDFM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
