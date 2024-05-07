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

template<typename T>
T Bound_f(const T &y, const double &a, const double &b){
  T x = a + (b-a)*(1/(1+exp(-y*2)));
  return x;
}

template<typename T>
T Bound_f_b(const T &y, const T &a, const T &b){
  T x = a + (b-a)*(1/(1+exp(-y*2)));
  return x;
}

template<typename T>
T fpsign(const T &x){
  T sign;
  if(x<0){
    sign= -1.0;
  } else {
    sign = 1.0;
  }
  return sign;
}

template<typename T>
T loggamma(T z) {
  double pi = 3.141592653589793238462643383280;
  // Coefficients used by the Lanczos approximation
  static const double g = 7;
  static const double lanczos[] = {
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
  };
  
  if (z < 0.5) {
    // Use the reflection formula
    return log(pi) - log(sin(pi * z)) - loggamma(1 - z);
  } else {
    z -= 1.0;
    
    T a = lanczos[0];
    T t = z + g + 0.5;
    for (int i = 1; i < g; i++) {
      a += lanczos[i] / (z + i);
    }
    
    return 0.5 * log(2 * pi) + (z + 0.5) * log(t) - t + log(a);
  }
}

// Calculate the coefficients used by Spouge's approximation (based on the C
// implemetation) taken from https://rosettacode.org/wiki/Gamma_function
RowVectorXd CalculateCoefficients(int numCoeff)
{
  RowVectorXd c(numCoeff);
  double k1_factrl = 1.0;
  double pi = 3.141592653589793238462643383280;
  
  c[0] = sqrt(2.0 * pi);
  for(size_t k=1; k < numCoeff; k++)
  {
    c[k] = exp(numCoeff-k) * pow(numCoeff-k, k-0.5) / k1_factrl;
    k1_factrl *= -(double)k;
  }
  return c;
}

// The Spouge approximation
template<typename T>
T Gamma(const RowVectorXd& coeffs, T x)
{
  const size_t numCoeff = coeffs.size();
  T accm = coeffs[0];
  for(size_t k=1; k < numCoeff; k++)
  {
    accm += coeffs[k] / ( x + k );
  }
  accm *= exp(-(x+numCoeff)) * pow(x+numCoeff, x+0.5);
  return accm/x;
}

template<typename T>
T Exp_Almon(const int &k,const int &K, const Trowvec<T> &theta_uc){
  Trowvec<T> theta(theta_uc.size());
  
  for (int i = 0; i<theta_uc.size(); i++){
    theta(i) = Bound_f(theta_uc(i),-0.1, 0.1);
  }
  T k2 = k;
  T a = theta(0)*k2+theta(1)*k2*k2;
  
  T c = 0;
  T b = 0;
  for (int Ki=0; Ki<K+1; Ki++){
    b = theta(0)*Ki+theta(1)*Ki*Ki;
    c = c + exp(b);
  }
  
  return exp(a)/c;
}

template<typename T>
T K_fun(const T &eta){
  double pi = 3.141592653589793238462643383280;
  RowVectorXd Coef = CalculateCoefficients(10);
  T K = loggamma((eta+1)/eta*0.5) - log(pi/eta)/2 - loggamma(1/(2*eta));
  
  return K;
}

template<typename T>
T gain(const T &pars, const Rcpp::String &model){
  
  T K = 0;
  
  if (model == "monthly") {
    K = exp(pars)/100;
  }
  
  if (model == "quarterly") {
    K = exp(pars)/100;
  }
  
  return K;
}

template<typename T>
T Skt_f(const T &v_t, const Trowvec<T> &parameters){
  T sigma = parameters[0];
  T nu1 = parameters[1];
  T alpha = parameters[3];
  
  T f_t;
  
  T eta = 1/nu1;
  f_t = exp(K_fun(eta) -0.5*log(sigma*sigma) - (1+eta)/(2*eta)*log(1 + eta*v_t*v_t/((1-fpsign(v_t)*alpha)*(1-fpsign(v_t)*alpha)*sigma*sigma)));
  
  if(f_t==f_t*2 or isnan(f_t) or f_t==0){
    f_t = 1e-100;
  }
  return f_t;
}

template<typename T>
T Skt_score(const T &v_t, const Trowvec<T> &parameters, const Rcpp::String &derivative){
  T sigma = parameters[0];
  T nu1 = parameters[1];
  T alpha = parameters[3];
  T score = 0.0;
  
  T eta = 1/nu1;
  T zeta = v_t/sigma;
  T w = (1+eta)/((1-fpsign(v_t)*alpha)*(1-fpsign(v_t)*alpha)+eta*zeta*zeta);
  if (derivative=="location"){
    score = score = 1/sigma*w*zeta;
  }
  
  if (derivative=="scale"){
    score = 1/(2*sigma*sigma)*(w*zeta*zeta-1);
  }
  
  if (derivative=="shape"){
    score = -fpsign(v_t)/(1-fpsign(v_t)*alpha)*w*zeta*zeta;
  }
  
  return score;
}

template<typename T>
T Skt_hess(const T &v_t, const Trowvec<T> &parameters, const Rcpp::String &derivative){
  T sigma = parameters[0];
  T nu1 = parameters[1];
  T alpha = parameters[3];
  
  T Hessian = 0.0;
    
  T eta = 1/nu1;
  
  if (derivative=="location"){
    Hessian = (1+eta)/(sigma*sigma*(1-alpha*alpha)*(1+3*eta));
  }
  
  if (derivative=="scale"){
    Hessian = 1/(2*(1+3*eta)*pow(sigma,4));
  }
  
  if (derivative=="shape"){
    Hessian = 3*(1+eta)/((1-alpha*alpha)*(1+3*eta));
  }
  
  return Hessian;
}

template<typename T>
Tmat<T> filter(const MatrixXd &y,
               const Tvec<T> &a_loc_init,
               const Tvec<T> &a_scale_init,
               const Tvec<T> &a_shape_init,
               const Tmat<T> &Z_loc,
               const Tmat<T> &Z_scale,
               const Tmat<T> &Z_shape,
               const Tmat<T> &T_loc,
               const Tmat<T> &T_scale,
               const Tmat<T> &T_shape,
               const Tmat<T> &K_loc,
               const Tmat<T> &K_scale,
               const Tmat<T> &K_shape,
               const Tmat<T> &parameters,
               const MatrixXi &NAs){
  
  Tmat<T> loglik(y.rows(),y.cols());
  loglik.setZero();
  
  Tmat<T> a_loc(T_loc.rows(), y.rows()+1);
  Tmat<T> a_scale(T_scale.rows(), y.rows()+1);
  Tmat<T> a_shape(T_shape.rows(), y.rows()+1);
  
  a_loc.col(0) = a_loc_init;
  a_scale.col(0) = a_scale_init;
  a_shape.col(0) = a_shape_init;
  
  Tvec<T> v(y.cols());
  int nb_not_NAs = 0;
  
  int id_cf_loc = y.cols()*5;
  int id_cf_scale = y.cols();
  int id_cf_shape = y.cols();
  
  T innov_loc;
  T innov_scale;
  T innov_shape;
  
  for (int t = 0; t < y.rows(); t++){
    
    for (int i = 0; i < y.cols(); i++) {
      v(i) = y(t,i) - Z_loc.row(i)*a_loc.col(t);
      a_loc(i*5+4, t+1) = a_loc(i*5+3, t);
      a_loc(i*5+3, t+1) = a_loc(i*5+2, t);
      a_loc(i*5+2, t+1) = a_loc(i*5+1, t);
      a_loc(i*5+1, t+1) = a_loc(i*5, t);
      a_loc(i*5, t+1) = a_loc(i*5, t);
      a_scale(i, t+1) = a_scale(i, t);
      a_shape(i, t+1) = a_shape(i, t);
    }
    
    nb_not_NAs = NAs.row(t).sum();
    
    Tvec<T> score_loc(nb_not_NAs);
    score_loc.setZero();
    Tvec<T> score_scale(nb_not_NAs);
    score_scale.setZero();
    Tvec<T> score_shape(nb_not_NAs);
    score_shape.setZero();
    
    innov_loc = 0;
    innov_scale = 0;
    innov_shape = 0;
    
    if(nb_not_NAs>0){
      int n_pars=0;
      for (int i = 0; i<y.cols(); i++){
        if(NAs(t,i)==1){
          Trowvec<T> parameters_i = parameters.row(i);
          
          T signal_scale = Z_scale.row(i)*a_scale.col(t);
          T signal_shape = Z_shape.row(i)*a_shape.col(t);
          
          parameters_i(0) = exp(signal_scale);
          parameters_i(3) = tanh(signal_shape);
          loglik(t,i) = log(Skt_f(v(i),parameters_i));
          
          score_loc(n_pars) = Skt_score(v(i), parameters_i, "location");
          score_scale(n_pars) = Skt_score(v(i), parameters_i, "scale")*2*exp(2*signal_scale);
          score_shape(n_pars) = Skt_score(v(i), parameters_i, "shape")*(1-tanh(signal_shape)*tanh(signal_shape));
                    
          a_loc(i*5+4, t+1) = a_loc(i*5+3, t);
          a_loc(i*5+3, t+1) = a_loc(i*5+2, t);
          a_loc(i*5+2, t+1) = a_loc(i*5+1, t);
          a_loc(i*5, t+1) = a_loc(i*5, t) + K_loc(i*5,i*5)*Z_loc(i,i*5)*score_loc(n_pars);
          a_loc(i*5+1, t+1) = a_loc(i*5, t+1);
          a_scale(i, t+1) = a_scale(i, t) + K_scale(i,i)*Z_scale(i,i)*score_scale(n_pars);
          a_shape(i, t+1) = a_shape(i, t) + K_shape(i,i)*Z_shape(i,i)*score_shape(n_pars);
          
          innov_loc += Z_loc(i,id_cf_loc)*score_loc(n_pars);
          innov_scale += Z_scale(i,id_cf_scale)*score_scale(n_pars);
          innov_shape += Z_shape(i,id_cf_shape)*score_shape(n_pars);
          
          n_pars += 1;
        }
      }
    }
    
    a_loc(id_cf_loc+4, t+1) = a_loc(id_cf_loc+3, t);
    a_loc(id_cf_loc+3, t+1) = a_loc(id_cf_loc+2, t);
    a_loc(id_cf_loc+2, t+1) = a_loc(id_cf_loc+1, t);
    a_loc(id_cf_loc, t) += K_loc(id_cf_loc,id_cf_loc)*innov_loc;
    a_loc(id_cf_loc, t+1) = T_loc(id_cf_loc,id_cf_loc)*a_loc(id_cf_loc, t) + T_loc(id_cf_loc,id_cf_loc+1)*a_loc(id_cf_loc+1, t);
    a_loc(id_cf_loc+1, t+1) = a_loc(id_cf_loc, t);
    
    a_scale(id_cf_scale+1, t+1) = a_scale(id_cf_scale, t);
    a_scale(id_cf_scale, t) += K_scale(id_cf_scale,id_cf_scale)*innov_scale;
    a_scale(id_cf_scale, t+1) = T_scale(id_cf_scale,id_cf_scale)*a_scale(id_cf_scale, t);
    
    a_shape(id_cf_shape+1, t+1) = a_shape(id_cf_shape, t);
    a_shape(id_cf_shape, t) +=  K_shape(id_cf_shape,id_cf_shape)*innov_shape;
    a_shape(id_cf_shape, t+1) = T_shape(id_cf_shape,id_cf_shape)*a_shape(id_cf_shape, t);
  }
  
  return loglik;
}

// [[Rcpp::export]]
Rcpp::List filter_list(const Rcpp::List &model){
  
  MatrixXd y = model["y"];
  Tvec<double> a_loc_init = model["a_loc_init"];
  Tvec<double> a_scale_init = model["a_scale_init"];
  Tvec<double> a_shape_init = model["a_shape_init"];
  Tmat<double> Z_loc = model["Z_loc"];
  Tmat<double> Z_scale = model["Z_scale"];
  Tmat<double> Z_shape = model["Z_shape"];
  Tmat<double> T_loc = model["T_loc"];
  Tmat<double> T_scale = model["T_scale"];
  Tmat<double> T_shape = model["T_shape"];
  Tmat<double> K_loc = model["K_loc"];
  Tmat<double> K_scale = model["K_scale"];
  Tmat<double> K_shape = model["K_shape"];
  Tmat<double> parameters = model["parameters"];
  MatrixXi NAs = model["NAs"];
  
  Tmat<double> loglik(y.rows(),y.cols());
  loglik.setZero();
  
  Tmat<double> a_loc(T_loc.rows(), y.rows()+1);
  Tmat<double> a_scale(T_scale.rows(), y.rows()+1);
  Tmat<double> a_shape(T_shape.rows(), y.rows()+1);
  
  Tmat<double> a_loc_f(T_loc.rows(), y.rows());
  Tmat<double> a_scale_f(T_scale.rows(), y.rows());
  Tmat<double> a_shape_f(T_shape.rows(), y.rows());
  
  Tmat<double> location_f(y.rows(), y.cols());
  Tmat<double> scale_f(y.rows(), y.cols());
  Tmat<double> shape_f(y.rows(), y.cols());
  
  a_loc.col(0) = a_loc_init;
  a_scale.col(0) = a_scale_init;
  a_shape.col(0) = a_shape_init;
  
  Tvec<double> v(y.cols());
  int nb_not_NAs = 0;
  
  for (int t = 0; t < y.rows(); t++){
    
    v = y.row(t) - Z_loc*a_loc.col(t);
    
    nb_not_NAs = NAs.row(t).sum();
    
    Tvec<double> score_loc(nb_not_NAs);
    score_loc.setZero();
    Tvec<double> score_scale(nb_not_NAs);
    score_scale.setZero();
    Tvec<double> score_shape(nb_not_NAs);
    score_shape.setZero();
    
    Tmat<double> Z_t(nb_not_NAs,T_loc.rows());
    Z_t.setZero();
    Tmat<double> Z_scale_t(nb_not_NAs,T_scale.rows());
    Z_scale_t.setZero();
    Tmat<double> Z_shape_t(nb_not_NAs,T_shape.rows());
    Z_shape_t.setZero();
    
    Tvec<double> signal_scale = Z_scale*a_scale.col(t);
    Tvec<double> signal_shape = Z_shape*a_shape.col(t);
    
    if(nb_not_NAs>0){
      int n_pars=0;
      for (int i = 0; i<y.cols(); i++){
        if(NAs(t,i)==1){
          Trowvec<double> parameters_i = parameters.row(i);
          
          parameters_i(0) = exp(signal_scale(i));
          parameters_i(3) = tanh(signal_shape(i));
          loglik(t,i) = log(Skt_f(v(i),parameters_i));
          
          score_loc(n_pars) = Skt_score(v(i), parameters_i, "location");
          score_scale(n_pars) = Skt_score(v(i), parameters_i, "scale")*2*exp(2*signal_scale[i]);
          score_shape(n_pars) = Skt_score(v(i), parameters_i, "shape")*(1-tanh(signal_shape[i])*tanh(signal_shape[i]));
           
          Z_t.row(n_pars) = Z_loc.row(i);
          Z_scale_t.row(n_pars) = Z_scale.row(i);
          Z_shape_t.row(n_pars) = Z_shape.row(i);
          
          n_pars += 1;
        }
      }
    }
    
    a_loc_f.col(t) = a_loc.col(t) + K_loc*Z_t.transpose()*score_loc;
    a_scale_f.col(t) = a_scale.col(t) + K_scale*Z_scale_t.transpose()*score_scale;
    a_shape_f.col(t) = a_shape.col(t) + K_shape*Z_shape_t.transpose()*score_shape;
    
    for (int i = 0; i < y.cols(); i++) {
      location_f(t, i) = Z_loc.row(i)*a_loc_f.col(t);
      scale_f(t, i) = exp(Z_scale.row(i)*a_scale_f.col(t));
      shape_f(t, i) = tanh(Z_shape.row(i)*a_shape_f.col(t));
    }
    
    a_loc.col(t+1) = T_loc*a_loc_f.col(t);
    a_scale.col(t+1) = T_scale*a_scale_f.col(t);
    a_shape.col(t+1) = T_shape*a_shape_f.col(t);
  }
  
  Rcpp::List output = Rcpp::List::create(Rcpp::Named("a") = a_loc,
                                         Rcpp::Named("a_scale") = a_scale,
                                         Rcpp::Named("a_shape") = a_shape,
                                         Rcpp::Named("a_f") = a_loc_f,
                                         Rcpp::Named("a_scale_f") = a_scale_f,
                                         Rcpp::Named("a_shape_f") = a_shape_f,
                                         Rcpp::Named("location_f") = location_f,
                                         Rcpp::Named("scale_f") = scale_f,
                                         Rcpp::Named("shape_f") = shape_f,
                                         Rcpp::Named("parameters") = parameters,
                                         Rcpp::Named("loglik") = loglik,
                                         Rcpp::Named("model") = model);
  
  return output;
}

template<typename T>
T loglik(const Tvec<T> &pars,
         const MatrixXd &y,
         Tmat<T> &T_loc,
         Tmat<T> &T_scale,
         Tmat<T> &T_shape,
         Tmat<T> &Z_loc,
         Tmat<T> &Z_scale,
         Tmat<T> &Z_shape,
         const bool &stoch_loc,
         const bool &stoch_vol,
         const bool &stoch_shape,
         const MatrixXi &NAs,
         const VectorXd &guess_loc,
         const VectorXd &guess_scale){
  
  int n_pars = 0;
  
  Tmat<T>  parameters(y.cols(),4);
  
  Tvec<T> a_loc_init(T_loc.rows());
  a_loc_init.setZero();
  Tvec<T> a_scale_init(T_scale.rows());
  a_scale_init.setZero();
  Tvec<T> a_shape_init(T_shape.rows());
  a_shape_init.setZero();
  
  Tmat<T> K_loc(T_loc.rows(),T_loc.rows());
  K_loc.setZero();
  Tmat<T> K_scale(T_scale.rows(),T_scale.rows());
  K_scale.setZero();
  Tmat<T> K_shape(T_shape.rows(),T_shape.rows());
  K_shape.setZero();
  
  T guess;
  T par;
  T factor_loading;
  
  Rcpp::String version;
  
  if (Z_loc(0,0)<1.0) {
    version = "quarterly";
  } else {
    version = "monthly";
  }
  
  //indexes of the common factors in the observation matrices
  int id_cf_loc = y.cols()*5;
  int id_cf_scale = y.cols();
  int id_cf_shape = y.cols();
  
  for (int i = 0; i<y.cols(); i++){
    
    parameters(i,1) = Bound_f(pars(n_pars++)-3,1.1,1e3);
    parameters(i,2) = parameters(i,1);
    
    guess = guess_loc(i) + pars(n_pars++);
    for (int j=0; j<5; j++){
      a_loc_init(i*5+j) = guess;
    }
    K_loc(i*5,i*5) = gain(pars(n_pars++), version);
    a_scale_init(i) = log(guess_scale(i)) + pars(n_pars++); 
  }
  
  Tvec<int> index_loc(y.cols() - 1);
  index_loc.setZero();
  // Common factor in location
  if (stoch_loc) {
    int j_loc = 0;
    
    for (int i = 0; i<y.cols(); i++){
      
      if(i>0){
        if(Z_loc(i,i*5)<1.0){
          //if series i is monthly
          factor_loading = pars(n_pars++);
          Z_loc(i,id_cf_loc) = factor_loading/3;
          Z_loc(i,id_cf_loc+1) = factor_loading*2/3;
          Z_loc(i,id_cf_loc+2) = factor_loading;
          Z_loc(i,id_cf_loc+3) = factor_loading*2/3;
          Z_loc(i,id_cf_loc+4) = factor_loading/3;
        } else {
          // if series i is quarterly
          Z_loc(i,id_cf_loc) = pars(n_pars++);
        }
        
        index_loc(j_loc++) = n_pars;
      }
    }
    K_loc(id_cf_loc,id_cf_loc) = gain(pars(n_pars++), version);
    T b = tanh(pars(n_pars++));
    T a = Bound_f_b(pars(n_pars++)+1,b-1,1-b);
    T_loc(id_cf_loc,id_cf_loc) = a;
    T_loc(id_cf_loc,id_cf_loc+1) = b;
  }
  
  Tvec<int> index_scale(y.cols() - 1);
  index_scale.setZero();
  // Common factor in scale
  if(stoch_vol){
    int j_scale = 0;
    
    for (int i = 0; i<y.cols(); i++){
      K_scale(i,i) = gain(pars(n_pars++), version);
      
      if (i > 0) {
        Z_scale(i,id_cf_scale) = pars(n_pars++);
        index_scale(j_scale++) = n_pars;
      }
    }
    
    K_scale(id_cf_scale,id_cf_scale) = gain(pars(n_pars++), version);
    T_scale(id_cf_scale,id_cf_scale) = tanh(pars(n_pars++)+1);
  }
  
  Tvec<int> index_shape(y.cols() - 1);
  index_shape.setZero();
  // Common factor in shape
  if(stoch_shape){
    int j_shape = 0;
    
    for (int i = 0; i<y.cols(); i++){
      a_shape_init(i) = pars(n_pars++);
      K_shape(i,i) = gain(pars(n_pars++), version);
      if (i > 0) {
        Z_shape(i,id_cf_shape) = pars(n_pars++);
        index_shape(j_shape++) = n_pars;
      }
    }
    
    K_shape(id_cf_shape,id_cf_shape) = gain(pars(n_pars++), version);
    T_shape(id_cf_shape,id_cf_shape) = tanh(pars(n_pars++)+1);
  }
  
  Tmat<T> loglik = filter(y,
                          a_loc_init, a_scale_init, a_shape_init,
                          Z_loc, Z_scale, Z_shape,
                          T_loc, T_scale, T_shape,
                          K_loc, K_scale, K_shape,
                          parameters,
                          NAs);
  
  T LL = -loglik.sum();
  
  return LL;
}

// [[Rcpp::export]]
double loglik_rcpp(const VectorXd &pars,
                   const MatrixXd &y,
                   MatrixXd &T_loc,
                   MatrixXd &T_scale,
                   MatrixXd &T_shape,
                   MatrixXd &Z_loc,
                   MatrixXd &Z_scale,
                   MatrixXd &Z_shape,
                   const bool &stoch_loc,
                   const bool &stoch_vol,
                   const bool &stoch_shape,
                   const MatrixXi &NAs,
                   const VectorXd &guess_loc,
                   const VectorXd &guess_scale)
{
  return loglik(pars, y,
                T_loc,
                T_scale,
                T_shape,
                Z_loc,
                Z_scale,
                Z_shape,
                stoch_loc, stoch_vol, stoch_shape,
                NAs,guess_loc,guess_scale);
}

// [[Rcpp::export]]
VectorXd loglik_deriv(const VectorXd &pars,
                      const MatrixXd &y,
                      MatrixXd &T_loc,
                      MatrixXd &T_scale,
                      MatrixXd &T_shape,
                      MatrixXd &Z_loc,
                      MatrixXd &Z_scale,
                      MatrixXd &Z_shape,
                      const bool &stoch_loc,
                      const bool &stoch_vol,
                      const bool &stoch_shape,
                      const MatrixXi &NAs,
                      const VectorXd &guess_loc,
                      const VectorXd &guess_scale)
{
  
  //Creating AD model matrices and vectors
  ////////////////////////////////////////////////////////
  Tvec<a_double> apars(pars.size());
  Tmat<a_double> aT_loc(T_loc.rows(),T_loc.cols());
  Tmat<a_double> aT_scale(T_scale.rows(),T_scale.cols());
  Tmat<a_double> aT_shape(T_shape.rows(),T_shape.cols());
  Tmat<a_double> aZ_loc(Z_loc.rows(),Z_loc.cols());
  Tmat<a_double> aZ_scale(Z_scale.rows(),Z_scale.cols());
  Tmat<a_double> aZ_shape(Z_shape.rows(),Z_shape.cols());
  
  for (int i = 0; i < pars.size(); ++i){
    apars(i) = pars(i);
  }
  for (int row = 0; row < aT_loc.rows(); ++row){
    for (int col = 0; col < aT_loc.cols(); ++col){
      aT_loc(row,col) = T_loc(row,col);
    }
  }
  for (int row = 0; row < aT_scale.rows(); ++row){
    for (int col = 0; col < aT_scale.cols(); ++col){
      aT_scale(row,col) = T_scale(row,col);
    }
  }
  for (int row = 0; row < aT_shape.rows(); ++row){
    for (int col = 0; col < aT_shape.cols(); ++col){
      aT_shape(row,col) = T_shape(row,col);
    }
  }
  for (int row = 0; row < aZ_loc.rows(); ++row){
    for (int col = 0; col < aZ_loc.cols(); ++col){
      aZ_loc(row,col) = Z_loc(row,col);
    }
  }
  for (int row = 0; row < aZ_scale.rows(); ++row){
    for (int col = 0; col < aZ_scale.cols(); ++col){
      aZ_scale(row,col) = Z_scale(row,col);
    }
  }
  for (int row = 0; row < aZ_shape.rows(); ++row){
    for (int col = 0; col < aZ_shape.cols(); ++col){
      aZ_shape(row,col) = Z_shape(row,col);
    }
  }
  ////////////////////////////////////////////////////////
  
  a_vector LogLik(1);
  ADFun<double> grad;
  Independent(apars);
  LogLik(0) = loglik(apars, y,
         aT_loc,
         aT_scale,
         aT_shape,
         aZ_loc,
         aZ_scale,
         aZ_shape,
         stoch_loc, stoch_vol, stoch_shape,
         NAs,guess_loc,guess_scale);
  grad.Dependent(apars, LogLik);
  return grad.Jacobian(pars);
}

// [[Rcpp::export]]
MatrixXd loglik_hess(const VectorXd &pars,
                     const MatrixXd &y,
                     MatrixXd &T_loc,
                     MatrixXd &T_scale,
                     MatrixXd &T_shape,
                     MatrixXd &Z_loc,
                     MatrixXd &Z_scale,
                     MatrixXd &Z_shape,
                     const bool &stoch_loc,
                     const bool &stoch_vol,
                     const bool &stoch_shape,
                     const MatrixXi &NAs,
                     const VectorXd &guess_loc,
                     const VectorXd &guess_scale)
{
  
  //Creating AD model matrices and vectors
  ////////////////////////////////////////////////////////
  Tvec<a_double> apars(pars.size());
  Tmat<a_double> aT_loc(T_loc.rows(),T_loc.cols());
  Tmat<a_double> aT_scale(T_scale.rows(),T_scale.cols());
  Tmat<a_double> aT_shape(T_shape.rows(),T_shape.cols());
  Tmat<a_double> aZ_loc(Z_loc.rows(),Z_loc.cols());
  Tmat<a_double> aZ_scale(Z_scale.rows(),Z_scale.cols());
  Tmat<a_double> aZ_shape(Z_shape.rows(),Z_shape.cols());
  
  for (int i = 0; i < pars.size(); ++i){
    apars(i) = pars(i);
  }
  for (int row = 0; row < aT_loc.rows(); ++row){
    for (int col = 0; col < aT_loc.cols(); ++col){
      aT_loc(row,col) = T_loc(row,col);
    }
  }
  for (int row = 0; row < aT_scale.rows(); ++row){
    for (int col = 0; col < aT_scale.cols(); ++col){
      aT_scale(row,col) = T_scale(row,col);
    }
  }
  for (int row = 0; row < aT_shape.rows(); ++row){
    for (int col = 0; col < aT_shape.cols(); ++col){
      aT_shape(row,col) = T_shape(row,col);
    }
  }
  for (int row = 0; row < aZ_loc.rows(); ++row){
    for (int col = 0; col < aZ_loc.cols(); ++col){
      aZ_loc(row,col) = Z_loc(row,col);
    }
  }
  for (int row = 0; row < aZ_scale.rows(); ++row){
    for (int col = 0; col < aZ_scale.cols(); ++col){
      aZ_scale(row,col) = Z_scale(row,col);
    }
  }
  for (int row = 0; row < aZ_shape.rows(); ++row){
    for (int col = 0; col < aZ_shape.cols(); ++col){
      aZ_shape(row,col) = Z_shape(row,col);
    }
  }
  ////////////////////////////////////////////////////////
  
  a_vector LogLik(1);
  ADFun<double> grad;
  Independent(apars);
  VectorXd vector_hessian(apars.size()*apars.size());
  
  LogLik(0) = loglik(apars, y,
         aT_loc,
         aT_scale,
         aT_shape,
         aZ_loc,
         aZ_scale,
         aZ_shape,
         stoch_loc,
         stoch_vol, stoch_shape,NAs,guess_loc,guess_scale);
  grad.Dependent(apars, LogLik);
  vector_hessian = grad.Hessian(pars,0);
  
  MatrixXd hessian(apars.size(),apars.size());
  int cursor = 0;
  for(int row = 0; row < apars.size(); row++)
  {
    for(int col = 0; col < apars.size(); col++)
    {
      hessian(row,col) = vector_hessian[cursor++]; 
    } 
  }
  
  return hessian;
}

// [[Rcpp::export]]
Rcpp::List loglik_list(const Tvec<double> &pars,
                       const Rcpp::List &model,
                       const bool &give_grad = false,
                       const bool &give_hess = false){
  
  MatrixXd y = model["y"];
  MatrixXd T_loc = model["T_loc"];
  MatrixXd T_scale = model["T_scale"];
  MatrixXd T_shape = model["T_shape"];
  MatrixXd Z_loc = model["Z_loc"];
  MatrixXd Z_scale = model["Z_scale"];
  MatrixXd Z_shape = model["Z_shape"];
  MatrixXi NAs = model["NAs"];
  Rcpp::List options = model["options"];
  bool stoch_loc = options["stoch_loc"];
  bool stoch_vol = options["stoch_vol"];
  bool stoch_shape = options["stoch_shape"];
  VectorXd guess_loc = model["guess_loc"];
  VectorXd guess_scale = model["guess_scale"];
  
  int n_pars = 0;
  
  Tmat<double>  parameters(y.cols(),4);
  
  Tvec<double> a_loc_init(T_loc.rows());
  a_loc_init.setZero();
  Tvec<double> a_scale_init(T_scale.rows());
  a_scale_init.setZero();
  Tvec<double> a_shape_init(T_shape.rows());
  a_shape_init.setZero();
  
  Tmat<double> K_loc(T_loc.rows(),T_loc.rows());
  K_loc.setZero();
  Tmat<double> K_scale(T_scale.rows(),T_scale.rows());
  K_scale.setZero();
  Tmat<double> K_shape(T_shape.rows(),T_shape.rows());
  K_shape.setZero();
  
  double guess;
  double par;
  double factor_loading;
  
  //indexes of the common factors in the observation matrices
  int id_cf_loc = y.cols()*5;
  int id_cf_scale = y.cols();
  int id_cf_shape = y.cols();
  
  Rcpp::String version;
  
  if (Z_loc(0,0)<1.0) {
    version = "quarterly";
  } else {
    version = "monthly";
  }
  
  for (int i = 0; i<y.cols(); i++){
    
    parameters(i,1) = Bound_f(pars(n_pars++)-3,1.1,1e3);
    parameters(i,2) = parameters(i,1);
    
    guess = guess_loc(i) + pars(n_pars++);
    for (int j=0; j<5; j++){
      a_loc_init(i*5+j) = guess;
    }
    K_loc(i*5,i*5) = gain(pars(n_pars++), version);
    a_scale_init(i) = log(guess_scale(i)) + pars(n_pars++); 
  }
  
  Tvec<int> index_loc(y.cols() - 1);
  index_loc.setZero();
  // Common factor in location
  if (stoch_loc) {
    int j_loc = 0;
    
    for (int i = 0; i<y.cols(); i++){
      
      if(i>0){
        if(Z_loc(i,i*5)<1.0){
          //if series i is monthly
          factor_loading = pars(n_pars++);
          Z_loc(i,id_cf_loc) = factor_loading/3;
          Z_loc(i,id_cf_loc+1) = factor_loading*2/3;
          Z_loc(i,id_cf_loc+2) = factor_loading;
          Z_loc(i,id_cf_loc+3) = factor_loading*2/3;
          Z_loc(i,id_cf_loc+4) = factor_loading/3;
        } else {
          // if series i is quarterly
          Z_loc(i,id_cf_loc) = pars(n_pars++);
        }
        
        index_loc(j_loc++) = n_pars;
      }
    }
    K_loc(id_cf_loc,id_cf_loc) = gain(pars(n_pars++), version);
    double b = tanh(pars(n_pars++));
    double a = Bound_f_b(pars(n_pars++)+1,b-1,1-b);
    T_loc(id_cf_loc,id_cf_loc) = a;
    T_loc(id_cf_loc,id_cf_loc+1) = b;
  }
  
  Tvec<int> index_scale(y.cols() - 1);
  index_scale.setZero();
  // Common factor in scale
  if(stoch_vol){
    int j_scale = 0;
    
    for (int i = 0; i<y.cols(); i++){
      K_scale(i,i) = gain(pars(n_pars++), version);
      
      if (i > 0) {
        Z_scale(i,id_cf_scale) = pars(n_pars++);
        index_scale(j_scale++) = n_pars;
      }
    }
    
    K_scale(id_cf_scale,id_cf_scale) = gain(pars(n_pars++), version);
    T_scale(id_cf_scale,id_cf_scale) = tanh(pars(n_pars++)+1);
  }
  
  Tvec<int> index_shape(y.cols() - 1);
  index_shape.setZero();
  // Common factor in shape
  if(stoch_shape){
    int j_shape = 0;
    
    for (int i = 0; i<y.cols(); i++){
      a_shape_init(i) = pars(n_pars++);
      K_shape(i,i) = gain(pars(n_pars++), version);
      if (i > 0) {
        Z_shape(i,id_cf_shape) = pars(n_pars++);
        index_shape(j_shape++) = n_pars;
      }
    }
    
    K_shape(id_cf_shape,id_cf_shape) = gain(pars(n_pars++), version);
    T_shape(id_cf_shape,id_cf_shape) = tanh(pars(n_pars++)+1);
  }
  
  Rcpp::List model_new = Rcpp::List::create(Rcpp::Named("y") = y,
                                            Rcpp::Named("a_loc_init") = a_loc_init,
                                            Rcpp::Named("a_scale_init") = a_scale_init,
                                            Rcpp::Named("a_shape_init") = a_shape_init,
                                            Rcpp::Named("Z_loc") = Z_loc,
                                            Rcpp::Named("Z_scale") = Z_scale,
                                            Rcpp::Named("Z_shape") = Z_shape,
                                            Rcpp::Named("T_loc") = T_loc,
                                            Rcpp::Named("T_scale") = T_scale,
                                            Rcpp::Named("T_shape") = T_shape,
                                            Rcpp::Named("K_loc") = K_loc,
                                            Rcpp::Named("K_scale") = K_scale,
                                            Rcpp::Named("K_shape") = K_shape,
                                            Rcpp::Named("parameters") = parameters,
                                            Rcpp::Named("index_loc") = index_loc,
                                            Rcpp::Named("index_scale") = index_scale,
                                            Rcpp::Named("index_shape") = index_shape,
                                            Rcpp::Named("NAs") = NAs,
                                            Rcpp::Named("n_pars") = n_pars);
  
  Rcpp::List output = filter_list(model_new);
  
  Tmat<double> loglik = filter(y,
                               a_loc_init, a_scale_init, a_shape_init,
                               Z_loc, Z_scale, Z_shape,
                               T_loc, T_scale, T_shape,
                               K_loc, K_scale, K_shape,
                               parameters,
                               NAs);
  
  VectorXd gradient(pars.size());
  gradient.setZero();
  if (give_grad) {
    gradient = loglik_deriv(pars,
                            y,
                            T_loc,
                            T_scale,
                            T_shape,
                            Z_loc,
                            Z_scale,
                            Z_shape,
                            stoch_loc,
                            stoch_vol,
                            stoch_shape,
                            NAs,
                            guess_loc,
                            guess_scale); 
  }
  
  MatrixXd hessian(pars.size(),pars.size());
  hessian.setZero();
  if (give_hess) {
    hessian = loglik_hess(pars,
                          y,
                          T_loc,
                          T_scale,
                          T_shape,
                          Z_loc,
                          Z_scale,
                          Z_shape,
                          stoch_loc,
                          stoch_vol,
                          stoch_shape,
                          NAs,
                          guess_loc,
                          guess_scale); 
  }
  
  return output;
}

// [[Rcpp::export]]
Rcpp::List hstep_filter(const MatrixXd &y,
                        const VectorXd &a_loc_0,
                        const VectorXd &a_scale_0,
                        const VectorXd &a_shape_0,
                        const MatrixXd &Z_loc,
                        const MatrixXd &Z_scale,
                        const MatrixXd &Z_shape,
                        const MatrixXd &T_loc,
                        const MatrixXd &T_scale,
                        const MatrixXd &T_shape,
                        const MatrixXd &K_loc,
                        const MatrixXd &K_scale,
                        const MatrixXd &K_shape,
                        const MatrixXd &parameters,
                        const MatrixXi &NAs){
  
  VectorXd v(y.cols());
  int nb_not_NAs = 0;
  
  VectorXd a_loc_t = a_loc_0;
  VectorXd a_scale_t = a_scale_0;
  VectorXd a_shape_t = a_shape_0;
  
  VectorXd a_loc_t_f = a_loc_0;
  VectorXd a_scale_t_f = a_scale_0;
  VectorXd a_shape_t_f = a_shape_0;
  
  for (int t = 0; t < y.rows(); t++){
    
    for (int i=0; i<y.cols();i++){
      v(i) = y(t,i) - Z_loc.row(i)*a_loc_t;
    }
    
    nb_not_NAs = NAs.row(t).sum();
    
    VectorXd scaled_score_loc(nb_not_NAs);
    scaled_score_loc.setZero();
    VectorXd scaled_score_scale(nb_not_NAs);
    scaled_score_scale.setZero();
    VectorXd scaled_score_shape(nb_not_NAs);
    scaled_score_shape.setZero();
    
    MatrixXd Z_t(nb_not_NAs,a_loc_t.size());
    Z_t.setZero();
    MatrixXd Z_scale_t(nb_not_NAs,a_scale_t.size());
    Z_scale_t.setZero();
    MatrixXd Z_shape_t(nb_not_NAs,a_shape_t.size());
    Z_shape_t.setZero();
    
    VectorXd signal_scale = Z_scale*a_scale_t;
    VectorXd signal_shape = Z_shape*a_shape_t;
    
    a_loc_t_f = a_loc_t;
    a_scale_t_f = a_scale_t;
    a_shape_t_f = a_shape_t;
    
    if(nb_not_NAs>0){
      int n_pars=0;
      for (int i = 0; i<y.cols(); i++){
        if(NAs(t,i)==1){
          Trowvec<double> parameters_i = parameters.row(i);
          
          parameters_i(0) = exp(signal_scale[i]);
          parameters_i(3) = tanh(signal_shape[i]);
          
          scaled_score_loc(n_pars) = Skt_score(v(i), parameters_i, "location");
          scaled_score_scale(n_pars) = Skt_score(v(i), parameters_i, "scale")*2*exp(2*signal_scale[i]);
          scaled_score_shape(n_pars) = Skt_score(v(i), parameters_i, "shape")*(1-parameters_i(3)*parameters_i(3));
          
          Z_t.row(n_pars) = Z_loc.row(i);
          Z_scale_t.row(n_pars) = Z_scale.row(i);
          Z_shape_t.row(n_pars) = Z_shape.row(i);
          
          n_pars += 1;
        }
      }
      a_loc_t_f = a_loc_t + K_loc*Z_t.transpose()*scaled_score_loc;
      a_scale_t_f = a_scale_t + K_scale*Z_scale_t.transpose()*scaled_score_scale;
      a_shape_t_f = a_shape_t + K_shape*Z_shape_t.transpose()*scaled_score_shape;
    }
    
    a_loc_t = T_loc*a_loc_t_f;
    a_scale_t = T_scale*a_scale_t_f;
    a_shape_t = T_shape*a_shape_t_f;
  }
  
  Rcpp::List output = Rcpp::List::create(Rcpp::Named("a_loc_t") = a_loc_t,
                                         Rcpp::Named("a_scale_t") = a_scale_t,
                                         Rcpp::Named("a_shape_t") = a_shape_t,
                                         Rcpp::Named("a_loc_t_f") = a_loc_t_f,
                                         Rcpp::Named("a_scale_t_f") = a_scale_t_f,
                                         Rcpp::Named("a_shape_t_f") = a_shape_t_f,
                                         Rcpp::Named("v") = v);
  
  return output;
}

template<typename T>
Tvec<T> filter_reg(const Tvec<double> &y,
                   const Tvec<T> &a_init,
                   const Tvec<T> &a_scale_init,
                   const Tvec<T> &a_shape_init,
                   const T &K,
                   const T &k_s,
                   const T &k_shape,
                   const Tvec<T> &X_loc,
                   const Tvec<T> &X_scale,
                   const Tvec<T> &X_shape,
                   const Trowvec<T> &parameters,
                   const Tvec<double> &NAs){
  
  Trowvec<T> parameters_b = parameters;
  
  Tvec<T> loglik(y.size());
  loglik.setZero();
  
  Tvec<T> a_t = a_init;
  Tvec<T> a_scale_t = a_scale_init;
  Tvec<T> a_shape_t = a_shape_init;
  
  T v;
  
  T S_scale_t;
  T S_loc_t;
  T S_shape_t;
  
  T score_loc;
  T score_scale;
  T score_shape;
  
  T location;
  T scale;
  T shape;
  
  bool there_is_nan;
  
  for (int t = 0; t < y.size(); t++){
    
    there_is_nan = false;
    if(NAs.size()>1){
      for (int j = 1; j<NAs.size(); j++){
        if (t==NAs[j] - 1){
          there_is_nan = true;
        }
      }
    }
    
    location = a_t(0) + X_loc(t)*a_t(1);
    v = y(t) - location;
    scale = a_scale_t(0) + X_scale(t)*a_scale_t(1);
    parameters_b(0) = exp(scale);
    shape = a_shape_t(0) + X_shape(t)*a_shape_t(1);
    parameters_b(3) = tanh(shape)*1;
    
    if(!there_is_nan){
      
      loglik(t) = log(Skt_f(v,parameters_b));
      
      S_loc_t =sqrt(Skt_hess(v, parameters_b, "location"));
      score_loc = Skt_score(v, parameters_b, "location")/S_loc_t;
      a_t(0) += K*score_loc;
      
      S_scale_t = sqrt(Skt_hess(v, parameters_b, "scale")*4*exp(4*scale));
      score_scale = Skt_score(v, parameters_b, "scale")*2*exp(2*scale)/S_scale_t;
      a_scale_t(0) += k_s*score_scale;
      
      S_shape_t = sqrt(Skt_hess(v, parameters_b, "shape")*pow(1-tanh(shape)*tanh(shape),2));
      score_shape = Skt_score(v, parameters_b, "shape")*(1-tanh(shape)*tanh(shape))/S_shape_t;
      a_shape_t(0) += k_shape*score_shape;
    }
    
  }
  
  return loglik;
  
}

// [[Rcpp::export]]
Rcpp::List filter_reg_list(const Rcpp::List &model){
  
  VectorXd y = model["y"];
  VectorXd a_init = model["a_init"];
  VectorXd a_scale_init = model["a_scale_init"];
  VectorXd a_shape_init = model["a_shape_init"];
  double K = model["K"];
  double k_s = model["k_s"];
  double k_shape = model["k_shape"];
  VectorXd X_loc = model["X_loc"];
  VectorXd X_scale = model["X_scale"];
  VectorXd X_shape = model["X_shape"];
  Trowvec<double> parameters = model["parameters"];
  VectorXd NAs = model["NAs"];
  
  Tvec<double> loglik(y.size());
  loglik.setZero();
  
  Tmat<double> a_loc(a_init.size(), y.size() + 1);
  Tmat<double> a_scale(a_scale_init.size(), y.size() + 1);
  Tmat<double> a_shape(a_shape_init.size(), y.size() + 1);
  
  a_loc.col(0) = a_init;
  a_scale.col(0) = a_scale_init;
  a_shape.col(0) = a_shape_init;
  
  double v;
  
  double S_scale_t;
  double S_loc_t;
  double S_shape_t;
  
  double score_loc;
  double score_scale;
  double score_shape;
  
  double scale_signal;
  double shape_signal;
  
  VectorXd location(y.size());
  location.setZero();
  VectorXd scale(y.size());
  scale.setZero();
  VectorXd shape(y.size());
  shape.setZero();
  
  bool there_is_nan;
  
  for (int t = 0; t < y.size(); t++){
    
    there_is_nan = false;
    if(NAs.size()>1){
      for (int j = 1; j<NAs.size(); j++){
        if (t==NAs[j] - 1){
          there_is_nan = true;
        }
      }
    }
    
    location(t) = a_loc(0, t) + X_loc(t)*a_loc(1, t);
    v = y(t) - location(t);
    scale_signal = a_scale(0, t) + X_scale(t)*a_scale(1, t);
    scale(t) = exp(scale_signal);
    parameters(0) = scale(t);
    shape_signal = a_shape(0, t) + X_shape(t)*a_shape(1, t);
    shape(t) = tanh(shape_signal)*1;
    parameters(3) = shape(t);
    
    if(!there_is_nan){
      
      loglik(t) = log(Skt_f(v,parameters));
      
      S_loc_t = sqrt(Skt_hess(v, parameters, "location"));
      score_loc = Skt_score(v, parameters, "location")/S_loc_t;
      a_loc(0, t+1) = a_loc(0, t) + K*score_loc;
      
      S_scale_t = sqrt(Skt_hess(v, parameters, "scale")*4*exp(4*scale_signal));
      score_scale = Skt_score(v, parameters, "scale")*2*exp(2*scale_signal)/S_scale_t;
      a_scale(0, t+1) = a_scale(0, t) + k_s*score_scale;
      
      S_shape_t = sqrt(Skt_hess(v, parameters, "shape")*pow(1-tanh(shape_signal)*tanh(shape_signal),2));
      score_shape = Skt_score(v, parameters, "shape")*(1-tanh(shape_signal)*tanh(shape_signal))/S_shape_t;
      a_shape(0, t+1) = a_shape(0, t) + k_shape*score_shape;
      
    } else {
      a_loc(0, t + 1) = a_loc(0, t);
      a_scale(0, t + 1) = a_scale(0, t);
      a_shape(0, t + 1) = a_shape(0, t);
    }
    
    a_loc(1, t + 1) = a_loc(1, t);
    a_scale(1, t + 1) = a_scale(1, t);
    a_shape(1, t + 1) = a_shape(1, t);
    
  }
  
  Rcpp::List output = Rcpp::List::create(Rcpp::Named("y") = y,
                                         Rcpp::Named("a_loc") = a_loc,
                                         Rcpp::Named("a_scale") = a_scale,
                                         Rcpp::Named("a_shape") = a_shape,
                                         Rcpp::Named("location") = location,
                                         Rcpp::Named("scale") = scale,
                                         Rcpp::Named("shape") = shape,
                                         Rcpp::Named("model") = model,
                                         Rcpp::Named("loglik") = loglik
  );
  
  return output;
  
}

template<typename T>
T loglik_reg(const Tvec<T> &pars,
             const Tvec<double> &y,
             const bool &stoch_loc,
             const bool &stoch_vol,
             const bool &stoch_shape,
             const Tvec<T> &X_loc,
             const Tvec<T> &X_scale,
             const Tvec<T> &X_shape,
             const Tvec<double> &NAs){
  
  Trowvec<T> parameters(4);
  int n_pars = 0;
  
  Tvec<T> a_init(2);
  a_init.setZero();
  
  Tvec<T> a_scale_init(2);
  a_scale_init.setZero();
  
  Tvec<T> a_shape_init(2);
  a_shape_init.setZero();
  
  Array<double,Eigen::Dynamic,1> y_array(y.size()-(NAs.size()-1));
  int count = 0;
  if(NAs.size()>1){
    bool there_is_no_nan;
    for (int i = 0; i<y.size(); i++){
      
      there_is_no_nan = true;
      for (int j = 1; j<NAs.size(); j++){
        if (i==NAs[j] - 1){
          there_is_no_nan = false;
        }
      }
      
      if(there_is_no_nan){
        y_array(count) = y(i);
        count = count + 1;
      }
    }
  } else {
    for (int i = 0; i<y.size(); i++){
      y_array(i) = y(i);
    }
    count = y.size();
  }
  
  T std_dev = sqrt((y_array - y_array.mean()).square().sum()/(count-1));
  
  parameters(1) = Bound_f(pars(n_pars++)-3,1.1,1e3);
  parameters(2) = parameters(1);
  
  a_init(0) = y_array.mean() + pars(n_pars++);
  a_scale_init(0) = std_dev + pars(n_pars++);
  
  if(stoch_shape){
    a_shape_init(0) = pars(n_pars++);
  }
  
  // Location
  T K = 0;
  
  if (stoch_loc) {
    K = exp(pars(n_pars++))/1e2;
  }
  
  if (X_loc.size() > 1) {
    a_init(1) = pars(n_pars++); 
  }
  
  // Scale
  T k_s = 0;
  if(stoch_vol){
    k_s = exp(pars(n_pars++))/1e2;
  }
  
  if (X_scale.size() > 1) {
    a_scale_init(1) = pars(n_pars++); 
  }
  
  // Shape
  T k_shape = 0;
  if(stoch_shape){
    k_shape = exp(pars(n_pars++))/1e2; 
  }
  
  if (X_shape.size() > 1) {
    a_shape_init(1) = pars(n_pars++);
  }
  
  // MIDAS
  // loc
  Trowvec<T> temporal_agg_loc(4+1);
  temporal_agg_loc.setZero();
  Trowvec<T> theta_loc(4+1);
  theta_loc.setZero();
  Tvec<T> WX_loc(y.rows());
  WX_loc.setZero();
  
  // scale
  Trowvec<T> temporal_agg_scale(4+1);
  temporal_agg_scale.setZero();
  Trowvec<T> theta_scale(4+1);
  theta_scale.setZero();
  Tvec<T> WX_scale(y.rows());
  WX_scale.setZero();
  
  // shape
  Trowvec<T> temporal_agg_shape(4+1);
  temporal_agg_shape.setZero();
  Trowvec<T> theta_shape(4+1);
  theta_shape.setZero();
  Tvec<T> WX_shape(y.rows());
  WX_shape.setZero();
  
  if (X_loc.size() > 1) {
    theta_loc(0) = pars(n_pars++); 
    theta_loc(1) = pars(n_pars++);
  }
  
  if(X_scale.size() > 1){
    theta_scale(0) = pars(n_pars++); 
    theta_scale(1) = pars(n_pars++);
  }
  
  if(X_shape.size() > 1){
    theta_shape(0) = pars(n_pars++); 
    theta_shape(1) = pars(n_pars++);
  }
  
  // agg
  for (int lag = 0; lag<4+1;lag++){
    temporal_agg_loc(lag) = Exp_Almon(lag,4,theta_loc);
    temporal_agg_scale(lag) = Exp_Almon(lag,4,theta_scale);
    temporal_agg_shape(lag) = Exp_Almon(lag,4,theta_shape);
  }
  
  int th = 4;
  for (int t = 0; t<y.rows(); t++){
    for (int lag = 0; lag<4+1;lag++){
      if(X_loc.size() > 1){
        WX_loc(t) += temporal_agg_loc(lag)*X_loc(th-lag);
      }
      if(X_scale.size() > 1){
        WX_scale(t) += temporal_agg_scale(lag)*X_scale(th-lag);
      }
      if(X_shape.size() > 1){
        WX_shape(t) += temporal_agg_shape(lag)*X_shape(th-lag);
      }
    }
    th += 3;
  }
  
  Tvec<T> loglik = filter_reg(y, a_init, a_scale_init, a_shape_init,
                              K, k_s, k_shape,
                              WX_loc, WX_scale, WX_shape,
                              parameters,NAs);
  
  
  T sum_loglik = loglik.sum();
  return  -sum_loglik;
}

// [[Rcpp::export]]
double loglik_reg_rcpp(const Tvec<double> &pars,
                       const Tvec<double> &y,
                       const bool &stoch_loc,
                       const bool &stoch_vol,
                       const bool &stoch_shape,
                       const Tvec<double> &X_loc,
                       const Tvec<double> &X_scale,
                       const Tvec<double> &X_shape,
                       const Tvec<double> &NAs){
  
  return loglik_reg(pars,
                    y,
                    stoch_loc,
                    stoch_vol,
                    stoch_shape,
                    X_loc,
                    X_scale,
                    X_shape,
                    NAs);
}

// [[Rcpp::export]]
Rcpp::List loglik_reg_list(const Tvec<double> &pars,
                           const Rcpp::List &model){
  
  VectorXd y = model["y"];
  Rcpp::List options = model["options"];
  bool stoch_loc = options["stoch_loc"];
  bool stoch_vol = options["stoch_vol"];
  bool stoch_shape = options["stoch_shape"];
  VectorXd X_loc = model["X_loc"];
  VectorXd X_scale = model["X_scale"];
  VectorXd X_shape = model["X_shape"];
  VectorXd NAs = model["NAs"];
  
  Trowvec<double> parameters(4);
  int n_pars = 0;
  
  VectorXd a_init(2);
  a_init.setZero();
  
  VectorXd a_scale_init(2);
  a_scale_init.setZero();
  
  VectorXd a_shape_init(2);
  a_shape_init.setZero();
  
  Array<double,Eigen::Dynamic,1> y_array(y.size()-(NAs.size()-1));
  int count = 0;
  if(NAs.size()>1){
    bool there_is_no_nan;
    for (int i = 0; i<y.size(); i++){
      
      there_is_no_nan = true;
      for (int j = 1; j<NAs.size(); j++){
        if (i==NAs[j] - 1){
          there_is_no_nan = false;
        }
      }
      
      if(there_is_no_nan){
        y_array(count) = y(i);
        count = count + 1;
      }
    }
  } else {
    for (int i = 0; i<y.size(); i++){
      y_array(i) = y(i);
    }
    count = y.size();
  }
  
  double std_dev = sqrt((y_array - y_array.mean()).square().sum()/(count-1));
  
  parameters(1) = Bound_f(pars(n_pars++) - 3, 1.1, 1e3);
  parameters(2) = parameters(1);
  
  a_init(0) = y_array.mean() + pars(n_pars++);
  a_scale_init(0) = std_dev + pars(n_pars++);
  
  if(stoch_shape){
    a_shape_init(0) = pars(n_pars++);
  }
  
  // Location
  double K = 0;
  
  if (stoch_loc) {
    K = exp(pars(n_pars++))/1e2;  
  }
  
  if (X_loc.size() > 1) {
    a_init(1) = pars(n_pars++); 
  }
  
  // Scale
  double k_s = 0;
  if(stoch_vol){
    k_s = exp(pars(n_pars++))/1e2; 
  }
  
  if (X_scale.size() > 1) {
    a_scale_init(1) = pars(n_pars++);
  }
  
  // Shape
  double k_shape = 0;
  if(stoch_shape){
    k_shape = exp(pars(n_pars++))/1e2; 
  }
  
  if (X_shape.size() > 1) {
    a_shape_init(1) = pars(n_pars++);
  }
  
  // MIDAS
  // loc
  Trowvec<double> temporal_agg_loc(4+1);
  temporal_agg_loc.setZero();
  Trowvec<double> theta_loc(4+1);
  theta_loc.setZero();
  VectorXd WX_loc(y.rows());
  WX_loc.setZero();
  
  // scale
  Trowvec<double> temporal_agg_scale(4+1);
  temporal_agg_scale.setZero();
  Trowvec<double> theta_scale(4+1);
  theta_scale.setZero();
  VectorXd WX_scale(y.rows());
  WX_scale.setZero();
  
  // shape
  Trowvec<double> temporal_agg_shape(4+1);
  temporal_agg_shape.setZero();
  Trowvec<double> theta_shape(4+1);
  theta_shape.setZero();
  VectorXd WX_shape(y.rows());
  WX_shape.setZero();
  
  if (X_loc.size() > 1) {
    theta_loc(0) = pars(n_pars++); 
    theta_loc(1) = pars(n_pars++);
  }
  
  if(X_scale.size() > 1){
    theta_scale(0) = pars(n_pars++); 
    theta_scale(1) = pars(n_pars++);
  }
  if(X_shape.size() > 1){
    theta_shape(0) = pars(n_pars++); 
    theta_shape(1) = pars(n_pars++);
  }
  // agg
  for (int lag = 0; lag<4+1;lag++){
    temporal_agg_loc(lag) = Exp_Almon(lag,4,theta_loc);
    temporal_agg_scale(lag) = Exp_Almon(lag,4,theta_scale);
    temporal_agg_shape(lag) = Exp_Almon(lag,4,theta_shape);
  }
  
  int th = 4;
  for (int t = 0; t < y.size(); t++){
    for (int lag = 0; lag<4+1;lag++){
      if(X_loc.size() > 1){
        WX_loc(t) += temporal_agg_loc(lag)*X_loc(th-lag);
      }
      if(X_scale.size() > 1){
        WX_scale(t) += temporal_agg_scale(lag)*X_scale(th-lag);
      }
      if(X_shape.size() > 1){
        WX_shape(t) += temporal_agg_shape(lag)*X_shape(th-lag);
      }
    }
    th += 3;
  }
  
  Rcpp::List model_new = Rcpp::List::create(Rcpp::Named("y") = y,
                                            Rcpp::Named("a_init") = a_init,
                                            Rcpp::Named("a_scale_init") = a_scale_init,
                                            Rcpp::Named("a_shape_init") = a_shape_init,
                                            Rcpp::Named("K") = K,
                                            Rcpp::Named("k_s") = k_s,
                                            Rcpp::Named("k_shape") = k_shape,
                                            Rcpp::Named("X_loc") = WX_loc,
                                            Rcpp::Named("X_scale") = WX_scale,
                                            Rcpp::Named("X_shape") = WX_shape,
                                            Rcpp::Named("temporal_agg_loc") = temporal_agg_loc,
                                            Rcpp::Named("temporal_agg_scale") = temporal_agg_scale,
                                            Rcpp::Named("temporal_agg_shape") = temporal_agg_shape,
                                            Rcpp::Named("parameters") = parameters,
                                            Rcpp::Named("NAs") = NAs,
                                            Rcpp::Named("n_pars") = n_pars);
  
  Rcpp::List output = filter_reg_list(model_new);
  
  return output;
}

// [[Rcpp::export]]
VectorXd loglik_reg_deriv(const VectorXd &pars,
                          const VectorXd &y,
                          const bool &stoch_loc,
                          const bool &stoch_vol,
                          const bool &stoch_shape,
                          const VectorXd &X_loc,
                          const VectorXd &X_scale,
                          const VectorXd &X_shape,
                          const VectorXd &NAs)
{
  
  
  //Creating AD model matrices and vectors
  ////////////////////////////////////////////////////////
  Tvec<a_double> apars(pars.rows());
  Tvec<a_double> aX_loc(X_loc.rows());
  Tvec<a_double> aX_scale(X_scale.rows());
  Tvec<a_double> aX_shape(X_shape.rows());
  
  for (int i = 0; i < pars.rows(); ++i){
    apars(i) = pars(i);
  }
  
  for (int row = 0; row < X_loc.rows(); ++row){
    aX_loc(row) = X_loc(row);
  }
  
  for (int row = 0; row < X_scale.rows(); ++row){
    aX_scale(row) = X_scale(row);
  }
  
  for (int row = 0; row < X_shape.rows(); ++row){
    aX_shape(row) = X_shape(row);
  }
    
  a_vector LogLik(1);
  ADFun<double> grad;
  Independent(apars);
  LogLik(0) = loglik_reg(apars, y, stoch_loc, stoch_vol, stoch_shape, aX_loc, aX_scale, aX_shape,NAs);
  grad.Dependent(apars, LogLik);
  return grad.Jacobian(pars);
}

// [[Rcpp::export]]
MatrixXd loglik_reg_hess(const VectorXd &pars,
                         const VectorXd &y,
                         const bool &stoch_loc,
                         const bool &stoch_vol,
                         const bool &stoch_shape,
                         const VectorXd &X_loc,
                         const VectorXd &X_scale,
                         const VectorXd &X_shape,
                         const VectorXd &NAs)
{
  
  
  //Creating AD model matrices and vectors
  ////////////////////////////////////////////////////////
  Tvec<a_double> apars(pars.rows());
  Tvec<a_double> aX_loc(X_loc.rows());
  Tvec<a_double> aX_scale(X_scale.rows());
  Tvec<a_double> aX_shape(X_shape.rows());
  
  for (int i = 0; i < pars.rows(); ++i){
    apars(i) = pars(i);
  }
  
  for (int row = 0; row < X_loc.rows(); ++row){
    aX_loc(row) = X_loc(row);
  }
  
  for (int row = 0; row < X_scale.rows(); ++row){
    aX_scale(row) = X_scale(row);
  }
  
  for (int row = 0; row < X_shape.rows(); ++row){
    aX_shape(row) = X_shape(row);
  }
  
  ////////////////////////////////////////////////////////
  
  a_vector LogLik(1);
  ADFun<double> grad;
  Independent(apars);
  VectorXd vector_hessian(apars.size()*apars.size());
  
  LogLik(0) = loglik_reg(apars, y, stoch_loc, stoch_vol, stoch_shape, aX_loc, aX_scale, aX_shape,NAs);
  grad.Dependent(apars, LogLik);
  vector_hessian = grad.Hessian(pars,0);
  
  MatrixXd hessian(apars.size(),apars.size());
  int cursor = 0;
  for(int row = 0; row < apars.size(); row++)
  {
    for(int col = 0; col < apars.size(); col++)
    {
      hessian(row,col) = vector_hessian[cursor++]; 
    } 
  }
  
  return hessian;
}

