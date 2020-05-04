#pragma once
#include <vector>
#include <variant>
#include <cmath>
#include <gsl/gsl_poly.h>
struct DataPoint {
  double E;
  double r;
  double cx;
  double cy;
  double ax;
  double ay;
  double beta;
  double amp;
};

// Fit types:
enum class fitType { CX, AX, AY, BETA };

//////////////////CX////////////////////

struct cFitResults {
  double a1, a2, a3;
  double b1, b2, b3;
  double chi2;
};

//////////////////////AX/Y////////////////////
struct aFitResults2D {
  double ke, kr;
  double a1, a2, a3;
  double b0, b1, b2, b3;
  double chi2;
};

struct aFitResults1D {
  double E;
  double k;
  double a, b, c, d;
  double chi2;
};

struct aFitResults {
  std::vector<aFitResults1D> low_e_fits;
  aFitResults2D high_e_fit;
};

////////////////////////BETA/////////////////////

struct betaFitResults2D {
  double a1, a2, a3;
  double b1, b2, b3;
  double chi2;
};

struct betaFitResults1D {
  double E;
  double a, b, c, d;
  double chi2;
};

struct betaFitResults {
  std::vector<betaFitResults1D> low_e_fits;
  betaFitResults2D high_e_fit;
};


struct MetaFitResults {
  double Ein;
  double T;
  cFitResults cx;
  aFitResults ax;
  aFitResults ay;
  betaFitResults beta;
};


cFitResults cx_fit(const std::vector<DataPoint>& data_points);
aFitResults ax_fit(const std::vector<DataPoint>& data_points, double crossover_point);
aFitResults ay_fit(const std::vector<DataPoint>& data_points, double crossover_point);
betaFitResults beta_fit(const std::vector<DataPoint>& data_points, double crossover_point);

template<fitType T, size_t dim>
struct fit_t_helper {
  typedef void type;
};

template<> struct fit_t_helper<fitType::CX,0>   { typedef cFitResults      type; };
template<> struct fit_t_helper<fitType::CX,2>   { typedef cFitResults      type; };
template<> struct fit_t_helper<fitType::AX,0>   { typedef aFitResults      type; };
template<> struct fit_t_helper<fitType::AX,1>   { typedef aFitResults1D    type; };
template<> struct fit_t_helper<fitType::AX,2>   { typedef aFitResults2D    type; };
template<> struct fit_t_helper<fitType::AY,0>   { typedef aFitResults      type; };
template<> struct fit_t_helper<fitType::AY,1>   { typedef aFitResults1D    type; };
template<> struct fit_t_helper<fitType::AY,2>   { typedef aFitResults2D    type; };
template<> struct fit_t_helper<fitType::BETA,0> { typedef betaFitResults   type; };
template<> struct fit_t_helper<fitType::BETA,1> { typedef betaFitResults1D type; };
template<> struct fit_t_helper<fitType::BETA,2> { typedef betaFitResults2D type; };

template<fitType T, size_t dim>
using fit_t_part = typename fit_t_helper<T,dim>::type;
template<fitType T>
using fit_t = typename fit_t_helper<T,0>::type;

// Helper templates
template<fitType T, size_t d>
constexpr bool good_fit_type() {
  if constexpr (d==2)
    return (T==fitType::CX || T==fitType::AX || T==fitType::AY || T==fitType::BETA);
  else if constexpr (d==1)
    return (T==fitType::AX || T==fitType::AY || T==fitType::BETA);
  else
    return false;
}

template<fitType T, size_t dim>
double eval(const fit_t_part<T, dim>& fit, const DataPoint& p) {
  // Evaluates the given fit at the given (E,r) point
  // and returns the result
  static_assert(good_fit_type<T,dim>());
  if constexpr (dim==2) {
    double pc_coefs[] = {1.0, fit.a1, fit.a2, fit.a3};
    double r_coefs[] = {0.0, fit.b1, fit.b2, fit.b3};
    double ke = 0.0, kr = 0.0;
    if constexpr(T==fitType::AX || T==fitType::AY) {
      r_coefs[0] = fit.b0;
      ke = fit.ke;
      kr = fit.kr;
    }
    return gsl_poly_eval(pc_coefs, 4, p.E)
      * gsl_poly_eval(r_coefs, 4, p.r)
      * std::exp(-(ke * p.E + kr * p.r));
  } else {
    double r_coefs[] = {fit.a, fit.b, fit.c, fit.d};
    if constexpr(T==fitType::BETA)
      return gsl_poly_eval(r_coefs, 4, p.r) * p.r;
    else
      return gsl_poly_eval(r_coefs, 4, p.r) * std::exp(-fit.k*p.r);
  }
}

