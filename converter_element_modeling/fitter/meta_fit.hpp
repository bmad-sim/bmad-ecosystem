#pragma once
#include <type_traits>
#include <vector>
#include <iostream>
#include <string>
#include <gsl/gsl_poly.h>
#include <cmath>

struct DataPoint {
  double E, r;
  double cx, ax, ay, beta, dxds_min, dxds_max, dyds_max;
  double amp;
  double chi2;
  size_t npts;
};

// Fit types:
enum class fitType { CX, AX, AY, BETA, DXDS_MIN, DXDS_MAX, DYDS_MAX };
std::ostream& operator<<(std::ostream&, fitType);
std::string fit_to_string(fitType);

struct FitResults1D {
  double E, a, b, c, d;
  double chi2;
};

struct FitResults2D {
  double a1, a2, a3;
  double b0, b1, b2, b3;
  double ke, kr, C;
  double chi2;
};

struct FitResults {
  std::vector<FitResults1D> low_e_fits;
  FitResults2D high_e_fit;
  FitResults() = default;
};


template <fitType T>
FitResults fit_routine(const std::vector<DataPoint>&, double);


struct MetaFitResults {
  double Ein;
  double T;
  FitResults cx, ax, ay, beta, dxds_min, dxds_max, dyds_max;
};


template<fitType T, size_t dim, typename F>
double eval(const F& fit, const DataPoint& p) {
  // Evaluates the given fit at the given (E,r) point
  // and returns the result
  if constexpr (dim==2) {
    static_assert(std::is_same_v<F, FitResults2D>);
    double pc_coefs[] = {1.0, fit.a1, fit.a2, fit.a3};
    double r_coefs[] = {fit.b0, fit.b1, fit.b2, fit.b3};
    double ke = fit.ke, kr = fit.kr, C = fit.C;
    return gsl_poly_eval(pc_coefs, 4, p.E)
      * gsl_poly_eval(r_coefs, 4, p.r)
      * std::exp(-(ke * p.E + kr * p.r)) + C;
  } else {
    static_assert(std::is_same_v<F, FitResults1D>);
    double r_coefs[] = {fit.a, fit.b, fit.c, fit.d};
    if constexpr(T==fitType::CX || T==fitType::BETA)
      return gsl_poly_eval(r_coefs, 4, p.r) * p.r;
    else
      return gsl_poly_eval(r_coefs, 4, p.r);
  }
}
