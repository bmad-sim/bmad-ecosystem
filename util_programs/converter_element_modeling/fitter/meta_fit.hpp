#pragma once
#include <type_traits>
#include <vector>
#include <iostream>
#include <string>
#include <gsl/gsl_poly.h>
#include <cmath>
#include "cauchy.hpp"
#include "read_data.hpp"

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
FitResults fit_routine(const std::vector<CauchyPoint>& cauchy, const TableData& er_table, double crossover_point);


struct MetaFitResults {
  double Ein;
  double T;
  FitResults cx, ax, ay, beta, dxds_min, dxds_max, dyds_max;
  TableData er_table, polz_table;

  MetaFitResults(double Ein, double T, FitResults&& cx, FitResults&& ax, FitResults&& ay,
      FitResults&& beta, FitResults&& dxds_min, FitResults&& dxds_max,
      FitResults&& dyds_max, TableData&& er_table, TableData&& polz_table);
};


template<fitType T, size_t dim, typename F>
double eval(const F& fit, const CauchyPoint& p) {
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

template<size_t n>
std::string TAB() {
  std::string tab = "  ";
  if constexpr (n==0) return "";
  else return tab + TAB<n-1>();
}

enum Comma_spec { COMMA_AFTER, NO_COMMA };
template<fitType T>
void output_bmad(const FitResults& fit, std::ostream& bmad_file, Comma_spec comma) {
  bmad_file << TAB<3>() << fit_to_string(T) << " = {\n";
  for (const auto& fit_1d : fit.low_e_fits) {
    bmad_file << TAB<4>() << "fit_1d_r = {pc_out = " << fit_1d.E
      << ", poly = [";
      if constexpr (T==fitType::CX || T==fitType::BETA)
        bmad_file << 0.0 << ", ";
      bmad_file << fit_1d.a << ", "
      << fit_1d.b << ", "
      << fit_1d.c << ", "
      << fit_1d.d << "]},\n";
  }
  const auto& fit_2d = fit.high_e_fit;
  if constexpr (T==fitType::DXDS_MIN)
    bmad_file << TAB<4>() << "C = " << fit_2d.C << ",\n";
  bmad_file << TAB<4>() << "fit_2d_pc = {k = " << fit_2d.ke << ", poly = [1.0, "
    << fit_2d.a1 << ", "
    << fit_2d.a2 << ", "
    << fit_2d.a3 << "]},\n";
  bmad_file << TAB<4>() << "fit_2d_r = {k = " << fit_2d.kr << ", poly = ["
    << fit_2d.b0 << ", "
    << fit_2d.b1 << ", "
    << fit_2d.b2 << ", "
    << fit_2d.b3 << "]}}";
  switch (comma) {
    case COMMA_AFTER: bmad_file << ",\n"; break;
    case NO_COMMA: bmad_file << "\n"; break;
  }
  return;
}
