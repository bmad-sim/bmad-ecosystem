#include <vector>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "cauchy.hpp"



const double PI = 4.0 * atan(1.0);

template<typename T>
inline constexpr T sqr(T x) { return x*x; }

inline CauchyPoint gsl_to_cauchy_point(const gsl_vector *fit_params) {
  return {0.0, 0.0, // E,r don't matter to fitting routine
    gsl_vector_get(fit_params, 0), // cx
    std::abs(gsl_vector_get(fit_params, 1)), // ax
    std::abs(gsl_vector_get(fit_params, 2)), // ay
    gsl_vector_get(fit_params, 3), // beta
    0.0, 0.0, 0.0, //dxds/dyds_min/max
    std::abs(gsl_vector_get(fit_params, 4)), // amp
    0.0, CauchyStatus::OK}; // chi2, npts, stat
}

int asym_cauchy_fit_function(const gsl_vector *fit_params, void *data_points,
    gsl_vector *fit_residuals) {
  // Computes fit function minus bin count for each
  // bin passed in through data_points, and stores the
  // results in fit_residuals
  std::vector<BinPoint>& bins = * (std::vector<BinPoint> *) data_points;

  auto p = gsl_to_cauchy_point(fit_params);

  size_t num_pts = bins.size();
  for (size_t i=0; i<num_pts; i++) {
    auto [xval, yval, binval] = bins[i];
    double fit_val = p.amp * (1 + p.beta * xval) / (1 + p.ax * sqr(xval-p.cx) + p.ay * sqr(yval));
    gsl_vector_set(fit_residuals, i, fit_val - binval);
  }

  return GSL_SUCCESS;
}

CauchyPoint asym_cauchy_fit(const XYBinnedData& data) {
  // Default guess
  CauchyPoint guess;
  guess.cx = 0.0;
  guess.ax = 1.0;
  guess.ay = 1.0;
  guess.beta = 1.0;
  guess.amp = 1.0;
  return asym_cauchy_fit(data, guess);
}

CauchyPoint asym_cauchy_fit(const XYBinnedData& data, const CauchyPoint& guess) {
  // This function fits a normal distribution to the binned
  // data supplied in bins using the GSL's nonlinear
  // fitting routine.  In case the fitting process fails,
  // all -1's are returned for the fit results

  // GSL boilerplate
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
  const size_t num_pts = data.bins.size();
  const size_t num_fit_params = 5; // c_x, a_x, a_y, beta, amp

  gsl_vector *fit_residuals;
  double initial_parameter_guess[] = { guess.cx, guess.ax, guess.ay, guess.beta, guess.amp };
  gsl_vector_view init_param_gsl = gsl_vector_view_array(initial_parameter_guess, num_fit_params);

  double chisq;
  int status, info;

  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;

  // Set up the fdf struct
  fdf.f = asym_cauchy_fit_function;
  fdf.df = nullptr;
  fdf.fvv = nullptr;
  fdf.n = num_pts;
  fdf.p = num_fit_params;
  fdf.params = (void *) &(data.bins);

  if (num_pts < num_fit_params) {
    std::cerr << "ERROR: cannot perform fit with less than " << num_fit_params << " data points!\n";
    std::cerr << "Problem: only " << num_pts << " bins for pc = "
      << data.E/1e6 << " MeV, r = " << data.r*1e2 << " cm\n";
    std::exit(1);
  }

  // Allocate fitting workspace
  w = gsl_multifit_nlinear_alloc(T, &fdf_params, num_pts, num_fit_params);

  // Initialize solver
  gsl_multifit_nlinear_init(&init_param_gsl.vector, &fdf, w);

  // Solve the system
  size_t max_iter = 1000;
  status = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol, nullptr, nullptr, &info, w);

  // Compute final residuals
  fit_residuals = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(fit_residuals, fit_residuals, &chisq);

  // Return results
  auto result = gsl_to_cauchy_point(w->x);
  switch (status) {
    case GSL_EMAXITER:
      //std::cerr << "Iteration limit reached for cauchy fit at pc_out = "
      //  << pc_out*1e-6 << " MeV, r = " << r*1e2 << " cm\n";
      result.stat = CauchyStatus::MAXITER;
      break;
    case GSL_ENOPROG:
      result.stat = CauchyStatus::NOPROG;
      break;
    default:
      result.stat = CauchyStatus::OK;
  }

  // Fill in the rest of the parameters to result
  result.ax = std::sqrt(result.ax);
  result.ay = std::sqrt(result.ay);
  result.E = data.E;
  result.r = data.r;
  // bins are sorted by increasing dxds, then increasing dyds
  result.dxds_min = data.bins.front().x;
  result.dxds_max = data.bins.back().x;
  result.dyds_max = data.bins.back().y;
  result.chi2 = sqrt(chisq)/(num_pts - num_fit_params);

  gsl_multifit_nlinear_free(w);
  return result;
}

