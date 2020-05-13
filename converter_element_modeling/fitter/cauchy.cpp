#include <vector>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "cauchy.hpp"



const double PI = 4.0 * atan(1.0);

template<typename T>
inline constexpr T sqr(T x) { return x*x; }

int asym_cauchy_fit_function(const gsl_vector *fit_params, void *data_points,
    gsl_vector *fit_residuals) {
  // Computes fit function minus bin count for each
  // bin passed in through data_points, and stores the
  // results in fit_residuals
  std::vector<BinPoint>& bins = * (std::vector<BinPoint> *) data_points;

  double c_x = gsl_vector_get(fit_params, 0);
  double c_y = gsl_vector_get(fit_params, 1);
  double b_x = gsl_vector_get(fit_params, 4);
  double a_x = std::abs(gsl_vector_get(fit_params, 2));
  double a_y = std::abs(gsl_vector_get(fit_params, 3));
  double amp = std::abs(gsl_vector_get(fit_params, 5));

  size_t num_pts = bins.size();
  for (size_t i=0; i<num_pts; i++) {
    auto [xval, yval, binval] = bins[i];
    double fit_val = amp * (1 + b_x * xval) / (1 + a_x * sqr(xval-c_x) + a_y * sqr(yval-c_y));
    gsl_vector_set(fit_residuals, i, fit_val - binval);
  }

  return GSL_SUCCESS;
}


cauchyFitResults asym_cauchy_fit(double pc_out, double r, const std::vector<BinPoint>& bins) {
  // This function fits a normal distribution to the binned
  // data supplied in bins using the GSL's nonlinear
  // fitting routine.  In case the fitting process fails,
  // all -1's are returned for the fit results

  // GSL boilerplate
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
  const size_t num_pts = bins.size();
  const size_t num_fit_params = 6; //mu_x, mu_y, s_x, s_y, amp

  gsl_vector *fit_residuals;
  //gsl_matrix *jacobian;
  //gsl_matrix *covar = gsl_matrix_alloc(num_fit_params, num_fit_params);
  double initial_parameter_guess[] = { 0.0, 0.0, 1.0, 1.0, 1.0, 100.0 };
  //double initial_parameter_guess[] = { 0.0, 1.0, 1.0, 1.0 };
  gsl_vector_view init_param_gsl = gsl_vector_view_array(initial_parameter_guess, num_fit_params);

  double chisq;//, chisq0;
  int status, info;

  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;

  // Set up the fdf struct
  fdf.f = asym_cauchy_fit_function;
  //fdf.df = fit_function_jac;
  fdf.df = nullptr;
  fdf.fvv = nullptr;
  fdf.n = num_pts;
  fdf.p = num_fit_params;
  fdf.params = (void *) &bins;

  // Compute weights

  // Allocate fitting workspace
  w = gsl_multifit_nlinear_alloc(T, &fdf_params, num_pts, num_fit_params);

  // Initialize solver
  gsl_multifit_nlinear_init(&init_param_gsl.vector, &fdf, w);

  // Initial residuals
  //gsl_blas_ddot(fit_residuals, fit_residuals, &chisq0);

  // Solve the system
  size_t max_iter = 100;
  status = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol, nullptr, nullptr, &info, w);

  // Compute final covariance
  //jacobian = gsl_multifit_nlinear_jac(w);
  //gsl_multifit_nlinear_covar(jacobian, 0.0, covar);

  // Compute final residuals
  fit_residuals = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(fit_residuals, fit_residuals, &chisq);

  // Return results
  if (status == GSL_EMAXITER)
    std::cerr << "Iteration limit reached for cauchy fit at pc_out = "
      << pc_out << " MeV, r = " << r << " cm\n";
  if (status == GSL_ENOPROG) return {-1, -1, -1, -1, -1, -1, -1};

  double fit_c_x = gsl_vector_get(w->x, 0);
  double fit_c_y = gsl_vector_get(w->x, 1);
  double fit_b_x = gsl_vector_get(w->x, 4);
  double fit_a_x = std::abs(gsl_vector_get(w->x, 2));
  double fit_a_y = std::abs(gsl_vector_get(w->x, 3));
  double fit_amp = std::abs(gsl_vector_get(w->x, 5));
  size_t dof = num_pts - num_fit_params;

  gsl_multifit_nlinear_free(w);
  return {fit_c_x, fit_c_y, fit_a_x, fit_a_y,
          fit_b_x, fit_amp, sqrt(chisq)/dof};
}

