#include <vector>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "normal_fit.hpp"

const double PI = 4.0 * atan(1.0);

template<typename T>
inline constexpr T sqr(T x) { return x*x; }

double normal2d(double x, double y, double mu_x, double mu_y, double s_x, double s_y) {
  s_x = std::abs(s_x);
  s_y = std::abs(s_y);

  double xterm = (x-mu_x)/s_x;
  double yterm = (y-mu_y)/s_y;
  return (1.0 / (2*PI*s_x*s_y)) * std::exp(-0.5 * (sqr(xterm) + sqr(yterm)));
}


int fit_function(const gsl_vector *fit_params, void *data_points,
    gsl_vector *fit_residuals) {
  // Computes fit function minus bin count for each
  // bin passed in through data_points, and stores the
  // results in fit_residuals
  std::vector<BinPoint>& bins = * (std::vector<BinPoint> *) data_points;

  double mu_x = gsl_vector_get(fit_params, 0);
  //double mu_y = gsl_vector_get(fit_params, 1);
  double mu_y = 0.0;
  double s_x = gsl_vector_get(fit_params, 1);
  double s_y = gsl_vector_get(fit_params, 2);
  double amp = gsl_vector_get(fit_params, 3);

  if (s_x==0.0 || s_y==0.0) return GSL_EDOM;

  size_t num_pts = bins.size();
  for (size_t i=0; i<num_pts; i++) {
    auto [xval, yval, binval] = bins[i];
    double fit_val = amp * normal2d(xval, yval, mu_x, mu_y, s_x, s_y);
    gsl_vector_set(fit_residuals, i, fit_val - binval);
  }

  return GSL_SUCCESS;
}


int cauchy_fit_function(const gsl_vector *fit_params, void *data_points,
    gsl_vector *fit_residuals) {
  // Computes fit function minus bin count for each
  // bin passed in through data_points, and stores the
  // results in fit_residuals
  std::vector<BinPoint>& bins = * (std::vector<BinPoint> *) data_points;

  double mu_x = gsl_vector_get(fit_params, 0);
  //double mu_y = gsl_vector_get(fit_params, 1);
  double mu_y = 0.0;
  double a_x = gsl_vector_get(fit_params, 1);
  double a_y = gsl_vector_get(fit_params, 2);
  double amp = gsl_vector_get(fit_params, 3);

  //if (s_x==0.0 || s_y==0.0) return GSL_EDOM;

  size_t num_pts = bins.size();
  for (size_t i=0; i<num_pts; i++) {
    auto [xval, yval, binval] = bins[i];
    double fit_val = amp / ( 1 + a_x * sqr(xval-mu_x) + a_y * sqr(yval-mu_y) );
    gsl_vector_set(fit_residuals, i, fit_val - binval);
  }

  return GSL_SUCCESS;
}


int fit_function_jac(const gsl_vector *fit_params, void *data_points,
    gsl_matrix *jac) {
  // Computes the jacobian for the normal distribution at the given
  // fit parameters for each data_point, and writes the results
  // to jac
  std::vector<BinPoint>& bins = * (std::vector<BinPoint> *) data_points;

  double mu_x = gsl_vector_get(fit_params, 0);
  //double mu_y = gsl_vector_get(fit_params, 1);
  double mu_y = 0.0;
  double s_x = gsl_vector_get(fit_params, 1);
  double s_y = gsl_vector_get(fit_params, 2);
  double amp = gsl_vector_get(fit_params, 3);

  if (s_x==0.0 || s_y==0.0) return GSL_EDOM;

  size_t num_pts = bins.size();
  for (size_t i=0; i<num_pts; i++) {
    auto [xval, yval, binval] = bins[i];
    double xterm = (xval - mu_x) / s_x;
    double yterm = (yval - mu_y) / s_y;

    gsl_matrix_set(jac, i, 0, amp * normal2d(xval, yval, mu_x, mu_y, s_x, s_y) * xterm / s_x);
    //gsl_matrix_set(jac, i, 1, amp * normal2d(xval, yval, mu_x, mu_y, s_x, s_y) * yterm / s_y);
    gsl_matrix_set(jac, i, 1, amp * normal2d(xval, yval, mu_x, mu_y, s_x, s_y) * sqr(xterm) / s_x);
    gsl_matrix_set(jac, i, 2, amp * normal2d(xval, yval, mu_x, mu_y, s_x, s_y) * sqr(yterm) / s_y);
    gsl_matrix_set(jac, i, 3, normal2d(xval, yval, mu_x, mu_y, s_x, s_y));
  }

  return GSL_SUCCESS;
}


void print_status(const size_t iter, void *,
    const gsl_multifit_nlinear_workspace *w) {
  // Prints the parameter values at each step of the iteration
  gsl_vector *fit_params = gsl_multifit_nlinear_position(w);
  gsl_vector *f = gsl_multifit_nlinear_residual(w);

  std::cout << "Iter: " << iter
    << "\tmu_x = " << gsl_vector_get(fit_params,0)
    //<< "\tmu_y = " << gsl_vector_get(fit_params,1)
    << "\ts_x = " << gsl_vector_get(fit_params,1) << '\n'
    << "\ts_y = " << gsl_vector_get(fit_params,2)
    << "\tamplitude = " << gsl_vector_get(fit_params,3)
    << "\tChi^2 = " << gsl_blas_dnrm2(f) << '\n';
  getchar();

  return;
}


FitResults normal_fit(const std::vector<BinPoint>& bins) {
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
  const size_t num_fit_params = 4; //mu_x, mu_y, s_x, s_y, amp

  gsl_vector *fit_residuals;
  //gsl_matrix *jacobian;
  //gsl_matrix *covar = gsl_matrix_alloc(num_fit_params, num_fit_params);
  std::vector<double> weights(num_pts);
  //double initial_parameter_guess[] = { 0.0, 0.0, 1.0, 1.0, 1.0 };
  double initial_parameter_guess[] = { 0.0, 1.0, 1.0, 100.0 };
  gsl_vector_view init_param_gsl = gsl_vector_view_array(initial_parameter_guess, num_fit_params);
  gsl_vector_view weights_gsl = gsl_vector_view_array(&weights[0], num_pts);

  double chisq;//, chisq0;
  int status, info;

  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;

  // Set up the fdf struct
  fdf.f = cauchy_fit_function;
  //fdf.df = fit_function_jac;
  fdf.df = nullptr;
  fdf.fvv = nullptr;
  fdf.n = num_pts;
  fdf.p = num_fit_params;
  fdf.params = (void *) &bins;

  // Compute weights
  size_t bin_ix;
  for (bin_ix=0; bin_ix<num_pts; bin_ix++) {
    weights[bin_ix] = 1/std::max(1.0, bins[bin_ix].count);
  }

  // Allocate fitting workspace
  w = gsl_multifit_nlinear_alloc(T, &fdf_params, num_pts, num_fit_params);

  // Initialize solver
  gsl_multifit_nlinear_winit(&init_param_gsl.vector, &weights_gsl.vector, &fdf, w);

  // Initial residuals
  //gsl_blas_ddot(fit_residuals, fit_residuals, &chisq0);

  // Solve the system
  size_t max_iter = 100;
  status = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol, print_status, nullptr, &info, w);

  // Compute final covariance
  //jacobian = gsl_multifit_nlinear_jac(w);
  //gsl_multifit_nlinear_covar(jacobian, 0.0, covar);

  // Compute final residuals
  fit_residuals = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(fit_residuals, fit_residuals, &chisq);

  // Return results
  if (status != GSL_SUCCESS) return {-1, -1, -1, -1, -1, -1};

  double fit_mu_x = gsl_vector_get(w->x, 0);
  //double fit_mu_y = gsl_vector_get(w->x, 1);
  double fit_mu_y = 0.0;
  double fit_s_x = gsl_vector_get(w->x, 1);
  double fit_s_y = gsl_vector_get(w->x, 2);
  double fit_amp = gsl_vector_get(w->x, 3);
  size_t dof = num_pts - num_fit_params;

  gsl_multifit_nlinear_free(w);
  return {fit_mu_x, fit_mu_y, fit_s_x, fit_s_y, fit_amp, chisq/dof};
}


