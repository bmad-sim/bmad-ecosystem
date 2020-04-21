#include <vector>
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <variant>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_poly.h>

#include "meta_fit.hpp"

#define ABS_PARAMS

const double PI = 4.0 * atan(1.0);

template<typename T>
inline constexpr T sqr(T x) { return x*x; }
template<typename T>
inline constexpr T cube(T x) { return x*x*x; }


//////////////////////////// 2D FIT FUNCTIONS ////////////////////////////////
template<fitType T>
int fit_function_2d(const gsl_vector *fit_params, void *data_points,
    gsl_vector *fit_residuals) {

  static_assert(T==fitType::CX || T==fitType::AX ||
                T==fitType::AY || T==fitType::BETA,
                "Error: unknown fit type selected");
  // Cast data_points to the appropriate type
  std::vector<DataPoint>& data = * (std::vector<DataPoint> *) data_points;

  // Extract fit coefficients from fit_params
  [[maybe_unused]] double ke, kr;
  std::array<double, 4> E_coefs, r_coefs;
  if constexpr (T==fitType::CX || T==fitType::BETA) {
    double a0 = gsl_vector_get(fit_params, 0);
    double a1 = gsl_vector_get(fit_params, 1);
    double a2 = gsl_vector_get(fit_params, 2);
    double b0 = gsl_vector_get(fit_params, 3);
    double b1 = gsl_vector_get(fit_params, 4);
    double b2 = gsl_vector_get(fit_params, 5);
    double b3 = gsl_vector_get(fit_params, 6);
    E_coefs = {a0, a1, a2, 1.0};
    r_coefs = {b0, b1, b2, b3};
  } else {
    ke = std::abs(gsl_vector_get(fit_params, 0));
    kr = std::abs(gsl_vector_get(fit_params, 1));
    double ae = gsl_vector_get(fit_params, 2);
    double be = gsl_vector_get(fit_params, 3);
    double ce = gsl_vector_get(fit_params, 4);
    //double de = gsl_vector_get(fit_params, 5);
    double ar = gsl_vector_get(fit_params, 5);
    double br = gsl_vector_get(fit_params, 6);
    double cr = gsl_vector_get(fit_params, 7);
    double dr = gsl_vector_get(fit_params, 8);
    E_coefs = {ae, be, ce, 1.0};
    r_coefs = {ar, br, cr, dr};
  }

  // Main loop: compute fit residual at each point
  size_t num_pts = data.size();
  for (size_t i=0; i<num_pts; i++) {
    auto E = data[i].E;
    auto r = data[i].r;
    // Select the correct variable to fit to
    double real_val;
    if constexpr (T==fitType::CX) real_val = data[i].cx;
    else if constexpr (T==fitType::AX) real_val = data[i].ax;
    else if constexpr (T==fitType::AY) real_val = data[i].ay;
    else if constexpr (T==fitType::BETA) real_val = data[i].beta;
    // Compute the correct fit value
    double fit_val = gsl_poly_eval(E_coefs.data(), 4, E)
      * gsl_poly_eval(r_coefs.data(), 4, r);
    if constexpr (T==fitType::AX || T==fitType::AY)
      fit_val *= std::exp(-(ke*E+kr*r));
    gsl_vector_set(fit_residuals, i, fit_val - real_val);
  }

  return GSL_SUCCESS;
}

// Specializations
//int cx_fit_function(const gsl_vector *fit_params, void *data_points,
//    gsl_vector *fit_residuals) {
//  return fit_function_2d<fitType::CX>(fit_params, data_points, fit_residuals);
//}
//int ax_fit_function_2d(const gsl_vector *fit_params, void *data_points,
//    gsl_vector *fit_residuals) {
//  return fit_function_2d<fitType::AX>(fit_params, data_points, fit_residuals);
//}
//int ay_fit_function_2d(const gsl_vector *fit_params, void *data_points,
//    gsl_vector *fit_residuals) {
//  return fit_function_2d<fitType::AY>(fit_params, data_points, fit_residuals);
//}
//int beta_fit_function_2d(const gsl_vector *fit_params, void *data_points,
//    gsl_vector *fit_residuals) {
//  return fit_function_2d<fitType::BETA>(fit_params, data_points, fit_residuals);
//}


/////////////////////////1D FIT FUNCTIONS//////////////////////////////////
template<fitType T>
int fit_function_1d(const gsl_vector *fit_params, void *data_points,
    gsl_vector *fit_residuals) {
  // Computes the fit residuals for the a_{x or y} (r) 1D fit
  // data_points should point to a vector of data points
  // which all have the same E value

  static_assert(T==fitType::AX || T==fitType::AY || T==fitType::BETA,
                "Error: 1D fit only used for ax, ay, beta");
  // Cast data_points to the appropriate type
  std::vector<DataPoint>& data = * (std::vector<DataPoint> *) data_points;

  std::array<double, 5> r_coefs;
  double k = gsl_vector_get(fit_params, 0);
  double a = gsl_vector_get(fit_params, 1);
  double b = gsl_vector_get(fit_params, 2);
  double c = gsl_vector_get(fit_params, 3);
  double d = gsl_vector_get(fit_params, 4);
  if constexpr (T==fitType::BETA) {
    r_coefs = {k, a, b, c, d};
  } else {
    r_coefs = {a, b, c, d, 0};
  }

  size_t num_pts = data.size();
  for (size_t i=0; i<num_pts; i++) {
    auto r = data[i].r;
    double real_val;
    if constexpr (T==fitType::AX) real_val = data[i].ax;
    else if constexpr (T==fitType::AY) real_val = data[i].ay;
    else if constexpr (T==fitType::BETA) real_val = data[i].beta;
    double fit_val = gsl_poly_eval(r_coefs.data(), 5, r);
    if constexpr (T!=fitType::BETA) fit_val *= std::exp(-k*r);
    gsl_vector_set(fit_residuals, i, fit_val - real_val);
  }

  return GSL_SUCCESS;
}

// Specializations
//int ax_fit_function_1d(const gsl_vector *fit_params, void *data_points,
//    gsl_vector *fit_residuals) {
//  return fit_function_1d<fitType::AX>(fit_params, data_points, fit_residuals);
//}
//int ay_fit_function_1d(const gsl_vector *fit_params, void *data_points,
//    gsl_vector *fit_residuals) {
//  return fit_function_1d<fitType::AY>(fit_params, data_points, fit_residuals);
//}
//int beta_fit_function_1d(const gsl_vector *fit_params, void *data_points,
//    gsl_vector *fit_residuals) {
//  return fit_function_1d<fitType::BETA>(fit_params, data_points, fit_residuals);
//}

/////////////////////////FITTING ROUTINES///////////////////////////
template<fitType T>
FitResults fit_routine(const std::vector<DataPoint>& data_points,
    double crossover_point) {
  // This function performs a hybrid 1D/2D fit to
  // ax/ay as functions of r or E and r, and returns the obtained fit
  // parameters.  If fitting fails, { {}, {-1,-1...}} is returned.
  // data_points should be sorted by E and then by r

  // Set up the return type
  FitResults result, junk;
  if constexpr (T==fitType::CX) {
    result = cFitResults{};
    junk = cFitResults{-1, -1, -1, -1, -1, -1, -1, -1};
  } else if constexpr(T==fitType::AX || T==fitType::AY) {
    result = aFitResults{};
    junk = aFitResults{{}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, }};
  } else if constexpr(T==fitType::BETA) {
    result = betaFitResults{};
    junk = betaFitResults{{}, {-1, -1, -1, -1, -1, -1, -1, -1}};
  }

  // Partition data_points into the 1D and 2D regions
  // for all but the cx fit
  // below crossover_point -> 1D
  // above crossover_point -> 2D
  // Also convert to MeV and cm here
  auto convert = [](DataPoint& p) { p.E *= 1e-6; p.r *= 1e2; };
  std::vector<std::vector<DataPoint>> data_points_1d{1};
  std::vector<DataPoint> data_points_2d;
  std::vector<double> low_E_vals;
  if constexpr (T!=fitType::CX) {
    low_E_vals.reserve(10);
    low_E_vals.push_back(data_points[0].E * 1e-6);
    size_t ix_1d=0;
    for (const auto& p : data_points) {
      if (p.E < crossover_point) { // goes in 1D lists
        if (p.E*1e-6 != low_E_vals.back()) { // new E value encountered
          low_E_vals.push_back(p.E * 1e-6);
          data_points_1d.push_back({});
          ix_1d++;
        }
        data_points_1d[ix_1d].push_back(p);
        convert(data_points_1d[ix_1d].back());
      } else { // goes in 2D list
        data_points_2d.push_back(p);
        convert(data_points_2d.back());
      }
    }
  } else {
    data_points_2d = data_points;
    for (auto& p : data_points_2d) convert(p);
  }

  // GSL boilerplate
  const gsl_multifit_nlinear_type *mnT = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_vector *fit_residuals;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
  constexpr size_t num_fit_params_1d = 5;
  size_t num_fit_params_2d;
  //double initial_parameter_guess_1d[num_fit_params_1d];
  //for (auto& p : initial_parameter_guess_1d) p = 1.0;
  //double initial_parameter_guess_2d[10];
  if constexpr (T==fitType::CX) {
    num_fit_params_2d = 7;
    //for (auto& p : initial_parameter_guess_2d) p = 0.0;
  } else if constexpr (T==fitType::AX || T==fitType::AY) {
    num_fit_params_2d = 9;
    //for (auto& p : initial_parameter_guess_2d) p = 1.0;
  } else {
    num_fit_params_2d = 7;
    //for (auto& p : initial_parameter_guess_2d) p = 0.0;
  }

  //gsl_vector_view init_param_gsl_1d = gsl_vector_view_array(initial_parameter_guess_1d, num_fit_params_1d);
  //gsl_vector_view init_param_gsl_2d = gsl_vector_view_array(initial_parameter_guess_2d, num_fit_params_2d);
  gsl_vector* init_param_gsl_1d = gsl_vector_alloc(num_fit_params_1d);
  gsl_vector_set_all(init_param_gsl_1d, 1.0);
  gsl_vector* init_param_gsl_2d = gsl_vector_alloc(num_fit_params_2d);
  gsl_vector_set_all(init_param_gsl_2d, 1.0);

  double chisq;
  int status, info;

  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;
  const size_t max_iter = 100;

  // First, do the 1D fits for everything except cx
  if constexpr(T!=fitType::CX) {
    fdf.f = fit_function_1d<T>;
    //if constexpr(T==fitType::AX)
    //  fdf.f = ax_fit_function_1d;
    //else if constexpr(T==fitType::AY)
    //  fdf.f = ay_fit_function_1d;
    //else if constexpr(T==fitType::BETA)
    //  fdf.f = beta_fit_function_1d;

    fdf.df = nullptr;
    fdf.fvv = nullptr;
    fdf.p = num_fit_params_1d;

    constexpr size_t ix = type_ix<T>();
    auto& fit_results_1d = std::get<ix>(result).low_e_fits;
    fit_results_1d.reserve(data_points_1d.size());

    size_t ix_1d=0;
    for (auto& v_1d : data_points_1d) {
      // more setup
      fdf.n = v_1d.size();
      fdf.params = (void *) &v_1d;
      w = gsl_multifit_nlinear_alloc(mnT, &fdf_params, v_1d.size(), num_fit_params_1d);
      // do the fit
      gsl_multifit_nlinear_init(init_param_gsl_1d, &fdf, w);
      status = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol, nullptr, nullptr, &info, w);
      // check for errors
      if (status == GSL_EMAXITER)
        std::cout << "Iteration limit reached\n";
      if (status == GSL_ENOPROG) {
        gsl_multifit_nlinear_free(w);
        gsl_vector_free(init_param_gsl_1d);
        gsl_vector_free(init_param_gsl_2d);
        return junk;
      }
      // compute chisq
      fit_residuals = gsl_multifit_nlinear_residual(w);
      gsl_blas_ddot(fit_residuals, fit_residuals, &chisq);
      size_t dof = v_1d.size() - num_fit_params_1d;
      chisq = sqrt(chisq/dof);
      // add results to fit_results_1d
      double fit_k = gsl_vector_get(w->x, 0);
      double fit_a = gsl_vector_get(w->x, 1);
      double fit_b = gsl_vector_get(w->x, 2);
      double fit_c = gsl_vector_get(w->x, 3);
      double fit_d = gsl_vector_get(w->x, 4);
      // rescale for correct units
      if constexpr (T==fitType::BETA) {
        fit_k *= 1e0; // used as a0
        fit_a *= 1e2; // used as a1
        fit_b *= 1e4; // used as a2
        fit_c *= 1e6; // used as a3
        fit_d *= 1e8; // used as a4
      } else {
        fit_k *= 1e2;
        fit_a *= 1e0;
        fit_b *= 1e2;
        fit_c *= 1e4;
        fit_d *= 1e6;
      }
      fit_results_1d.push_back({low_E_vals[ix_1d], fit_k, fit_a, fit_b, fit_c, fit_d, chisq});
      gsl_multifit_nlinear_free(w);
      ix_1d++;
    }
  }


  // Now do 2D fits
  fdf.f = fit_function_2d<T>;
  //if constexpr(T==fitType::CX)
  //  fdf.f = cx_fit_function;
  //else if constexpr(T==fitType::AX)
  //  fdf.f = ax_fit_function_2d;
  //else if constexpr(T==fitType::AY)
  //  fdf.f = ay_fit_function_2d;
  //else if constexpr(T==fitType::BETA)
  //  fdf.f = beta_fit_function_2d;
  fdf.df = nullptr;
  fdf.p = num_fit_params_2d;
  fdf.n = data_points_2d.size();
  fdf.params = (void *) &data_points_2d;
  w = gsl_multifit_nlinear_alloc(mnT, &fdf_params, data_points_2d.size(), num_fit_params_2d);
  // do the fit
  gsl_multifit_nlinear_init(init_param_gsl_2d, &fdf, w);
  status = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol, nullptr, nullptr, &info, w);
  // check for errors
  if (status == GSL_EMAXITER)
    std::cout << "Iteration limit reached\n";
  if (status == GSL_ENOPROG) {
    gsl_multifit_nlinear_free(w);
    gsl_vector_free(init_param_gsl_1d);
    gsl_vector_free(init_param_gsl_2d);
    return junk;
  }
  // compute chisq
  fit_residuals = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(fit_residuals, fit_residuals, &chisq);
  size_t dof = data_points_2d.size() - num_fit_params_2d;
  chisq = sqrt(chisq/dof);
  // add results to fit_results_1d
  constexpr size_t ix = type_ix<T>();
  if constexpr(T==fitType::CX || T==fitType::BETA) {
    double a0 = 1e+18*gsl_vector_get(w->x, 0);
    double a1 = 1e+12*gsl_vector_get(w->x, 1);
    double a2 = 1e+06*gsl_vector_get(w->x, 2);
    double b0 = 1e-18*gsl_vector_get(w->x, 3);
    double b1 = 1e-16*gsl_vector_get(w->x, 4);
    double b2 = 1e-14*gsl_vector_get(w->x, 5);
    double b3 = 1e-12*gsl_vector_get(w->x, 6);
    if constexpr (T==fitType::CX) {
      auto& fit_result = std::get<ix>(result);
      fit_result = {a0, a1, a2, b0, b1, b2, b3, chisq};
    } else {
      auto& fit_result = std::get<ix>(result).high_e_fit;
      fit_result = {a0, a1, a2, b0, b1, b2, b3, chisq};
    }
  } else {
    double ke = 1e-06 * std::abs(gsl_vector_get(w->x, 0));
    double kr = 1e+02 * std::abs(gsl_vector_get(w->x, 1));
    double ae = 1e+18 * gsl_vector_get(w->x, 2);
    double be = 1e+12 * gsl_vector_get(w->x, 3);
    double ce = 1e+06 * gsl_vector_get(w->x, 4);
    //double de = gsl_vector_get(w->x, 5);
    double ar = 1e-18 * gsl_vector_get(w->x, 5);
    double br = 1e-16 * gsl_vector_get(w->x, 6);
    double cr = 1e-14 * gsl_vector_get(w->x, 7);
    double dr = 1e-12 * gsl_vector_get(w->x, 8);
    auto& fit_result = std::get<ix>(result).high_e_fit;
    fit_result = {ke, kr, ae, be, ce, ar, br, cr, dr, chisq};
  }
  gsl_multifit_nlinear_free(w);
  gsl_vector_free(init_param_gsl_1d);
  gsl_vector_free(init_param_gsl_2d);

  return result;
}

// Specializations
cFitResults cx_fit(const std::vector<DataPoint>& data_points) {
  auto result = fit_routine<fitType::CX>(data_points, 0);
  return std::get<cFitResults>(result);
}
aFitResults ax_fit(const std::vector<DataPoint>& data_points, double xpt) {
  auto result = fit_routine<fitType::AX>(data_points, xpt);
  return std::get<aFitResults>(result);
}
aFitResults ay_fit(const std::vector<DataPoint>& data_points, double xpt) {
  auto result = fit_routine<fitType::AY>(data_points, xpt);
  return std::get<aFitResults>(result);
}
betaFitResults beta_fit(const std::vector<DataPoint>& data_points, double xpt) {
  auto result = fit_routine<fitType::BETA>(data_points, xpt);
  return std::get<betaFitResults>(result);
}
