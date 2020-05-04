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

// Helper templates
template<fitType T>
constexpr fit_t<T> make_junk() {
  static_assert(good_fit_type<T,2>());
  if constexpr (T==fitType::CX) {
    return fit_t<T>{-1, -1, -1, -1, -1, -1, -1};
  } else if constexpr (T==fitType::AX || T==fitType::AY) {
    return fit_t<T>{{}, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, }};
  } else {
    return fit_t<T>{{}, {-1, -1, -1, -1, -1, -1, -1}};
  }
}

template<fitType T, size_t dim>
void rescale(fit_t_part<T,dim>& fit) {
  // Rescales the given fit from MeV, cm to eV, m
  static_assert(good_fit_type<T,dim>());
  if constexpr (dim==1) {
    if constexpr (T==fitType::BETA) {
      fit.a *= 1e2; fit.b *= 1e4; fit.c *= 1e6; fit.d *= 1e8;
    } else {
      fit.k *= 1e2;
      fit.a *= 1e0; fit.b *= 1e2; fit.c *= 1e4; fit.d *= 1e6;
    }
  } else {
    fit.a1 *= 1e-6; fit.a2 *= 1e-12; fit.a3 *= 1e-18;
    fit.b1 *= 1e2; fit.b2 *= 1e4; fit.b3 *= 1e6;
    if constexpr (T==fitType::AX || T==fitType::AY) {
      fit.ke *= 1e-6; fit.kr *= 1e2;
    }
  }
}


template<fitType T, size_t dim>
inline fit_t_part<T,dim> gsl_to_fit(const gsl_vector *fit_params) {
  // Pulls the fit coefficients out of fit_params and
  // returns them in an appropriately-typed fit struct
  static_assert(good_fit_type<T,dim>());
  if constexpr (dim==2) {
    if constexpr(T==fitType::CX || T==fitType::BETA) {
      return fit_t_part<T,dim>{gsl_vector_get(fit_params, 0),
              gsl_vector_get(fit_params, 1),
              gsl_vector_get(fit_params, 2),
              gsl_vector_get(fit_params, 3),
              gsl_vector_get(fit_params, 4),
              gsl_vector_get(fit_params, 5), 0.0};
    } else {
      return fit_t_part<T,dim>{std::abs(gsl_vector_get(fit_params, 0)),
              std::abs(gsl_vector_get(fit_params, 1)),
              gsl_vector_get(fit_params, 2),
              gsl_vector_get(fit_params, 3),
              gsl_vector_get(fit_params, 4),
              gsl_vector_get(fit_params, 5),
              gsl_vector_get(fit_params, 6),
              gsl_vector_get(fit_params, 7),
              gsl_vector_get(fit_params, 8), 0.0};
    }
  } else {
    if constexpr(T==fitType::BETA) {
      return fit_t_part<T,dim>{0.0, gsl_vector_get(fit_params, 0),
              gsl_vector_get(fit_params, 1),
              gsl_vector_get(fit_params, 2),
              gsl_vector_get(fit_params, 3), 0.0};
    } else {
      return fit_t_part<T,dim>{0.0, std::abs(gsl_vector_get(fit_params, 0)),
              gsl_vector_get(fit_params, 1),
              gsl_vector_get(fit_params, 2),
              gsl_vector_get(fit_params, 3),
              gsl_vector_get(fit_params, 4), 0.0};
    }
  }
}


template<fitType T>
inline double get_real_val(const DataPoint& p) {
  // Returns the value from p corresponding to the given fit type
  static_assert(good_fit_type<T,2>());
  if      constexpr (T==fitType::CX) return p.cx;
  else if constexpr (T==fitType::AX) return p.ax;
  else if constexpr (T==fitType::AY) return p.ay;
  else                               return p.beta;
}



//////////////////////////// FIT FUNCTION //////////////////////////////
template<fitType T, size_t dim>
int fit_function(const gsl_vector *fit_params, void *data_points,
    gsl_vector *fit_residuals) {
  // Computes the residuals for each fit point in the fit,
  // taylored to the calling convention required by gsl_multifit_nlinear
  static_assert(good_fit_type<T,dim>());
  std::vector<DataPoint> data = * (std::vector<DataPoint> *) data_points;

  fit_t_part<T, dim> fit = gsl_to_fit<T,dim>(fit_params);
  size_t ix = 0;
  for (const auto& p : data)
    gsl_vector_set(fit_residuals, ix++, eval<T,dim>(fit, p) - get_real_val<T>(p));

  return GSL_SUCCESS;
}


///////////////////////// FITTING ROUTINE ///////////////////////////
template<fitType T>
fit_t<T> fit_routine(const std::vector<DataPoint>& data_points,
    double crossover_point) {
  // This function performs a hybrid 1D/2D fit to
  // ax/ay as functions of r or E and r, and returns the obtained fit
  // parameters.  If fitting fails, { {}, {-1,-1...}} is returned.
  // data_points should be sorted by E and then by r

  // Set up the return value
  fit_t<T> result = make_junk<T>();


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
  size_t num_fit_params_1d, num_fit_params_2d;
  if constexpr (T==fitType::CX) {
    num_fit_params_1d = 0;
    num_fit_params_2d = 6;
  } else if constexpr (T==fitType::AX || T==fitType::AY) {
    num_fit_params_1d = 5;
    num_fit_params_2d = 9;
  } else {
    num_fit_params_1d = 4;
    num_fit_params_2d = 6;
  }

  gsl_vector* init_param_gsl_1d;
  if constexpr (T!=fitType::CX) {
    init_param_gsl_1d = gsl_vector_alloc(num_fit_params_1d);
    gsl_vector_set_all(init_param_gsl_1d, 1.0);
  }
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
    fdf.f = fit_function<T,1>;
    fdf.df = nullptr;
    fdf.fvv = nullptr;
    fdf.p = num_fit_params_1d;

    //constexpr size_t ix = type_ix<T>();
    auto& fit_results_1d = result.low_e_fits;
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
        return make_junk<T>();
      }
      // compute chisq
      fit_residuals = gsl_multifit_nlinear_residual(w);
      gsl_blas_ddot(fit_residuals, fit_residuals, &chisq);
      size_t dof = v_1d.size() - num_fit_params_1d;
      chisq = sqrt(chisq/dof);
      // add results to fit_results_1d
      auto fit_result = gsl_to_fit<T,1>(w->x);
      rescale<T,1>(fit_result);
      fit_result.E = 1e6 * low_E_vals[ix_1d];
      fit_result.chi2 = chisq;
      fit_results_1d.push_back(fit_result);
      gsl_multifit_nlinear_free(w);
      ix_1d++;
    }
  }


  // Now do 2D fits
  fdf.f = fit_function<T,2>;
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
    if constexpr (T!=fitType::CX) gsl_vector_free(init_param_gsl_1d);
    gsl_vector_free(init_param_gsl_2d);
    return make_junk<T>();
  }
  // compute chisq
  fit_residuals = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(fit_residuals, fit_residuals, &chisq);
  size_t dof = data_points_2d.size() - num_fit_params_2d;
  chisq = sqrt(chisq/dof);
  // add results to fit_results_1d
  auto fit_result = gsl_to_fit<T,2>(w->x);
  rescale<T,2>(fit_result);
  fit_result.chi2 = chisq;
  if constexpr (T==fitType::CX) result = fit_result;
  else result.high_e_fit = fit_result;

  gsl_multifit_nlinear_free(w);
  if constexpr (T!=fitType::CX) gsl_vector_free(init_param_gsl_1d);
  gsl_vector_free(init_param_gsl_2d);

  return result;
}

// Specializations
cFitResults cx_fit(const std::vector<DataPoint>& data_points) { return fit_routine<fitType::CX>(data_points, 0); }
aFitResults ax_fit(const std::vector<DataPoint>& data_points, double xpt) { return fit_routine<fitType::AX>(data_points, xpt); }
aFitResults ay_fit(const std::vector<DataPoint>& data_points, double xpt) { return fit_routine<fitType::AY>(data_points, xpt); }
betaFitResults beta_fit(const std::vector<DataPoint>& data_points, double xpt) { return fit_routine<fitType::BETA>(data_points, xpt); }
