#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_poly.h>

#include "meta_fit.hpp"

// Make fitType printable
using namespace std::string_literals;
std::string fit_to_string(fitType T) {
  switch (T) {
    case fitType::CX:       return "c_x"s;
    case fitType::AX:       return "alpha_x"s;
    case fitType::AY:       return "alpha_y"s;
    case fitType::BETA:     return "beta"s;
    case fitType::DXDS_MIN: return "dxds_min"s;
    case fitType::DXDS_MAX: return "dxds_max"s;
    case fitType::DYDS_MAX: return "dyds_max"s;
    default:                return ""s;
  }
}
std::ostream& operator<<(std::ostream& out, fitType T) { return (out << fit_to_string(T)); }

// ~move constructor for MetaFitResults
MetaFitResults:: MetaFitResults(double p_Ein, double p_T, FitResults&& p_cx, FitResults&& p_ax,
    FitResults&& p_ay, FitResults&& p_beta, FitResults&& p_dxds_min, FitResults&& p_dxds_max,
    FitResults&& p_dyds_max, TableData&& p_ertable, TableData&& p_polztable)
  : Ein{p_Ein}, T{p_T}, cx{p_cx}, ax{p_ax}, ay{p_ay}, beta{p_beta},
    dxds_min{p_dxds_min}, dxds_max{p_dxds_max}, dyds_max{p_dyds_max},
    er_table{p_ertable}, polz_table{p_polztable} {}

////////////////////// Fit Routine implementation ///////////////////
// Helper templates
inline FitResults make_junk() {
  FitResults f;
  f.high_e_fit.chi2 = -1;
  return f;
}


template<fitType T, typename F>
void rescale(F& fit) {
  // Rescales the given fit from MeV, cm to eV, m
  if constexpr (std::is_same_v<F, FitResults1D>) {
    if constexpr (T==fitType::CX || T==fitType::BETA) {
      fit.a *= 1e2; fit.b *= 1e4; fit.c *= 1e6; fit.d *= 1e8;
    } else {
      fit.a *= 1e0; fit.b *= 1e2; fit.c *= 1e4; fit.d *= 1e6;
    }
  } else {
    fit.a1 *= 1e-6; fit.a2 *= 1e-12; fit.a3 *= 1e-18;
    fit.b1 *= 1e2; fit.b2 *= 1e4; fit.b3 *= 1e6;
    fit.ke *= 1e-6; fit.kr *= 1e2;
  }
}


template<fitType T, typename F>
inline F gsl_to_fit(const gsl_vector *fit_params) {
  // Pulls the fit coefficients out of fit_params and
  // returns them in an appropriately-typed fit struct
  if constexpr (std::is_same_v<F, FitResults2D>) {
    if constexpr(T==fitType::DXDS_MIN) {
      return F{gsl_vector_get(fit_params, 0),
              gsl_vector_get(fit_params, 1),
              gsl_vector_get(fit_params, 2),
              gsl_vector_get(fit_params, 3),
              gsl_vector_get(fit_params, 4),
              gsl_vector_get(fit_params, 5),
              gsl_vector_get(fit_params, 6),
              std::abs(gsl_vector_get(fit_params, 7)),
              std::abs(gsl_vector_get(fit_params, 8)),
              gsl_vector_get(fit_params, 9), 0.0};
    } else if constexpr(T==fitType::CX || T==fitType::BETA) {
      return F{gsl_vector_get(fit_params, 0),
              gsl_vector_get(fit_params, 1),
              gsl_vector_get(fit_params, 2),
              0.0, // b0 set to 0
              gsl_vector_get(fit_params, 3),
              gsl_vector_get(fit_params, 4),
              gsl_vector_get(fit_params, 5),
              std::abs(gsl_vector_get(fit_params, 6)),
              std::abs(gsl_vector_get(fit_params, 7)), 0.0, 0.0};
    } else {
      return F{gsl_vector_get(fit_params, 0),
              gsl_vector_get(fit_params, 1),
              gsl_vector_get(fit_params, 2),
              gsl_vector_get(fit_params, 3),
              gsl_vector_get(fit_params, 4),
              gsl_vector_get(fit_params, 5),
              gsl_vector_get(fit_params, 6),
              std::abs(gsl_vector_get(fit_params, 7)),
              std::abs(gsl_vector_get(fit_params, 8)), 0.0, 0.0};
    }
  } else {
    return F{0.0, gsl_vector_get(fit_params, 0),
            gsl_vector_get(fit_params, 1),
            gsl_vector_get(fit_params, 2),
            gsl_vector_get(fit_params, 3), 0.0};
  }
}


template<fitType T>
inline double get_real_val(const CauchyPoint& p) {
  // Returns the value from p corresponding to the given fit type
  if      constexpr (T==fitType::CX)       return p.cx;
  else if constexpr (T==fitType::AX)       return p.ax;
  else if constexpr (T==fitType::AY)       return p.ay;
  else if constexpr (T==fitType::BETA)     return p.beta;
  else if constexpr (T==fitType::DXDS_MIN) return p.dxds_min;
  else if constexpr (T==fitType::DXDS_MAX) return p.dxds_max;
  else                                     return p.dyds_max;
}



//////////////////////////// FIT FUNCTION //////////////////////////////
template<fitType T, typename F>
int fit_function(const gsl_vector *fit_params, void *data_points,
    gsl_vector *fit_residuals) {
  // Computes the residuals for each fit point in the fit,
  // taylored to the calling convention required by gsl_multifit_nlinear
  std::vector<CauchyPoint> data = * (std::vector<CauchyPoint> *) data_points;

  F fit = gsl_to_fit<T,F>(fit_params);
  size_t ix = 0;
  constexpr size_t dim = (std::is_same_v<F, FitResults1D>)?1:2;
  for (const auto& p : data)
    gsl_vector_set(fit_residuals, ix++, eval<T,dim>(fit, p) - get_real_val<T>(p));

  return GSL_SUCCESS;
}


///////////////////////// FITTING ROUTINE ///////////////////////////
template<fitType T>
FitResults fit_routine(const std::vector<CauchyPoint>& data_points,
    const TableData& er_table, double crossover_point) {
  // This function performs a hybrid 1D/2D fit to
  // ax/ay as functions of r or E and r, and returns the obtained fit
  // parameters.  If fitting fails, { {}, {-1,-1...}} is returned.
  // data_points should be sorted by E and then by r
  // The fits are weighted using the er_table of amplitudes, so that
  // bins which are more likely are weighted higher

  FitResults result;

  // Partition data_points into the 1D and 2D regions
  // for all but the cx fit
  // below crossover_point -> 1D
  // above crossover_point -> 2D
  // Also convert to MeV and cm here
  auto convert = [](CauchyPoint& p) { p.E *= 1e-6; p.r *= 1e2; };
  std::vector<std::vector<CauchyPoint>> data_points_1d{1};
  std::vector<CauchyPoint> data_points_2d;
  std::vector<double> low_E_vals;
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
  // If all points fell into data_poinst_2d, clear low_e_vals
  // and data_points_1d
  if (data_points_1d[0].size() == 0) {
    data_points_1d.clear();
    low_E_vals.clear();
  }

  // GSL boilerplate
  const gsl_multifit_nlinear_type *mnT = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_vector *fit_residuals;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
  size_t num_fit_params_1d = 4, num_fit_params_2d;
  if constexpr (T==fitType::DXDS_MIN)
    num_fit_params_2d = 10;
  else if constexpr (T==fitType::CX || T==fitType::BETA)
    num_fit_params_2d = 8;
  else
    num_fit_params_2d = 9;

  gsl_vector* init_param_gsl_1d;
  init_param_gsl_1d = gsl_vector_alloc(num_fit_params_1d);
  gsl_vector_set_all(init_param_gsl_1d, 1.0);
  gsl_vector* init_param_gsl_2d = gsl_vector_alloc(num_fit_params_2d);
  gsl_vector_set_all(init_param_gsl_2d, 1.0);

  double chisq;
  int status, info;

  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;
  const size_t max_iter = 1000;

  // First, do the 1D fits
  fdf.f = fit_function<T,FitResults1D>;
  fdf.df = nullptr;
  fdf.fvv = nullptr;
  fdf.p = num_fit_params_1d;

  //constexpr size_t ix = type_ix<T>();
  auto& fit_results_1d = result.low_e_fits;
  fit_results_1d.reserve(data_points_1d.size());

  ix_1d=0;
  for (const auto& v_1d : data_points_1d) {
    // more setup
    fdf.n = v_1d.size();
    // Warn and continue if not enough points
    if (fdf.p > fdf.n) {
      std::cerr << "WARNING: not enough points available for 1D fit to "
        << T << " at pc_out = " << low_E_vals[ix_1d++] << " MeV\n";
      continue;
    }
    fdf.params = (void *) &v_1d;
    gsl_vector *weights = gsl_vector_alloc(v_1d.size());
    gsl_vector_set_all(weights, 0);

    for (size_t w_ix = 0; w_ix < v_1d.size(); w_ix++)
      gsl_vector_set(weights, w_ix, er_table.data[w_ix + ix_1d*v_1d.size()]);
    w = gsl_multifit_nlinear_alloc(mnT, &fdf_params, v_1d.size(), num_fit_params_1d);
    // do the fit
    gsl_multifit_nlinear_winit(init_param_gsl_1d, weights, &fdf, w);
    status = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol, nullptr, nullptr, &info, w);
    // check for errors
    //if (status == GSL_EMAXITER)
    //  std::cout << "Iteration limit reached for 1D fit to " <<
    //    T << " at pc_out = " << low_E_vals[ix_1d] << " MeV\n";
    if (status == GSL_ENOPROG) {
      gsl_multifit_nlinear_free(w);
      gsl_vector_free(init_param_gsl_1d);
      gsl_vector_free(init_param_gsl_2d);
      gsl_vector_free(weights);
      return make_junk();
    }
    // compute chisq
    fit_residuals = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(fit_residuals, fit_residuals, &chisq);
    size_t dof = v_1d.size() - num_fit_params_1d;
    chisq = sqrt(chisq/dof);
    // add results to fit_results_1d
    auto fit_result = gsl_to_fit<T,FitResults1D>(w->x);
    rescale<T,FitResults1D>(fit_result);
    fit_result.E = 1e6 * low_E_vals[ix_1d];
    fit_result.chi2 = chisq;
    fit_results_1d.push_back(fit_result);
    gsl_multifit_nlinear_free(w);
    gsl_vector_free(weights);
    ix_1d++;
  }


  // Now do 2D fits
  fdf.f = fit_function<T,FitResults2D>;
  fdf.p = num_fit_params_2d;
  fdf.n = data_points_2d.size();
  if (fdf.p > fdf.n) {
    std::cerr << "WARNING: not enough points available for 2D fit to "
      << T << '\n';
    gsl_vector_free(init_param_gsl_1d);
    gsl_vector_free(init_param_gsl_2d);
    return make_junk();
  }
  if constexpr (T==fitType::DXDS_MIN) // Initial guess of all 1.0 does not work well for dxds_min
    gsl_vector_set_all(init_param_gsl_2d, 0.0);
  gsl_vector *weights = gsl_vector_alloc(fdf.n);

  size_t total_1d_size = data_points_1d.size() * data_points_1d[0].size();
  for (size_t w_ix = 0; w_ix < fdf.n; w_ix++)
    gsl_vector_set(weights, w_ix, er_table.data[total_1d_size + w_ix]);
  fdf.params = (void *) &data_points_2d;
  w = gsl_multifit_nlinear_alloc(mnT, &fdf_params, data_points_2d.size(), num_fit_params_2d);
  // do the fit
  gsl_multifit_nlinear_winit(init_param_gsl_2d, weights, &fdf, w);
  status = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol, nullptr, nullptr, &info, w);
  // check for errors
  //if (status == GSL_EMAXITER)
  //  std::cout << "Iteration limit reached for 2D fit to " << T << '\n';
  if (status == GSL_ENOPROG) {
    gsl_multifit_nlinear_free(w);
    gsl_vector_free(init_param_gsl_1d);
    gsl_vector_free(init_param_gsl_2d);
    gsl_vector_free(weights);
    return make_junk();
  }
  // compute chisq
  fit_residuals = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(fit_residuals, fit_residuals, &chisq);
  size_t dof = data_points_2d.size() - num_fit_params_2d;
  chisq = sqrt(chisq/dof);
  // add results to fit_results_1d
  auto fit_result = gsl_to_fit<T,FitResults2D>(w->x);
  rescale<T,FitResults2D>(fit_result);
  fit_result.chi2 = chisq;
  result.high_e_fit = fit_result;

  gsl_multifit_nlinear_free(w);
  gsl_vector_free(init_param_gsl_1d);
  gsl_vector_free(init_param_gsl_2d);
  gsl_vector_free(weights);

  return result;
}


// Instantiate fit_routine
template FitResults fit_routine<fitType::CX>      (const std::vector<CauchyPoint>&, const TableData&, double);
template FitResults fit_routine<fitType::AX>      (const std::vector<CauchyPoint>&, const TableData&, double);
template FitResults fit_routine<fitType::AY>      (const std::vector<CauchyPoint>&, const TableData&, double);
template FitResults fit_routine<fitType::BETA>    (const std::vector<CauchyPoint>&, const TableData&, double);
template FitResults fit_routine<fitType::DXDS_MIN>(const std::vector<CauchyPoint>&, const TableData&, double);
template FitResults fit_routine<fitType::DXDS_MAX>(const std::vector<CauchyPoint>&, const TableData&, double);
template FitResults fit_routine<fitType::DYDS_MAX>(const std::vector<CauchyPoint>&, const TableData&, double);
