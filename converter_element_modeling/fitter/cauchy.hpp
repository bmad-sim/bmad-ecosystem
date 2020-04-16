#pragma once
#include <vector>
struct BinPoint {
  double x;
  double y;
  double count;
};

//struct FitResults {
//  double mu_x;
//  double mu_y;
//  double sigma_x;
//  double sigma_y;
//  double amp;
//  double reduced_chi_square;
//};

struct cauchyFitResults {
  double mu_x;
  double mu_y;
  double sigma_x;
  double sigma_y;
  double beta_x;
  double amp;
  double reduced_chi_square;
};

//struct BigFitResults {
//  double mu_x;
//  double mu_y;
//  double s_x;
//  double s_y;
//  double n_amp;
//  double c_x;
//  double c_y;
//  double a_x;
//  double a_y;
//  double c_amp;
//  double reduced_chi_square;
//};
//
//struct asymBigFitResults {
//  double c_x1;
//  double c_y1;
//  double a_x1;
//  double a_y1;
//  double beta1;
//  double amp1;
//  double c_x2;
//  double c_y2;
//  double a_x2;
//  double a_y2;
//  double beta2;
//  double amp2;
//  double reduced_chi_square;
//};

//double normal2d(double x, double y, double mu_x, double mu_y, double s_x, double s_y);
//FitResults normal_fit(const std::vector<BinPoint>& bins);
//FitResults unweighted_cauchy_fit(const std::vector<BinPoint>& bins);
//FitResults weighted_cauchy_fit(const std::vector<BinPoint>& bins);
//FitResults reverse_weighted_cauchy_fit(const std::vector<BinPoint>& bins);
cauchyFitResults asym_cauchy_fit(const std::vector<BinPoint>& bins);
//aFitResults rw_asym_cauchy_fit(const std::vector<BinPoint>& bins);
//BigFitResults mixed_fit(const std::vector<BinPoint>& bins);
//asymBigFitResults asym_double_fit(const std::vector<BinPoint>& bins);
