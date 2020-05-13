#pragma once
#include <vector>
struct BinPoint {
  double x;
  double y;
  double count;
};

struct cauchyFitResults {
  double mu_x;
  double mu_y;
  double sigma_x;
  double sigma_y;
  double beta_x;
  double amp;
  double reduced_chi_square;
};

cauchyFitResults asym_cauchy_fit(double pc_out, double r, const std::vector<BinPoint>& bins);
