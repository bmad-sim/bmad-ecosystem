#pragma once
#include <vector>
struct BinPoint {
  double x;
  double y;
  double count;
};

struct FitResults {
  double mu_x;
  double mu_y;
  double sigma_x;
  double sigma_y;
  double amp;
  double reduced_chi_square;
};

double normal2d(double x, double y, double mu_x, double mu_y, double s_x, double s_y);
FitResults normal_fit(const std::vector<BinPoint>& bins);
