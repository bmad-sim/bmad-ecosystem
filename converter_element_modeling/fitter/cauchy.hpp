#pragma once
#include <vector>
struct BinPoint {
  double x;
  double y;
  double density;
};

enum class CauchyStatus { OK, NOPROG, MAXITER };

struct CauchyPoint {
  double E, r;
  double cx, ax, ay, beta, dxds_min, dxds_max, dyds_max;
  double amp;
  double chi2;
  CauchyStatus stat;
};

CauchyPoint asym_cauchy_fit(double pc_out, double r, const std::vector<BinPoint>& bins);
CauchyPoint asym_cauchy_fit(double pc_out, double r, const std::vector<BinPoint>& bins, const CauchyPoint& guess);
