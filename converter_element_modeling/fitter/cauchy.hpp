#pragma once
#include <vector>
struct BinPoint {
  double x;
  double y;
  double density;
};

struct XYBinnedData {
  double E;
  double r;
  std::vector<BinPoint> bins;
};

enum class CauchyStatus { OK, NOPROG, MAXITER };

struct CauchyPoint {
  double E, r;
  double cx, ax, ay, beta, dxds_min, dxds_max, dyds_max;
  double amp;
  double chi2;
  CauchyStatus stat;
};

CauchyPoint asym_cauchy_fit(const XYBinnedData& data);
CauchyPoint asym_cauchy_fit(const XYBinnedData& data, const CauchyPoint& guess);
