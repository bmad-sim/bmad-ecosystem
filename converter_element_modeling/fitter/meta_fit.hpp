#pragma once
#include <vector>
#include <variant>
struct DataPoint {
  double E;
  double r;
  double cx;
  double cy;
  double ax;
  double ay;
  double beta;
  double amp;
};

// Fit types:
enum class fitType { CX, AX, AY, BETA };

//////////////////CX////////////////////

struct cFitResults {
  double a0, a1, a2;
  double b0, b1, b2, b3;
  double chi2;
};

//////////////////////AX/Y////////////////////
struct aFitResults2D {
  double ke;
  double kr;
  double ae;
  double be;
  double ce;
  //double de;
  double ar;
  double br;
  double cr;
  double dr;
  double chi2;
};

struct aFitResults1D {
  double E;
  double k;
  double a;
  double b;
  double c;
  double d;
  double chi2;
};

struct aFitResults {
  std::vector<aFitResults1D> low_e_fits;
  aFitResults2D high_e_fit;
};

////////////////////////BETA/////////////////////

struct betaFitResults2D {
  double a0, a1, a2;
  double b0, b1, b2, b3;
  double chi2;
};

struct betaFitResults1D {
  double E;
  double a0;
  double a1;
  double a2;
  double a3;
  double a4;
  double chi2;
};

struct betaFitResults {
  std::vector<betaFitResults1D> low_e_fits;
  betaFitResults2D high_e_fit;
};

using FitResults = std::variant<cFitResults, aFitResults, betaFitResults>;

template<fitType T>
constexpr size_t type_ix() {
  if constexpr(T==fitType::CX) return 0;
  else if constexpr(T==fitType::AX || T==fitType::AY) return 1;
  else return 2;
}


struct MetaFitResults {
  double Ein;
  double T;
  cFitResults cx;
  aFitResults ax;
  aFitResults ay;
  betaFitResults beta;
};


cFitResults cx_fit(const std::vector<DataPoint>& data_points);

aFitResults ax_fit(const std::vector<DataPoint>& data_points, double crossover_point);
aFitResults ay_fit(const std::vector<DataPoint>& data_points, double crossover_point);

betaFitResults beta_fit(const std::vector<DataPoint>& data_points, double crossover_point);
