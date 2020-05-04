#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include "meta_fit.hpp"

using Dvec = std::vector<DataPoint>;

template<fitType T>
void write_poly_coefs(const fit_t_part<T, 2>& fit, std::ostream& out) {
  // Pulls any polynomial coefficients from the fit struct and
  // writes them on separate lines in the given ostream
  out << "a1 = " << fit.a1 << '\n';
  out << "a2 = " << fit.a2 << '\n';
  out << "a3 = " << fit.a3 << '\n';
  if constexpr (T==fitType::AX||T==fitType::AY)
    out << "b0 = " << fit.b0 << '\n';
  else
    out << "b0 = " << 0.0 << '\n';
  out << "b1 = " << fit.b1 << '\n';
  out << "b2 = " << fit.b2 << '\n';
  out << "b3 = " << fit.b3 << '\n';
  return;
}


template<fitType T>
void write_1d_files(const fit_t<T>& fit, const Dvec& data, const char* dir) {
  // Writes the 1D fit files for each of the low_e_fits in fit
  // DATA MUST BE SORTED BY E AND THEN R
  static_assert(T!=fitType::CX);
  char prefix[10];
  if constexpr (T==fitType::AX) strcpy(prefix, "ax");
  else if constexpr (T==fitType::AY) strcpy(prefix, "ay");
  else strcpy(prefix, "beta");
  std::ofstream gpfile;
  char filename[100];
  size_t ix = 0;
  size_t slice_num = 0;
  double cur_E = data[0].E;
  for (const auto& f : fit.low_e_fits) {
    sprintf(filename, "%s/%s_1d_%lu.dat", dir, prefix, slice_num);
    gpfile.open(filename);
    if (!gpfile.good()) return;
    while (data[ix].E == cur_E) {
      gpfile << data[ix].E << '\t'
        << data[ix].r << '\t'
        << eval<T,1>(f, data[ix]) << '\n';
      ix++;
    }
    cur_E = data[ix].E;
    gpfile.close();
    slice_num++;
  }

  sprintf(filename, "%s/%s_1d_master.gp", dir, prefix);
  gpfile.open(filename);
  gpfile << "splot for [i=0:" << slice_num-1 <<
    "] sprintf('" << prefix << "_1d_%d.dat', i) w lp title '1D slice fit'\n";
  gpfile.close();
  return;
}


template<fitType T>
void output_gp(const fit_t<T>& fit, const Dvec& data, const char* dir, double xpt) {
  std::ofstream gpfile;
  char prefix[10];
  if constexpr (T==fitType::CX) strcpy(prefix, "cx");
  else if constexpr (T==fitType::AX) strcpy(prefix, "ax");
  else if constexpr (T==fitType::AY) strcpy(prefix, "ay");
  else strcpy(prefix, "beta");
  char filename[100];
  sprintf(filename, "%s/%s_2d.gp", dir, prefix);
  gpfile.open(filename);
  if (!gpfile.good()) return;

  if constexpr (T!=fitType::CX) {
    write_1d_files<T>(fit, data, dir);
    write_poly_coefs<T>(fit.high_e_fit, gpfile);
  } else {
    write_poly_coefs<T>(fit, gpfile);
  }
  if constexpr (T==fitType::AX||T==fitType::AY) {
    gpfile << "ke = " << fit.high_e_fit.ke << '\n';
    gpfile << "kr = " << fit.high_e_fit.kr << '\n';
  }

  gpfile << prefix << "_fit(E,r) = (1 + a1*E + a2*E**2 + a3*E**3)"
    << " * (b0 + b1*r + b2*r**2 + b3*r**3)";
  if constexpr (T==fitType::AX || T==fitType::AY)
    gpfile << " * exp(-(ke*E + kr*r))";
  gpfile << '\n';

  gpfile << prefix << "_fit_trimmed(E,r) = (E>" << xpt
    << ")?" << prefix << "_fit(E,r):(1/0)";
  if constexpr (T==fitType::AX || T==fitType::AY)
    gpfile << " * exp(-(ke*E + kr*r))";
  gpfile << '\n';

  gpfile.close();

  sprintf(filename, "%s/%s_master.gp", dir, prefix);
  gpfile.open(filename);
  gpfile << "set title '" << prefix << " Fit Results'\n";
  gpfile << "set xlabel 'pc out (eV)'\n";
  gpfile << "set ylabel 'r out (m)'\n";
  gpfile << "set zlabel '" << prefix << "'\n";
  gpfile << "call '" << prefix << "_2d.gp'\n";
  if constexpr (T!=fitType::CX) {
    gpfile << "call '" << prefix << "_1d_master.gp'\n";
    gpfile << "rep 'coef.dat' u 1:2:";
  } else {
    gpfile << "splot 'coef.dat' u 1:2:";
  }
  if      constexpr (T==fitType::CX) gpfile << '3';
  else if constexpr (T==fitType::AX) gpfile << '5';
  else if constexpr (T==fitType::AY) gpfile << '6';
  else                               gpfile << '7';
  gpfile << " w p title 'Cauchy data'\n";
  gpfile << "rep " << prefix << "_fit";
  if constexpr (T!=fitType::CX) gpfile << "_trimmed";
  gpfile << "(x,y) title '2D fit'";
  gpfile.close();
  return;
}

// Instantiate output_gp
void output_gp_cx(const fit_t<fitType::CX>& fit, const Dvec& data, const char* dir, double xpt) { output_gp<fitType::CX>(fit, data, dir, xpt); }
void output_gp_ax(const fit_t<fitType::AX>& fit, const Dvec& data, const char* dir, double xpt) { output_gp<fitType::AX>(fit, data, dir, xpt); }
void output_gp_ay(const fit_t<fitType::AY>& fit, const Dvec& data, const char* dir, double xpt) { output_gp<fitType::AY>(fit, data, dir, xpt); }
void output_gp_beta(const fit_t<fitType::BETA>& fit, const Dvec& data, const char* dir, double xpt) { output_gp<fitType::BETA>(fit, data, dir, xpt); }
