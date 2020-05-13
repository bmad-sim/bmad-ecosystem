#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include "meta_fit.hpp"
#include "gnuplot.hpp"

void output_cauchy_impl(std::ostream& gpfile, const DataPoint& p, double pc, double r) {
  gpfile << "set isosamples 20\n"; // TODO: Use num_bins from config file instead
  gpfile << "set xlabel 'dx/ds'\n";
  gpfile << "set ylabel 'dy/ds'\n";
  gpfile << "set zlabel 'Bin count'\n";
  gpfile << "set title 'pc out = " << pc << " MeV, r = " << r << " cm'\n";
  // Add a label with the fit parameters
  gpfile << "unset label\n";
  gpfile << "set label 'cx = "   << p.cx  << "' at graph 0.6, graph 0.5, graph 1.0\n";
  gpfile << "set label 'ax = "   << p.ax  << "' at graph 0.6, graph 0.5, graph 0.9\n";
  gpfile << "set label 'ay = "   << p.ay  << "' at graph 0.6, graph 0.5, graph 0.8\n";
  gpfile << "set label 'beta = " << p.beta << "' at graph 0.6, graph 0.5, graph 0.7\n";
  // Add the plots
  char bincmd[200];
  sprintf(bincmd, "splot 'E%0.2lf_r%0.3lf_bin.dat' u 1:2:3 w p title 'Binned data'\n", pc, r);
  gpfile << bincmd;
  gpfile << "cx = " << p.cx << '\n';
  gpfile << "ax = " << p.ax << '\n';
  gpfile << "ay = " << p.ay << '\n';
  gpfile << "beta = " << p.beta << '\n';
  gpfile << "amp = " << p.amp << '\n';
  gpfile << "rep amp * (1 + beta*x)/(1 + ax**2*(x-cx)**2 + ay**2*y**2) title 'Cauchy fit'\n";
}

void output_cauchy_gp(const char* foldername, const DataPoint& p) {
  // Writes a gnuplot file to plot the cauchy distribution
  // obtained from the direct fits alongside the observed data
  // Rescale E and r
  double pc_out = p.E/1e6;
  double r_out = p.r*1e2;
  std::ofstream cauchy_gp;
  char cauchy_gp_name[100];
  sprintf(cauchy_gp_name, "%s/cauchy_E%0.2lf_r%0.3lf.gp", foldername, pc_out, r_out);
  cauchy_gp.open(cauchy_gp_name);

  output_cauchy_impl(cauchy_gp, p, pc_out, r_out);

  cauchy_gp.close();
  return;
}

template<size_t dim>
void output_cauchy_meta_gp(const char* foldername, const MetaFitResults& mf, double pc_out, double r) {
  // Writes a gnuplot file to plot the cauchy distribution
  // obtained from the meta fits alongside the observed data
  std::ofstream cauchy_gp;
  char cauchy_gp_name[100];
  sprintf(cauchy_gp_name, "%s/meta_E%0.2lf_r%0.3lf.gp", foldername, pc_out, r);
  cauchy_gp.open(cauchy_gp_name);

  DataPoint p;
  p.E = pc_out;
  p.r = r;
  if constexpr (dim==1) {
    size_t i = 0;
    for (const auto& f : mf.cx.low_e_fits) {
      if (f.E - pc_out < 1e6) break;
      i++;
    }
    p.cx = eval<fitType::CX, dim>(mf.cx.low_e_fits[i], p);
    p.ax = eval<fitType::AX, dim>(mf.ax.low_e_fits[i], p);
    p.ay = eval<fitType::AY, dim>(mf.ay.low_e_fits[i], p);
    p.beta = eval<fitType::BETA, dim>(mf.beta.low_e_fits[i], p);
  } else {
    p.cx = eval<fitType::CX, dim>(mf.cx.high_e_fit, p);
    p.ax = eval<fitType::AX, dim>(mf.ax.high_e_fit, p);
    p.ay = eval<fitType::AY, dim>(mf.ay.high_e_fit, p);
    p.beta = eval<fitType::BETA, dim>(mf.beta.high_e_fit, p);
  }


  output_cauchy_impl(cauchy_gp, p, pc_out, r);

  cauchy_gp.close();
  return;
}

template<fitType T>
void write_poly_coefs(const FitResults2D& fit, std::ostream& out) {
  // Pulls any polynomial coefficients from the fit struct and
  // writes them on separate lines in the given ostream
  out << "a1 = " << fit.a1 << '\n';
  out << "a2 = " << fit.a2 << '\n';
  out << "a3 = " << fit.a3 << '\n';
  out << "b0 = " << fit.b0 << '\n';
  out << "b1 = " << fit.b1 << '\n';
  out << "b2 = " << fit.b2 << '\n';
  out << "b3 = " << fit.b3 << '\n';
  out << "ke = " << fit.ke << '\n';
  out << "kr = " << fit.kr << '\n';
  out << "C = " << fit.C << '\n';
  return;
}


template<fitType T>
void write_1d_files(const FitResults& fit, const std::vector<DataPoint>& data, const char* dir) {
  // Writes the 1D fit files for each of the low_e_fits in fit
  // DATA MUST BE SORTED BY E AND THEN R
  std::ofstream gpfile;
  char filename[100];
  std::string fit_name = fit_to_string(T); // l-value created on the stack so it lives until the end of this function
  const char *prefix = fit_name.c_str();
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
void output_metafit_gp(const FitResults& fit, const std::vector<DataPoint>& data, const char* dir, double xpt) {
  std::ofstream gpfile;
  std::string fit_name = fit_to_string(T); // l-value created on the stack so it lives until the end of this function
  const char *prefix = fit_name.c_str();
  char filename[100];
  sprintf(filename, "%s/%s_2d.gp", dir, prefix);
  gpfile.open(filename);
  if (!gpfile.good()) {
    std::cerr << "what";
    return;
  }

  write_1d_files<T>(fit, data, dir);
  write_poly_coefs<T>(fit.high_e_fit, gpfile);

  gpfile << prefix << "_fit(E,r) = (1 + a1*E + a2*E**2 + a3*E**3)*exp(-ke*E)"
     << " * (b0 + b1*r + b2*r**2 + b3*r**3) * exp(-kr*r) + C\n";

  gpfile << prefix << "_fit_trimmed(E,r) = (E>" << xpt
    << ")?" << prefix << "_fit(E,r):(1/0)\n";

  gpfile.close();

  sprintf(filename, "%s/%s_master.gp", dir, prefix);
  gpfile.open(filename);
  gpfile << "set title '" << prefix << " Fit Results'\n";
  gpfile << "unset label\n";
  gpfile << "set xlabel 'pc out (eV)'\n";
  gpfile << "set ylabel 'r out (m)'\n";
  gpfile << "set zlabel '" << prefix << "'\n";
  gpfile << "call '" << prefix << "_2d.gp'\n";
  gpfile << "call '" << prefix << "_1d_master.gp'\n";
  gpfile << "rep 'coef.dat' u 'E':'r':'" << fit_name << "' w p title 'Cauchy data'\n";
  gpfile << "rep " << prefix << "_fit_trimmed(x,y) title '2D fit'";
  gpfile.close();
  return;
}

// Instantiate output_gp
template void output_metafit_gp<fitType::CX>(const FitResults&, const std::vector<DataPoint>&, const char*, double);
template void output_metafit_gp<fitType::AX>(const FitResults&, const std::vector<DataPoint>&, const char*, double);
template void output_metafit_gp<fitType::AY>(const FitResults&, const std::vector<DataPoint>&, const char*, double);
template void output_metafit_gp<fitType::BETA>(const FitResults&, const std::vector<DataPoint>&, const char*, double);
template void output_metafit_gp<fitType::DXDS_MIN>(const FitResults&, const std::vector<DataPoint>&, const char*, double);
template void output_metafit_gp<fitType::DXDS_MAX>(const FitResults&, const std::vector<DataPoint>&, const char*, double);
template void output_metafit_gp<fitType::DYDS_MAX>(const FitResults&, const std::vector<DataPoint>&, const char*, double);
