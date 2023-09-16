#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include "meta_fit.hpp"
#include "gnuplot.hpp"

void output_cauchy_impl(std::ostream& gpfile, const CauchyPoint& p, double pc, double r) {
  gpfile << "set isosamples 20\n"; // TODO: Use num_bins from config file instead
  gpfile << "set xlabel 'dx/ds'\n";
  gpfile << "set ylabel 'dy/ds'\n";
  gpfile << "set zlabel 'Bin density' rotate by 90\n";
  gpfile << "set title 'pc out = " << pc << " MeV, r = " << r << " cm'\n";
  // Add a label with the fit parameters
  gpfile << "unset label\n";
  gpfile << "set label 'c_x = "   << p.cx  << "' at graph 0.6, graph 0.5, graph 1.0\n";
  gpfile << "set label 'alpha_x = "   << p.ax  << "' at graph 0.6, graph 0.5, graph 0.9\n";
  gpfile << "set label 'alpha_y = "   << p.ay  << "' at graph 0.6, graph 0.5, graph 0.8\n";
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

void output_cauchy_gp(const char* foldername, const CauchyPoint& p) {
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

// helper integration method for output_meta_cauchy_gp
inline double sqr(double x) { return x*x; }
double integrate(const CauchyPoint& p) {
  // basic Riemann integration routine
  auto f = [&p](double x, double y) { return (1.0 + p.beta*x) / (1.0 + sqr(p.ax*(x - p.cx)) + sqr(p.ay*y)); };
  double sum, Iold, Inew=0.0;
  double x, y, dx, dy;
  unsigned n=1;
  do {
    sum = 0.0;
    n *= 2;
    Iold = Inew;
    dx = (p.dxds_max - p.dxds_min) / (n-1);
    dy = 2 * p.dyds_max / (n-1);
    for (unsigned i=0; i<n; i++) {
      x = p.dxds_min + i * dx;
      for (unsigned j=0; j<n; j++) {
        y = -p.dyds_max + j*dy;
        sum += f(x,y);
      }
    }
    Inew = sum * dx * dy;
  } while ((std::abs(Inew - Iold) > 1e-4) || n < 8);
  return Inew;
}


void output_meta_cauchy_gp(const char* foldername, const MetaFitResults& mf, double pc_out, double r) {
  // Writes a gnuplot file to plot the cauchy distribution
  // obtained from the meta fits alongside the observed data
  std::ofstream cauchy_gp;
  char cauchy_gp_name[100];
  sprintf(cauchy_gp_name, "%s/meta_E%0.2lf_r%0.3lf.gp", foldername, pc_out, r);
  cauchy_gp.open(cauchy_gp_name);

  CauchyPoint p;
  p.E = pc_out;
  p.r = r;
  // Detect if pc_out is in the 1D region or 2D region
  const auto& cfits = mf.cx.low_e_fits; // used for checking pc outs for all low_e_fits
  if (pc_out < mf.cx.low_e_fits.back().E) {
    // 1D -> find bounding 1D fits and average their parameters
    unsigned low_ix = 0, high_ix = 1;
    for (const auto& fit : cfits) {
      if (fit.E > pc_out) break;
      low_ix++, high_ix++;
    }
    double low_weight = (pc_out - cfits[low_ix].E) / (cfits[high_ix].E - cfits[low_ix].E);
    double high_weight = 1 - low_weight;
    p.cx   = low_weight * eval<fitType::CX, 1>(mf.cx.low_e_fits[low_ix], p)     + high_weight * eval<fitType::CX, 1>(mf.cx.low_e_fits[high_ix], p);
    p.ax   = low_weight * eval<fitType::AX, 1>(mf.ax.low_e_fits[low_ix], p)     + high_weight * eval<fitType::AX, 1>(mf.ax.low_e_fits[high_ix], p);
    p.ay   = low_weight * eval<fitType::AY, 1>(mf.ay.low_e_fits[low_ix], p)     + high_weight * eval<fitType::AY, 1>(mf.ay.low_e_fits[high_ix], p);
    p.beta = low_weight * eval<fitType::BETA, 1>(mf.beta.low_e_fits[low_ix], p) + high_weight * eval<fitType::BETA, 1>(mf.beta.low_e_fits[high_ix], p);
    p.dxds_min = low_weight * eval<fitType::DXDS_MIN, 1>(mf.dxds_min.low_e_fits[low_ix], p) + high_weight * eval<fitType::DXDS_MIN, 1>(mf.dxds_min.low_e_fits[high_ix], p);
    p.dxds_max = low_weight * eval<fitType::DXDS_MAX, 1>(mf.dxds_max.low_e_fits[low_ix], p) + high_weight * eval<fitType::DXDS_MAX, 1>(mf.dxds_max.low_e_fits[high_ix], p);
    p.dyds_max = low_weight * eval<fitType::DYDS_MAX, 1>(mf.dyds_max.low_e_fits[low_ix], p) + high_weight * eval<fitType::DYDS_MAX, 1>(mf.dyds_max.low_e_fits[high_ix], p);
  } else {
    // 2D -> just evaluate
    p.cx = eval<fitType::CX, 2>(mf.cx.high_e_fit, p);
    p.ax = eval<fitType::AX, 2>(mf.ax.high_e_fit, p);
    p.ay = eval<fitType::AY, 2>(mf.ay.high_e_fit, p);
    p.beta = eval<fitType::BETA, 2>(mf.beta.high_e_fit, p);
    p.dxds_min = eval<fitType::DXDS_MIN, 2>(mf.dxds_min.high_e_fit, p);
    p.dxds_max = eval<fitType::DXDS_MAX, 2>(mf.dxds_max.high_e_fit, p);
    p.dyds_max = eval<fitType::DYDS_MAX, 2>(mf.dyds_max.high_e_fit, p);
  }

  // Compute amp so that distribution is normalized
  p.amp = 1.0 / integrate(p);

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
void write_1d_files(const FitResults& fit, const std::vector<CauchyPoint>& data, const char* dir) {
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
  if (slice_num)
    gpfile << "splot for [i=0:" << slice_num-1 <<
      "] sprintf('" << prefix << "_1d_%d.dat', i) w lp title '1D slice fit'\n";
  else
    gpfile << "set label 'No 1D slice fits' at graph 0.75, graph 0.75, graph 0.9\n";
  gpfile.close();
  return;
}


template<fitType T>
void output_metafit_gp(const FitResults& fit, const std::vector<CauchyPoint>& data, const char* dir, double xpt) {
  std::ofstream gpfile;
  std::string fit_name = fit_to_string(T); // l-value created on the stack so it lives until the end of this function
  const char *prefix = fit_name.c_str();
  char filename[100];
  sprintf(filename, "%s/%s_2d.gp", dir, prefix);
  gpfile.open(filename);
  if (!gpfile.good()) {
    std::cerr << "ERROR: could not open " << filename << '\n';
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
  gpfile << ((fit.low_e_fits.size()) ? "rep" : "splot")
    << " 'coef.dat' u 'E':'r':'" << fit_name << "' w p title 'Cauchy data'\n";
  gpfile << "rep " << prefix << "_fit_trimmed(x,y) title '2D fit'";
  gpfile.close();
  return;
}

void output_master_gp(const char* dir) {
  char name[100];
  sprintf(name, "%s/master.gp", dir);
  std::ofstream gpfile;
  gpfile.open(name);
  if (gpfile.fail()) return;

  gpfile <<
R"x(call 'c_x_master.gp'
system('read')
call 'alpha_x_master.gp'
system('read')
call 'alpha_y_master.gp'
system('read')
call 'beta_master.gp'
system('read')
call 'dxds_min_master.gp'
system('read')
call 'dxds_max_master.gp'
system('read')
call 'dyds_max_master.gp')x";
  gpfile.close();
  return;
}


// Instantiate output_gp
template void output_metafit_gp<fitType::CX>(const FitResults&, const std::vector<CauchyPoint>&, const char*, double);
template void output_metafit_gp<fitType::AX>(const FitResults&, const std::vector<CauchyPoint>&, const char*, double);
template void output_metafit_gp<fitType::AY>(const FitResults&, const std::vector<CauchyPoint>&, const char*, double);
template void output_metafit_gp<fitType::BETA>(const FitResults&, const std::vector<CauchyPoint>&, const char*, double);
template void output_metafit_gp<fitType::DXDS_MIN>(const FitResults&, const std::vector<CauchyPoint>&, const char*, double);
template void output_metafit_gp<fitType::DXDS_MAX>(const FitResults&, const std::vector<CauchyPoint>&, const char*, double);
template void output_metafit_gp<fitType::DYDS_MAX>(const FitResults&, const std::vector<CauchyPoint>&, const char*, double);
