#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <filesystem>

#include "parser.hpp"
#include "cauchy.hpp"
#include "read_data.hpp"
#include "meta_fit.hpp"
template<typename T>
inline constexpr T sqr(T x) { return x*x; }
template<typename T>
inline constexpr T cube(T x) { return x*x*x; }

std::string tab = "  ";
template<size_t n>
std::string TAB() { if constexpr (n==0) return ""; else return tab + TAB<n-1>(); }

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
  // Given many dxds/dyds bin files of the form "E%lf_r%lf_bin.dat",
  // first fits cauchy distributions to each file containing a
  // sufficient number of data points, then fits various
  // functional forms to the obtained cauchy fit parameters

  // Step 0: determine where the data directory is
  std::string data_dir;
  if (argc==1) {
    // No arguments -> read data directory from config.txt
    bool success = parse_config_file("config.txt", data_dir);
    if (!success) {
      return 1;
    }
  } else {
    // 2+ arguments; assume argv[1] is the name of the data directory
    data_dir = argv[1];
  }

  std::ofstream bmad_file;
  bmad_file.open(data_dir + "/converter.bmad");
  std::vector<MetaFitResults> metafits;


  // Loop over each folder in dir_dat, running the fitting process on each
  std::string folder_name_format = data_dir + "/dir_dat/" + "E%lf_T%lf";
  for (auto& E_T_folder : fs::directory_iterator(data_dir + "/dir_dat")) {
    double Ein, T, Eout, r;
    // Parse Ein, T out of the directory name
    const char *E_T_folder_name = (E_T_folder.path()).c_str();
    sscanf(E_T_folder_name, folder_name_format.c_str(), &Ein, &T);
    Ein *= 1e6; // convert to eV
    T *= 1e-2; // convert to meters
    // Vector for holding fitting results
    std::vector<DataPoint> data_points;
    data_points.reserve(200);

    // Loop over each file in the folder
    std::string file_name_format = std::string(E_T_folder_name) + "/E%lf_r%lf_bin.dat";
    for (auto& bin_file : fs::directory_iterator(E_T_folder.path())) {
      // Parse Eout, r from the file name
      const char *bin_file_name = (bin_file.path()).c_str();
      sscanf(bin_file_name, file_name_format.c_str(), &Eout, &r);
      Eout *= 1e6; // convert to eV
      r *= 1e-2; // convert to meters

      // Read in the data from the file
      std::vector<BinPoint> bins;
      read_list_data(bin_file_name, bins);

      // Skip this file if not enough data is present
      const size_t per_bin_threshold = 5.0;
      size_t data_count = std::accumulate(bins.begin(), bins.end(), 0,
          [](size_t count, const BinPoint& bin) { return count + bin.count; });
      if (data_count / bins.size() < per_bin_threshold)
        continue;

      // Do the cauchy fit
      auto [c_x, c_y, a_x, a_y, beta, amp, chi2] = asym_cauchy_fit(bins);
      // Check for fit success
      if (chi2 == -1) {
        std::cout << "Asymmetric cauchy fit failed for " << bin_file_name << '\n';
      } else {
        // Record the obtained coefficients on success
        data_points.push_back({Eout, r, c_x, c_y, a_x, a_y, beta, amp});
      }
    }

    // Now that a cauchy fit has been performed for each file in this
    // directory, do the meta-fits
    // TODO: read this from config file
    double crossover_point = 1e7; // in eV
    // sort by E and then r
    std::sort(data_points.begin(), data_points.end(),
        [](const DataPoint& p1, const DataPoint& p2) {
          return p1.r < p2.r;
        });
    std::stable_sort(data_points.begin(), data_points.end(),
        [](const DataPoint& p1, const DataPoint& p2) {
          return p1.E < p2.E;
        });

    // Do the meta fits
    auto cx_fit_results = cx_fit(data_points);
    auto ax_fit_results = ax_fit(data_points, crossover_point);
    auto ay_fit_results = ay_fit(data_points, crossover_point);
    auto beta_fit_results = beta_fit(data_points, crossover_point);

    // Store the results of the fit
    metafits.push_back({Ein, T, cx_fit_results,
        ax_fit_results, ay_fit_results, beta_fit_results});
  }
  // Sort fits by T and then by Ein
  std::sort(metafits.begin(), metafits.end(),
      [](const auto& f1, const auto& f2) {
          return f1.Ein < f2.Ein;
      });
  std::stable_sort(metafits.begin(), metafits.end(),
      [](const auto& f1, const auto& f2) {
          return f1.T < f2.T;
      });

  // Now that all fits have been performed and recorded, write the bmad file
  //size_t mf_ix=0;

  double T = metafits[0].T;
  bmad_file << "distribution = {\n";
  bmad_file << TAB<1>() << "thickness = " << T << ",\n";
  // TODO: USE ACTUAL dxy_ds_max
  bmad_file << TAB<1>() << "dxy_ds_max = " << 10;
  for (const auto& mf : metafits) {
    if (mf.T != T) { // start a new distribution
      T = mf.T;
      bmad_file << "},\ndistribution = {\n";
      bmad_file << TAB<1>() << "thickness = " << mf.T << ",\n";
      // TODO: USE ACTUAL dxy_ds_max
      bmad_file << TAB<1>() << "dxy_ds_max = " << 10;
    }
    bmad_file << ",\n" << TAB<1>() << "sub_distribution = {\n";
    bmad_file << TAB<2>() << "pc_in = " << mf.Ein << ",\n";
    bmad_file << TAB<2>() << "prob_pc_r = {\n";
    bmad_file << TAB<3>() << "r_values = [";
    // Open the ER file to write the E/r probability table
    std::ifstream er_file;
    char er_file_name[100];
    std::string er_name_format = data_dir + "/E%0.0lf_T%0.3lf_er.dat";
    sprintf(er_file_name, er_name_format.c_str(), mf.Ein*1e-6, mf.T*1e2);
    er_file.open(er_file_name);
    // Read r values from first row
    size_t num_cols;
    er_file >> num_cols;
    std::vector<double> r_values(num_cols);
    for (auto& r : r_values) er_file >> r;
    for (auto& r : r_values) r *= 1e-2; // convert to meters
    // Write r values to table (no comma after last one)
    for (size_t i=0; i<num_cols-1; i++) bmad_file << r_values[i] << ", ";
    bmad_file << r_values[num_cols-1] << "]";
    // Copy over each row of the table
    double Eout;
    std::vector<double> probs(num_cols);
    for (;;) {
      er_file >> Eout;
      if (er_file.fail()) break;
      Eout *= 1e6; // convert to eV
      for (auto& p : probs) er_file >> p;
      bmad_file << ",\n" << TAB<3>() << "row = {pc_out = " << Eout;
      bmad_file << ", prob = [";
      for (size_t i=0; i<num_cols-1; i++) bmad_file << probs[i] << ", ";
      bmad_file << probs[num_cols-1] << "]}";
    }
    bmad_file << "\n" << TAB<2>() << "},\n";
    er_file.close();

    // Output the fit parameters
    bmad_file << TAB<2>() << "direction_out = {\n";
    // beta
    bmad_file << TAB<3>() << "beta = {\n";
    for (const auto& beta1d : mf.beta.low_e_fits) {
      bmad_file << TAB<4>() << "fit_1d_r = {pc_out = " << beta1d.E
        << ", poly = ["
        << beta1d.a0 << ", "
        << beta1d.a1 << ", "
        << beta1d.a2 << ", "
        << beta1d.a3 << ", "
        << beta1d.a4 << "]},\n";
    }
    const auto& fit2b = mf.beta.high_e_fit;
    bmad_file << TAB<4>() << "poly_pc = ["
      << fit2b.a0 << ", "
      << fit2b.a1 << ", "
      << fit2b.a2 << "],\n";
    bmad_file << "        poly_r = ["
      << fit2b.b0 << ", "
      << fit2b.b1 << ", "
      << fit2b.b2 << ", "
      << fit2b.b3 << "]},\n";

    // cx
    bmad_file << TAB<3>() << "c_x = {\n";
    bmad_file << TAB<4>() << "poly_pc = ["
      << mf.cx.a0 << ", "
      << mf.cx.a1 << ", "
      << mf.cx.a2 << "],\n";
    bmad_file << TAB<4>() << "poly_r = ["
      << mf.cx.b0 << ", "
      << mf.cx.b1 << ", "
      << mf.cx.b2 << ", "
      << mf.cx.b3 << "]},\n";

    // ax
    bmad_file << TAB<3>() << "alpha_x = {\n";
    for (const auto& ax1d : mf.ax.low_e_fits) {
      bmad_file << TAB<4>() << "fit_1d_r = {pc_out = " << ax1d.E
        << ", k = " << ax1d.k << ", " << "poly = ["
        << ax1d.a << ", "
        << ax1d.b << ", "
        << ax1d.c << ", "
        << ax1d.d << "]},\n";
    }
    const auto& fit2ax = mf.ax.high_e_fit;
    bmad_file << TAB<4>() << "fit_2d_pc = {k = " << fit2ax.ke << ", poly = ["
      << fit2ax.ae << ", "
      << fit2ax.be << ", "
      << fit2ax.ce << ", "
      << fit2ax.de << "]},\n";
    bmad_file << TAB<4>() << "fit_2d_r = {k = " << fit2ax.kr << ", poly = ["
      << fit2ax.ar << ", "
      << fit2ax.br << ", "
      << fit2ax.cr << ", "
      << fit2ax.dr << "]}},\n";

    // ay
    bmad_file << TAB<3>() << "alpha_y = {\n";
    for (const auto& ax1d : mf.ay.low_e_fits) {
      bmad_file << TAB<4>() << "fit_1d_r = {pc_out = " << ax1d.E
        << ", k = " << ax1d.k << ", " << "poly = ["
        << ax1d.a << ", "
        << ax1d.b << ", "
        << ax1d.c << ", "
        << ax1d.d << "]},\n";
    }
    const auto& fit2ay = mf.ay.high_e_fit;
    bmad_file << TAB<4>() << "fit_2d_pc = {k = " << fit2ay.ke << ", poly = ["
      << fit2ay.ae << ", "
      << fit2ay.be << ", "
      << fit2ay.ce << ", "
      << fit2ay.de << "]},\n";
    bmad_file << TAB<4>() << "fit_2d_r = {k = " << fit2ax.kr << ", poly = ["
      << fit2ay.ar << ", "
      << fit2ay.br << ", "
      << fit2ay.cr << ", "
      << fit2ay.dr << "]}}\n";

    bmad_file << TAB<2>() << "}\n" << TAB<1>() << "}";
  }
  bmad_file << "\n}\n";




















  return 0;
}
