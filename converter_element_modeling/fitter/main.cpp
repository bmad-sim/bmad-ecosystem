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
#include "gnuplot.hpp"
template<typename T>
inline constexpr T sqr(T x) { return x*x; }
template<typename T>
inline constexpr T cube(T x) { return x*x*x; }

std::string tab = "  ";
template<size_t n>
std::string TAB() { if constexpr (n==0) return ""; else return tab + TAB<n-1>(); }

namespace fs = std::filesystem;

double dp_diff(fitType T, const DataPoint& p1, const DataPoint& p2) {
  if      (T==fitType::CX) return std::abs(p1.cx - p2.cx);
  else if (T==fitType::AX) return std::abs(p1.ax - p2.ax);
  else if (T==fitType::AY) return std::abs(p1.ay - p2.ay);
  else                     return std::abs(p1.beta - p2.beta);
}

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

  char coef_file_name[100];
  std::ofstream coef_file;
  for (auto& E_T_folder : fs::directory_iterator(data_dir + "/dir_dat")) {
    double Ein, T, Eout, r;
    // Parse Ein, T out of the directory name
    const char *E_T_folder_name = (E_T_folder.path()).c_str();
    auto ret = sscanf(E_T_folder_name, folder_name_format.c_str(), &Ein, &T);
    if (ret != 2) continue; // parsing failure
    Ein *= 1e6; // convert to eV
    T *= 1e-2; // convert to meters
    // Vector for holding fitting results
    std::vector<DataPoint> data_points;
    data_points.reserve(200);

    sprintf(coef_file_name, "%s/coef.dat", E_T_folder_name);
    coef_file.open(coef_file_name);
    coef_file << "E\tr\tcx\tcy\tax\tay\tbeta\tamp\tchi2\n";

    // Loop over each file in the folder
    std::string file_name_format = std::string(E_T_folder_name) + "/E%lf_r%lf_bin.dat";
    for (auto& bin_file : fs::directory_iterator(E_T_folder.path())) {
      // Parse Eout, r from the file name
      const char *bin_file_name = (bin_file.path()).c_str();
      ret = sscanf(bin_file_name, file_name_format.c_str(), &Eout, &r);
      if (ret!=2) continue; // parsing failure
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


    // Remove outliers
    const double MAX_DIFF = 20.0;
    bool removed=true;
    while (removed) {
      removed = false;
      for (int fit_i=0; fit_i<4; fit_i++) {
        size_t data_len = data_points.size();
        std::vector<bool> good(data_points.size());
        std::fill(good.begin(), good.end(), true);
        fitType fitT = (fitType) fit_i;
        double new_diff, old_diff = dp_diff(fitT, data_points[1], data_points[0]);
        for (size_t i=2; i<data_len; i++) {
          new_diff = dp_diff(fitT, data_points[i], data_points[i-1]);
          if (((new_diff/old_diff > MAX_DIFF) || !good[i-1]) &&
              (data_points[i].E == data_points[i-1].E) &&
              (data_points[i].E > crossover_point)) {
            good[i] = false;
            removed = true;
          }
          old_diff = new_diff;
          //std::cout << "Checks: "
          //  << (new_diff/old_diff > MAX_DIFF || !good[i-1]) << ", "
          //  << (data_points[i].E == data_points[i-1].E) << ", "
          //  << (data_points[i].E > crossover_point) << '\n';
        }
        auto checker = [&good,i=0](const DataPoint& p) mutable { return !good[i++]; };
        data_points.erase(std::remove_if(data_points.begin(), data_points.end(), checker), data_points.end());
      }
    }
    data_points.erase(data_points.end()-1, data_points.end());
    //auto& pend1 = data_points[data_points.size()-1];
    //auto& pend2 = data_points[data_points.size()-2];
    //auto& pend3 = data_points[data_points.size()-3];
    //std::cout << "Checks at end: "
    //  << (dp_diff(fitType::AX, pend1, pend2)/dp_diff(fitType::AX, pend2, pend3) > MAX_DIFF ) << ", "
    //  << (pend1.E == pend2.E) << ", "
    //  << (pend1.E > crossover_point) << '\n';
    //std::cout << "Diff at end: " << dp_diff(fitType::AX, data_points[data_points.size()-1], data_points[data_points.size()-2]) << '\n';
    //getchar();

    // Do the meta fits
    auto cx_fit_results = cx_fit(data_points);
    auto ax_fit_results = ax_fit(data_points, crossover_point);
    auto ay_fit_results = ay_fit(data_points, crossover_point);
    auto beta_fit_results = beta_fit(data_points, crossover_point);

    // Write gnuplot files
    output_gp_cx(cx_fit_results, data_points, E_T_folder_name, crossover_point);
    output_gp_ax(ax_fit_results, data_points, E_T_folder_name, crossover_point);
    output_gp_ay(ay_fit_results, data_points, E_T_folder_name, crossover_point);
    output_gp_beta(beta_fit_results, data_points, E_T_folder_name, crossover_point);

    for (const auto& p : data_points) {
      coef_file << p.E << '\t' << p.r << '\t'
        << p.cx << '\t'       << p.cy << '\t'
        << p.ax << '\t'       << p.ay << '\t'
        << p.beta << '\t'      << p.amp << '\n';
    }
    coef_file.close();

    // Store the results of the fit
    metafits.push_back({Ein, T, cx_fit_results,
        ax_fit_results, ay_fit_results, beta_fit_results});
  }
  if (!metafits.size()) {
    std::cout << "Error: could not find data to fit to\n";
    return 1;
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
  // TODO: read material from config file
  bmad_file << TAB<1>() << "material = tungsten" << ",\n";
  bmad_file << TAB<1>() << "species_out = positron" << ",\n";
  bmad_file << TAB<1>() << "thickness = " << T << ",\n";
  // TODO: USE ACTUAL dxy_ds_max
  bmad_file << TAB<1>() << "dxy_ds_max = " << 1.8;
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
    bmad_file << TAB<3>() << "r_values = [0.0, "; // fix r=0
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
      bmad_file << ", prob = [0.0, "; // r=0 -> prob = 0
      for (size_t i=0; i<num_cols-1; i++) bmad_file << 1e-4 * probs[i] << ", "; // convert to prob / (eV * m)
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
        << 0.0 << ", "
        << beta1d.a << ", "
        << beta1d.b << ", "
        << beta1d.c << ", "
        << beta1d.d << "]},\n";
    }
    const auto& fit2b = mf.beta.high_e_fit;
    bmad_file << TAB<4>() << "poly_pc = [1.0, "
      << fit2b.a1 << ", "
      << fit2b.a2 << ", "
      << fit2b.a3 << "],\n";
    bmad_file << "        poly_r = ["
      //<< fit2b.b0 << ", "
      << 0.0 << ", "
      << fit2b.b1 << ", "
      << fit2b.b2 << ", "
      << fit2b.b3 << "]},\n";

    // cx
    bmad_file << TAB<3>() << "c_x = {\n";
    bmad_file << TAB<4>() << "poly_pc = [1.0, "
      << mf.cx.a1 << ", "
      << mf.cx.a2 << ", "
      << mf.cx.a3 << "],\n";
    bmad_file << TAB<4>() << "poly_r = ["
      //<< mf.cx.b0 << ", "
      << 0.0 << ", "
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
    bmad_file << TAB<4>() << "fit_2d_pc = {k = " << fit2ax.ke << ", poly = [1.0, "
      << fit2ax.a1 << ", "
      << fit2ax.a2 << ", "
      << fit2ax.a3 << "]},\n";
    bmad_file << TAB<4>() << "fit_2d_r = {k = " << fit2ax.kr << ", poly = ["
      << fit2ax.b0 << ", "
      << fit2ax.b1 << ", "
      << fit2ax.b2 << ", "
      << fit2ax.b3 << "]}},\n";

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
    bmad_file << TAB<4>() << "fit_2d_pc = {k = " << fit2ay.ke << ", poly = [1.0, "
      << fit2ay.a1 << ", "
      << fit2ay.a2 << ", "
      << fit2ay.a3 << "]},\n";
    bmad_file << TAB<4>() << "fit_2d_r = {k = " << fit2ax.kr << ", poly = ["
      << fit2ay.b0 << ", "
      << fit2ay.b1 << ", "
      << fit2ay.b2 << ", "
      << fit2ay.b3 << "]}}\n";

    bmad_file << TAB<2>() << "}\n" << TAB<1>() << "}";
  }
  bmad_file << "\n}\n";




















  return 0;
}
