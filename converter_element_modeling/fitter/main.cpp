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
namespace fs = std::filesystem;

inline constexpr double dp_diff(const DataPoint& p1, const DataPoint& p2) {
  return std::abs(p1.chi2/p1.amp - p2.chi2/p2.amp);
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
    coef_file << "E\tr\tc_x\ta_x\ta_y\tbeta\tamp\tdxds_min\tdxds_max\tdyds_max\tnpts\tchisq\n";

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

      // Record dxds_min, dxds_max, dyds_max
      double dxds_min = bins.front().x;
      double dxds_max = bins.back().x;
      double dyds_max = bins.back().y;
      // Do the cauchy fit
      auto [c_x, c_y, a_x, a_y, beta, amp, chi2] = asym_cauchy_fit(Eout/1e6, r/1e-2, bins);
      a_x = std::sqrt(a_x);
      a_y = std::sqrt(a_y);
      // Check for fit success
      if (chi2 == -1) {
        std::cout << "Asymmetric cauchy fit failed for " << bin_file_name << '\n';
      } else {
        // Record the obtained coefficients on success
        data_points.push_back({Eout, r, c_x, a_x, a_y, beta, dxds_min, dxds_max, dyds_max, amp, chi2, data_count});
        output_cauchy_gp(E_T_folder_name, data_points.back());
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
    //const double MAX_DIFF = 20.0;
    //    size_t data_len = data_points.size();
    //    std::vector<bool> good(data_points.size());
    //    std::fill(good.begin(), good.end(), true);
    //    double new_diff, old_diff = dp_diff(data_points[1], data_points[0]);
    //    size_t old_ix = 1; // last ix that was good
    //    for (size_t i=2; i<data_len; i++) {
    //      new_diff = dp_diff(data_points[i], data_points[old_ix]);
    //      if ((new_diff/old_diff > MAX_DIFF) && (data_points[i].E == data_points[old_ix].E)) {
    //        good[i] = false;
    //      } else {
    //        old_diff = new_diff;
    //        old_ix++;
    //      }
    //    }
    //    auto checker = [&good,i=0](const DataPoint& p) mutable { return !good[i++]; };
    //    data_points.erase(std::remove_if(data_points.begin(), data_points.end(), checker), data_points.end());


    // Do the meta fits
    auto cx_fit_results = fit_routine<fitType::CX>(data_points, crossover_point);
    auto ax_fit_results = fit_routine<fitType::AX>(data_points, crossover_point);
    auto ay_fit_results = fit_routine<fitType::AY>(data_points, crossover_point);
    auto beta_fit_results = fit_routine<fitType::BETA>(data_points, crossover_point);
    auto xmin_fit_results = fit_routine<fitType::DXDS_MIN>(data_points, crossover_point);
    auto xmax_fit_results = fit_routine<fitType::DXDS_MAX>(data_points, crossover_point);
    auto ymax_fit_results = fit_routine<fitType::DYDS_MAX>(data_points, crossover_point);

    // Write gnuplot files
    output_metafit_gp<fitType::CX>(cx_fit_results, data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::AX>(ax_fit_results, data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::AY>(ay_fit_results, data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::BETA>(beta_fit_results, data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::DXDS_MIN>(xmin_fit_results, data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::DXDS_MAX>(xmax_fit_results, data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::DYDS_MAX>(ymax_fit_results, data_points, E_T_folder_name, crossover_point);

    for (const auto& p : data_points) {
      coef_file << p.E << '\t' << p.r << '\t' << p.cx << '\t'
        << p.ax << '\t'       << p.ay << '\t'
        << p.beta << '\t'      << p.amp << '\t'
        << p.dxds_min << '\t' << p.dxds_max << '\t' << p.dyds_max << '\t'
        << p.npts << '\t'       << p.chi2 << '\n';
    }
    coef_file.close();

    // Store the results of the fit
    metafits.push_back({Ein, T, cx_fit_results,
        ax_fit_results,
        ay_fit_results,
        beta_fit_results,
        xmin_fit_results,
        xmax_fit_results,
        ymax_fit_results});
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
  size_t mf_ix=0;

  double T = metafits[0].T;
  bmad_file << "distribution = {\n";
  // TODO: read material from config file
  bmad_file << TAB<1>() << "material = tungsten" << ",\n";
  bmad_file << TAB<1>() << "species_out = positron" << ",\n";
  bmad_file << TAB<1>() << "thickness = " << T;
  for (const auto& mf : metafits) {
    if (mf.T != T) { // start a new distribution
      T = mf.T;
      bmad_file << "},\ndistribution = {\n";
      bmad_file << TAB<1>() << "thickness = " << mf.T << ",\n";
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
    output_bmad<fitType::CX>(mf.cx, bmad_file);
    output_bmad<fitType::AX>(mf.ax, bmad_file);
    output_bmad<fitType::AY>(mf.ay, bmad_file);
    output_bmad<fitType::BETA>(mf.beta, bmad_file);
    output_bmad<fitType::DXDS_MIN>(mf.dxds_min, bmad_file);
    output_bmad<fitType::DXDS_MAX>(mf.dxds_max, bmad_file);
    output_bmad<fitType::DYDS_MAX>(mf.dyds_max, bmad_file);

    bmad_file << TAB<2>() << "}\n" << TAB<1>() << "}";
  }
  bmad_file << "\n}\n";
  bmad_file.close();

  return 0;
}
