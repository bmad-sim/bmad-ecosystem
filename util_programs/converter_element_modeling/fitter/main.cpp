#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <filesystem>

#include "parser.hpp"
#include "cauchy.hpp"
#include "read_data.hpp"
#include "meta_fit.hpp"
#include "gnuplot.hpp"
#include "bmad.hpp"
#include "chisq.hpp"
namespace fs = std::filesystem;


int main(int argc, char* argv[]) {
  // Given many dxds/dyds bin files of the form "E%lf_r%lf_bin.dat",
  // first fits cauchy distributions to each file containing a
  // sufficient number of data points, then fits various
  // functional forms to the obtained cauchy fit parameters

  // Step 0: parse config file
  const auto settings = parse_config_file("config.txt");
  if (settings.valid()) {
    std::cout << settings;
  } else {
    std::cout << "Please correct the above issues in config.txt\n";
    return 1;
  }
  const auto& data_dir = settings.output_directory;

  std::vector<MetaFitResults> metafits;

  std::vector<XYBinnedData> dxds_dyds_data;

  // Loop over each folder in dir_dat, running the fitting process on each
  std::string folder_name_format = data_dir + "/dir_dat/" + "E%lf_T%lf";

  char coef_file_name[100];
  std::ofstream coef_file;
  for (auto& E_T_folder : fs::directory_iterator(data_dir + "/dir_dat")) {
    double Ein, T;//, Eout, r;
    // Parse Ein, T out of the directory name
    const char *E_T_folder_name = (E_T_folder.path()).c_str();
    auto ret = sscanf(E_T_folder_name, folder_name_format.c_str(), &Ein, &T);
    if (ret != 2) continue; // parsing failure
    Ein *= 1e6; // convert to eV
    T *= 1e-2; // convert to meters
    // Read table of E/r probabilities
    TableData er_table;
    read_table_data(data_dir.c_str(), Ein, T, "er", er_table);
    // Read table of z polarizations
    TableData polz_table;
    read_table_data(data_dir.c_str(), Ein, T, "polz", polz_table);
    // Vector for holding fitting results
    std::vector<CauchyPoint> data_points;
    data_points.reserve(er_table.data.size());

    sprintf(coef_file_name, "%s/coef.dat", E_T_folder_name);
    coef_file.open(coef_file_name);
    coef_file << "E\tr\tc_x\talpha_x\talpha_y\tbeta\tamp\tdxds_min\tdxds_max\tdyds_max\tnpts\tchisq\n";

    // Loop over each file in the folder
    std::cout << "Running fits for pc_in = " << Ein*1e-6 << " MeV, T = " << T*1e2 << " cm...\n";
    std::ifstream xy_binlist;
    char list_filename[200];
    sprintf(list_filename, "%s/E%0.0lf_T%0.3lf_xy_bins.txt", data_dir.c_str(), Ein*1e-6, T*1e2);
    xy_binlist.open(list_filename);
    char bin_filename[300];
    std::string file_name_format = std::string(E_T_folder_name) + "/E%lf_r%lf_bin.dat";
    for(;;) {
      xy_binlist.getline(bin_filename, 300);
      if (!xy_binlist) break;
      if (!strlen(bin_filename)) continue;
      double Eout, r;
      ret = sscanf(bin_filename, file_name_format.c_str(), &Eout, &r);
      if (ret!=2) continue; // parsing failure
      Eout *= 1e6; // convert to eV
      r *= 1e-2; // convert to meters
      XYBinnedData xydata{};
      xydata.E = Eout;
      xydata.r = r;
      read_list_data(bin_filename, xydata);
      if (!xydata.bins.size()) continue;

      CauchyPoint result;
      if (data_points.size()) // use previous fit results as initial param guess
        result = asym_cauchy_fit(xydata, data_points.back());
      else
        result = asym_cauchy_fit(xydata);
      // Check for fit success
      if (result.stat == CauchyStatus::NOPROG) {
        std::cerr << "Asymmetric cauchy fit failed for pc_out = " << Eout*1e6 << " MeV, r = " << r*1e2 << " cm\n";
      } else {
        // Record the obtained coefficients on success
        data_points.push_back(result);
        output_cauchy_gp(E_T_folder_name, result);
        // Save the bin data for use in goodness-of-fit analysis later
        dxds_dyds_data.push_back(xydata);
      }
    }



    // Now that a cauchy fit has been performed for each file in this
    // directory, do the meta-fits
    double crossover_point = settings.fit_crossover * 1e6; // convert to eV
    // sort by E and then r
    std::sort(data_points.begin(), data_points.end(),
        [](const CauchyPoint& p1, const CauchyPoint& p2) {
          return p1.r < p2.r;
        });
    std::stable_sort(data_points.begin(), data_points.end(),
        [](const CauchyPoint& p1, const CauchyPoint& p2) {
          return p1.E < p2.E;
        });

    // Do the meta fits
    auto cx_fit_results   = fit_routine<fitType::CX>      (data_points, er_table, crossover_point);
    auto ax_fit_results   = fit_routine<fitType::AX>      (data_points, er_table, crossover_point);
    auto ay_fit_results   = fit_routine<fitType::AY>      (data_points, er_table, crossover_point);
    auto beta_fit_results = fit_routine<fitType::BETA>    (data_points, er_table, crossover_point);
    auto xmin_fit_results = fit_routine<fitType::DXDS_MIN>(data_points, er_table, crossover_point);
    auto xmax_fit_results = fit_routine<fitType::DXDS_MAX>(data_points, er_table, crossover_point);
    auto ymax_fit_results = fit_routine<fitType::DYDS_MAX>(data_points, er_table, crossover_point);

    // Write gnuplot files
    output_metafit_gp<fitType::CX>      (cx_fit_results,   data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::AX>      (ax_fit_results,   data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::AY>      (ay_fit_results,   data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::BETA>    (beta_fit_results, data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::DXDS_MIN>(xmin_fit_results, data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::DXDS_MAX>(xmax_fit_results, data_points, E_T_folder_name, crossover_point);
    output_metafit_gp<fitType::DYDS_MAX>(ymax_fit_results, data_points, E_T_folder_name, crossover_point);

    output_master_gp(E_T_folder_name);


    for (const auto& p : data_points) {
      coef_file << p.E << '\t' << p.r << '\t' << p.cx << '\t'
        << p.ax << '\t'       << p.ay << '\t'
        << p.beta << '\t'      << p.amp << '\t'
        << p.dxds_min << '\t' << p.dxds_max << '\t'
        << p.dyds_max << '\t' << p.chi2 << '\n';
    }
    coef_file.close();

    // Store the results of the fit
    metafits.emplace_back(Ein, T,
        std::move(cx_fit_results),
        std::move(ax_fit_results),
        std::move(ay_fit_results),
        std::move(beta_fit_results),
        std::move(xmin_fit_results),
        std::move(xmax_fit_results),
        std::move(ymax_fit_results),
        std::move(er_table),
        std::move(polz_table));

    // Chi-square analysis
    output_chisq(metafits.back(), data_points, dxds_dyds_data, settings);

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
  write_bmad_file(metafits, settings);

  return 0;
}
