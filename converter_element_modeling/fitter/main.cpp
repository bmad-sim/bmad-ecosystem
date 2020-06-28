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
  auto& data_dir = settings.output_directory;
  //std::string data_dir;
  //if (argc==1) {
  //  // No arguments -> read data directory from config.txt
  //  bool success = parse_config_file("config.txt", data_dir);
  //  if (!success) {
  //    return 1;
  //  }
  //} else {
  //  // 2+ arguments; assume argv[1] is the name of the data directory
  //  data_dir = argv[1];
  //}

  std::ofstream bmad_file;
  bmad_file.open(data_dir + "/converter.bmad");
  std::vector<MetaFitResults> metafits;

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
    ER_table er_table;
    read_er_data(data_dir.c_str(), Ein, T, er_table);
    // Vector for holding fitting results
    std::vector<CauchyPoint> data_points;
    data_points.reserve(er_table.probs.size());

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
      std::vector<BinPoint> bins;
      read_list_data(bin_filename, bins);
      if (!bins.size()) continue;

      CauchyPoint result;
      if (data_points.size()) // use previous fit results as initial param guess
        result = asym_cauchy_fit(Eout, r, bins, data_points.back());
      else
        result = asym_cauchy_fit(Eout, r, bins);
      // Check for fit success
      if (result.stat == CauchyStatus::NOPROG) {
        std::cerr << "Asymmetric cauchy fit failed for pc_out = " << Eout*1e6 << " MeV, r = " << r*1e2 << " cm\n";
      } else {
        // Record the obtained coefficients on success
        data_points.push_back(result);
        output_cauchy_gp(E_T_folder_name, result);
      }
    }






    //for (auto Eout : er_table.pc_vals) {
    //  //Eout *= 1e6; // convert to eV
    //  for (auto r : er_table.r_vals) {
    //    // Parse Eout, r from the file name
    //    //const char *bin_file_name = (bin_file.path()).c_str();
    //    //ret = sscanf(bin_file_name, file_name_format.c_str(), &Eout, &r);
    //    //if (ret!=2) continue; // parsing failure
    //    //Eout *= 1e6; // convert to eV
    //    //r *= 1e-2; // convert to meters

    //    // Read in the data from the file
    //    //char bin_file_name[200];
    //    //sprintf(bin_file_name, "%s/E%0.2lf_r%0.3lf_bin.dat", E_T_folder_name, Eout * 1e-6, r * 1e2);
    //    std::vector<BinPoint> bins;
    //    read_list_data(E_T_folder_name, Eout, r, bins);
    //    if (!bins.size()) {
    //      std::cerr << "WARNING: no data read for pc_out = " << Eout*1e-6 << " MeV, r = " << r * 1e2 << " cm\n";
    //      continue;
    //    }

    //    // Do the cauchy fit
    //    CauchyPoint result;
    //    if (data_points.size()) // use previous fit results as initial param guess
    //      result = asym_cauchy_fit(Eout, r, bins, data_points.back());
    //    else
    //      result = asym_cauchy_fit(Eout, r, bins);
    //    // Check for fit success
    //    if (result.stat == CauchyStatus::NOPROG) {
    //      std::cerr << "Asymmetric cauchy fit failed for pc_out = " << Eout*1e6 << " MeV, r = " << r*1e2 << " cm\n";
    //    } else {
    //      // Record the obtained coefficients on success
    //      data_points.push_back(result);
    //      output_cauchy_gp(E_T_folder_name, result);
    //    }
    //  }
    //}


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

    // With data_points sorted, assign amp variable to each CauchyPoint
    // from er_table for use with weighting meta fits
    size_t ix=0;
    for (auto& p : data_points) p.amp = er_table.probs[ix++];

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
        std::move(er_table));

    // Write cauchy gp filess with parameters from metafits
    //.auto print_status = [](const char * message, unsigned n, unsigned N) {
    //.  std::cout << "\033[K" << message << "(" << 100.0 * n / N << "%)\n\033[1A";
    //.};
    //.unsigned prog = 0;
    //.for (const auto& p : data_points) {
    //.  output_meta_cauchy_gp(E_T_folder_name, metafits.back(), p.E, p.r);
    //.  print_status("Outputing meta-fit cauchy files...", prog++, data_points.size());
    //.}
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
  bmad_file << TAB<1>() << "material = " << settings.target_material << ",\n";
  bmad_file << TAB<1>() << "species_out = positron" << ",\n";
  bmad_file << TAB<1>() << "thickness = " << T;
  for (const auto& mf : metafits) {
    if (mf.T != T) { // start a new distribution
      T = mf.T;
      bmad_file << "},\ndistribution = {\n";
      bmad_file << TAB<1>() << "material = tungsten" << ",\n";
      bmad_file << TAB<1>() << "species_out = positron" << ",\n";
      bmad_file << TAB<1>() << "thickness = " << mf.T;
    }
    bmad_file << ",\n" << TAB<1>() << "sub_distribution = {\n";
    bmad_file << TAB<2>() << "pc_in = " << mf.Ein << ",\n";
    bmad_file << TAB<2>() << "prob_pc_r = {\n";
    bmad_file << TAB<3>() << "r_values = [0.0, "; // fix r=0
    // Write r values to table (no comma after last one)
    auto num_rows = mf.table.pc_vals.size();
    auto num_cols = mf.table.r_vals.size();
    for (unsigned i=0; i<num_cols-1; i++) bmad_file << mf.table.r_vals[i] << ", ";
    bmad_file << mf.table.r_vals[num_cols-1] << "]";
    // Copy over each row of the table
    for (unsigned i=0; i<num_rows; i++) {
      bmad_file << ",\n" << TAB<3>() << "row = {pc_out = " << mf.table.pc_vals[i];
      bmad_file << ", prob = [0.0, "; // r=0 -> prob = 0
      for (unsigned j=0; j<num_cols-1; j++) bmad_file << 1e-4 * mf.table.probs[i*num_cols + j] << ", "; // convert to prob / (eV * m)
      bmad_file << 1e-4 * mf.table.probs[(i+1)*num_cols-1] << "]}";
    }
    bmad_file << "\n" << TAB<2>() << "},\n";

    // Output the fit parameters
    bmad_file << TAB<2>() << "direction_out = {\n";
    output_bmad<fitType::CX>(mf.cx, bmad_file, COMMA_AFTER);
    output_bmad<fitType::AX>(mf.ax, bmad_file, COMMA_AFTER);
    output_bmad<fitType::AY>(mf.ay, bmad_file, COMMA_AFTER);
    output_bmad<fitType::BETA>(mf.beta, bmad_file, COMMA_AFTER);
    output_bmad<fitType::DXDS_MIN>(mf.dxds_min, bmad_file, COMMA_AFTER);
    output_bmad<fitType::DXDS_MAX>(mf.dxds_max, bmad_file, COMMA_AFTER);
    output_bmad<fitType::DYDS_MAX>(mf.dyds_max, bmad_file, NO_COMMA);

    bmad_file << TAB<2>() << "}\n" << TAB<1>() << "}";
  }
  bmad_file << "\n}\n";
  bmad_file.close();

  return 0;
}
