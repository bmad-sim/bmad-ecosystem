#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cctype>
#include <csignal>
#include <cstdlib>
#include <filesystem>

#include "GeantMain.hpp"
#include "parser.hpp"
#include "binner.hpp"
#include "cal_binner.hpp"

namespace fs = std::filesystem;

// Comment this out to use manual binning
#define AUTO_BIN

// Global state needed to support response to
// CTRL-C interupt signal
bool HALT_SIGNAL = false;

void signal_handler(int) {
  HALT_SIGNAL = true;
}


int main() {
  // First: Read in run settings from config file
  const auto settings = parse_config_file("config.txt");

  if (settings.valid()) {
    std::cout << settings;
  } else {
    std::cout << "Please correct the above issues in config.txt\n";
    return 1;
  }


  // Now, make sure output directory exists and is empty
  auto& od = settings.output_directory; // short var name
  if (!fs::is_directory(od)) {
    fs::create_directory(od);
  } else if (!fs::is_empty(od)) {
    std::cout << "Output directory \"" << od << "\" is not empty.\n"
      << "Remove contents before continuing? [yN] ";
    std::string response;
    std::cin >> response;
    // The following makes the response be lower case
    for (auto& c: response) c = std::tolower(c);

    if (response == "y" || response == "yes") {
      for (auto& p: fs::directory_iterator(od))
        fs::remove_all(p);
    } else {
      std::cout << "Please move/remove files from " << od
        << "before continuing\n";
      return 1;
    }
  }
  // Set up directory structure for dxds/dyds data
  for (auto E: settings.pc_in) {
    for (auto T: settings.target_thickness) {
      char dirname[100];
      sprintf(dirname, "%s/dir_dat/E%0.0lf_T%0.3lf", od.c_str(), E, T);
      fs::create_directories(dirname);
    }
  }

  // Next: Run simulation for each energy, thickness combination
  char outfile_name[100];
  char outfile_name_lowpc[100];

  std::ofstream eff_file;
  eff_file.open(od + "/yields.dat");
  eff_file << "E\tT\tY\n";
  auto runManager = Initialize_Geant();
  for (auto E: settings.pc_in) {
    for (auto T: settings.target_thickness) {
      std::cout << "Simulating target thickness T = "
        << T << "cm, incoming electron energy " << E << " MeV\n";
      std::ofstream outfile;
      std::ofstream outfile_lowpc;
      sprintf(outfile_name, "%s/E%0.0lf_T%0.3lf_er.dat", od.c_str(), E, T);
      sprintf(outfile_name_lowpc, "E%0.0lf_T%0.3lf_er_lowpc.dat", E, T);

      CalibrationBinner calbin(runManager, settings, E, T);
      auto [num_E_bins, num_r_bins] = calbin.calibrate();

      // Construct the real binner using settings obtained during calibration
      Binner my_binner(std::move(calbin));

      // Set up the handling of C-c interupt
      std::signal(SIGINT, signal_handler);
      ////////////

      // Main data collection loop
      int num_runs=0;
      constexpr size_t RUN_LENGTH = 10000;
      while (!my_binner.has_enough_data() && !HALT_SIGNAL) {
        run_simulation(runManager, settings.target_material, E, T, &my_binner, 10000);
        num_runs++;
#ifndef CONTINUOUS_PRINTOUT
      }
#endif
      size_t num_elec_in = num_runs * RUN_LENGTH;
      outfile.open(outfile_name);
      outfile_lowpc.open(outfile_name_lowpc);
      int pos_tot = my_binner.get_total();
      // set up first line with r values
      outfile << num_r_bins;
      for (size_t j=0; j<num_r_bins; j++)
        outfile << '\t' << my_binner.get_r_val(j);
      outfile << '\n';
      // each line should be one E value
      for (size_t i=0; i<num_E_bins; i++) {
        outfile << my_binner.get_E_val(i);
        for (size_t j=0; j<num_r_bins; j++) {
          double bin_E = my_binner.get_E_val(i);
          double bin_r = my_binner.get_r_val(j);
          const auto& bin = my_binner.get_bin(i,j);
          outfile << '\t' << bin.density(num_elec_in);
          bin.bin_momenta(od.c_str(), E, T, bin_E, bin_r);
        }
        outfile << '\n';
      }

      outfile_lowpc << "r\tpc_min\n";
      for (size_t j=0; j<num_r_bins; j++) {
        double bin_r = my_binner.get_r_val(j);
        auto low_pc = my_binner.lowest_pc_vals[j];
        outfile_lowpc << bin_r << '\t' << low_pc << '\n';
      }

#ifdef CONTINUOUS_PRINTOUT
      outfile.close();
      outfile_lowpc.close();
    }
#endif
      std::cout << "Efficiency: " << ((double) pos_tot) / (10000*num_runs) << '\n';
      eff_file << E << '\t' << T << '\t' << ((double) pos_tot) / (10000*num_runs) << '\n';
      if (HALT_SIGNAL) goto END;
    }
  }

END:
  eff_file.close();
  delete runManager;
  return 0;

}

