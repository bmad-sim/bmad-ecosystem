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
#include "bin.hpp"
#include "cal_binner.hpp"
#include "point_cache.hpp"

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
  auto runManager = Initialize_Geant();
  // Make a PointCache for runManager to use
  PointCache point_cache;
  for (auto E: settings.pc_in) {
    for (auto T: settings.target_thickness) {
      std::cout << "Simulating target thickness T = "
        << T << "cm, incoming electron energy " << E << " MeV\n";

      CalibrationBinner calbin(runManager, &point_cache, settings, E, T);
      calbin.calibrate();


      // Construct the real binner using settings obtained during calibration
      //Binner my_binner(std::move(calbin));
      //my_binner.set_point_cache(point_cache);
      auto& my_binner = calbin;

      // Set up the handling of C-c interupt
      std::signal(SIGINT, signal_handler);
      ////////////

      // Main data collection loop
      constexpr size_t RUN_LENGTH = 4000;
      while (!my_binner.has_enough_data() && !HALT_SIGNAL) {
        calbin.run(RUN_LENGTH);
      }
      calbin.write_data();

      if (HALT_SIGNAL) goto END;
    }
  }

END:
  delete runManager;
  return 0;

}

