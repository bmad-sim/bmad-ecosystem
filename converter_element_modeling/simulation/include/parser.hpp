#pragma once
#include <vector>
#include <string>
#include <iostream>
struct SimSettings {
  double out_energy_min;
  double out_energy_max;
  double out_r_max;
  double dxy_ds_max;
  std::string target_material;
  std::string output_directory;
  std::vector<double> in_energies;
  std::vector<double> target_thicknesses;

  SimSettings() : out_energy_min{0},
                   out_energy_max{100},
                   out_r_max{0.2},
                   dxy_ds_max{10},
                   target_material{"tungsten"},
                   output_directory{"sim_data"},
                   in_energies{},
                   target_thicknesses{} {}
};

// Make SimSettings printable
std::ostream& operator<<(std::ostream& out, const SimSettings& s);

bool parse_config_file(const char* config_file_name,
    SimSettings& settings);
