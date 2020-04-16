#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cctype>

#include "parser.hpp"

bool parse_config_file(const char* config_file_name,
    SimSettings& settings) {
  // Parses the supplied config file for
  // Target material,
  // E- values, and
  // Target thicknesses which should be used in the simulation
  // Mandatory settings
  bool material_set = false, energy_set = false, thickness_set = false, out_energy_max_set = false, out_r_max_set = false;
  // Optional settings
  bool output_dir_set = false, out_energy_min_set = false, dxy_ds_max_set = false;
  std::ifstream config_file;
  config_file.open(config_file_name);

  std::string line;
  if (!config_file.is_open()) {
    std::cout << "Could not open " << config_file_name << '\n';
    return false;
  }

  // Parse the config file line by line:
  while (!config_file.eof()) {
    std::getline(config_file, line);
    // Remove leading whitespace
    while (isspace(line[0])) line.erase(0,1);
    // Remove comments (started with !)
    auto comment_start = std::find(line.begin(), line.end(), '!');
    if (comment_start != line.end()) line.erase(comment_start, line.end());
    if (!line.size()) continue;

    std::stringstream stream_line(line);
    std::vector<std::string> words;
    while (stream_line.good()) {
      words.push_back("");
      stream_line >> *(words.end() - 1);
    }

    if (words.size() < 3) continue;
    if (words[0] == "out_energy_min") {
      try {
        double Emin = stod(words[2]);
        settings.out_energy_min = Emin;
        out_energy_min_set = true;
      } catch (...) {
        std::cout << "Error parsing out_energy_min from config file:\n";
        for (const auto& x : words) std::cout << x << ' ';
        std::cout << '\n';
      }
    } else if (words[0] == "out_energy_max") {
      try {
        double Emax = stod(words[2]);
        settings.out_energy_max = Emax;
        out_energy_max_set = true;
      } catch (...) {
        std::cout << "Error parsing out_energy_max from config file:\n";
        for (const auto& x : words) std::cout << x << ' ';
        std::cout << '\n';
      }
    } else if (words[0] == "dxy_ds_max") {
      try {
        double dxy_ds_max = stod(words[2]);
        settings.dxy_ds_max = dxy_ds_max;
        dxy_ds_max_set = true;
      } catch (...) {
        std::cout << "Error parsing dxy_ds_max from config file:\n";
        for (const auto& x : words) std::cout << x << ' ';
        std::cout << '\n';
      }
    } else if (words[0] == "out_r_max") {
      try {
        double out_r_max = stod(words[2]);
        settings.out_r_max = out_r_max;
        out_r_max_set = true;
      } catch (...) {
        std::cout << "Error parsing out_r_max from config file:\n";
        for (const auto& x : words) std::cout << x << ' ';
        std::cout << '\n';
      }
    } else if (words[0] == "output_directory") {
      settings.output_directory = words[2];
      output_dir_set = true;
      if (words.size() > 3) {
        std::cout << "Warning: output directory name trimmed from \n\"";
        for (size_t i=2; i<words.size(); i++)
          std::cout << words[i] << ' ';
        std::cout << "\"\n\tto " << words[2] << '\n';
      }
    } else if (words[0] == "target_material") {
      settings.target_material = words[2];
      material_set = true;
    } else if (words[0] == "in_energies") {
      double E;
      std::for_each(words.begin()+2, words.end(), [&](std::string w) {
          if (w.length() > 0) {
            try {
              E = stod(w);
              settings.in_energies.push_back(E);
              energy_set = true;
            } catch(...) {
              std::cout << "Error parsing in_energies from config file:\n";
              for (const auto& x : words) std::cout << x << ' ';
              std::cout << '\n';
            }
          }
          return;
        } );
    } else if (words[0] == "target_thicknesses") {
      double T;
      std::for_each(words.begin()+2, words.end(), [&](std::string w) {
          if (w.length() > 0) {
            try {
              T = stod(w);
              settings.target_thicknesses.push_back(T);
              thickness_set = true;
            } catch(...) {
              std::cout << "Error parsing in_energies from config file:\n";
              for (const auto& x : words) std::cout << x << ' ';
              std::cout << '\n';
            }
          }
          return;
        } );
    }
  }
  // Check that all settings have been specified
  if (!out_energy_min_set) {
    std::cout << "Minimum outgoing positron energy not specified in " << config_file_name << '\n';
    std::cout << "Defaulting to " << settings.out_energy_min << '\n';
  }
  if (!out_energy_max_set) {
    std::cout << "Maximum outgoing positron energy not specified in " << config_file_name << '\n';
  }
  if (!out_r_max_set) {
    std::cout << "Maximum positron radial displacement not specified in " << config_file_name << '\n';
  }
  if (!dxy_ds_max_set) {
    std::cout << "Maximum dx'/ds, dy'/ds magnitude not specified in " << config_file_name << '\n';
    std::cout << "Defaulting to " << settings.dxy_ds_max << '\n';
  }
  if (!output_dir_set) {
    std::cout << "Output directory not specified in " << config_file_name << '\n';
    std::cout << "Defaulting to " << settings.output_directory << '\n';
  }
  if (!material_set) std::cout << "Target material not specified in " << config_file_name << '\n';
  if (!energy_set) std::cout << "Incoming electron energy not specified in " << config_file_name << '\n';
  if (!thickness_set) std::cout << "Target thicknesses not specified in " << config_file_name << '\n';

  return (material_set && energy_set && thickness_set && out_energy_max_set && out_r_max_set);
}

std::ostream& operator<<(std::ostream& out, const SimSettings& s) {
  out << "Minimum outgoing electron energy: " << s.out_energy_min << " MeV\n";
  out << "Maximum outgoing electron energy: " << s.out_energy_max << " MeV\n";
  out << "Maximum dx'/ds, dy'/ds magnitude: " << s.dxy_ds_max << '\n';
  out << "Target material: " << s.target_material << '\n';
  out << "Incoming electron energies: {";
  for (auto E : s.in_energies) out << E << ", ";
  out << "} MeV\n";
  out << "Target thicknesses: {";
  for (auto T : s.target_thicknesses) out << T << ", ";
  out << "} cm\n";
  out << "Output directory: {" << s.output_directory << '\n';
  return out;
}
