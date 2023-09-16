#pragma once
#include <vector>
#include <string>
#include <iostream>

enum class SetType {
  energy, length, count, name, generic
};

struct SettingsSpec {
  bool is_vec;
  SetType type;
  void * ele;
};

enum class Unit {
  eV, keV, MeV, GeV, TeV, mm, cm, m, none
};

struct StrNum {
  Unit unit;
  double value;
  std::string value_str;
  std::string unit_str;

  StrNum() : unit{Unit::none}, value{0.0}, value_str{}, unit_str{} {}

  bool parse();
};

bool operator==(const StrNum& s1, const StrNum& s2);

struct SimSettings {
  size_t /* num_bins,*/ num_pc_bins, num_r_bins;
  double out_pc_min, out_pc_max;
  double dxy_ds_max;
  double fit_crossover;
  std::string target_material;
  std::string output_directory;
  std::vector<double> pc_in;
  std::vector<double> target_thickness;
  std::vector<double> polarization_in;
  std::vector<double> pc_bin_points;
  std::vector<double> r_bin_points;
  std::vector<double> pc_bin_widths;
  std::vector<double> r_bin_widths;

  SimSettings() = default;
  //SimSettings() :  num_pc_bins{0},
  //                 num_r_bins{0},
  //                 out_pc_min{0.0},
  //                 out_pc_max{0.0},
  //                 dxy_ds_max{0.0},
  //                 fit_crossover{0.0},
  //                 target_material{},
  //                 output_directory{},
  //                 pc_in{},
  //                 target_thickness{},
  //                 polarization_in{} {}

  bool valid() const;
  SettingsSpec lookup(const std::string& name);
};



// Make SimSettings printable
std::ostream& operator<<(std::ostream& out, const SimSettings& s);

SimSettings parse_config_file(const char* config_file_name);
