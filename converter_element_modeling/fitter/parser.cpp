#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cctype>

#include "parser.hpp"
using namespace std::string_literals;

SettingsSpec SimSettings::lookup(const std::string& name) {
  //if (name=="num_bins"s) {
  //  return { false, SetType::count, (void *) &(this->num_bins) };
  //} else if (name=="num_pc_bins"s) {
  if (name=="num_pc_bins"s) {
    return { false, SetType::count, (void *) &(this->num_pc_bins) };
  } else if (name=="num_r_bins"s) {
    return { false, SetType::count, (void *) &(this->num_r_bins) };
  } else if (name=="out_pc_min"s) {
    return { false, SetType::energy, (void *) &(this->out_pc_min) };
  } else if (name=="out_pc_max"s) {
    return { false, SetType::energy, (void *) &(this->out_pc_max) };
  } else if (name=="dxy_ds_max"s) {
    return { false, SetType::generic, (void *) &(this->dxy_ds_max) };
  } else if (name=="fit_crossover"s) {
    return { false, SetType::energy, (void *) &(this->fit_crossover) };
  } else if (name=="material"s) {
    return { false, SetType::name, (void *) &(this->target_material) };
  } else if (name=="output_directory"s) {
    return { false, SetType::name, (void *) &(this->output_directory) };
  } else if (name=="pc_in"s) {
    return { true, SetType::energy, (void *) &(this->pc_in) };
  } else if (name=="thicknesses"s) {
    return { true, SetType::length, (void *) &(this->target_thickness) };
  } else return {false, SetType::generic, nullptr};
}


bool SimSettings::valid() const {
  // Checks whether or not this is in a valid state

  bool bins_set = /* (num_bins > 0) || */ (num_pc_bins > 0 && num_r_bins > 0);
  bool out_pc_set = (out_pc_min >= 0) && (out_pc_min < out_pc_max);
  bool dxy_ds_set = dxy_ds_max > 0;
  bool fit_crossover_set = (fit_crossover >= 0);
  bool material_set = (target_material != ""s);
  bool output_dir_set = (output_directory != ""s);
  auto is_pos = [](const double x) { return x>0; };
  bool pc_in_set = (pc_in.size() != 0) &&
    std::all_of(pc_in.cbegin(), pc_in.cend(), is_pos);
  bool thickness_set = (target_thickness.size() != 0) &&
    std::all_of(target_thickness.cbegin(), target_thickness.cend(), is_pos);

  if (!bins_set)
    std::cerr << "ERROR: num_pc_bins and/or num_r_bins missing or invalid\n";
  if (!out_pc_set)
    std::cerr << "ERROR: out_pc_min and/or out_pc_max missing or invalid\n";
  if (!dxy_ds_set)
    std::cerr << "ERROR: dxy_ds_max missing or invalid\n";
  if (!fit_crossover_set)
    std::cerr << "ERROR: fit_crossover missing or invalid\n";
  if (!material_set)
    std::cerr << "ERROR: material missing or invalid\n";
  if (!output_dir_set)
    std::cerr << "ERROR: output_dir missing or invalid\n";
  if (!pc_in_set)
    std::cerr << "ERROR: pc_in missing or invalid\n";
  if (!thickness_set)
    std::cerr << "ERROR: thicknesses missing or invalid\n";

  return bins_set && out_pc_set && dxy_ds_set && fit_crossover_set
    && material_set && output_dir_set && pc_in_set && thickness_set;
}

bool operator==(const StrNum& s1, const StrNum& s2) {
  return (s1.unit==s2.unit) &&
    (s1.value==s2.value) &&
    (s1.value_str==s2.value_str) &&
    (s1.unit_str==s2.unit_str);
}


bool StrNum::parse() {
  // Attempts to parse unit_str and value_str
  // If successful, returns true and sets unit and
  // value accordingly.  Otherwise, returns false

  // Parse value
  try {
    double x = stod(value_str);
    value = x;
  } catch(...) {
    return false;
  }

  // Parse unit
  if (unit_str.size() == 0) return true;
  for (auto& c : unit_str) c = std::tolower(c);
  if (unit_str.back() == ',') unit_str.pop_back();
  if (unit_str == "ev") { unit = Unit::eV; }
  else if (unit_str == "kev") unit = Unit::keV;
  else if (unit_str == "mev") unit = Unit::MeV;
  else if (unit_str == "gev") unit = Unit::GeV;
  else if (unit_str == "tev") unit = Unit::TeV;
  else if (unit_str == "mm") unit = Unit::mm;
  else if (unit_str == "cm") unit = Unit::cm;
  else if (unit_str == "m") unit = Unit::m;
  else return false;
  return true;
}


double convert(const StrNum& v, const SetType T) {
  // converts v to MeV or cm depending on T
  switch (T) {
    case SetType::energy:
      switch (v.unit) {
        case Unit::none: // default to eV
        case Unit::eV: return v.value * 1e-6;
        case Unit::keV: return v.value * 1e-3;
        case Unit::MeV: return v.value;
        case Unit::GeV: return v.value * 1e3;
        default: return -1; //illegal unit
      }
    case SetType::length:
      switch (v.unit) {
        case Unit::none: // default to eV
        case Unit::m:  return v.value * 1e2;
        case Unit::cm: return v.value;
        case Unit::mm: return v.value * 1e-1;
        default: return -1; // illegal unit
      }
    default: // default to no conversion
      return v.value;
  }
}


SimSettings parse_config_file(const char* config_file_name) {
  // Parses the supplied config file and returns a
  // SimSettings object with the values set in the config file
  SimSettings settings;

  std::ifstream config_file(config_file_name);
  if (!config_file.is_open()) {
    std::cout << "Could not open " << config_file_name << '\n';
    return settings;
  }

  // Parse the config file line by line:
  std::string line;
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
      stream_line >> words.back();
    }

    auto print_error_message = [&words]() {
      std::cout << "Error parsing " << words[0] << " from config file:\n";
      for (const auto& x : words) std::cout << x << ' ';
      std::cout << '\n';
    };

    if (words.size() < 3) continue;
    auto spec = settings.lookup(words[0]);
    if (!spec.ele) continue;
    if (spec.is_vec) {
      auto& vec = * (std::vector<double> *) spec.ele;
      std::vector<StrNum> values(1);
      size_t words_used = 0;
      for (size_t i=2; i<words.size(); i++) {
        switch (words_used) {
          case 0: values.back().value_str = words[i]; break;
          case 1: values.back().unit_str = words[i]; break;
        }
        if (words[i].back() == ',') {
          values.push_back(StrNum());
          words_used = 0;
        } else {
          words_used++;
        }
      }
      // Remove the last element of values in case the last
      // character on the line was a comma
      if (values.back() == StrNum()) values.pop_back();
      for (auto& v : values) {
        auto success = v.parse();
        if (!success) print_error_message();
      }
      // Convert values in units to MeV/cm and add to settings
      for (const auto& v : values)
        vec.push_back(convert(v, spec.type));
    } else {
      switch (spec.type) {
        case SetType::name:
          * (std::string *) spec.ele = words[2];
          break;
        case SetType::count:
          try {
            * (size_t *) spec.ele = std::stoi(words[2]);
          } catch(...) {
            print_error_message();
          }
          break;
        default:
          StrNum value;
          value.value_str = words[2];
          if (spec.type != SetType::generic && words.size() > 3)
            value.unit_str = words[3];
          bool success = value.parse();
          if (!success)
            print_error_message();
          else
            * (double *) spec.ele = convert(value, spec.type);
      }
    }
  }

  return settings;
}

std::ostream& operator<<(std::ostream& out, const SimSettings& s) {
  out << "Minimum outgoing electron pc: " << s.out_pc_min << " MeV\n";
  out << "Maximum outgoing electron pc: " << s.out_pc_max << " MeV\n";
  out << "Maximum dx/ds, dy/ds magnitude: " << s.dxy_ds_max << '\n';
  out << "Target material: " << s.target_material << '\n';
  out << "Incoming electron pc";
  if (s.pc_in.size() > 1) {
    out<< "'s: {";
    for (unsigned i=0; i<s.pc_in.size()-1; i++) out << s.pc_in[i] << ", ";
    out << s.pc_in.back() << '}';
  } else {
    out << ": " << s.pc_in.front();
  }
  out << " MeV\n";
  out << "Converter thickness";
  if (s.target_thickness.size() > 1) {
    out << "es: {";
    for (unsigned i=0; i<s.target_thickness.size()-1; i++)
      out << s.target_thickness[i] << ", ";
    out << s.target_thickness.back() << '}';
  } else {
    out << ": " << s.target_thickness.back();
  }
  out << " cm\n";
  out << "Output directory: " << s.output_directory << '\n';
  return out;
}
