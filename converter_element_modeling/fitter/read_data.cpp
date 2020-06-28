#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "read_data.hpp"

bool fuzzy_file_open(std::ifstream& file, const char* folder, double pc_out, double r) {
  // Fuzzes the floating point values pc_out and r to
  // try and find the correct file to open
  char name[200];
  double pc_offsets[] = {0.0, -0.01, 0.01};
  double r_offsets[] = {0.0, -0.001, 0.001};
  bool exact = true;
  for (auto dp : pc_offsets) {
    for (auto dr : r_offsets) {
      sprintf(name, "%s/E%0.2lf_r%0.3lf_bin.dat", folder, pc_out*1e-6+dp, r*1e2+dr);
      file.open(name);
      if (!file.fail()) {
        if (!exact)
          std::cout << "Note: reading file " << name
            << " for pc_out = " << pc_out*1e-6 << " MeV, r = "
            << r * 1e2 << " cm\n";
        return true;
      }
      file.close();
      file.clear();
      exact = false;
    }
  }
  return false;
}

void read_list_data(const char * folder, double pc_out, double r, std::vector<BinPoint>& bins) {
  std::ifstream datafile;
  if (!fuzzy_file_open(datafile, folder, pc_out, r)) {
    std::cerr << "WARNING: could not open data file in directory " << folder
      << " for pc_out = " << pc_out << " MeV, r = " << r << " cm\n";
    return;
  }
  double x, y, count;
  for (;;) {
    datafile >> x;
    datafile >> y;
    datafile >> count;
    if (datafile.fail()) break;
    bins.push_back({x, y, count});
  }
  datafile.close();
  return;
}

void read_list_data(const char* filename, std::vector<BinPoint>& bins) {
  std::ifstream datafile;
  datafile.open(filename);
  if (!datafile) {
    std::cerr << "WARNING: could not open " << filename << "\n";
    return;
  }
  double x, y, count;
  for (;;) {
    datafile >> x;
    datafile >> y;
    datafile >> count;
    if (datafile.fail()) break;
    bins.push_back({x, y, count});
  }
  datafile.close();
  return;
}

void read_er_data(const char* folder, double pc_in, double T, ER_table& table) {
  // Reads the data in the E-r file for the given pc_in and target thickness,
  // and stores the data in the ER_table struct
  // probs is a 1D vector to avoid extraneous new/delete
  // when resizing.  Values are stored in probs in the order in which
  // they are read from the file (left to right, top to bottom)

  // Setup
  table.pc_vals.clear();
  table.r_vals.clear();
  table.probs.clear();

  // Open the file
  char filename[200];
  sprintf(filename, "%s/E%0.0lf_T%0.3lf_er.dat", folder, pc_in*1e-6, T*1e2);
  std::ifstream er_file;
  er_file.open(filename);
  if (er_file.fail()) return;

  // Read the r values from the first row
  unsigned num_cols;
  er_file >> num_cols;
  table.r_vals.resize(num_cols);
  for (auto& r : table.r_vals) er_file >> r;
  for (auto& r : table.r_vals) r *= 1e-2; // convert to m

  table.pc_vals.reserve(num_cols);
  table.probs.reserve(num_cols*num_cols);

  // Read a row at a time from the file
  double pc_out;
  for (;;) {
    er_file >> pc_out;
    if (er_file.eof() || er_file.fail()) break;
    pc_out *= 1e6; // convert to eV
    table.pc_vals.push_back(pc_out);
    for (const auto& r : table.r_vals) { // run loop for exactly one row
      table.probs.emplace_back();
      er_file >> table.probs.back();
    }
  }

  er_file.close();
  return;
}







