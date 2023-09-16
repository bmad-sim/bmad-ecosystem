#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
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

void read_list_data_impl(std::ifstream& datafile, XYBinnedData& data) {
  // Check first line to see how many columns there are
  std::string line1;
  std::getline(datafile, line1);
  line1.push_back('\n'); // keep the string stream state as "good" after extracting the last value
  std::istringstream line1_stream(line1);
  double junk;
  int num_cols = 0;
  for (;;) {
    line1_stream >> junk;
    if (!line1_stream.good()) break;
    num_cols++;
  }

  // Read the subsequent lines
  std::vector<double> row;
  row.resize(num_cols);
  // Don't forget the line that was extracted above
  datafile.seekg(0);
  for (;;) {
    for (int i=0; i<num_cols; i++)
      datafile >> row[i];
    if (datafile.fail()) break;
    data.bins.push_back({row[0], row[1], row[2]});
  }
  datafile.close();
}

void read_list_data(const char * folder, double pc_out, double r, XYBinnedData& data) {
  std::ifstream datafile;
  if (!fuzzy_file_open(datafile, folder, pc_out, r)) {
    std::cerr << "WARNING: could not open data file in directory " << folder
      << " for pc_out = " << pc_out << " MeV, r = " << r << " cm\n";
    return;
  }
  read_list_data_impl(datafile, data);
  return;
}

void read_list_data(const char* filename, XYBinnedData& data) {
  std::ifstream datafile;
  datafile.open(filename);
  if (!datafile) {
    std::cerr << "WARNING: could not open " << filename << "\n";
    return;
  }
  read_list_data_impl(datafile, data);
  return;
}

void read_table_data(const char* folder, double pc_in, double T, const char *param, TableData& table) {
  // Reads the data in the specified table file for the given pc_in and target thickness,
  // and stores the data in the TableData struct
  // data is a 1D vector to avoid extraneous new/delete
  // when resizing.  Values are stored in data in the order in which
  // they are read from the file (left to right, top to bottom)

  // Setup
  table.pc_vals.clear();
  table.r_vals.clear();
  table.data.clear();

  // Open the file
  char filename[200];
  sprintf(filename, "%s/E%0.0lf_T%0.3lf_%s.dat", folder, pc_in*1e-6, T*1e2, param);
  std::ifstream table_file;
  table_file.open(filename);
  if (table_file.fail()) return;

  // Read the r values from the first row
  unsigned num_cols;
  table_file >> num_cols;
  table.r_vals.resize(num_cols);
  for (auto& r : table.r_vals) table_file >> r;
  for (auto& r : table.r_vals) r *= 1e-2; // convert to m

  table.pc_vals.reserve(num_cols);
  table.data.reserve(num_cols*num_cols);

  // Read a row at a time from the file
  double pc_out;
  for (;;) {
    table_file >> pc_out;
    if (table_file.eof() || table_file.fail()) break;
    pc_out *= 1e6; // convert to eV
    table.pc_vals.push_back(pc_out);
    for (const auto& r : table.r_vals) { // run loop for exactly one row
      table.data.emplace_back();
      table_file >> table.data.back();
    }
  }

  return;
}







