#pragma once
#include <vector>

#include "cauchy.hpp"

struct TableData {
  std::vector<double> pc_vals, r_vals, data;
};

void read_list_data(const char * folder, double pc_out, double r, XYBinnedData& data);
void read_list_data(const char* filename, XYBinnedData& data);
void read_table_data(const char *folder, double pc_in, double T, const char *param, TableData& table);
