#pragma once
#include <vector>

#include "cauchy.hpp"

struct ER_table {
  std::vector<double> pc_vals, r_vals, probs;
};

void read_list_data(const char * folder, double pc_out, double r, std::vector<BinPoint>& bins);
void read_er_data(const char *folder, double pc_in, double T, ER_table& table);
