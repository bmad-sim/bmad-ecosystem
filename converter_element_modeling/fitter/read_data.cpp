#include <iostream>
#include <fstream>
#include <vector>

#include "read_data.hpp"

void read_list_data(const char * filename, std::vector<BinPoint>& bins) {
  std::ifstream datafile;
  datafile.open(filename);
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

