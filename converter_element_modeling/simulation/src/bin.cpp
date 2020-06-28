#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <cstdio>

#include "bin.hpp"


Bin::Bin()
  : count(0),
  area(0),
  p_list{},
  binned_ptr(nullptr) {}

Bin::~Bin() { if (binned_ptr) delete binned_ptr; }

Bin::Bin(const Bin& o) : count(o.count), area(o.area), p_list(o.p_list) {
  if (o.binned_ptr) {
    binned_ptr = new std::vector<BinPoint>(*(o.binned_ptr));
  }
}

Bin& Bin::operator=(const Bin& o) {
  count = o.count;
  area = o.area;
  p_list = o.p_list;
  if (binned_ptr) delete binned_ptr;
  if (o.binned_ptr)
    binned_ptr = new std::vector<BinPoint>(*(o.binned_ptr));
  return *this;
}

void Bin::add_point(double dxds, double dyds) {
  p_list.push_back({dxds, dyds});
  count++;
}
double Bin::density(size_t num_elec_in) const {
  return count/(num_elec_in * area);
}

void Bin::write_momenta(const char * output_dir,
    double E_elec, double T, double E, double r, FILE *binlist) const {
  if (!binned_ptr) bin_momenta();
  if (!binned_ptr) return;
  auto& bins = *binned_ptr;
  // Write binned data to a file
  std::ofstream binfile;
  char binfile_name[100];
  sprintf(binfile_name, "%s/dir_dat/E%0.0lf_T%0.3lf/E%0.2lf_r%0.3lf_bin.dat",
      output_dir, E_elec, T, E, r);
  binfile.open(binfile_name);
  for (auto [bin_x, bin_y, bin_count, bin_area] : bins) {
    binfile << bin_x << '\t'
      << bin_y << '\t'
      << bin_count / (p_list.size() * bin_area) << '\n'; // normalized to 1
  }
  // Write the filename to binlist
  fprintf(binlist, "%s\n", binfile_name);
}

void Bin::binary_write_momenta(FILE *binary_file) const {
  if (!binned_ptr) bin_momenta();
  if (!binned_ptr) return;

  std::fwrite(binned_ptr->data(), sizeof(BinPoint), binned_ptr->size(), binary_file);
}

void Bin::bin_momenta() const {
  // Bins the (dxds, dyds) data and writes the
  // binned data to a file, named using the
  // values of E and r
  if (p_list.size() == 0) return;

  size_t num_bins = 20;

  if (binned_ptr) delete binned_ptr;
  binned_ptr = new std::vector<BinPoint>(num_bins*num_bins, {0,0,0,0});
  auto& bins = *binned_ptr;

  // Determine cutoff points for binning
  size_t ix_min = p_list.size() * 0.03;
  size_t ix_max = p_list.size() * 0.97;
  double dxds_min, dxds_max, dyds_min, dyds_max;

  std::sort(p_list.begin(), p_list.end(),
      [](xy_point p1, xy_point p2) { return p1.x < p2.x; });
  dxds_min = p_list[ix_min].x;
  dxds_max = p_list[ix_max].x;

  std::sort(p_list.begin(), p_list.end(),
      [](xy_point p1, xy_point p2) { return p1.y < p2.y; });
  dyds_min = p_list[ix_min].y;
  dyds_max = p_list[ix_max].y;

  // Remove out of range points
  auto outside_range = [dxds_min, dxds_max, dyds_min, dyds_max](xy_point p) -> bool {
    return !(p.x >= dxds_min && p.x<dxds_max && p.y>=dyds_min && p.y<dyds_max);
  };
  p_list.erase(std::remove_if(p_list.begin(), p_list.end(), outside_range), p_list.end());

  // Make sure we have different min/max values before continuing
  //if (p_list.size() <= 4) return 0;
  if (dxds_min == dxds_max) return;
  if (dyds_min == dyds_max) return;

  // Bin the (dxds,dyds) points
  double dxds_width = (dxds_max - dxds_min)/num_bins;
  double dyds_width = (dyds_max - dyds_min)/num_bins;
  for (const auto& p : p_list) {
    size_t x_ix = (p.x - dxds_min) / dxds_width;
    size_t y_ix = (p.y - dyds_min) / dyds_width;
    size_t bin_ix = x_ix + num_bins * y_ix;
    bins[bin_ix].count++;
  }

  // Fill in bin x/y values
  for (size_t y_ix=0; y_ix<num_bins; y_ix++) {
    for (size_t x_ix=0; x_ix<num_bins; x_ix++) {
      double xval = dxds_min + (x_ix + 0.5) * dxds_width;
      double yval = dyds_min + (y_ix + 0.5) * dyds_width;
      size_t bin_ix = x_ix + num_bins * y_ix;
      bins[bin_ix].x = xval;
      bins[bin_ix].y = yval;
      bins[bin_ix].area = dxds_width * dyds_width;
    }
  }

  return;
}



