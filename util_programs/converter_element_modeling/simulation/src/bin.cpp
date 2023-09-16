#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iostream>

#include "bin.hpp"


//ERBin::ERBin()
//  : count{0},
//  area{0},
//  p_list{},
//  xy_bins{} {}


void ERBin::add_point(const GeantParticle& p) {
  p_list.push_back(p);
  count++;
}


double ERBin::density(size_t num_elec_in) const {
  return count/(num_elec_in * area);
}


void ERBin::write_momenta(const char * output_dir,
    double E_elec, double T, double E, double r, std::ofstream& binlist) const {
  // Writes the dxds, dyds bins to a file, and records
  // the name of the written file to binlist
  if (!xy_bins) {
    std::cerr << "Warning: dxds, dyds have not been binned!\n";
    return;
  }
  auto& bins = *xy_bins;
  // Write binned data to a file
  std::ofstream binfile;
  char binfile_name[100];
  sprintf(binfile_name, "%s/dir_dat/E%0.0lf_T%0.3lf/E%0.2lf_r%0.3lf_bin.dat",
      output_dir, E_elec, T, E, r);
  binfile.open(binfile_name);
  for (auto& b : bins) {
    binfile << b.dxds << '\t'
      << b.dyds << '\t'
      << b.count / (p_list.size() * b.area) << '\t' // normalized to 1
      << b.avg_polx << '\t'
      << b.avg_poly << '\t'
      << b.avg_polz << '\n';
  }
  // Write the filename to binlist
  binlist << binfile_name << '\n';
}


void ERBin::bin_momenta() {
  // Bins the (dxds, dyds) data and writes the
  // binned data to a file, named using the
  // values of E and r
  if (p_list.size() == 0) return;

  size_t num_bins = 20;

  // XYBin{} is a zero-initialized XYBin
  xy_bins.emplace(std::vector<XYBin>(num_bins*num_bins, XYBin{}));
  auto& bins = *xy_bins;

  // Determine cutoff points for binning
  size_t ix_min = p_list.size() * 0.03;
  size_t ix_max = p_list.size() * 0.97;
  double dxds_min, dxds_max, dyds_min, dyds_max;

  std::sort(p_list.begin(), p_list.end(),
      [](const auto& p1, const auto& p2) { return p1.dxds < p2.dxds; });
  dxds_min = p_list[ix_min].dxds;
  dxds_max = p_list[ix_max].dxds;

  std::sort(p_list.begin(), p_list.end(),
      [](const auto& p1, const auto& p2) { return p1.dyds < p2.dyds; });
  dyds_min = p_list[ix_min].dyds;
  dyds_max = p_list[ix_max].dyds;

  // Remove out of range points
  auto outside_range = [dxds_min, dxds_max, dyds_min, dyds_max](const auto& p) -> bool {
    return !(p.dxds >= dxds_min && p.dxds < dxds_max
          && p.dyds >= dyds_min && p.dyds < dyds_max);
  };
  p_list.erase(std::remove_if(p_list.begin(), p_list.end(), outside_range), p_list.end());

  // Make sure we have different min/max values before continuing
  if (dxds_min == dxds_max) return;
  if (dyds_min == dyds_max) return;

  // Bin the (dxds,dyds) points
  double dxds_width = (dxds_max - dxds_min)/num_bins;
  double dyds_width = (dyds_max - dyds_min)/num_bins;
  for (const auto& p : p_list) {
    size_t x_ix = (p.dxds - dxds_min) / dxds_width;
    size_t y_ix = (p.dyds - dyds_min) / dyds_width;
    size_t bin_ix = x_ix + num_bins * y_ix;
    bins[bin_ix].count++;
    bins[bin_ix].avg_polx += p.polx;
    bins[bin_ix].avg_poly += p.poly;
    bins[bin_ix].avg_polz += p.polz;
  }

  // Divide avg_pol(x,y,z) by count to make them actually be the average
  for (auto & b : bins) {
    if (b.count) {
      b.avg_polx /= b.count;
      b.avg_poly /= b.count;
      b.avg_polz /= b.count;
    }
  }

  // Fill in bin dxds/dyds values
  for (size_t y_ix=0; y_ix<num_bins; y_ix++) {
    for (size_t x_ix=0; x_ix<num_bins; x_ix++) {
      double xval = dxds_min + (x_ix + 0.5) * dxds_width;
      double yval = dyds_min + (y_ix + 0.5) * dyds_width;
      size_t bin_ix = x_ix + num_bins * y_ix;
      bins[bin_ix].dxds = xval;
      bins[bin_ix].dyds = yval;
      bins[bin_ix].area = dxds_width * dyds_width;
    }
  }

  return;
}


double ERBin::get_average(const std::function<double(const GeantParticle&)>& property) const {
  // Takes the average of property(p) over all p in p_list
  // property is a function that takes a GeantParticle by const reference
  // and returns a double; for example,
  // [](const GeantParticle& p) { return p.polx; }, or
  //
  // [](const GeantParticle& p) { return p.polx * cos(p.theta); }, or
  if (p_list.size())
    return std::accumulate(p_list.begin(), p_list.end(), 0.0,
        [&](double sum, const auto& p) { return sum + property(p); }) / p_list.size();
  else
    return 0.0;
}

