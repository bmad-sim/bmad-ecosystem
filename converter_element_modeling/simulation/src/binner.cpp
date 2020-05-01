#include <vector>
#include <utility>
#include <tuple>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>
#include "binner.hpp"
#include "GeantMain.hpp"

const double PI = 4.0 * atan(1.0);

inline double sqr(double x) { return x*x; }

Bin::Bin()
  : count(0),
  total_count(nullptr),
  area(0),
  dxds_tot(0),
  dyds_tot(0),
  dxds_square_tot(0),
  dyds_square_tot(0),
  p_list{} {}

Bin::~Bin() {}

void Bin::add_point(double dxds, double dyds) {
  p_list.push_back({dxds, dyds});
  dxds_tot += dxds;
  dyds_tot += dyds;
  dxds_square_tot += dxds*dxds;
  dyds_square_tot += dyds*dyds;
  count++;
}

void Bin::bin_momenta(const char * output_dir,
    double E_elec, double T, double E, double r) const {
  // Bins the (dxds, dyds) data and writes the
  // binned data to a file, named using the
  // values of E and r
  if (p_list.size() == 0) return;

  size_t num_bins;
  size_t num_bins_min = 7;
  size_t num_bins_max = 20;
  double per_bin_target = 5.0;
  size_t target_num_bins = (size_t) std::sqrt(count / per_bin_target);
  if (target_num_bins < num_bins_min) num_bins = num_bins_min;
  else if (target_num_bins > num_bins_max) num_bins = num_bins_max;
  else num_bins = target_num_bins;

  std::vector<BinPoint> bins(num_bins*num_bins, {0,0,0});

  // Determine cutoff points for binning
  size_t ix_min = p_list.size() * 0.025;
  size_t ix_max = p_list.size() * 0.975;
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
    }
  }

  // Write binned data to a file
  std::ofstream binfile;
  char binfile_name[100];
  sprintf(binfile_name, "%s/dir_dat/E%0.0lf_T%0.3lf/E%0.2lf_r%0.3lf_bin.dat",
      output_dir, E_elec, T, E, r);
  binfile.open(binfile_name);
  for (auto [bin_x, bin_y, bin_count] : bins) {
    binfile << bin_x << '\t'
      << bin_y << '\t'
      << bin_count << '\n';
  }
  binfile.close();

  return;
}


double Bin::mu_dxds() const {
  return dxds_tot/count;
}

double Bin::mu_dyds() const {
  return dyds_tot/count;
}

double Bin::sd_dxds() const {
  return sqrt(dxds_square_tot/count - (mu_dxds())*(mu_dxds()));
}

double Bin::sd_dyds() const {
  return sqrt(dyds_square_tot/count - (mu_dyds())*(mu_dyds()));
}

double Bin::density(size_t num_elec_in) const {
  //return count/(*total_count * area);
  return count/(num_elec_in * area);
}

void Bin::clear() {
  // Returns bin to its default (empty) state
  count=0;
  dxds_tot=0;
  dyds_tot=0;
  dxds_square_tot=0;
  dyds_square_tot=0;
}

// Pure virtual destructors must have a body
BinnerBase::~BinnerBase() {}

Binner::Binner() : Binner(1,1, 30, 0.2) {}

Binner::Binner(int N_E_bins, int N_r_bins) : Binner(N_E_bins, N_r_bins, 30, 0.2) {}

Binner::Binner(int N_E_bins, int N_r_bins, double p_E_max, double p_r_max)
  : Binner(N_E_bins, N_r_bins, 0.0, 0.0, p_E_max, p_r_max) {}

Binner::Binner(int N_E_bins, int N_r_bins, double p_E_min, double p_r_min, double p_E_max, double p_r_max)
  : E_min(p_E_min), r_min(p_r_min),
  E_max(p_E_max), r_max(p_r_max),
  dxds_bound(10), dyds_bound(10),
  E_edges(N_E_bins+1), r_edges(N_r_bins+1),
  bin_counts(N_E_bins, std::vector<Bin>(N_r_bins)),
  lowest_pc_vals(N_r_bins, 100),
  total_count(0) {
    std::iota(E_edges.begin(), E_edges.end(), 0.0);
    std::iota(r_edges.begin(), r_edges.end(), 0.0);
    for (auto& edge : E_edges) edge *= (E_max-E_min)/(N_E_bins+1);
    for (auto& edge : r_edges) edge *= (r_max-r_min)/(N_r_bins+1);
    size_t E_ix, r_ix;
    E_ix = 0;
    for (auto& row : bin_counts) {
      r_ix = 0;
      for (auto& bin : row) {
        bin.total_count = &total_count;
        bin.area = (E_edges[E_ix+1] - E_edges[E_ix])
          * (r_edges[r_ix+1] - r_edges[r_ix]);
        r_ix++;
      }
      E_ix++;
    }
}

Binner::~Binner() {}


std::pair<int, int> Binner::get_bin_num(double E, double r) {
  auto E_it = std::find_if(E_edges.begin(), E_edges.end(),
      [E](double edge) { return edge >= E; });
  size_t n_E = std::distance(E_edges.begin(), E_it)-1;
  auto r_it = std::find_if(r_edges.begin(), r_edges.end(),
      [r](double edge) { return edge >= r; });
  size_t n_r = std::distance(r_edges.begin(), r_it)-1;
  if (n_E >= bin_counts.size()) n_E = bin_counts.size()-1;
  if (n_r >= bin_counts[0].size()) n_r = bin_counts[0].size()-1;
  return std::make_pair(n_E, n_r);
}

bool Binner::in_range(DataPoint p) const {
  return (p.E > E_min) && (p.E < E_max)
    && (p.r > r_min) && (p.r < r_max);
    // && (std::abs(p.dxds) < dxds_bound)
    // && (std::abs(p.dyds) < dyds_bound);
}

void Binner::add_point(DataPoint p) {
  if (p.r > r_edges.front() && p.r < r_edges.back()) {
    auto [junk, n_r] = get_bin_num(E_edges.front(), p.r);
    lowest_pc_vals[n_r] = std::min(p.E, lowest_pc_vals[n_r]);
  }
  if (in_range(p)) {
    auto [n_E, n_r] = get_bin_num(p.E, p.r);
    bin_counts[n_E][n_r].add_point(p.dxds, p.dyds);
    total_count++;
  }
  return;
}

double Binner::get_E_val(int n_E) const {
  //return E_min + n_E*E_bin_width + 0.5*E_bin_width;
  return (E_edges[n_E] + E_edges[n_E+1])/2;
}


double Binner::get_r_val(int n_r) const {
  //return r_min + n_r*r_bin_width + 0.5*r_bin_width;
  return (r_edges[n_r] + r_edges[n_r+1])/2;
}

const Bin& Binner::get_bin(int n_E, int n_r) const {
  return bin_counts[n_E][n_r];
}

bool Binner::has_empty_bins() const {
  int empty_bins = 0;
  for (const auto& row : bin_counts) {
    for (const auto& bin : row) {
      if (bin.count == 0) empty_bins++;
    }
  }
  std::cout << "Number of empty bins: " << empty_bins
    << '/' << (E_edges.size()-1)*(r_edges.size()-1) << '\n';
  return (empty_bins != 0);
}

bool Binner::has_enough_data() const {

  const char* CLEAR_LINE = "\033[K";
  const char* CURSOR_UP = "\033[4A";
  size_t num_bins = bin_counts.size() * bin_counts[0].size();
  std::cout << CLEAR_LINE << "Simulation progress:\n";
  std::cout << CLEAR_LINE << "\tAverage bin count: " << (double) total_count/num_bins << '\n';

  size_t empty_bins = 0;
  for (const auto& row : bin_counts) {
    for (const auto& bin : row) {
      if (bin.count == 0) empty_bins++;
    }
  }
  std::cout << CLEAR_LINE << "\tNumber of empty bins: " << empty_bins
    << "/" << num_bins << '\n';
  std::cout << CLEAR_LINE << "\t\"Average\" bin count variance: "
    << std::sqrt((double) total_count / num_bins) << '\n';

  bool has_enough = total_count > num_bins*1e5;
  // If !has_enough, we will be rerunning this function next, so we need to move
  // the cursor back up to write over our output
  if (!has_enough)
    std::cout << CURSOR_UP;
  return has_enough;
}

