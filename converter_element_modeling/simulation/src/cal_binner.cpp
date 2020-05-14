#include <vector>
#include <utility>
#include <tuple>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

#include "cal_binner.hpp"
#include "bin.hpp"
#include "GeantMain.hpp"


// Constructor
CalibrationBinner::CalibrationBinner(RunManager_t* p_run_manager, PointCache* p_cache, const SimSettings& p_s, double p_in_energy, double p_target_thickness) :
  cal_run({}),
  E_edges({}), r_edges({}),
  bins(),
  total_count(0), elec_in(0),
  runManager(p_run_manager), point_cache(p_cache),
  target_material(p_s.target_material),
  out_dir(p_s.output_directory),
  in_energy(p_in_energy), target_thickness(p_target_thickness) {
    if (p_s.num_pc_bins==0)
      E_edges.resize(p_s.num_bins + 1);
    else
      E_edges.resize(p_s.num_pc_bins + 1);
    E_edges.front() = p_s.out_pc_min;
    E_edges.back() = p_s.out_pc_max;
    bins.resize(E_edges.size()-1);

    if (p_s.num_r_bins==0)
      r_edges.resize(p_s.num_bins + 1);
    else
      r_edges.resize(p_s.num_r_bins + 1);
    r_edges.front() = 0;
    r_edges.back() = 1e6; // set to a ridiculous number initially, fix later
}


// Main methods
void CalibrationBinner::calibrate() {
  // Short run
  size_t calibration_length = 10000;
  run_simulation(runManager, target_material, in_energy, target_thickness, point_cache, calibration_length);

  // Read data points from point_cache
  point_cache->Lock();
  point_cache->DumpData(&cal_run);
  point_cache->Unlock();


  // Sorting functions
  auto esort = [](const auto& p1, const auto& p2) {
    return p1.E < p2.E;
  };
  auto rsort = [](const auto& p1, const auto& p2) {
    return p1.r < p2.r;
  };

  // Select E_max, r_max to capture 99% of data
  size_t ix99 = 0.99 * cal_run.size();
  std::sort(cal_run.begin(), cal_run.end(), rsort);
  r_edges.back() = cal_run[ix99].r;
  std::sort(cal_run.begin(), cal_run.end(), esort);
  E_edges.back() = std::min(cal_run[ix99].E, E_edges.back());

  // -> remove out of range points
  cal_run.erase(std::remove_if(cal_run.begin(), cal_run.end(),
        [this](DataPoint p) { return !in_range(p); }),
      cal_run.end());

  // Set E,r edges
  size_t stride = (cal_run.size()-1) / (E_edges.size() - 1);
  // E edges
  std::sort(cal_run.begin(), cal_run.end(), esort);
  size_t c_ix = 0;
  auto gen_e = [&,c_ix]() mutable {
    c_ix += stride;
    return cal_run[c_ix].E;
  };
  std::generate(E_edges.begin()+1, E_edges.end()-1, gen_e);
  // r edges
  std::sort(cal_run.begin(), cal_run.end(), rsort);
  c_ix = 0;
  auto gen_r = [&,c_ix]() mutable {
    c_ix += stride;
    return cal_run[c_ix].r;
  };
  std::generate(r_edges.begin()+1, r_edges.end()-1, gen_r);

  // Set bin areas
  unsigned E_ix=0, r_ix;
  for (auto & row : bins) {
    r_ix = 0;
    row.resize(r_edges.size()-1);
    for (auto& b : row) {
      b.area = (E_edges[E_ix+1] - E_edges[E_ix])
        * (r_edges[r_ix+1] - r_edges[r_ix]);
      r_ix++;
    }
    E_ix++;
  }


  // No point in wasting data -> add cal run to bins
  for (auto& p : cal_run) add_point(p);
}


void CalibrationBinner::add_point(DataPoint p) {
  if (!in_range(p)) return;
  auto [n_E, n_r] = get_bin_num(p.E, p.r);
  bins[n_E][n_r].add_point(p.dxds, p.dyds);
  total_count++;
}


void CalibrationBinner::run(size_t run_length) {
  point_cache->Lock();
  point_cache->Clear();
  point_cache->Unlock();
  run_simulation(runManager, target_material, in_energy, target_thickness,
      point_cache, run_length);
  for (auto& p : *point_cache) add_point(p);
}


void CalibrationBinner::write_data() {
  std::ofstream outfile;
  char filename[200];
  sprintf(filename, "%s/E%0.0lf_T%0.3lf_er.dat",
      out_dir.c_str(), in_energy, target_thickness);
  outfile.open(filename);

  // set up first line with r values
  auto num_E_bins = E_edges.size()-1;
  auto num_r_bins = r_edges.size()-1;
  outfile << num_r_bins;
  for (size_t j=0; j<num_r_bins; j++)
    outfile << '\t' << get_r_val(j);
  outfile << '\n';
  // each line should be one E value
  for (size_t i=0; i<num_E_bins; i++) {
    outfile << get_E_val(i);
    for (size_t j=0; j<num_r_bins; j++) {
      double bin_E = get_E_val(i);
      double bin_r = get_r_val(j);
      auto& bin = get_bin(i,j);
      outfile << '\t' << bin.density(elec_in);
      bin.bin_momenta(out_dir.c_str(), in_energy, target_thickness, bin_E, bin_r);
    }
    outfile << '\n';
  }
  outfile.close();
}


// Support methods
double CalibrationBinner::get_E_val(int n_E) const {
  return 0.5*(E_edges[n_E] + E_edges[n_E+1]);
}
double CalibrationBinner::get_r_val(int n_r) const {
  return 0.5*(r_edges[n_r] + r_edges[n_r+1]);
}
bool CalibrationBinner::in_range(DataPoint p) const {
  return (p.E>E_edges.front()) && (p.E<E_edges.back()) && (p.r<r_edges.back());
}
Bin& CalibrationBinner::get_bin(int n_E, int n_r) {
  return bins[n_E][n_r];
}
std::pair<unsigned, unsigned> CalibrationBinner::get_bin_num(double E, double r) const {
  auto E_it = std::find_if(E_edges.begin(), E_edges.end(),
      [E](double edge) { return edge >= E; });
  size_t n_E = std::distance(E_edges.begin(), E_it)-1;
  auto r_it = std::find_if(r_edges.begin(), r_edges.end(),
      [r](double edge) { return edge >= r; });
  size_t n_r = std::distance(r_edges.begin(), r_it)-1;
  if (n_E >= bins.size()) n_E = bins.size()-1;
  if (n_r >= bins[0].size()) n_r = bins[0].size()-1;
  return std::make_pair(n_E, n_r);
}
bool CalibrationBinner::has_enough_data() const {
  const char* CLEAR_LINE = "\033[K";
  const char* CURSOR_UP = "\033[4A";
  size_t num_bins = bins.size() * bins[0].size();
  std::cout << CLEAR_LINE << "Simulation progress:\n";
  std::cout << CLEAR_LINE << "\tAverage bin count: " << (double) total_count/num_bins << '\n';

  size_t empty_bins = 0;
  for (const auto& row : bins) {
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
