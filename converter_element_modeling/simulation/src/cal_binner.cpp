#include <vector>
#include <utility>
#include <tuple>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "cal_binner.hpp"
#include "GeantMain.hpp"
#include "binner.hpp"

CalibrationBinner::CalibrationBinner(G4RunManager* p_run_manager, const SimSettings& p_s, double p_in_energy, double p_target_thickness) :
  E_min(p_s.out_pc_min),
  E_max(p_s.out_pc_max),
  r_max(0),
  dxds_bound(p_s.dxy_ds_max), dyds_bound(p_s.dxy_ds_max),
  cal_run({}),
  E_edges({}), r_edges({}),
  total_count(0),
  runManager(p_run_manager), target_material(p_s.target_material),
  in_energy(p_in_energy), target_thickness(p_target_thickness) {
    if (p_s.num_pc_bins==0)
      E_edges.resize(p_s.num_bins + 1);
    else
      E_edges.resize(p_s.num_pc_bins + 1);

    if (p_s.num_r_bins==0)
      r_edges.resize(p_s.num_bins + 1);
    else
      r_edges.resize(p_s.num_r_bins + 1);
}


// CalibrationBinner to Binner constructor
template <typename CBinT>
Binner::Binner(CBinT&& cb) :
  E_min(cb.E_min), r_min(0),
  E_max(cb.E_max), r_max(cb.r_max),
  dxds_bound(cb.dxds_bound), dyds_bound(cb.dyds_bound),
  E_edges(std::move(cb.E_edges)),
  r_edges(std::move(cb.r_edges)),
  bin_counts(E_edges.size()-1,
      std::vector<Bin>(r_edges.size()-1)),
  lowest_pc_vals(r_edges.size()-1, 100),
  total_count(0) {
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
    std::cout << "Binner configured from CalibrationBinnner\n";
    std::cout << "E_edges: {";
    for (auto E : E_edges) { std::cout << E << ", "; }
    std::cout << "}\n";
    std::cout << "r_edges: {";
    for (auto r : r_edges) { std::cout << r << ", "; }
    std::cout << "}\n";
}

template Binner::Binner(CalibrationBinner&&);

std::pair<size_t, size_t> CalibrationBinner::calibrate() {
  // Short run
  cal_run.reserve(calibration_length);
  run_simulation(runManager, target_material, in_energy, target_thickness, this, calibration_length);

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
  r_max = cal_run[ix99].r;
  std::sort(cal_run.begin(), cal_run.end(), esort);
  if (E_max == 0) E_max = cal_run[ix99].E;

  // -> remove out of range points
  cal_run.erase(std::remove_if(cal_run.begin(), cal_run.end(),
        [this](DataPoint p) { return !in_range(p); }),
      cal_run.end());

  // Set E,r edges
  // E edges
  std::sort(cal_run.begin(), cal_run.end(), esort);
  size_t e_stride = (cal_run.size()-1) / (E_edges.size() - 1);
  size_t c_ix = -e_stride;
  auto gen_e = [&,c_ix]() mutable {
    c_ix += e_stride;
    return cal_run[c_ix].E;
  };
  std::generate(E_edges.begin(), E_edges.end(), gen_e);
  // r edges
  std::sort(cal_run.begin(), cal_run.end(), rsort);
  size_t r_stride = (cal_run.size()-1) / (E_edges.size() - 1);
  c_ix = -r_stride;
  auto gen_r = [&,c_ix]() mutable {
    c_ix += r_stride;
    return cal_run[c_ix].r;
  };
  std::generate(r_edges.begin(), r_edges.end(), gen_r);

  // return the number of bins used
  return {E_edges.size()-1, r_edges.size()-1};
}

