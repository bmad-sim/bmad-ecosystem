#include <vector>
#include <utility>
#include <tuple>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <numeric>
#include "cal_binner.hpp"
#include "gnuplot.hpp"
#include "GeantMain.hpp"
#include "binner.hpp"

#define PRINT_DIAGNOSTICS

CalibrationBinner::CalibrationBinner(G4RunManager* p_run_manager, const std::string& p_target_material, double p_in_energy, double p_target_thickness) :
  E_min(0),
  E_max(0), r_max(0),
  dxds_bound(10), dyds_bound(10),
  cal_run({}), //E_sorted_cal_run({}), r_sorted_cal_run({}),
  /*bins({}), */ E_edges({}), r_edges({}),
  row_totals({}), col_totals({}),
  fine_row_totals({}), fine_col_totals({}),
  total_count(0),
  runManager(p_run_manager), target_material(p_target_material),
  in_energy(p_in_energy), target_thickness(p_target_thickness) {}


CalibrationBinner::CalibrationBinner(G4RunManager* p_run_manager, const SimSettings& p_s, double p_in_energy, double p_target_thickness) :
  E_min(p_s.out_energy_min),
  E_max(p_s.out_energy_max),
  r_max(p_s.out_r_max),
  dxds_bound(p_s.dxy_ds_max), dyds_bound(p_s.dxy_ds_max),
  cal_run({}), //E_sorted_cal_run({}), r_sorted_cal_run({}),
  /*bins({}), */ E_edges({}), r_edges({}),
  row_totals({}), col_totals({}),
  fine_row_totals({}), fine_col_totals({}),
  total_count(0),
  runManager(p_run_manager), target_material(p_s.target_material),
  in_energy(p_in_energy), target_thickness(p_target_thickness) {}


std::pair<int,int> CalibrationBinner::calibrate() {
  // Short run
  cal_run.reserve(calibration_length);
  run_simulation(runManager, target_material, in_energy, target_thickness, this, calibration_length);

  // Determine max E and r values
  constexpr int max_ptile = 95;
  static_assert(max_ptile < 100 && max_ptile > 0);
  size_t ix_max = cal_run.size() * max_ptile / 100;

  std::sort(cal_run.begin(), cal_run.end(),
      [](DataPoint p1, DataPoint p2) { return p1.E < p2.E; });
  E_max = cal_run[ix_max].E;

  std::sort(cal_run.begin(), cal_run.end(),
      [](DataPoint p1, DataPoint p2) { return p1.r < p2.r; });
  r_max = cal_run[ix_max].r;

  // Remove out-of range datapoints from cal_run
  cal_run.erase(std::remove_if(cal_run.begin(), cal_run.end(),
        [this](DataPoint p) { return !in_range(p); }),
      cal_run.end());

  double threshold = 0.10;
  // Set initial E,r edges
  constexpr size_t initial_edge_count=2;
  static_assert(initial_edge_count >= 2);

  E_edges.resize(initial_edge_count);
  std::iota(E_edges.begin(), E_edges.end(), 0.0);
  for (auto& E : E_edges) E *= E_max/(initial_edge_count-1);

  r_edges.resize(initial_edge_count);
  std::iota(r_edges.begin(), r_edges.end(), 0.0);
  for (auto& r : r_edges) r *= r_max/(initial_edge_count-1);

  do {
    add_bins(); // handles initial bins, edges setup
    rebin();
  } while (find_bad_bins(threshold) > threshold);

  // return the number of bins used
  return std::make_pair(E_edges.size()-1, r_edges.size()-1);
}


void CalibrationBinner::rebin() {
  // Bins the data points in cal_run
  // this->bins should already be resized and zeroed
  // to be in-sync with E_edges and r_edges,
  // and each bin should have its current area assigned

  // Zero the counts
  for (auto& x : row_totals) x=0;
  for (auto& x : col_totals) x=0;
  for (auto& x : fine_row_totals) x={0,0,0,0};
  for (auto& x : fine_col_totals) x={0,0,0,0};

  // bin the data
  size_t E_ix, r_ix;
  size_t fine_E_ix, fine_r_ix;
  double E_width, r_width;
  for (auto p : cal_run) {
    auto left_E_edge = std::adjacent_find(
        E_edges.begin(), E_edges.end(),
        [p](double low, double high) {
          return ((low<=p.E)&&(high>=p.E));
        });
    E_ix = std::distance(E_edges.begin(), left_E_edge);
    E_width = E_edges[E_ix+1] - E_edges[E_ix];
    fine_E_ix = (p.E - E_edges[E_ix])*4 / E_width;

    auto left_r_edge = std::adjacent_find(
        r_edges.begin(), r_edges.end(),
        [p](double low, double high) {
          return ((low<p.r)&&(high>p.r));
        });
    r_ix = std::distance(r_edges.begin(), left_r_edge);
    r_width = r_edges[r_ix+1] - r_edges[r_ix];
    fine_r_ix = (p.r - r_edges[r_ix])*4 / r_width;

    //if (!in_range(p)) {
    //  std::cout << "Out-of-range datapoint present in calibration run\n";
    //  getchar();
    //}
    //if (left_E_edge == E_edges.end()) {
    //  std::cout << "Walked off the edge of E_edges\n";
    //  getchar();
    //}
    //if (left_r_edge == r_edges.end()) {
    //  std::cout << "Walked off the edge of r_edges\n";
    //  getchar();
    //}

    row_totals[E_ix]++;
    col_totals[r_ix]++;
    fine_row_totals[E_ix][fine_E_ix]++;
    fine_col_totals[r_ix][fine_r_ix]++;
  }
  return;
}


double CalibrationBinner::find_bad_bins(double threshold) {
  // This function returns the largest difference
  // present between two adjacent bins as a fraction of
  // the maximum probability density in any bin, and inserts
  // new edges into E/r edges to fix the bins whose split
  // difference is greater than the threshold

  int num_E_bins = row_totals.size();
  int num_r_bins = col_totals.size();
  double split_diff;
  double max_diff = 0;
  std::vector<size_t> bad_E_indices, bad_r_indices;
  //std::set<int> bad_E_bins, bad_r_bins;
  //size_t count_threshold = total_count * 1e-3;
  bool at_least_one_candidate = false;

  double max_row_density=0, max_col_density=0;

  // Check adjacent rows (neighboring E values)
  for (int i=0; i<num_E_bins; i++) {
    max_row_density = std::max(max_row_density, row_density(i));
  }

#ifdef PRINT_DIAGNOSTICS
  std::cout << "Max row density: " << max_row_density << '\n';
#endif

  for (int i=0; i<num_E_bins; i++) {
    if (fine_row_total_diff(i) > 5*row_density_margin_of_error(i)) {
      at_least_one_candidate = true;
      split_diff = 4*fine_row_total_diff(i)/
        ((E_edges[i+1] - E_edges[i]) * max_row_density);
      if (split_diff > max_diff) {
        max_diff = split_diff;
#ifdef PRINT_DIAGNOSTICS
        std::cout << "--->New max diff found for row "
          << i << " with counts "
          << fine_row_totals[i][0] << ", "
          << fine_row_totals[i][1] << ", "
          << fine_row_totals[i][2] << ", "
          << fine_row_totals[i][3] << '\n'
          << "    New max diff is " << max_diff << '\n';
#endif
      }
      if (split_diff > threshold) {
        std::cout << "Adding new E edge between "
          << E_edges[i] << " and " << E_edges[i+1] << '\n';
        bad_E_indices.push_back(i);
      }
    }
  }

  // Check adjacent differences in columns (neighboring E values)
  for (int j=0; j<num_r_bins; j++) {
    max_col_density = std::max(max_col_density, col_density(j));
  }

#ifdef PRINT_DIAGNOSTICS
  std::cout << "Max col density: " << max_col_density << '\n';
#endif

  for (int i=0; i<num_r_bins; i++) {
    if (fine_col_total_diff(i) > 5*col_density_margin_of_error(i)) {
      at_least_one_candidate = true;
      split_diff = 4*fine_col_total_diff(i)/
        ((r_edges[i+1] - r_edges[i]) * max_col_density);
      if (split_diff > max_diff) {
        max_diff = split_diff;
#ifdef PRINT_DIAGNOSTICS
        std::cout << "--->New max diff found for col "
          << i << " with counts "
          << fine_col_totals[i][0] << ", "
          << fine_col_totals[i][1] << ", "
          << fine_col_totals[i][2] << ", "
          << fine_col_totals[i][3] << '\n'
          << "    New max diff is " << max_diff << '\n';
#endif
      }
      if (split_diff > threshold) {
        std::cout << "Adding new r edge between "
          << r_edges[i] << " and " << r_edges[i+1] << '\n';
        bad_r_indices.push_back(i);
      }
    }
  }

  // Add new edges to fix bad bins
  double new_edge;
  for (auto E_ix : bad_E_indices) {
    new_edge = (E_edges[E_ix] + E_edges[E_ix+1])/2;
    E_edges.push_back(new_edge);
  }
  for (auto r_ix : bad_r_indices) {
    new_edge = (r_edges[r_ix] + r_edges[r_ix+1])/2;
    r_edges.push_back(new_edge);
  }
  // Since new edges are being added onto the end of E_edges and r_edges,
  // the indices stored in bad_E_bins and bad_r_bins will not
  // be shifted around in the above process

  // Resort the edges
  std::sort(E_edges.begin(), E_edges.end());
  std::sort(r_edges.begin(), r_edges.end());

  if (at_least_one_candidate) {
#ifdef PRINT_DIAGNOSTICS
    std::cout << "Max diff: " << max_diff << '\n';
    getchar();
#endif
    return max_diff;
  } else {
#ifdef PRINT_DIAGNOSTICS
    std::cout << "Adding more points to calibration run...\n";
#endif
    add_more_points();
    rebin();
    return find_bad_bins(threshold);
  }
}

void CalibrationBinner::add_more_points() {
  // Runs the simulation some more to add more
  // points to the calibration run, then removes any that are
  // out of range
  run_simulation(runManager, target_material, in_energy, target_thickness, this, calibration_length);

  // remove out-of-range points
  cal_run.erase(std::remove_if(cal_run.begin(), cal_run.end(),
        [this](DataPoint p) { return !in_range(p); }),
      cal_run.end());

  // update total_count
  total_count = cal_run.size();
  std::cout << "total count is now " << total_count << '\n';
  // rebin
  //add_bins();
  rebin();
  return;
}

void CalibrationBinner::add_bins() {
  // Increases the length of E_edges and r_edges
  // and sets their contents so that each row
  // and column contains the same number of points
  // Also resizes bins and row/col totals
  //size_t inc = 30;
  //E_edges.resize(E_edges.size()+inc);
  //r_edges.resize(r_edges.size()+inc);

  //E_edges[0]=0;
  //E_edges[E_edges.size()-1]=E_max;
  //r_edges[0]=0;
  //r_edges[r_edges.size()-1]=r_max;

  //size_t cal_run_stride = E_sorted_cal_run.size()/E_edges.size();

  //for (size_t E_ix=1; E_ix<E_edges.size()-1; E_ix++)
  //  E_edges[E_ix] = E_sorted_cal_run[E_ix * cal_run_stride].E;

  //for (size_t r_ix=1; r_ix<r_edges.size()-1; r_ix++)
  //  r_edges[r_ix] = r_sorted_cal_run[r_ix * cal_run_stride].r;

  //bins.resize(E_edges.size()-1);
  //int E_ix, r_ix;
  //E_ix = 0;
  //for (auto& row : bins) {
  //  row.resize(r_edges.size()-1);
  //  r_ix = 0;
  //  for (auto& bin : row) {
  //    bin.count = 0;
  //    bin.total_count = &total_count;
  //    bin.area = get_bin_area(E_ix, r_ix);
  //    r_ix++;
  //  }
  //  E_ix++;
  //}
  // Resize and zero row and column totals
  row_totals.resize(E_edges.size()-1);
  col_totals.resize(r_edges.size()-1);
  fine_row_totals.resize(E_edges.size()-1);
  fine_col_totals.resize(r_edges.size()-1);
  return;
}


double CalibrationBinner::get_bin_area(int E_ix, int r_ix) const {
  return (E_edges[E_ix+1] - E_edges[E_ix]) * (r_edges[r_ix+1] - r_edges[r_ix]);
}

double CalibrationBinner::row_density(int E_ix) const {
  return row_totals[E_ix]/(E_edges[E_ix+1]-E_edges[E_ix]);
}


double CalibrationBinner::col_density(int r_ix) const {
  return col_totals[r_ix]/(r_edges[r_ix+1]-r_edges[r_ix]);
}

size_t CalibrationBinner::fine_row_total_diff(int E_ix) const {
  // Returns the difference in counts for the two fine row bins
  // for bin number E_ix which are the furthest apart
  return *std::max_element(fine_row_totals[E_ix].begin(), fine_row_totals[E_ix].end())
    - *std::min_element(fine_row_totals[E_ix].begin(), fine_row_totals[E_ix].end());
}

size_t CalibrationBinner::fine_col_total_diff(int r_ix) const {
  // Returns the difference in counts for the two fine column bins
  // for bin number r_ix which are the furthest apart
  return *std::max_element(fine_col_totals[r_ix].begin(), fine_col_totals[r_ix].end())
    - *std::min_element(fine_col_totals[r_ix].begin(), fine_col_totals[r_ix].end());
}

double CalibrationBinner::row_density_margin_of_error(int E_ix) const {
  return std::sqrt(*std::max_element(fine_row_totals[E_ix].begin(), fine_row_totals[E_ix].end()))
    + std::sqrt(*std::min_element(fine_row_totals[E_ix].begin(), fine_row_totals[E_ix].end()));
}

double CalibrationBinner::col_density_margin_of_error(int r_ix) const {
  return std::sqrt(*std::max_element(fine_col_totals[r_ix].begin(), fine_col_totals[r_ix].end()))
    + std::sqrt(*std::min_element(fine_col_totals[r_ix].begin(), fine_col_totals[r_ix].end()));
}

bool CalibrationBinner::significant_row_diff(int E_ix) const {
  return std::abs<int>(row_totals[E_ix] - row_totals[E_ix+1])
    > std::sqrt(row_totals[E_ix]) + std::sqrt(row_totals[E_ix+1]);
}

bool CalibrationBinner::significant_col_diff(int r_ix) const {
  return std::abs<int>(col_totals[r_ix] - col_totals[r_ix+1])
    > std::sqrt(col_totals[r_ix]) + std::sqrt(col_totals[r_ix+1]);
}


// CalibrationBinner to Binner constructor
template <typename CBinT>
Binner::Binner(CBinT&& cb) :
  E_min(0), r_min(0),
  E_max(cb.E_max), r_max(cb.r_max),
  dxds_bound(10), dyds_bound(10),
  E_edges(std::move(cb.E_edges)),
  r_edges(std::move(cb.r_edges)),
  bin_counts(E_edges.size()-1,
      std::vector<Bin>(r_edges.size()-1)),
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

std::pair<int,int> CalibrationBinner::new_calibrate() {
  // Short run
  cal_run.reserve(calibration_length);
  run_simulation(runManager, target_material, in_energy, target_thickness, this, calibration_length);

  // E_min, E_max, r_max are already set by constructor
  // -> remove out of range points
  cal_run.erase(std::remove_if(cal_run.begin(), cal_run.end(),
        [this](DataPoint p) { return !in_range(p); }),
      cal_run.end());

  // Set E,r edges
  constexpr size_t edge_count=16;
  static_assert(edge_count >= 2);
  E_edges.resize(edge_count);
  r_edges.resize(edge_count);
  // E edges
  std::sort(cal_run.begin(), cal_run.end(),
      [](const auto& p1, const auto& p2) {
        return p1.E < p2.E;
      });
  for (size_t i=0; i<edge_count; i++) {
    size_t c_ix = i * cal_run.size() / (edge_count-1);
    if (i == edge_count-1) c_ix--;
    E_edges[i] = cal_run[c_ix].E;
  }
  // r edges
  std::sort(cal_run.begin(), cal_run.end(),
      [](const auto& p1, const auto& p2) {
        return p1.r < p2.r;
      });
  for (size_t i=0; i<edge_count; i++) {
    size_t c_ix = i * cal_run.size() / (edge_count-1);
    if (i == edge_count-1) c_ix--;
    r_edges[i] = cal_run[c_ix].r;
  }

  // return the number of bins used
  return std::make_pair(E_edges.size()-1, r_edges.size()-1);
}

//void pretty_vec_print(const std::vector<std::vector<CalBin>>& mat) {
//  for (const auto& row: mat) {
//    for (const auto& bin : row) std::cout << bin.count << '\t';
//    std::cout << '\n';
//  }
//  return;
//}
