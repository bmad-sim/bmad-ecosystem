#include <vector>
#include <deque>
#include <utility>
#include <tuple>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <chrono>

#include "cal_binner.hpp"
#include "bin.hpp"
#include "GeantMain.hpp"

using namespace std::string_literals;

CalibrationBinner::CalibrationBinner(G4RunManager* p_run_manager, PointCache* p_cache, const SimSettings& p_s, double p_in_energy, double p_target_thickness) :
  cal_run(),
  E_edges(), r_edges(),
  bins(),
  spin_bins(),
  total_count(0), elec_in(0),
  eff_approx(0.0),
  runManager(p_run_manager), point_cache(p_cache),
  target_material(p_s.target_material),
  out_dir(p_s.output_directory),
  in_energy(p_in_energy), target_thickness(p_target_thickness) {
    if (p_s.polarization_in.size() != 0) {
      in_polarization = G4ThreeVector(p_s.polarization_in[0],
                                      p_s.polarization_in[1],
                                      p_s.polarization_in[2]);
    } else {
      in_polarization = G4ThreeVector(0, 0, 0);
    }
    E_edges.resize(p_s.num_pc_bins + 1);
    E_edges.front() = p_s.out_pc_min;
    E_edges.back() = p_s.out_pc_max;
    bins.resize(E_edges.size()-1);
    spin_bins.resize(bins.size());

    r_edges.resize(p_s.num_r_bins + 1);
    r_edges.front() = 0;
    r_edges.back() = 1e6; // set to a ridiculous number initially, fix later
}


// Main methods
void CalibrationBinner::calibrate() {
  // Gather some data for calibration
  unsigned n_runs = 0;
  size_t calibration_length = 10000;
  do {
    // Short run
    run_simulation(runManager, target_material, in_energy, target_thickness, in_polarization, point_cache, calibration_length);

    // Read data points from point_cache
    point_cache->Lock();
    point_cache->DumpData(&cal_run);
    point_cache->Unlock();
    n_runs++;
  } while (cal_run.size() < 3000);

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
        [this](GeantParticle p) { return !in_range(p); }),
      cal_run.end());

  // Compute approximate efficiency
  eff_approx = cal_run.size() / (n_runs * calibration_length);

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
  stride = (cal_run.size()-1) / (r_edges.size() - 1);
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

  for (auto& row : spin_bins)
    row.resize(r_edges.size()-1);


  // No point in wasting data -> add cal run to bins
  for (auto& p : cal_run) add_point(p);
}


void CalibrationBinner::add_point(const GeantParticle& p) {
  if (!in_range(p)) return;
  auto [n_E, n_r] = get_bin_num(p.E, p.r);
  bins[n_E][n_r].add_point(p);
  spin_bins[n_E][n_r].add_point(p);
  total_count++;
}


void CalibrationBinner::run(size_t num_positrons) {
  point_cache->Lock();
  point_cache->Clear();
  point_cache->Unlock();
  size_t run_length;
  if (eff_approx > 0) run_length = num_positrons / eff_approx;
  else run_length = num_positrons;
  run_simulation(runManager, target_material, in_energy, target_thickness,
      in_polarization, point_cache, run_length);
  for (auto& p : *point_cache) add_point(p);
  elec_in += run_length;
}


void CalibrationBinner::write_data() {
  // Performs the following:
  // * Writes the E,r histogram to E..._T..._er.dat
  // * Writes the average x polarizations to E..._T..._polx.dat
  // * Writes the average y polarizations to E..._T..._poly.dat
  // * Writes the average z polarizations to E..._T..._polz.dat
  // * For each E,r bin, writes the dxds,dyds histogram
  //   to dir_dat/E..._T.../E..._r..._bin.dat, and records
  //   these filenames to E..._T..._xy_bins.txt

  // Open list for dxds, dyds filenames
  char list_filename[200];
  sprintf(list_filename, "%s/E%0.0lf_T%0.3lf_xy_bins.txt",
      out_dir.c_str(), in_energy, target_thickness);
  std::ofstream binlist;
  binlist.open(list_filename);
  if (!binlist) {
    std::cerr << "ERROR: Could not open " << list_filename << " for writing\n";
    return;
  }
  char spinlist_filename[200];
  sprintf(spinlist_filename, "%s/E%0.0lf_T%0.3lf_spin_bins.txt",
      out_dir.c_str(), in_energy, target_thickness);
  std::ofstream sbinlist;
  sbinlist.open(spinlist_filename);
  if (!sbinlist) {
    std::cerr << "ERROR: Could not open " << spinlist_filename << " for writing\n";
    return;
  }

  // Write the dxds, dyds files
  auto num_E_bins = E_edges.size()-1;
  auto num_r_bins = r_edges.size()-1;
  for (size_t i=0; i<num_E_bins; i++) {
    for (size_t j=0; j<num_r_bins; j++) {
      double bin_E = get_E_val(i);
      double bin_r = get_r_val(j);
      auto& bin = get_bin(i,j);
      auto& sbin = get_sbin(i,j);
      bin.bin_momenta();
      bin.write_momenta(out_dir.c_str(), in_energy, target_thickness, bin_E, bin_r, binlist);
      sbin.bin_spins();
      sbin.write_spins(out_dir.c_str(), in_energy, target_thickness, bin_E, bin_r, sbinlist);
    }
  }
  binlist.close();

  // Write the E,r histogram, and polarization files
  write_table("er"s);
  write_table("polx"s);
  write_table("poly"s);
  write_table("polz"s);
  write_table("polxcos"s);
  write_table("polycos"s);
  write_table("polzcos"s);
  write_table("polxsin"s);
  write_table("polysin"s);
  write_table("polzsin"s);
}

void CalibrationBinner::write_table(const std::string& param) const {
  // Averaging functions
  constexpr auto polx = [](const GeantParticle& p) { return p.polx; };
  constexpr auto poly = [](const GeantParticle& p) { return p.poly; };
  constexpr auto polz = [](const GeantParticle& p) { return p.polz; };

  constexpr auto polxcos = [](const GeantParticle& p) { return p.polx * cos(p.theta); };
  constexpr auto polycos = [](const GeantParticle& p) { return p.poly * cos(p.theta); };
  constexpr auto polzcos = [](const GeantParticle& p) { return p.polz * cos(p.theta); };

  constexpr auto polxsin = [](const GeantParticle& p) { return p.polx * sin(p.theta); };
  constexpr auto polysin = [](const GeantParticle& p) { return p.poly * sin(p.theta); };
  constexpr auto polzsin = [](const GeantParticle& p) { return p.polz * sin(p.theta); };

  // Initialize file io
  std::ofstream table_file;
  char filename[200];
  sprintf(filename, "%s/E%0.0lf_T%0.3lf_%s.dat",
      out_dir.c_str(), in_energy, target_thickness, param.c_str());
  table_file.open(filename);
  if (!table_file) {
    std::cerr << "ERROR: could not open " << filename << '\n';
    return;
  }

  // set up first line with r values
  auto num_E_bins = E_edges.size()-1;
  auto num_r_bins = r_edges.size()-1;
  table_file << num_r_bins;
  for (size_t j=0; j<num_r_bins; j++)
    table_file << '\t' << get_r_val(j);
  table_file << '\n';
  // each line should be one E value
  for (size_t i=0; i<num_E_bins; i++) {
    table_file << get_E_val(i);
    for (size_t j=0; j<num_r_bins; j++) {
      double bin_E = get_E_val(i);
      double bin_r = get_r_val(j);
      if (param != "er"s && j == 0) bin_r = 0; // r=0 bin for spin data
      auto& bin = get_bin(i,j);
      if (param == "er"s) table_file << '\t' << bin.density(elec_in);
      else if (param == "polx"s) table_file << '\t' << bin.get_average(polx);
      else if (param == "poly"s) table_file << '\t' << bin.get_average(poly);
      else if (param == "polz"s) table_file << '\t' << bin.get_average(polz);
      else if (param == "polxcos"s) table_file << '\t' << bin.get_average(polxcos);
      else if (param == "polycos"s) table_file << '\t' << bin.get_average(polycos);
      else if (param == "polzcos"s) table_file << '\t' << bin.get_average(polzcos);
      else if (param == "polxsin"s) table_file << '\t' << bin.get_average(polxsin);
      else if (param == "polysin"s) table_file << '\t' << bin.get_average(polysin);
      else if (param == "polzsin"s) table_file << '\t' << bin.get_average(polzsin);
    }
    table_file << '\n';
  }
  table_file.close();
}



// Support methods
double CalibrationBinner::get_E_val(int n_E) const {
  return 0.5*(E_edges[n_E] + E_edges[n_E+1]);
}
double CalibrationBinner::get_r_val(int n_r) const {
  return 0.5*(r_edges[n_r] + r_edges[n_r+1]);
}
bool CalibrationBinner::in_range(const GeantParticle& p) const {
  return (p.E>E_edges.front()) && (p.E<E_edges.back()) && (p.r<r_edges.back());
}
ERBin& CalibrationBinner::get_bin(int n_E, int n_r) {
  return bins[n_E][n_r];
}
SpinBinner& CalibrationBinner::get_sbin(int n_E, int n_r) {
  return spin_bins[n_E][n_r];
}
const ERBin& CalibrationBinner::get_bin(int n_E, int n_r) const {
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

// Helper class for CalibrationBinner::has_enough_data()
struct BatchStats {
  size_t num_elec, num_pos;
  std::chrono::duration<double> time;
};
BatchStats operator+(const BatchStats& a, const BatchStats& b) {
  return { a.num_elec + b.num_elec, a.num_pos + b.num_pos, a.time + b.time };
}

bool CalibrationBinner::has_enough_data() const {
  constexpr double per_bin_target = 1e5;
  // Record of most recent (past minute) incoming/outgoing throughput
  static std::deque<BatchStats> recent_batches{};
  constexpr double time_to_keep = 60; // seconds

  // Timing
  static auto prev_t = std::chrono::steady_clock::now();
  auto current = std::chrono::steady_clock::now();
  std::chrono::duration<double> delta_t = current - prev_t;

  static size_t prev_pos_count = 0;
  static size_t prev_elec_count = 0;

  auto delta_elec = (elec_in - prev_elec_count);
  auto delta_pos = (total_count - prev_pos_count);
  recent_batches.push_back({delta_elec, delta_pos, delta_t});

  auto recent_total = std::accumulate(recent_batches.begin(), recent_batches.end(),
      BatchStats{}); // default constructed BatchStats{} has all members == 0

  double elec_throughput = recent_total.num_elec / recent_total.time.count();
  double pos_throughput = recent_total.num_pos / recent_total.time.count();

  size_t num_bins = bins.size() * bins[0].size();
  double count_per_bin = static_cast<double>(total_count) / num_bins;

  auto bin_is_empty = [](const ERBin& bin) { return bin.count == 0; };
  size_t empty_bins = std::accumulate(bins.begin(), bins.end(), 0ull,
      [&bin_is_empty](const auto& total, const auto& row) {
        return total + std::count_if(row.begin(), row.end(), bin_is_empty);
      });

  double efficiency = static_cast<double>(total_count) / elec_in;

  const char* CLEAR_LINE = "\033[K";
  const char* CURSOR_UP = "\033[9F";
  printf(CLEAR_LINE);
  printf("Simulation progress:\n");

  printf(CLEAR_LINE);
  printf("\tAverage bin count: %0.03e\n", count_per_bin);

  printf(CLEAR_LINE);
  printf("\tNumber of empty bins: %lu/%lu\n", empty_bins, num_bins);

  printf(CLEAR_LINE);
  printf("\tAverage bin count s.d.: %0.03e\n", std::sqrt(count_per_bin));

  printf(CLEAR_LINE);
  printf("\tNumber of incoming particles tracked: %0.03e / %0.02e\n",
      static_cast<double>(elec_in),
      num_bins * per_bin_target / efficiency);

  printf(CLEAR_LINE);
  printf("\tNumber of outgoing particles tracked: %0.03e / %0.02e\n",
      static_cast<double>(total_count),
      static_cast<double>(num_bins * per_bin_target));

  if (delta_t.count() > 1e-6) {
    printf(CLEAR_LINE);
    printf("\tIncoming particle throughput: %0.02e/sec\n", elec_throughput);

    printf(CLEAR_LINE);
    printf("\tOutgoing particle throughput: %0.02e/sec\n", pos_throughput);

    auto time_left = std::max(0.0, (per_bin_target * num_bins - total_count) / pos_throughput);
    unsigned hr = time_left / 3600;
    unsigned mn = (time_left - 3600*hr) / 60;
    unsigned sc = (time_left - 3600*hr - 60*mn);
    printf(CLEAR_LINE);
    printf("\tEstimated time to completion: %u:%02u:%02u\n", hr, mn, sc);
  } else {
    printf("\n\n\n");
  }

  // Update static variables
  prev_elec_count = elec_in;
  prev_pos_count = total_count;
  prev_t = current;
  if (recent_total.time.count() > time_to_keep) // only keep 60 seconds of data
    recent_batches.pop_front();


  bool has_enough = count_per_bin > per_bin_target;
  // If !has_enough, we will be rerunning this function next, so we need to move
  // the cursor back up to write over our output
  if (!has_enough) {
    printf(CURSOR_UP);
  } else {
    // Clear static variables; the next time this function is called,
    // it will be for a different pc_in/T combination
    prev_elec_count = 0;
    prev_pos_count = 0;
    recent_batches.clear();
  }
  return has_enough;
}
