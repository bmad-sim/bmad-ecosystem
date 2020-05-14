#pragma once
#include <vector>
#include <array>
#include <utility>
#include <tuple>
#include <string>
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "parser.hpp"
#include "bin.hpp"
#include "point_cache.hpp"

#ifdef G4MULTITHREADED
typedef G4RunManager RunManager_t;
#else
typedef G4RunManager RunManager_t;
#endif


class CalibrationBinner {
  private:
    std::vector<DataPoint> cal_run;
    std::vector<double> E_edges, r_edges;
    std::vector<std::vector<Bin>> bins;
    size_t total_count, elec_in;
    RunManager_t* runManager;
    PointCache* point_cache;
    std::string target_material, out_dir;
    const double in_energy, target_thickness;

    std::pair<unsigned, unsigned> get_bin_num(double E, double r) const;

  public:
    CalibrationBinner() = delete;
    CalibrationBinner(RunManager_t*, PointCache*, const SimSettings&, double, double);
    CalibrationBinner(const CalibrationBinner&) = default;
    CalibrationBinner& operator=(const CalibrationBinner&) = default;

    void calibrate();
    void run(size_t run_length);
    void write_data();

    void add_point(DataPoint p);
    double get_E_val(int n_E) const;
    double get_r_val(int n_r) const;
    Bin& get_bin(int n_E, int n_r);
    bool in_range(DataPoint p) const;
    bool has_empty_bins() const;
    bool has_enough_data() const;
    size_t get_total() const { return total_count; }

};
