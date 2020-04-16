#pragma once
#include <vector>
#include <array>
#include <utility>
#include <tuple>
#include <string>
#include "G4RunManager.hh"
#include "binner.hpp"
#include "parser.hpp"

class CalibrationBinner : public BinnerBase {
  private:
    //double E_min=0, r_min=0;
    //double E_max, r_max;
    //size_t num_E_bins, num_r_bins;
    double E_min, E_max, r_max;
    double dxds_bound, dyds_bound;
    std::vector<DataPoint> cal_run;
    //  E_sorted_cal_run, r_sorted_cal_run;
    //std::vector<std::vector<CalBin>> bins;
    std::vector<double> E_edges, r_edges;
    std::vector<size_t> row_totals, col_totals;
    std::vector<std::array<size_t, 4>> fine_row_totals, fine_col_totals;
    size_t total_count;
    G4RunManager* runManager;
    std::string target_material;
    const double in_energy, target_thickness;
    int calibration_length = 10000;

  public:
    CalibrationBinner() = delete;
    CalibrationBinner(G4RunManager*, const std::string&, double, double);
    CalibrationBinner(G4RunManager*, const SimSettings&, double, double);
    CalibrationBinner(const CalibrationBinner&) = default;
    CalibrationBinner& operator=(const CalibrationBinner&) = default;
    void add_point(DataPoint p) { cal_run.push_back(p); }
    bool in_range(DataPoint p) const {
      return (p.E>E_min) && (p.E<E_max) && (p.r<r_max)
        && (std::abs(p.dxds) < dxds_bound)
        && (std::abs(p.dyds) < dyds_bound);
    }

    std::pair<int, int> calibrate();
    std::pair<int, int> new_calibrate();

    void add_more_points();
    void rebin();
    void add_bins();
    double find_bad_bins(double threshold);
    double get_bin_area(int E_ix, int r_ix) const;
    double row_density(int E_ix) const;
    double col_density(int r_ix) const;
    size_t fine_row_total_diff(int E_ix) const;
    size_t fine_col_total_diff(int r_ix) const;
    double row_density_margin_of_error(int E_ix) const;
    double col_density_margin_of_error(int r_ix) const;
    bool significant_row_diff(int E_ix) const;
    bool significant_col_diff(int r_ix) const;
    friend class Binner;

};
