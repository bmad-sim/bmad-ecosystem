#pragma once
#include <vector>
#include <utility>
//#include "GeantMain.hpp"
#include "G4RunManager.hh"

struct DataPoint {
  double E;
  double r;
  double dxds;
  double dyds;
};

struct xy_point { double x; double y; };

struct Bin {
  public:
    size_t count;
    const size_t* total_count;
    double area;
    double dxds_tot;
    double dyds_tot;
    double dxds_square_tot;
    double dyds_square_tot;
  private:
    mutable std::vector<xy_point> p_list;

  public:
    Bin();
    Bin(const Bin&) = default;
    Bin& operator=(const Bin&) = default;
    ~Bin();
    void add_point(double dxds, double dyds);
    double normal_chi2(double E, double r) const;
    void bin_momenta(const char * output_dir,
        double E_elec, double T,
        double E, double r) const;
    double mu_dxds() const;
    double mu_dyds() const;
    double sd_dxds() const;
    double sd_dyds() const;
    double density(size_t num_elec_in) const;
    void clear();
};


class BinnerBase {
  public:
    virtual ~BinnerBase() = 0;
    virtual void add_point(DataPoint p) = 0;
};


class Binner : public BinnerBase {
  protected:
    double E_min, r_min;
    double E_max, r_max;
    double dxds_bound, dyds_bound;
    //int num_E_bins, num_r_bins;
    //double E_bin_width, r_bin_width;
    std::vector<double> E_edges, r_edges;
    std::vector<std::vector<Bin>> bin_counts;
    size_t total_count;

    std::pair<int, int> get_bin_num(double E, double r);
    //void resize_bins();
    bool in_range(DataPoint p);
    //friend std::pair<int, int> calibrate_binner(G4RunManager *runManager, const std::string& target_material, double in_energy, double target_thickness, Binner *binner);

  public:
    Binner();
    Binner(int N_E_bins, int N_r_bins);
    Binner(int N_E_bins, int N_r_bins, double p_E_max, double p_r_max);
    Binner(int N_E_bins, int N_r_bins, double p_E_min, double p_r_min, double p_E_max, double p_r_max);
    template <typename CBinT> Binner(CBinT&&);
    ~Binner();
    void add_point(DataPoint p);
    double get_E_val(int n_E) const;
    double get_r_val(int n_r) const;
    const Bin& get_bin(int n_E, int n_r) const;
    bool has_empty_bins() const;
    bool has_enough_data() const;
    int get_total() const { return total_count; }
    //double bin_area() const;
};

//class CalibrationBinner : public Binner {
//  protected:
//    std::vector<DataPoint> calibration_set;
//    std::vector<std::vector<Bin>> alt_bin_counts;
//    std::vector<std::vector<Bin>> *coarse_bins, *fine_bins;
//    double growth_factor = 1.25;
//    int calibration_length = 10000;
//    int fine_num_E_bins;
//    int fine_num_r_bins;
//    double fine_E_bin_width;
//    double fine_r_bin_width;
//    std::pair<int, int> get_fine_bin_num(DataPoint p);
//
//  public:
//    CalibrationBinner();
//    void calibrate(G4RunManager *runManager, const std::string& target_material, double in_energy, double target_thickness);
//    void add_point(DataPoint p);
//    void resize_bins();
//    double bin_diff();
//};
//
