#pragma once
#include <vector>

struct BinPoint {
  double x, y;
  size_t count;
};

struct DataPoint { double E, r, dxds, dyds; };

struct xy_point { double x, y; };

struct Bin {
  public:
    size_t count;
    double area;
  private:
    mutable std::vector<xy_point> p_list;

  public:
    Bin();
    Bin(const Bin&) = default;
    Bin& operator=(const Bin&) = default;
    ~Bin() = default;
    void add_point(double dxds, double dyds);
    double density(size_t num_elec_in) const;
    void bin_momenta(const char * output_dir,
        double E_elec, double T,
        double E, double r) const;
};

