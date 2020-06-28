#pragma once
#include <vector>
#include <cstdio>

struct BinPoint {
  double x, y;
  size_t count;
  double area;
};

struct DataPoint { double E, r, dxds, dyds; };

struct xy_point { double x, y; };

struct Bin {
  public:
    size_t count;
    double area;
  private:
    mutable std::vector<xy_point> p_list;
    mutable std::vector<BinPoint> *binned_ptr;

  public:
    Bin();
    Bin(const Bin&);
    Bin& operator=(const Bin&);
    ~Bin();
    void add_point(double dxds, double dyds);
    double density(size_t num_elec_in) const;
    void write_momenta(const char * output_dir,
        double E_elec, double T,
        double E, double r, FILE *binlist) const;
    void binary_write_momenta(FILE *binary_file) const;
    void bin_momenta() const;
};

