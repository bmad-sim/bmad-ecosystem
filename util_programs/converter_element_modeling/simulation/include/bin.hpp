#pragma once
#include <iostream>
#include <vector>
#include <optional>
#include <functional>

struct XYBin {
  double dxds, dyds;
  double avg_polx, avg_poly, avg_polz;
  size_t count;
  double area;
};

struct GeantParticle { double E, r, theta, dxds, dyds, polx, poly, polz; };

struct ERBin {
  public:
    size_t count;
    double area;
  private:
    std::vector<GeantParticle> p_list;
    std::optional<std::vector<XYBin>> xy_bins;

  public:
    ERBin() = default;
    void add_point(const GeantParticle& p);
    double density(size_t num_elec_in) const;
    void write_momenta(const char * output_dir,
        double E_elec, double T,
        double E, double r, std::ofstream& binlist) const;
    void bin_momenta();
    double get_average(const std::function<double(const GeantParticle&)>& property) const;
};

