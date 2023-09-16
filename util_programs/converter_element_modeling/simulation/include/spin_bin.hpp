#pragma once
#include <vector>
#include <optional>
#include <iostream>
#include "bin.hpp"

struct SpinBin {
  double sx_val, sy_val, sz_val;
  unsigned sx_count, sy_count, sz_count;
};

class SpinBinner {
  private:
    std::vector<GeantParticle> p_list;
    std::optional<std::vector<SpinBin>> spin_bins;

  public:
    SpinBinner() = default;
    inline void add_point(const GeantParticle& p) { p_list.push_back(p); }
    void bin_spins();
    void write_spins(const char * output_dir,
        double E_elec, double T,
        double E, double r, std::ofstream& binlist) const;
};
