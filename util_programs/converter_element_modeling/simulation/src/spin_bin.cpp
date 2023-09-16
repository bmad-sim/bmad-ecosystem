#include <fstream>
#include <numeric>
#include <algorithm>
#include <cmath>
#include "spin_bin.hpp"

void SpinBinner::bin_spins() {
  //constexpr double min_spin = -1;
  //constexpr double max_spin = 1;
  constexpr unsigned num_bins = 200;
  //constexpr double ds = (max_spin - min_spin) / num_bins;

  spin_bins.emplace();
  spin_bins->resize(num_bins);

  // Determine bounds for sx, sy, sz
  auto sortx = [](const GeantParticle& p1, const GeantParticle& p2) {
    return p1.polx < p2.polx;
  };
  auto sorty = [](const GeantParticle& p1, const GeantParticle& p2) {
    return p1.poly < p2.poly;
  };
  auto sortz = [](const GeantParticle& p1, const GeantParticle& p2) {
    return p1.polz < p2.polz;
  };

  auto pnum = p_list.size();
  if (pnum < 20) return;

  std::sort(p_list.begin(), p_list.end(), sortx);
  const double min_x = p_list[pnum*0.05].polx;
  const double max_x = p_list[pnum*0.95].polx;
  const double dsx = (max_x - min_x) / num_bins;

  std::sort(p_list.begin(), p_list.end(), sorty);
  const double min_y = p_list[pnum*0.05].poly;
  const double max_y = p_list[pnum*0.95].poly;
  const double dsy = (max_y - min_y) / num_bins;

  std::sort(p_list.begin(), p_list.end(), sortz);
  const double min_z = p_list[pnum*0.05].polz;
  const double max_z = p_list[pnum*0.95].polz;
  const double dsz = (max_z - min_z) / num_bins;

  // Initialize bin values
  unsigned ix = 0;
  for (auto& sb : *spin_bins) {
    sb.sx_val = min_x + ix * dsx + 0.5*dsx;
    sb.sy_val = min_y + ix * dsy + 0.5*dsy;
    sb.sz_val = min_z + ix * dsz + 0.5*dsz;
    sb.sx_count = 0;
    sb.sy_count = 0;
    sb.sz_count = 0;
    ix++;
  }

  // Bin particles
  auto in_sx_range = [=](const GeantParticle& p) { return (p.polx > min_x && p.polx < max_x); };
  auto in_sy_range = [=](const GeantParticle& p) { return (p.poly > min_y && p.poly < max_y); };
  auto in_sz_range = [=](const GeantParticle& p) { return (p.polz > min_z && p.polz < max_z); };
  auto get_sx_bin_num = [=](const GeantParticle& p) {
    return std::min(static_cast<unsigned>((p.polx - min_x)/dsx), num_bins-1);
  };
  auto get_sy_bin_num = [=](const GeantParticle& p) {
    return std::min(static_cast<unsigned>((p.poly - min_y)/dsy), num_bins-1);
  };
  auto get_sz_bin_num = [=](const GeantParticle& p) {
    return std::min(static_cast<unsigned>((p.polz - min_z)/dsz), num_bins-1);
  };
  for (const auto& p : p_list) {
    if (in_sx_range(p)) (*spin_bins)[get_sx_bin_num(p)].sx_count++;
    if (in_sy_range(p)) (*spin_bins)[get_sy_bin_num(p)].sy_count++;
    if (in_sz_range(p)) (*spin_bins)[get_sz_bin_num(p)].sz_count++;
  }
}


void SpinBinner::write_spins(const char * output_dir,
    double E_elec, double T,
    double E, double r, std::ofstream& binlist) const {
  if (!spin_bins) {
    std::cerr << "Warning: sx,sy,sz have not been binned!\n";
    return;
  }
  auto& bins = *spin_bins;
  // Write binned data to a file
  std::ofstream binfile;
  char binfile_name[100];
  sprintf(binfile_name, "%s/spin_dat/E%0.0lf_T%0.3lf/E%0.2lf_r%0.3lf_bin.dat",
      output_dir, E_elec, T, E, r);
  binfile.open(binfile_name);
  for (auto& b : bins) {
    binfile << b.sx_val << '\t'
      << b.sx_count << '\t'
      << b.sy_val << '\t'
      << b.sy_count << '\t'
      << b.sz_val << '\t'
      << b.sz_count << '\n';
  }
  // Write the filename to binlist
  binlist << binfile_name << '\n';
}



