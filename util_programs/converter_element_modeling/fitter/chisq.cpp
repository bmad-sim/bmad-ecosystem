#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <fstream>
#include "chisq.hpp"

template<typename T>
constexpr T sqr(T x) { return x*x; }

double cauchy_eval(const CauchyPoint& c, double x, double y) {
  // Evaluates the given cauchy fit at the specified (x,y) point
  return c.amp * ( 1 + c.beta * x) / ( 1 + sqr(c.ax*(x-c.cx)) + sqr(c.ay*y));
}

struct ChisqResult {
  double pc, r, chisq;
};

std::ostream& operator<<(std::ostream& os, const ChisqResult& chi) {
  return (os << "pc = " << chi.pc/1e6
             << " MeV, r = " << chi.r*1e2
             << " cm, chisq = " << chi.chisq);
}

void output_chisq(const MetaFitResults& mf,
    const std::vector<CauchyPoint>& cauchy_results,
    const std::vector<XYBinnedData>& xydata,
    const SimSettings& settings) {
  // Performs goodness-of-fit analysis for the cauchy fits and meta fits,
  // and writes output to E..._T..._chisq.txt


  // Initial setup needed to get er bin areas
  std::vector<double> e_edges = {settings.out_pc_min}, r_edges = {0.0};
  for (auto e_val : mf.er_table.pc_vals)
    e_edges.push_back(e_edges.back() + 2*(e_val/1e6 - e_edges.back()));
  for (auto r_val : mf.er_table.r_vals)
    r_edges.push_back(r_edges.back() + 2*(r_val*1e2 - r_edges.back()));

  std::vector<double> er_areas;
  er_areas.resize(mf.er_table.data.size());
  auto e_len = mf.er_table.pc_vals.size();
  auto r_len = mf.er_table.r_vals.size();
  for (unsigned i=0; i<e_len; i++)
    for (unsigned j=0; j<r_len; j++)
      er_areas[i*r_len + j] = (e_edges[i+1]-e_edges[i])*(r_edges[j+1]-r_edges[j]);

  // Lambda for computing chi-square
  auto compute_chisq = [](const CauchyPoint& c, const std::vector<BinPoint>& xy_bins) {
    double min_nonzero = (std::min_element(xy_bins.begin(), xy_bins.end(),
        [](const auto& bin1, const auto& bin2) {
          return (bin1.density < bin2.density) && (bin1.density != 0);
        }))->density;

    double chisq = std::accumulate(xy_bins.begin(), xy_bins.end(), 0.0,
        [&c,&min_nonzero](double sum, const auto& bin) {
          return sqr(cauchy_eval(c, bin.x, bin.y) - bin.density) / std::max(bin.density, min_nonzero);
        });
    chisq /= xy_bins.size();
    return chisq;
  };


  // Goodness of fit for cauchy fits vs binned data from geant
  unsigned ix = 0; // xydata, mf.er_table, and cauchy_results are all sorted the same,
                   // so a single index can be used to access corresponding elements
  std::vector<ChisqResult> cauchy_chisq;
  for (const auto& c : cauchy_results) {
    const auto& xy_bins = xydata[ix].bins;

    cauchy_chisq.push_back({c.E, c.r, compute_chisq(c, xy_bins)});
    ix++;
  }


  // Goodness of fit for meta fits vs data from geant
  std::vector<ChisqResult> meta_chisq;
  ix=0;
  for (const auto& c : cauchy_results) {
    const auto& xy_bins = xydata[ix].bins;
    auto meta_c = c; // take amp from c

    // Determine which 1D/2D fit should be used
    const auto& low_e = mf.cx.low_e_fits; // for convenience
    auto fit1d_it = std::find_if(low_e.begin(), low_e.end(),
        [&c](const auto& fit1d) { return fit1d.E >= c.E; });
    // ^ this works because low_e is sorted in order of increasing E
    if (fit1d_it == low_e.end()) { // 2D
      meta_c.cx   = eval<fitType::CX,2>   (mf.cx.high_e_fit,    c);
      meta_c.ax   = eval<fitType::AX,2>   (mf.ax.high_e_fit,    c);
      meta_c.ay   = eval<fitType::AY,2>   (mf.ay.high_e_fit,    c);
      meta_c.beta = eval<fitType::BETA,2> (mf.beta.high_e_fit,  c);
    } else { // 1D
      meta_c.cx   = eval<fitType::CX,1>   (*fit1d_it, c);
      meta_c.ax   = eval<fitType::AX,1>   (*fit1d_it, c);
      meta_c.ay   = eval<fitType::AY,1>   (*fit1d_it, c);
      meta_c.beta = eval<fitType::BETA,1> (*fit1d_it, c);
    }

    meta_chisq.push_back({c.E, c.r, compute_chisq(meta_c, xy_bins)});
    ix++;
  }


  // Output results
  std::ofstream chisq_file;
  char filename[200];
  sprintf(filename, "/E%0.0lf_T%0.3lf_chisq.txt", mf.Ein/1e6, mf.T*1e2);
  chisq_file.open(settings.output_directory + filename);

  auto chisq_output = [&chisq_file](std::vector<ChisqResult>& chisq_vec) {
    auto sort_func = [](const auto& cr1, const auto& cr2) { return cr1.chisq < cr2.chisq; };

    auto [min_chisq, max_chisq] = std::minmax_element(chisq_vec.begin(), chisq_vec.end(), sort_func);
    chisq_file << "\tBest Chi-square: " << *min_chisq << '\n';
    chisq_file << "\tWorst Chi-square: " << *max_chisq << '\n';

    double avg_chi2 = std::accumulate(chisq_vec.begin(), chisq_vec.end(), 0.0,
          [](double sum, const auto& cr) { return sum + cr.chisq; })/chisq_vec.size();
    chisq_file << "\tAverage Chi-square: " << avg_chi2 << '\n';

    chisq_file << "\tBins with Chi-square more than five times the average:\n";
    chisq_vec.erase(std::remove_if(chisq_vec.begin(), chisq_vec.end(),
                          [avg_chi2](const auto& cr) { return cr.chisq < 5*avg_chi2; }),
                      chisq_vec.end());
    std::sort(chisq_vec.begin(), chisq_vec.end(), sort_func);
    for (const auto& cr : chisq_vec)
      chisq_file << "\t\t" << cr << '\n';
  };


  chisq_file << "Goodness of fit for cauchy fits:\n";

  chisq_output(cauchy_chisq);

  const char* SEPARATOR = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
  chisq_file << SEPARATOR;

  chisq_file << "Goodness of fit for meta fits:\n";

  chisq_output(meta_chisq);




  return;
}
