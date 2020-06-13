#pragma once
#include <vector>
#include "cauchy.hpp"
#include "meta_fit.hpp"

void output_cauchy_gp(const char* foldername, const CauchyPoint& p);
void output_meta_cauchy_gp(const char* foldername, const MetaFitResults& mf, double E, double r);
template<fitType T>
void output_metafit_gp(const FitResults& fit, const std::vector<CauchyPoint>& data, const char* dir, double xpt);
void output_master_gp(const char* dir);
