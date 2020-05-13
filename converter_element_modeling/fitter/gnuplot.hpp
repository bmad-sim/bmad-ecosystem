#pragma once
#include <vector>
#include "meta_fit.hpp"

void output_cauchy_gp(const char* foldername, const DataPoint& p);
void output_meta_cauchy_gp(const char* foldername, const MetaFitResults& mf, double E, double r);
template<fitType T>
void output_metafit_gp(const FitResults& fit, const std::vector<DataPoint>& data, const char* dir, double xpt);
