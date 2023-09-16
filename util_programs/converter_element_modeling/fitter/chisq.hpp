#pragma once
#include <vector>
#include "cauchy.hpp"
#include "meta_fit.hpp"
#include "parser.hpp"
void output_chisq(const MetaFitResults& mf,
    const std::vector<CauchyPoint>& cauchy_results,
    const std::vector<XYBinnedData>& xydata,
    const SimSettings& settings);
