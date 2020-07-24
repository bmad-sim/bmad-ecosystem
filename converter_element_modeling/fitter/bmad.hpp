#pragma once
#include <vector>
#include "meta_fit.hpp"
#include "parser.hpp"

void write_bmad_file(const std::vector<MetaFitResults>& metafits, const SimSettings& settings);
