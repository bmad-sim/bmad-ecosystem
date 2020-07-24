#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "bmad.hpp"


void write_bmad_file(const std::vector<MetaFitResults>& metafits, const SimSettings& settings) {
  // Writes the converter file for the given meta fits
  std::string bmad_file_name = settings.output_directory + "/converter.bmad";
  std::ofstream bmad_file;
  bmad_file.open(bmad_file_name);
  if (!bmad_file) {
    std::cerr << "ERROR: could not open \"" << bmad_file_name << "\"\n";
    return;
  }

  size_t mf_ix=0;

  double T = metafits[0].T;
  bmad_file << "distribution = {\n";
  bmad_file << TAB<1>() << "material = " << settings.target_material << ",\n";
  bmad_file << TAB<1>() << "species_out = positron" << ",\n";
  bmad_file << TAB<1>() << "thickness = " << T;
  for (const auto& mf : metafits) {
    if (mf.T != T) { // start a new distribution
      T = mf.T;
      bmad_file << "},\ndistribution = {\n";
      bmad_file << TAB<1>() << "material = tungsten" << ",\n";
      bmad_file << TAB<1>() << "species_out = positron" << ",\n";
      bmad_file << TAB<1>() << "thickness = " << mf.T;
    }
    bmad_file << ",\n" << TAB<1>() << "sub_distribution = {\n";
    bmad_file << TAB<2>() << "pc_in = " << mf.Ein << ",\n";
    bmad_file << TAB<2>() << "spin_in = [";
    if (settings.polarization_in.size()) {
      bmad_file << settings.polarization_in[0] << ", ";
      bmad_file << settings.polarization_in[1] << ", ";
      bmad_file << settings.polarization_in[2] << "],\n";
    } else {
      bmad_file << "0, 0, 0],\n";
    }
    bmad_file << TAB<2>() << "prob_pc_r = {\n";
    bmad_file << TAB<3>() << "r_values = [0.0, "; // fix r=0
    // Write r values to table (no comma after last one)
    auto num_rows = mf.er_table.pc_vals.size();
    auto num_cols = mf.er_table.r_vals.size();
    for (unsigned i=0; i<num_cols-1; i++) bmad_file << mf.er_table.r_vals[i] << ", ";
    bmad_file << mf.er_table.r_vals[num_cols-1] << "]";
    // Copy over each row of the table
    for (unsigned i=0; i<num_rows; i++) {
      bmad_file << ",\n" << TAB<3>() << "row = {pc_out = " << mf.er_table.pc_vals[i];
      bmad_file << ", prob = [0.0, "; // r=0 -> prob = 0
      for (unsigned j=0; j<num_cols-1; j++) bmad_file << 1e-4 * mf.er_table.data[i*num_cols + j] << ", "; // convert to prob / (eV * m)
      bmad_file << 1e-4 * mf.er_table.data[(i+1)*num_cols-1] << "]}";
    }
    bmad_file << "\n" << TAB<2>() << "},\n";


    // Polarization table
    bmad_file << TAB<2>() << "spin_z_out = {\n";
    bmad_file << TAB<3>() << "r_values = [";
    // Write r values to table (no comma after last one)
    num_rows = mf.polz_table.pc_vals.size();
    num_cols = mf.polz_table.r_vals.size();
    for (unsigned i=0; i<num_cols-1; i++) bmad_file << mf.polz_table.r_vals[i] << ", ";
    bmad_file << mf.polz_table.r_vals[num_cols-1] << "]";
    // Copy over each row of the table
    for (unsigned i=0; i<num_rows; i++) {
      bmad_file << ",\n" << TAB<3>() << "row = {pc_out = " << mf.polz_table.pc_vals[i];
      bmad_file << ", spin_z = ["; // r=0 -> prob = 0
      for (unsigned j=0; j<num_cols-1; j++) bmad_file << mf.polz_table.data[i*num_cols + j] << ", ";
      bmad_file << mf.polz_table.data[(i+1)*num_cols-1] << "]}";
    }
    bmad_file << "\n" << TAB<2>() << "},\n";

    // Output the fit parameters
    bmad_file << TAB<2>() << "direction_out = {\n";
    output_bmad<fitType::CX>(mf.cx, bmad_file, COMMA_AFTER);
    output_bmad<fitType::AX>(mf.ax, bmad_file, COMMA_AFTER);
    output_bmad<fitType::AY>(mf.ay, bmad_file, COMMA_AFTER);
    output_bmad<fitType::BETA>(mf.beta, bmad_file, COMMA_AFTER);
    output_bmad<fitType::DXDS_MIN>(mf.dxds_min, bmad_file, COMMA_AFTER);
    output_bmad<fitType::DXDS_MAX>(mf.dxds_max, bmad_file, COMMA_AFTER);
    output_bmad<fitType::DYDS_MAX>(mf.dyds_max, bmad_file, NO_COMMA);

    bmad_file << TAB<2>() << "}\n" << TAB<1>() << "}";
  }
  bmad_file << "\n}\n";

  return;
}

