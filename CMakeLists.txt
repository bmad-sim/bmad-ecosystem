cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set (EXE_SPECS
  ## cmake.converter_simulation
  ## cmake.converter_fitter
  cmake.beam_file_translate_format
  cmake.controller_function_plot
  cmake.ptc_flat_file_to_bmad
  cmake.make_a_matching_knob
  cmake.wake_fit
  cmake.bmad_to_astra
  cmake.bmad_to_autocad
  cmake.bmad_to_blender
  cmake.bmad_to_csrtrack
  cmake.bmad_to_gpt
  cmake.bmad_to_mad_and_sad
  cmake.bmad_to_merlin
  cmake.compare_tracking_methods_plot
  cmake.compare_tracking_methods_text
  cmake.element_attributes
  cmake.f77_to_f90
  cmake.lattice_cleaner
  cmake.sad_to_bmad_postprocess
  cmake.cartesian_map_fit
)

set (TEST_EXE_SPECS
  cmake.converter_simulation
  cmake.converter_fitter
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
