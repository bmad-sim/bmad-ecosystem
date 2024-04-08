cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set (EXE_SPECS
  ## cmake_files/cmake.converter_simulation
  ## cmake_files/cmake.converter_fitter
  cmake_files/cmake.beam_file_translate_format
  cmake_files/cmake.controller_response_plot
  cmake_files/cmake.ptc_flat_file_to_bmad
  cmake_files/cmake.make_a_matching_knob
  cmake_files/cmake.wake_fit
  cmake_files/cmake.bmad_to_astra
  cmake_files/cmake.bmad_to_blender
  cmake_files/cmake.bmad_to_csrtrack
  cmake_files/cmake.bmad_to_gpt
  cmake_files/cmake.bmad_to_julia
  cmake_files/cmake.bmad_to_mad_sad_elegant
  cmake_files/cmake.bmad_to_merlin
  cmake_files/cmake.bmad_to_slicktrack
  cmake_files/cmake.compare_tracking_methods_plot
  cmake_files/cmake.compare_tracking_methods_text
  cmake_files/cmake.element_attributes
  cmake_files/cmake.f77_to_f90
  cmake_files/cmake.generalized_gradient_fit
  cmake_files/cmake.lattice_cleaner
  cmake_files/cmake.errors_mad_to_bmad
  cmake_files/cmake.mais_ripken
  cmake_files/cmake.photon_init_plot
  cmake_files/cmake.photon_surface_data
  cmake_files/cmake.ptc_profiler
  cmake_files/cmake.sad_to_bmad_postprocess
  cmake_files/cmake.cartesian_map_fit
)

set (TEST_EXE_SPECS
  cmake_files/cmake.converter_simulation
  cmake_files/cmake.converter_fitter
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
