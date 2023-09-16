cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set(INC_DIRS 
)

set (SRC_DIRS
)

#OpenMP Flags
#set (FFLAGS "-openmp")
#set (LINK_FLAGS "-openmp")

set (EXE_SPECS 
  ## cmake_files/cmake.parallel_track_example  ## This causes problems in distribution builds.
  ## cmake_files/cmake.coarray_example         ## For gfortran need to modify cmake_files/cmake.coarray_example to link in coarray library.
  cmake_files/cmake.spin_amplitude_dependent_tune
  cmake_files/cmake.ptc_layout_example
  cmake_files/cmake.construct_taylor_map
  cmake_files/cmake.dispersion_simulation
  cmake_files/cmake.lapack_examples         ## This is broken with lapack version upgrade.
  cmake_files/cmake.lattice_geometry_example
  cmake_files/cmake.multi_turn_tracking_example
  cmake_files/cmake.beam_track_example
  cmake_files/cmake.bmad_to_opal_example
  cmake_files/cmake.em_field_query_example
  cmake_files/cmake.csr_example
  cmake_files/cmake.particle_track_example
  cmake_files/cmake.simple_bmad_program
  cmake_files/cmake.plot_example
  cmake_files/cmake.tune_tracker
  cmake_files/cmake.ptc_spin_orbital_normal_form
  cmake_files/cmake.spin_matching
)

# Only build mpi_mp if both ACC_ENABLE_OPENMP and ACC_ENABLE_MPI are set to "Y".

if ("$ENV{ACC_ENABLE_OPENMP}" AND "$ENV{ACC_ENABLE_MPI}")
  list (APPEND EXE_SPECS "cmake_files/cmake.mpi_mp")
endif ()


include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
