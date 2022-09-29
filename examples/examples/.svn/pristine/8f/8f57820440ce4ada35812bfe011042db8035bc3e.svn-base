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
  ## cmake.parallel_track_example  ## This causes problems in distribution builds.
  ## cmake.spacecharge_example # Depends on open_spacecharge, under development
  ## cmake.coarray_example         ## For gfortran need to modify cmake.coarray_example to link in coarray library.
  cmake.spin_amplitude_dependent_tune
  cmake.ptc_layout_example
  cmake.construct_taylor_map
  cmake.dispersion_simulation
  cmake.lapack_examples         ## This is broken with lapack version upgrade.
  cmake.lattice_geometry_example
  cmake.multi_turn_tracking_example
  cmake.beam_track_example
  cmake.bmad_to_opal_example
  cmake.em_field_query_example
  cmake.csr_example
  cmake.particle_track_example
  cmake.simple_bmad_program
  cmake.plot_example
  cmake.tune_tracker
  cmake.ptc_spin_orbital_normal_form
  cmake.spin_matching
)

# Only build mpi_mp if both ACC_ENABLE_OPENMP and ACC_ENABLE_MPI are set to "Y".

if ("$ENV{ACC_ENABLE_OPENMP}" AND "$ENV{ACC_ENABLE_MPI}")
  list (APPEND EXE_SPECS "cmake.mpi_mp")
endif ()


include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
