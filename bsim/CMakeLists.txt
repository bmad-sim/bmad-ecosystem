if (DEFINED ENV{BUILD_TYPE})
  set (LIBNAME bsim)
  cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
  project(ACC)

  IF ($ENV{ACC_CONDA_BUILD})
    SET (INC_DIRS
      $ENV{ACC_CONDA_PATH}/include
      $ENV{ACC_CONDA_PATH}/include/xraylib/
    )
  ELSE ()
    SET (INC_DIRS
    )
  ENDIF ()

  set (SRC_DIRS
    code
    code_synrad
    code_synrad3d
    synrad3d/custom
    modules
  )

  set (EXE_SPECS
    #cmake_files/cmake.ansga2
    #cmake_files/cmake.aspea2

    cmake_files/cmake.e_cooling
    cmake_files/cmake.long_term_tracking
    cmake_files/cmake.tune_scan
    cmake_files/cmake.dark_current_tracker
    cmake_files/cmake.envelope_ibs
    cmake_files/cmake.frequency_map
    cmake_files/cmake.bbu
    cmake_files/cmake.dynamic_aperture
    cmake_files/cmake.ibs_linac
    cmake_files/cmake.ibs_ring
    cmake_files/cmake.sodom2
    cmake_files/cmake.spin_stroboscope
    cmake_files/cmake.synrad
    cmake_files/cmake.synrad3d
    cmake_files/cmake.synrad_aperture_to_wall
    cmake_files/cmake.tune_plane_res_plot
    cmake_files/cmake.wall_generator
  )

  #set (FFLAGS "-qopenmp")
  #set (LINK_FLAGS "-qopenmp")

  if($ENV{ACC_ENABLE_MPI})
    list(APPEND EXE_SPECS "cmake_files/cmake.moga")
    list(APPEND EXE_SPECS "cmake_files/cmake.aperture_by_tracking")
    list(APPEND EXE_SPECS "cmake_files/cmake.touschek_background")
    list(APPEND EXE_SPECS "cmake_files/cmake.long_term_tracking_mpi")
    list(APPEND EXE_SPECS "cmake_files/cmake.tune_scan_mpi")
  endif()

  include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
else()
  file(GLOB src
    "code/*.f90"
    "code_synrad/*.f90"
    "code_synrad3d/*.f90"
    "synrad3d/custom/*.f90"
    "modules/*.f90"
  )
  add_library(bsim SHARED ${src})
  set_target_properties(bsim PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  )
  target_link_libraries(bsim bmad)
  target_include_directories(bsim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
    $<INSTALL_INTERFACE:modules>
  )

  add_executable(e_cooling e_cooling/e_cooling.f90 e_cooling/e_cooling_mod.f90 e_cooling/track1_bunch_hook.f90)
  target_link_libraries(e_cooling bmad)

  add_executable(long_term_tracking long_term_tracking/long_term_tracking.f90)
  target_link_libraries(long_term_tracking bsim)

  add_executable(tune_scan tune_scan/tune_scan.f90)
  target_link_libraries(tune_scan bsim)

  add_executable(dark_current_tracker dark_current_tracker/dark_current_mod.f90 dark_current_tracker/dark_current_tracker.f90)
  target_link_libraries(dark_current_tracker bmad)

  add_executable(envelope_ibs envelope_ibs/envelope_ibs.f90)
  target_link_libraries(envelope_ibs bmad)

  add_executable(frequency_map frequency_map/frequency_map.f90)
  target_link_libraries(frequency_map bsim)

  add_executable(bbu bbu/bbu_program.f90)
  target_link_libraries(bbu bsim)

  add_executable(dynamic_aperture
    dynamic_aperture/da_program_mod.f90
    dynamic_aperture/dynamic_aperture.f90
    dynamic_aperture/lat_make_mat6_hook.f90
    dynamic_aperture/track_many_hook.f90
  )
  target_link_libraries(dynamic_aperture bsim)

  add_executable(ibs_linac ibs_linac/ibs_linac.f90)
  target_link_libraries(ibs_linac bmad)

  add_executable(ibs_ring ibs_ring/ibs_ring.f90)
  target_link_libraries(ibs_ring bmad)

  add_executable(sodom2 sodom2/sodom2.f90 sodom2/sodom2_mod.f90)
  target_link_libraries(sodom2 bmad)

  add_executable(spin_stroboscope spin_stroboscope/spin_stroboscope.f90 spin_stroboscope/track1_custom.f90)
  target_link_libraries(spin_stroboscope bmad)

  add_executable(synrad synrad/synrad.f90 synrad/synrad_plot_mod.f90)
  target_link_libraries(synrad bsim)

  add_executable(synrad3d synrad3d/synrad3d.f90)
  target_link_libraries(synrad3d bsim)

  add_executable(synrad_aperture_to_wall synrad_aperture_to_wall/synrad_aperture_to_wall.f90)
  target_link_libraries(synrad_aperture_to_wall bmad)

  add_executable(tune_plane_res_plot tune_plane_res_plot/tune_plane_res_mod.f90 tune_plane_res_plot/tune_plane_res_plot.f90)
  target_link_libraries(tune_plane_res_plot sim_utils)

  add_executable(wall_generator wall_generator/wall_generator.f90)
  target_link_libraries(wall_generator bmad)

  install(TARGETS bsim EXPORT Bmad)
  install(TARGETS
    e_cooling
    long_term_tracking
    tune_scan
    dark_current_tracker
    envelope_ibs
    frequency_map
    bbu
    dynamic_aperture
    ibs_linac
    ibs_ring
    sodom2
    spin_stroboscope
    synrad
    synrad3d
    synrad_aperture_to_wall
    tune_plane_res_plot
    wall_generator)
  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/" DESTINATION modules FILES_MATCHING PATTERN "*.mod")
  foreach(d code code_synrad code_synrad3d synrad3d/custom modules)
    install(DIRECTORY ${d} DESTINATION src/bsim FILES_MATCHING PATTERN "*.f90")
  endforeach()
endif()
