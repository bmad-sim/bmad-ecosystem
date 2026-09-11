if (DEFINED ENV{ACC_SET_GMAKE_JOBS} AND DEFINED ENV{GMAKE})
  cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
  project(ACC)

  SET (LIBNAME sim_utils)

  SET (CREATE_SHARED true)

  SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

  # If this is a conda build, fgsl is contained
  # one level deeper in the include folder
  IF ($ENV{ACC_CONDA_BUILD})
    SET (INC_DIRS
      $ENV{ACC_CONDA_PATH}/include/
      $ENV{ACC_CONDA_PATH}/include/fgsl/
    )
  ELSE ()
    SET (INC_DIRS
    )
  ENDIF ()

  SET (SRC_DIRS
    file 
    interfaces 
    io 
    math 
    matrix 
    misc
    optimizers
    plot 
    string
    geodesic_lm
    special_functions
  )

  SET (SHARED_DEPS
    ${PLOT_LINK_LIBS}
    $ENV{ACC_BMAD_LINK_LIBS}
    X11
    readline
  )

  include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
else()
  foreach(d file interfaces io math matrix misc optimizers plot string geodesic_lm special_functions)
    file(GLOB s "${d}/*.f")
    set(src ${src} ${s})
    file(GLOB s "${d}/*.f90")
    set(src ${src} ${s})
    file(GLOB s "${d}/*.c")
    set(src ${src} ${s})
    file(GLOB s "${d}/*.cc")
    set(src ${src} ${s})
    file(GLOB s "${d}/*.cpp")
    set(src ${src} ${s})
  endforeach()

  find_library(lapack95_lib NAMES lapack95)
  find_path(lapack95_inc NAMES f95_lapack.mod HINTS "${CMAKE_PREFIX_PATH}")
  if (lapack95_lib AND lapack95_inc)
    add_library(lapack95::lapack95 UNKNOWN IMPORTED)
    set_target_properties(lapack95::lapack95 PROPERTIES IMPORTED_LOCATION ${lapack95_lib})
    target_include_directories(lapack95::lapack95 INTERFACE "${lapack95_inc}")
  else()
    message(FATAL_ERROR "Could not find lapack95 library")
  endif()

  pkg_check_modules(plplot_fortran REQUIRED IMPORTED_TARGET plplot-fortran)
  pkg_check_modules(fgsl REQUIRED IMPORTED_TARGET fgsl)
  find_package(X11 REQUIRED)
  pkg_check_modules(fftw3 REQUIRED IMPORTED_TARGET fftw3)
  pkg_check_modules(readline REQUIRED IMPORTED_TARGET readline)

  add_library(sim_utils SHARED ${src})
  set_target_properties(sim_utils PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  )
  set_source_files_properties(
    plot/quick_plot.f90
    plot/plplot_interface.f90
    plot/pgsubs.f90
    plot/noplot_interface.f90
    plot/pgplot_interface.f90
    file/filename_mod.f90
    PROPERTIES Fortran_PREPROCESS ON)
  set_source_files_properties(
    plot/quick_plot.f90
    plot/plplot_interface.f90
    PROPERTIES COMPILE_DEFINITIONS CESR_PLPLOT)
  target_link_libraries(sim_utils PUBLIC
    bmad_options
    lapack95::lapack95
    PkgConfig::fgsl
    PkgConfig::plplot_fortran
    X11::X11
    PkgConfig::readline
    PkgConfig::fftw3)
  if (ENABLE_OPENMP)
    find_library(FFTW3_OMP_LIB NAMES "fftw3_omp" REQUIRED)
    add_library(fftw3_omp UNKNOWN IMPORTED)
    set_target_properties(fftw3_omp PROPERTIES IMPORTED_LOCATION ${FFTW3_OMP_LIB})
    target_link_libraries(sim_utils PUBLIC fftw3_omp OpenMP::OpenMP_Fortran)
  endif()
  target_include_directories(sim_utils PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:modules>
  )
  install(TARGETS sim_utils EXPORT Bmad)
  foreach(d file interfaces io math matrix misc optimizers plot string geodesic_lm special_functions)
    install(DIRECTORY ${d} DESTINATION src/sim_utils
      FILES_MATCHING
      PATTERN "*.f"
      PATTERN "*.f90"
      PATTERN "*.c"
      PATTERN "*.cc"
      PATTERN "*.cpp")
  endforeach()
  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/" DESTINATION modules FILES_MATCHING PATTERN "*.mod")
endif()
