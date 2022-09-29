cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(ACC)

set (LIBNAME mad_tpsa)

set(INC_DIRS
 code
)

set (SRC_DIRS
  code
  code/sse
  interface
)

set (EXE_SPECS
  cmake.gtpsa_ex1_c
  cmake.gtpsa_ex2_c
  cmake.gtpsa_ex3_c
  cmake.gtpsa_ex4_c
  cmake.gtpsa_ex5_c
  cmake.gtpsa_ex6_c
  cmake.gtpsa_scan_c
  cmake.gtpsa_ex1_f90
  cmake.gtpsa_ex2_f90
  cmake.gtpsa_ex3_f90
  cmake.gtpsa_ex4_f90
  cmake.gtpsa_ex5_f90
  cmake.gtpsa_ex6_f90
  cmake.tpsa_test
)

#set (FFLAGS "-pg")  # For profiling
#set (CFLAGS "-pg")

set (CFLAGS "-ffast-math -ftree-vectorize")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
