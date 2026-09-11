if (DEFINED ENV{BUILD_TYPE})
  cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
  project(ACC)

  SET (LIBNAME bmad)

  SET (CREATE_SHARED true)

  SET (PLOT_LINK_LIBS $ENV{PLOT_LINK_LIBS})

  # If this is a conda build, xraylib is contained
  # one level deeper in the include folder

  IF ($ENV{ACC_CONDA_BUILD})
    SET (INC_DIRS
      $ENV{ACC_CONDA_PATH}/include
      $ENV{ACC_CONDA_PATH}/include/xraylib/
      $ENV{ACC_CONDA_PATH}/include/fgsl/
    )
  ELSE ()
    SET (INC_DIRS
    )
  ENDIF ()

  SET (SRC_DIRS 
    code 
    ### custom  # Using function pointers now
    geometry
    interface
    modules
    multiparticle
    output
    ptc
    photon
    parsing
    low_level
    hdf5
    space_charge
    spin
  )

  SET (SHARED_DEPS
    sim_utils
    $ENV{ACC_BMAD_LINK_LIBS}
    X11
    readline
  )

  include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
else()
  foreach(d code geometry interface modules multiparticle output ptc photon parsing low_level hdf5 space_charge spin)
    file(GLOB s "${d}/*.f90")
    set(src ${src} ${s})
  endforeach()

  # Have I mentioned that HDF5 is a pain?
  find_package(HDF5 COMPONENTS Fortran HL)
  if (HDF5_FOUND)
    if(ENABLE_SHARED OR ENABLE_SHARED_ONLY)
      if(TARGET hdf5_hl_fortran-shared)
        set(HDF5_TARGET hdf5_hl_fortran-shared)
      elseif(TARGET hdf5_hl_fortran-static)
        set(HDF5_TARGET hdf5_hl_fortran-static)
      endif()
    else()
      if(TARGET hdf5_hl_fortran-static)
        set(HDF5_TARGET hdf5_hl_fortran-static)
      elseif(TARGET hdf5_hl_fortran-shared)
        set(HDF5_TARGET hdf5_hl_fortran-shared)
      endif()
    endif()
    if (NOT DEFINED HDF5_TARGET)
      if(TARGET hdf5::hdf5_hl_fortran AND TARGET hdf5::hdf5_fortran)
        # Sometimes the hdf5_hl_fortran target gets screwed up
        get_property(hl_loc TARGET hdf5::hdf5_hl_fortran PROPERTY IMPORTED_LOCATION SET)
        if (hl_loc)
          get_property(hl_loc TARGET hdf5::hdf5_hl_fortran PROPERTY IMPORTED_LOCATION)
          set(prop IMPORTED_LOCATION)
        endif()
        if (NOT hl_loc)
          get_property(hl_loc TARGET hdf5::hdf5_hl_fortran PROPERTY IMPORTED_CONFIGURATIONS SET)
          if (hl_loc)
            set(hl_loc OFF)
            get_property(configs TARGET hdf5::hdf5_hl_fortran PROPERTY IMPORTED_CONFIGURATIONS)
            foreach(c ${configs})
              get_property(hl_loc_t TARGET hdf5::hdf5_hl_fortran PROPERTY IMPORTED_LOCATION_${c} SET)
              if (hl_loc_t)
                get_property(hl_loc_t TARGET hdf5::hdf5_hl_fortran PROPERTY IMPORTED_LOCATION_${c})
                set(prop IMPORTED_LOCATION_${c})
              endif()
              if (hl_loc_t)
                set(hl_loc ${hl_loc_t})
              endif()
            endforeach()
          endif()
        endif()
        if (NOT hl_loc AND prop)
          # The hdf5_hl_fortran target gets screwed up
          foreach(lib ${HDF5_Fortran_HL_LIBRARIES})
            if (NOT hl_loc AND lib AND lib MATCHES hl_fortran)
              set(hl_loc ${lib})
              set_target_properties(hdf5::hdf5_hl_fortran PROPERTIES ${prop} "${lib}")
            endif()
          endforeach()
        endif()
        if(NOT hl_loc)
          message(FATAL_ERROR "Problem with HDF5 hl_fortran Target")
        endif()
        set(HDF5_TARGET hdf5::hdf5_hl_fortran hdf5::hdf5_fortran)
      else()
        # Old CMake, we need to roll our own targets...
        if (NOT TARGET hdf5::hdf5_fortran)
          add_library(hdf5::hdf5_fortran UNKNOWN IMPORTED)
          foreach(lib ${HDF5_Fortran_LIBRARIES})
            if(lib MATCHES hdf5_fortran)
              list(REMOVE_ITEM HDF5_Fortran_LIBRARIES ${lib})
              set_target_properties(hdf5::hdf5_fortran PROPERTIES IMPORTED_LOCATION ${lib})
              foreach(d ${CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES})
	        list(REMOVE_ITEM HDF5_Fortran_INCLUDE_DIRS ${d})
              endforeach()
              target_include_directories(hdf5::hdf5_fortran INTERFACE ${HDF5_Fortran_INCLUDE_DIRS})
              target_link_libraries(hdf5::hdf5_fortran INTERFACE ${HDF5_Fortran_LIBRARIES})
              break()
            endif()
          endforeach()
        endif()
        if (NOT TARGET hdf5::hdf5_hl_fortran)
          add_library(hdf5::hdf5_hl_fortran UNKNOWN IMPORTED)
          foreach(lib ${HDF5_Fortran_HL_LIBRARIES})
            if(lib MATCHES hl_fortran)
              list(REMOVE_ITEM HDF5_Fortran_HL_LIBRARIES ${lib})
              set_target_properties(hdf5::hdf5_hl_fortran PROPERTIES IMPORTED_LOCATION ${lib})
              foreach(d ${CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES})
	        list(REMOVE_ITEM HDF5_Fortran_HL_INCLUDE_DIRS ${d})
              endforeach()
              target_include_directories(hdf5::hdf5_hl_fortran INTERFACE ${HDF5_Fortran_HL_INCLUDE_DIRS})
              target_link_libraries(hdf5::hdf5_hl_fortran INTERFACE ${HDF5_Fortran_HL_LIBRARIES} hdf5::hdf5_fortran)
              break()
            endif()
          endforeach()
        endif()
        set(HDF5_TARGET hdf5::hdf5_hl_fortran)
      endif()
    endif()
  else()
    message(FATAL_ERROR "HDF5 Target not found")
  endif()

  pkg_check_modules(xraylib REQUIRED IMPORTED_TARGET libxrlf03)

  add_library(bmad SHARED ${src})
  set_target_properties(bmad PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  )
  set_source_files_properties(space_charge/fft_interface_mod.f90 multiparticle/envelope_mod.f90 PROPERTIES Fortran_PREPROCESS ON)

  target_link_libraries(bmad PUBLIC bmad_options sim_utils forest ${HDF5_TARGET} PkgConfig::xraylib)
  target_include_directories(bmad PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
    $<INSTALL_INTERFACE:include>
    $<INSTALL_INTERFACE:modules>
  )

  install(TARGETS bmad EXPORT Bmad)
  install(DIRECTORY code DESTINATION src/bmad
    FILES_MATCHING
    PATTERN "*.f90")
  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/" DESTINATION modules FILES_MATCHING PATTERN "*.mod")
endif()
