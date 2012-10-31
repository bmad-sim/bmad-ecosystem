set (EXENAME interface_test)

set (SRC_FILES 
  interface_test
	util
)

set (LINK_LIBS
  bmad
  sim_utils
  recipes_f-90_LEPP 
  xsif 
  forest
)

# This is so CMake will not be confused and know that the main program is in Fortran.
set (LINKER_LANGUAGE_PROP Fortran)
