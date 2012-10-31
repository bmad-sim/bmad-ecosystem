set (EXENAME interface_test)
set (SRC_FILES 
  program/interface_test.f90
)

set (LINK_LIBS
  bmad
  sim_utils
  recipes_f-90_LEPP 
  xsif 
  forest
)

