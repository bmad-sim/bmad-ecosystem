set (EXENAME slice_test)
set (SRC_FILES
  slice_test/slice_test.f90
)

set (INC_DIRS
  ../include
  include
)

set (LINK_LIBS
  tao 
  bmad 
  sim_utils
  recipes_f-90_LEPP 
  xsif 
  pgplot 
  forest 
)
