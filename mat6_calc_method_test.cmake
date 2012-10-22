set (EXENAME mat6_calc_method_test)
set (SRC_FILES
  mat6_calc_method_test/mat6_calc_method_test.f90
)

set (INC_DIRS
  ../include
  include
)

set (LINK_LIBS
  bmad 
  sim_utils
  recipes_f-90_LEPP 
  forest
  xsif 
)