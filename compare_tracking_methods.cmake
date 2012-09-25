set (EXENAME compare_tracking_methods)
set (SRC_FILES
  compare_tracking_methods/compare_tracking_methods.f90
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
  pgplot 
)