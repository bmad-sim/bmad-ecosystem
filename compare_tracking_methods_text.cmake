set (EXENAME compare_tracking_methods_text)
set (SRC_FILES
  compare_tracking_methods_text/compare_tracking_methods_text.f90
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