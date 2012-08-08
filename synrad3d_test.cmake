set (EXENAME synrad3d_test)
set (SRC_FILES
  synrad3d_test/synrad3d_test.f90
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
  readline
  curses
  termcap
)
