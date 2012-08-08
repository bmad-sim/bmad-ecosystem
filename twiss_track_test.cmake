set (EXENAME twiss_track_test)
set (SRC_FILES
  twiss_track_test/twiss_track_test.f90
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
