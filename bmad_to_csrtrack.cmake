set (EXENAME bmad_to_csrtrack)
set (SRC_FILES
  bmad_to_csrtrack/bmad_to_csrtrack.f90
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
)