set (EXENAME bmad_to_mad_and_xsif)
set (SRC_FILES
  bmad_to_mad_and_xsif/bmad_to_mad_and_xsif.f90
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