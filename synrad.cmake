set (EXENAME synrad)
set (SRC_FILES
  synrad/synrad.f90
  synrad/synrad_write_power_mod.f90
)

set (INC_DIRS
  include
)

set (LINK_LIBS
  bsim
  bmadz
  cesr_utils
  bmad
  sim_utils
  recipes_f-90_LEPP
  forest
  pgplot
  xsif
)