set (EXENAME dynamic_aperture)
set (SRC_FILES
  dynamic_aperture/dynamic_aperture_test.f90
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