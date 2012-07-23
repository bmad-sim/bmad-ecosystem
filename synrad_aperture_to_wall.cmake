set (EXENAME synrad_aperture_to_wall)
set (SRC_FILES
  synrad_aperture_to_wall/synrad_aperture_to_wall.f90
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
  -Bdynamic /usr/lib64/libX11.so
  readline
)