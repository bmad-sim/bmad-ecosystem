set (EXENAME synrad3d)
set (SRC_FILES
  synrad3d/synrad3d.f90
  synrad3d/synrad3d_output_mod.f90
  synrad3d/synrad3d_plot_mod.f90
  synrad3d/synrad3d_struct.f90
  synrad3d/synrad3d_test_mod.f90
  synrad3d/synrad3d_track_mod.f90
  synrad3d/synrad3d_utils.f90
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