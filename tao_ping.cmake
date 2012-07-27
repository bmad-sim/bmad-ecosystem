set (EXENAME tao_ping)
set (SRC_FILES
  tao_ping/tao_ping_program.f90
  tao_ping/tao_hook_init_read_lattice_info.f90
  tao_ping/tao_ping_struct.f90
  tao_ping/tao_ping_utils.f90
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