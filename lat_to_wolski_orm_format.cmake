set (EXENAME lat_to_wolski_orm_format)
set (SRC_FILES
  lat_to_wolski_orm_format/lat_to_wolski_orm_format.f90
)

set (INC_DIRS
  ../include
  include
)

set (LINK_LIBS
  cesr_utils
  bmad 
  sim_utils
  recipes_f-90_LEPP 
  forest
)