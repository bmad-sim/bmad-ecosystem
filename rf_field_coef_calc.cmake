set (EXENAME rf_field_coef_calc)
set (SRC_FILES
  rf_field_coef_calc/rf_field_coef_calc.f90
  rf_field_coef_calc/rf_field_coef_calc_mod.f90
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
  stdc++
)