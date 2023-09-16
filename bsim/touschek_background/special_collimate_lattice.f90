MODULE special_collimate_lattice_mod

USE bmad
USE superimpose_mod

TYPE collimator_struct
  REAL(rp) location
  REAL(rp) radius
END TYPE collimator_struct

CONTAINS

SUBROUTINE special_collimate_lattice(lat,collimators)

IMPLICIT NONE

TYPE(lat_struct), INTENT(INOUT) :: lat
TYPE(collimator_struct), INTENT(IN) :: collimators(:) 
TYPE(ele_struct) :: collimator

INTEGER i
logical err_flag

!

CALL init_ele(collimator)

collimator%name = 'ECOLLIMATOR01'
collimator%type = 'ecollimator'
collimator%key = ecollimator$
collimator%value(l$) = 0.0_rp

CALL make_mat6(collimator,lat%param)

DO i=1, size(collimators)
  collimator%s = collimators(i)%location
  collimator%value(x1_limit$) = collimators(i)%radius
  collimator%value(x2_limit$) = collimators(i)%radius
  CALL add_superimpose(lat,collimator,0, err_flag)
ENDDO

CALL lat_make_mat6(lat)
CALL twiss_propagate_all(lat)

END SUBROUTINE special_collimate_lattice

END MODULE special_collimate_lattice_mod






