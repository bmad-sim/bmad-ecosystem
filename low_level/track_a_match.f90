!+
! Subroutine track_a_match (orbit, ele, param, err_flag, mat6, make_matrix)
!
! Bmad_standard tracking through a match element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Match element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_match (orbit, ele, param, err_flag, mat6, make_matrix)

use track1_mod, except_dummy => track_a_match

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)

logical, optional :: err_flag, make_matrix
logical err

!

if (is_true(ele%value(match_end_orbit$))) then
  ele%value(x0$)  = orbit%vec(1)
  ele%value(px0$) = orbit%vec(2)
  ele%value(y0$)  = orbit%vec(3)
  ele%value(py0$) = orbit%vec(4)
  ele%value(z0$)  = orbit%vec(5)
  ele%value(pz0$) = orbit%vec(6)
  orbit%vec = [ele%value(x1$), ele%value(px1$), ele%value(y1$), ele%value(py1$), ele%value(z1$), ele%value(pz1$)]
  return
endif

! Until match_end = False, use unit matrix.

call match_ele_to_mat6 (ele, err)
if (err) then
  ! Since there are cases where this error may be raised many 
  ! times, do not print an error message.
  if (present(err_flag)) err_flag = .true.
  orbit%state = lost$
  return
endif

if (present(err_flag)) err_flag = .false.
orbit%vec = matmul (ele%mat6, orbit%vec) + ele%vec0

end subroutine
