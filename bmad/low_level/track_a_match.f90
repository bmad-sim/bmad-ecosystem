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

use equal_mod, except_dummy => track_a_match

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp) xmat(6,6), vec0(6), z0, beta0, dt
real(rp), optional :: mat6(6,6)

logical, optional :: make_matrix
logical, optional :: err_flag
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

call match_ele_to_mat6 (ele, orbit, xmat, vec0, err, include_delta_time = .false.)
if (err) then
  ! Since there are cases where this error may be raised many times, do not print an error message.
  if (present(err_flag)) err_flag = .true.
  orbit%state = lost$
  return
endif

if (is_false(ele%value(phase_trombone_input$)) .and. is_false(ele%value(match_end_input$)) .and. is_false(ele%value(match_end_orbit_input$))) then
  ele%mat6 = xmat
  ele%vec0 = vec0
endif

!

if (logic_option(.false., make_matrix)) then
  if (ele%value(delta_time$) == 0) then
    mat6 = xmat
  else
    call match_ele_to_mat6 (ele, orbit, mat6, vec0, err, include_delta_time = .true.)
  endif
endif

!

if (ele%value(delta_time$) /= 0) then
  dt = orbit%time_dir * ele%value(delta_time$)
  orbit%t = orbit%t + dt
  orbit%vec(5) = orbit%vec(5) - orbit%beta * c_light * dt
endif

beta0 = orbit%beta
z0 = orbit%vec(5)
if (orbit%time_dir == 1) then
  orbit%vec = matmul (xmat, orbit%vec) + vec0
else
  call mat_inverse (xmat, xmat) 
  orbit%vec = matmul (xmat, orbit%vec - vec0)
endif

call convert_pc_to ((1+orbit%vec(6))*orbit%p0c, orbit%species, beta = orbit%beta)
orbit%t = orbit%t + (z0 / beta0 - orbit%vec(5) / orbit%beta) / c_light

end subroutine
