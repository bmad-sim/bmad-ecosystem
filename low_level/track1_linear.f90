!+
! Subroutine track1_linear (orbit, ele, param)
!
! Particle tracking through a single element assuming linearity.
! That is, just using ele%mat6.
!
! Input:
!   orbit      -- Coord_struct: Starting position
!   ele        -- Ele_struct: Element
!   param      -- lat_param_struct:
!
! Output:
!   orbit     -- Coord_struct: End position
!   param     -- lat_param_struct:
!-

subroutine track1_linear (orbit, ele, param)

use bmad_interface, except_dummy => track1_linear

implicit none

type (coord_struct) :: start_orb
type (coord_struct) :: orbit
type (ele_struct) :: ele
type (lat_param_struct) :: param

real(rp) dtime_ref, mat6(6,6), vec(6), v0(6)
character(*), parameter :: r_name = 'track1_linear'

! If ele%mat6 is out-of-date must recompute this first.
! Note: A match element does not have a static_linear_map attrib.

if (ele%bookkeeping_state%mat6 /= ok$ .and. (is_false(ele%value(static_linear_map$)) .or. ele%key == match$)) then
  if (ele%mat6_calc_method == tracking$) then
    call out_io(s_error$, r_name, 'MAT6_CALC_METHOD = TRACKING INCOMPATIBLE WITH TRACKING METHOD = LINEAR.', &
                                  'FOR ELEMENT: ' // ele_full_name(ele))
  else
    call make_mat6(ele, param)
  endif
endif

! Note: ele%mat6 holds the matrix for forward tracking (orbit%direction == 1) independent
! of whether the element is reversed (ele%orientation = -1) or not.

if (orbit%time_dir == -1) then
  call mat_inverse(ele%mat6, mat6)
  v0 = -matmul(mat6, ele%vec0)
else
  mat6 = ele%mat6
  v0 = ele%vec0
endif

start_orb = orbit
orbit%p0c = ele%value(p0c$)

if (orbit%direction == 1) then
  orbit%vec = matmul (mat6, orbit%vec) + v0

else
  mat6 = mat_symp_conj(mat6)
  orbit%vec(2) = -orbit%vec(2)
  orbit%vec(4) = -orbit%vec(4)
  orbit%vec = matmul(mat6, orbit%vec)
  orbit%vec(2) = -orbit%vec(2)
  orbit%vec(4) = -orbit%vec(4)
  orbit%vec(5) = start_orb%vec(5) - (orbit%vec(5) - start_orb%vec(5))

  vec = matmul(mat6, v0)
  vec = [vec(1), -vec(2), vec(3), -vec(4), vec(5), vec(6)]
  orbit%vec = orbit%vec - vec
endif

! If delta_ref_time has not been set then just assume that the particle has constant velocity.

dtime_ref = ele%value(delta_ref_time$)
if (dtime_ref == 0) dtime_ref = ele%value(l$) / (orbit%beta * c_light)
dtime_ref = dtime_ref * orbit%direction * orbit%time_dir 

call convert_pc_to (ele%value(p0c$) * (1 + orbit%vec(6)), orbit%species, beta = orbit%beta)

orbit%t = orbit%t + dtime_ref + (start_orb%vec(5) / start_orb%beta - orbit%vec(5) / orbit%beta) / c_light

!

if (orbit%direction*orbit%time_dir == 1) then
  orbit%s = ele%s
else
  orbit%s = ele%s_start
endif

end subroutine
