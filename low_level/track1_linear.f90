!+
! Subroutine track1_linear (start_orb, ele, param, end_orb)
!
! Particle tracking through a single element assuming linearity.
! That is, just using ele%mat6.
!
! Input:
!   start_orb  -- Coord_struct: Starting position
!   ele        -- Ele_struct: Element
!   param      -- lat_param_struct:
!
! Output:
!   end_orb   -- Coord_struct: End position
!   param     -- lat_param_struct:
!-

subroutine track1_linear (start_orb, ele, param, end_orb)

use bmad_interface, except_dummy => track1_linear

implicit none

type (coord_struct) :: start_orb, start2_orb
type (coord_struct) :: end_orb
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

! Note: ele%mat6 holds the matrix for forward tracking (start_orb%direction == 1) independent
! of whether the element is reversed (ele%orientation = -1) or not.

if (start_orb%time_dir == -1) then
  call mat_inverse(ele%mat6, mat6)
  v0 = -matmul(mat6, ele%vec0)
else
  mat6 = ele%mat6
  v0 = ele%vec0
endif

start2_orb = start_orb
end_orb = start_orb
end_orb%p0c = ele%value(p0c$)

if (start_orb%direction == 1) then
  end_orb%vec = matmul (mat6, start_orb%vec) + v0

else
  mat6 = mat_symp_conj(mat6)
  end_orb%vec(2) = -end_orb%vec(2)
  end_orb%vec(4) = -end_orb%vec(4)
  end_orb%vec = matmul(mat6, end_orb%vec)
  end_orb%vec(2) = -end_orb%vec(2)
  end_orb%vec(4) = -end_orb%vec(4)
  end_orb%vec(5) = start_orb%vec(5) - (end_orb%vec(5) - start_orb%vec(5))

  vec = matmul(mat6, v0)
  vec = [vec(1), -vec(2), vec(3), -vec(4), vec(5), vec(6)]
  end_orb%vec = end_orb%vec - vec
endif

! If delta_ref_time has not been set then just assume that the particle has constant velocity.

dtime_ref = ele%value(delta_ref_time$)
if (dtime_ref == 0) dtime_ref = ele%value(l$) / (end_orb%beta * c_light)
dtime_ref = dtime_ref * start_orb%direction*start_orb%time_dir 

call convert_pc_to (ele%value(p0c$) * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)

if (ele%value(p0c$) == ele%value(p0c_start$)) then
  end_orb%t = start2_orb%t + dtime_ref + (start2_orb%vec(5) - end_orb%vec(5)) / (end_orb%beta * c_light)
else
  end_orb%t = start2_orb%t + dtime_ref + start2_orb%vec(5) / (start2_orb%beta * c_light) - end_orb%vec(5) / (end_orb%beta * c_light)
endif

!

if (end_orb%direction*end_orb%time_dir == 1) then
  end_orb%s = ele%s
else
  end_orb%s = ele%s_start
endif

end subroutine
