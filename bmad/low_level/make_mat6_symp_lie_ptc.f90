!+
! Subroutine make_mat6_symp_lie_ptc (ele, param, start_orb, end_orb)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Input:
!   ele        -- Ele_struct: Element with transfer matrix
!   param      -- lat_param_struct: Parameters are needed for some elements.
!   start_orb  -- Coord_struct: Coordinates at the beginning of element. 
!
! Output:
!   ele        -- Ele_struct: Element with transfer matrix.
!     %vec0      -- 0th order map component
!     %mat6      -- 6x6 transfer matrix.
!   end_orb    -- Coord_struct: Coordinates at end of element.
!-

subroutine make_mat6_symp_lie_ptc (ele, param, start_orb, end_orb)

use ptc_interface_mod, except_dummy => make_mat6_symp_lie_ptc

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: start_orb, end_orb, s_orb
type (lat_param_struct)  param
type (taylor_struct) bmad_taylor(6), spin_taylor(0:3)

real(rp) dtime_ref, quat(0:3)

!

s_orb = start_orb

call ele_to_taylor(ele, param, start_orb, .true., orbital_taylor = bmad_taylor, spin_taylor = spin_taylor)
call taylor_to_mat6 (bmad_taylor, start_orb%vec, ele%vec0, ele%mat6)

end_orb = start_orb

if (bmad_com%spin_tracking_on) then
  quat = track_taylor (start_orb%vec, spin_taylor, bmad_taylor%ref)
  end_orb%spin = quat_rotate(quat/norm2(quat), start_orb%spin)
  call kill_taylor (spin_taylor)
endif

end_orb%vec = track_taylor (end_orb%vec, bmad_taylor)
call kill_taylor (bmad_taylor)

! Time change of particle

dtime_ref = ele%value(delta_ref_time$)

if (s_orb%vec(6) == end_orb%vec(6) .and. ele%value(p0c$) == ele%value(p0c_start$)) then
  end_orb%t = s_orb%t + dtime_ref + (s_orb%vec(5) - end_orb%vec(5)) / (end_orb%beta * c_light)
else
  call convert_pc_to (ele%value(p0c$) * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)
  end_orb%t = s_orb%t + dtime_ref + s_orb%vec(5) / (s_orb%beta * c_light) - end_orb%vec(5) / (end_orb%beta * c_light)
endif

! Set s-position

if (end_orb%direction == 1) then
  end_orb%s = ele%s
else
  end_orb%s = ele%s_start
endif

end subroutine

