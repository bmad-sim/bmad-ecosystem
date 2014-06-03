!+
! Subroutine track1_linear (start%orb, ele, param, end_orb)
!
! Particle tracking through a single element assuming linearity.
! That is, just using ele%mat6.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start%orb  -- Coord_struct: Starting position
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
real(rp) dtime_ref

! 

start2_orb = start_orb
end_orb = start_orb
end_orb%vec = matmul (ele%mat6, start_orb%vec) + ele%vec0

end_orb%s = ele%s
end_orb%p0c = ele%value(p0c$)

! If delta_ref_time has not been set then just assume that the particle has constant velocity.

dtime_ref = ele%value(delta_ref_time$)
if (dtime_ref == 0) dtime_ref = ele%value(l$) / (end_orb%beta * c_light)

if (ele%value(p0c$) == ele%value(p0c_start$)) then
  end_orb%t = start2_orb%t + dtime_ref + (start2_orb%vec(5) - end_orb%vec(5)) / (end_orb%beta * c_light)
else
  call convert_pc_to (ele%value(p0c$) * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)
  end_orb%t = start2_orb%t + dtime_ref + start2_orb%vec(5) / (start2_orb%beta * c_light) - end_orb%vec(5) / (end_orb%beta * c_light)
endif

end subroutine
