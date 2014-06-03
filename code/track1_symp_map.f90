!+
! Subroutine track1_symp_map (start_orb, ele, param, end_orb)
!
! Particle tracking through a single element using a partially inverted taylor
! map (In PTC/FPP this is called a genfield). 
!
! Modules Needed:
!   use bmad
!
! Input:
!   start_orb  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- lat_param_struct:
!
! Output:
!   end_orb   -- Coord_struct: End position
!-

subroutine track1_symp_map (start_orb, ele, param, end_orb)

use ptc_interface_mod, except_dummy => track1_symp_map
use tpsalie_analysis, only: assignment(=), operator(*), lnv

implicit none

type (coord_struct) :: start_orb, start2_orb
type (coord_struct) :: end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param

real(dp) re(lnv), beta0, dtime_ref

! Put in offsets if needed.

start2_orb = start_orb
end_orb = start_orb

if (ele%taylor_map_includes_offsets) then  ! simple case
  call track1_this_body

else
  call offset_particle (ele, param, set$, end_orb, .false., .false., .false.)
  call track1_this_body
  call offset_particle (ele, param, unset$, end_orb, .false., .false., .false.)
endif


end_orb%s = ele%s
end_orb%p0c = ele%value(p0c$)

! If delta_ref_time has not been set then just assume that the particle has constant velocity.

dtime_ref = ele%value(delta_ref_time$)
if (dtime_ref == 0) dtime_ref = ele%value(l$) / (end_orb%beta * c_light)

if (ele%value(p0c$) == ele%value(p0c_start$)) then
  end_orb%t = start2_orb%t + dtime_ref + (start2_orb%vec(5) - end_orb%vec(5)) / (end_orb%beta * c_light)
else
  call convert_pc_to (ele%value(p0c$) * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)
  end_orb%t = start2_orb%t + dtime_ref + &
                            start2_orb%vec(5) / (start2_orb%beta * c_light) - end_orb%vec(5) / (end_orb%beta * c_light)
endif

!---------------------------------------------------------------------
contains

subroutine track1_this_body

! Make the genfield map if needed.

if (.not. (associated(ele%ptc_genfield) .and. associated(ele%taylor(1)%term))) then
  if (.not. associated(ele%taylor(1)%term)) call ele_to_taylor(ele, param, end_orb)
  call kill_ptc_genfield (ele%ptc_genfield)  ! clean up if necessary
  allocate (ele%ptc_genfield)
  call taylor_to_genfield (ele%taylor, ele%ptc_genfield, ele%gen0)
endif

! track and add the constant term back in

beta0 = ele%value(p0c_start$) / ele%value(e_tot_start$)
call vec_bmad_to_ptc (end_orb%vec, beta0, re(1:6))

re = ele%ptc_genfield * re

beta0 = ele%value(p0c$) / ele%value(e_tot$)
call vec_ptc_to_bmad (re(1:6), beta0, end_orb%vec)
end_orb%vec = end_orb%vec + ele%gen0

end subroutine

end subroutine
