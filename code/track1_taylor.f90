!+
! Subroutine track1_taylor (start_orb, ele, param, end_orb)
!
! Subroutine to track through an element using the element's taylor map.
! If the taylor map does not exist, one will be created using the old
! reference (ele%taylor%ref) trajectory.
!
! Moudules needed:
!   use bmad
!
! Input:
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- lat_param_struct: Beam parameters.
!     %enegy     -- Energy in GeV
!     %particle  -- Particle type [positron$, or electron$]
!
! Output:
!   end_orb    -- Coord_struct: Ending coords.
!-

subroutine track1_taylor (start_orb, ele, param, end_orb)

use ptc_interface_mod, except_dummy => track1_taylor

implicit none

type (coord_struct) :: start_orb, end_orb, start2_orb
type (coord_struct) :: orb0
type (lat_param_struct) :: param
type (ele_struct) :: ele
real(rp) dtime_ref
character(*), parameter :: r_name = 'track1_taylor'


!

if (.not. associated(ele%taylor(1)%term)) then
  if (global_com%type_out) then
    ! 'WARNING: TAYLOR SERIES NOT PRESENT FOR: ' // ele%name
    ! 'I WILL MAKE A TAYLOR SERIES AROUND THE ZERO ORBIT...'
  endif
  orb0%vec = ele%taylor%ref
  call ele_to_taylor(ele, param)
endif

if (abs(relative_tracking_charge(start_orb, ele, param) - param%default_rel_tracking_charge) > 1e-10) then
  call out_io (s_fatal$, r_name, 'DEFAULT_REL_TRACKING_CHARGE DOES NOT AGREE WITH CHARGE OF TRACKED PARTICLE', &
                                 'FOR TRACKING OF ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

! If the Taylor map does not have the offsets included then do the appropriate
! tracking.

start2_orb = start_orb
end_orb = start_orb

if (ele%taylor_map_includes_offsets) then  ! simple case
  call track_taylor (end_orb%vec, ele%taylor, end_orb%vec)

else
  call offset_particle (ele, param, set$, end_orb, set_multipoles = .false., set_hvkicks = .false.)
  call track_taylor (end_orb%vec, ele%taylor, end_orb%vec)
  call offset_particle (ele, param, unset$, end_orb, set_multipoles = .false., set_hvkicks = .false.)
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
  end_orb%t = start2_orb%t + dtime_ref + start2_orb%vec(5) / (start2_orb%beta * c_light) - end_orb%vec(5) / (end_orb%beta * c_light)
endif

end subroutine
