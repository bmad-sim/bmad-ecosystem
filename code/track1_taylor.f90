!+
! Subroutine track1_taylor (start_orb, ele, param, end_orb, taylor_in)
!
! Subroutine to track through an element using the element's taylor map.
! If the taylor map does not exist, one will be created using the old
! reference (ele%taylor%ref) trajectory.
!
! Moudules needed:
!   use bmad
!
! Input:
!   start_orb     -- Coord_struct: Starting coords.
!   ele           -- Ele_struct: Element to track through.
!   param         -- lat_param_struct: Beam parameters.
!     %enegy        -- Energy in GeV
!     %particle     -- Particle type [positron$, or electron$]
!   taylor_in(6)  -- taylor_struct, optional: Alternative map to use instead of ele%taylor. 

! Output:
!   end_orb    -- Coord_struct: Ending coords.
!-

subroutine track1_taylor (start_orb, ele, param, end_orb, taylor_in)

use ptc_interface_mod, except_dummy => track1_taylor

implicit none

type (coord_struct) :: start_orb, end_orb, start2_orb
type (coord_struct) :: orb0
type (taylor_struct), optional, target :: taylor_in(6)
type (lat_param_struct) :: param
type (ele_struct), target :: ele
type (taylor_struct), pointer :: taylor(:)

real(rp) dtime_ref
integer dir
character(*), parameter :: r_name = 'track1_taylor'


! Which map to use?

if (present(taylor_in)) then
  taylor => taylor_in
else
  taylor => ele%taylor
endif

! Err checking

if (.not. associated(taylor(1)%term)) then
  ! call out_io (s_warn$, r_name, &
  !       'TAYLOR SERIES NOT PRESENT FOR: ' // ele%name, &
  !       'I WILL MAKE A TAYLOR SERIES AROUND THE ZERO ORBIT...')
  orb0%vec = taylor%ref
  call ele_to_taylor(ele, param, taylor)
endif

!if (abs(rel_tracking_charge_to_mass(start_orb, param) - param%default_rel_tracking_charge) > 1d-10) then
!  call out_io (s_fatal$, r_name, 'DEFAULT_REL_TRACKING_CHARGE DOES NOT AGREE WITH CHARGE OF TRACKED PARTICLE', &
!                                 'FOR TRACKING OF ELEMENT: ' // ele%name)
!  if (global_com%exit_on_error) call err_exit
!endif

! If tracking backwards then need to invert the Taylor map

if (start_orb%direction /= 1) then
  if (ele%key == rfcavity$ .or. ele%key == lcavity$) then
    call out_io (s_fatal$, r_name, 'CANNOT INVERT A TAYLOR MAP FOR BACKWARDS TRACKING IN AN ELEMENT WITH RF FIELDS: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif
  allocate (taylor(6))
  call taylor_inverse (ele%taylor, taylor)
endif

! If the Taylor map does not have the offsets included then do the appropriate
! tracking.

start2_orb = start_orb
end_orb = start_orb

if (.not. ele%taylor_map_includes_offsets) then  ! simple case
  call offset_particle (ele, param, set$, end_orb, set_multipoles = .false., set_hvkicks = .false.)
endif

if (start_orb%direction == 1) then
  call track_taylor (end_orb%vec, taylor, end_orb%vec)
else
  end_orb%vec(2) = -end_orb%vec(2)
  end_orb%vec(4) = -end_orb%vec(4)
  end_orb%vec(5) = -end_orb%vec(5)
  call track_taylor (end_orb%vec, taylor, end_orb%vec)
  end_orb%vec(2) = -end_orb%vec(2)
  end_orb%vec(4) = -end_orb%vec(4)
  end_orb%vec(5) = -end_orb%vec(5)
endif

if (.not. ele%taylor_map_includes_offsets) then  ! simple case
  call offset_particle (ele, param, unset$, end_orb, set_multipoles = .false., set_hvkicks = .false.)
endif

end_orb%s = ele%s
end_orb%p0c = ele%value(p0c$)

! Time change of particle

dtime_ref = ele%value(delta_ref_time$)

if (ele%value(p0c$) == ele%value(p0c_start$)) then
  end_orb%t = start2_orb%t + dtime_ref + (start2_orb%vec(5) - end_orb%vec(5)) / (end_orb%beta * c_light)
else
  call convert_pc_to (ele%value(p0c$) * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)
  end_orb%t = start2_orb%t + dtime_ref + start2_orb%vec(5) / (start2_orb%beta * c_light) - end_orb%vec(5) / (end_orb%beta * c_light)
endif

! Cleanup

if (start_orb%direction /= 1) deallocate(taylor)

end subroutine
