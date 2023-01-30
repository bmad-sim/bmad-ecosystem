!+
! Subroutine track1_taylor (start_orb, ele, param, end_orb, taylor, mat6, make_matrix)
!
! Subroutine to track through an element using the element's taylor map.
! If the taylor map does not exist, one will be created using the old
! reference (ele%taylor%ref) trajectory.
!
! Input:
!   start_orb     -- Coord_struct: Starting coords.
!   ele           -- Ele_struct: Element to track through.
!   param         -- lat_param_struct: Beam parameters.
!   make_matrix   -- logical, optional: Propagate the transfer matrix? Default is false.
!   taylor        -- taylor_struct, optional: Alternative map to use instead of ele%taylor. 
!
! Output:
!   end_orb    -- Coord_struct: Ending coords.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track1_taylor (start_orb, ele, param, end_orb, taylor, mat6, make_matrix)

use ptc_interface_mod, except_dummy => track1_taylor

implicit none

type (coord_struct) :: start_orb, end_orb, start2_orb
type (coord_struct) :: orb0
type (lat_param_struct) :: param
type (ele_struct), target :: ele
type (taylor_struct), target :: taylor2(6)
type (taylor_struct), optional, target :: taylor(6)
type (taylor_struct), pointer :: taylor_ptr(:)

real(rp), optional :: mat6(6,6)
real(rp) dtime_ref, z0

logical, optional :: make_matrix

character(*), parameter :: r_name = 'track1_taylor'

! Some init

start2_orb = start_orb
end_orb = start_orb
end_orb%p0c = ele%value(p0c$)

! Which map to use?

if (present(taylor)) then
  taylor_ptr => taylor
else
  taylor_ptr => ele%taylor
endif

! Err check

if (.not. associated(taylor_ptr(1)%term)) then
  if (present(taylor)) then
    call out_io (s_warn$, r_name, 'BAD TAYLOR MAP ARGUMENT!')
    if (global_com%exit_on_error) call err_exit
    return
  endif
  ! Else create a Taylor map around the zero orbit.
  call ele_to_taylor(ele, param)
endif

! Note: ele%mat6 holds the matrix for forward tracking (start_orb%direction == 1) independent
! of whether the element is reversed (ele%orientation = -1) or not.
! If tracking backwards then need to invert the Taylor map.

if (start_orb%direction == -1 .and. (ele%key == rfcavity$ .or. ele%key == lcavity$ .or. ele%key == crab_cavity$)) then
  call out_io (s_fatal$, r_name, 'CANNOT INVERT A TAYLOR MAP FOR BACKWARDS TRACKING IN AN ELEMENT WITH RF FIELDS: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  return
endif

if (start_orb%direction*start_orb%time_dir == -1) then
  call taylor_inverse (taylor_ptr, taylor2)
  taylor_ptr => taylor2
endif

! If the Taylor map does not have the offsets included then do the appropriate
! tracking.

if (.not. ele%taylor_map_includes_offsets) then  ! simple case
  call offset_particle (ele, set$, end_orb, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)
endif

!

if (logic_option(.false., make_matrix)) call taylor_to_mat6 (taylor_ptr, end_orb%vec, ele%vec0, mat6)

!

if (start_orb%direction == 1) then
  end_orb%vec = track_taylor (end_orb%vec, taylor_ptr)

else
  end_orb%vec(2) = -end_orb%vec(2)
  end_orb%vec(4) = -end_orb%vec(4)
  z0 = end_orb%vec(5)

  end_orb%vec = track_taylor (end_orb%vec, taylor_ptr)

  end_orb%vec(2) = -end_orb%vec(2)
  end_orb%vec(4) = -end_orb%vec(4)
  end_orb%vec(5) = z0 - (end_orb%vec(5) - z0)

  call kill_taylor(taylor2)
endif

!

if (.not. ele%taylor_map_includes_offsets) then  ! simple case
  call offset_particle (ele, unset$, end_orb, set_hvkicks = .false.)
endif

! Time change of particle

dtime_ref = ele%value(delta_ref_time$) * end_orb%direction * end_orb%time_dir

if (start2_orb%vec(6) == end_orb%vec(6) .and. ele%value(p0c$) == ele%value(p0c_start$)) then
  end_orb%t = start2_orb%t + dtime_ref + (start2_orb%vec(5) - end_orb%vec(5)) / (end_orb%beta * c_light)
else
  call convert_pc_to (ele%value(p0c$) * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)
  end_orb%t = start2_orb%t + dtime_ref + start2_orb%vec(5) / (start2_orb%beta * c_light) - end_orb%vec(5) / (end_orb%beta * c_light)
endif

!

if (end_orb%direction*end_orb%time_dir == 1) then
  end_orb%s = ele%s
else
  end_orb%s = ele%s_start
endif

end subroutine
