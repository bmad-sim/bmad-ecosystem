!+
! Subroutine make_mat6_taylor (ele, param, orb_in)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element with transfer matrix
!   param  -- lat_param_struct: Parameters are needed for some elements.
!   orb_in -- Coord_struct: Coordinates at the beginning of element. 
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %vec0  -- 0th order map component
!     %mat6  -- 1st order map component (6x6 transfer matrix).
!-

subroutine make_mat6_taylor (ele, param, orb_in)

use ptc_interface_mod, except_dummy => make_mat6_taylor
use make_mat6_mod, except_dummy2 => make_mat6_taylor

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: orb0, orb_in, orb_out
type (lat_param_struct)  param

! If ele%taylor_map_includes_offsets = False then the Taylor map does not have
! any offsets in it and we must put them in explicitly using offset_particle.

if (.not. associated(ele%taylor(1)%term)) call ele_to_taylor(ele, param, ele%taylor, orb_in)

if (ele%taylor_map_includes_offsets) then
  call taylor_to_mat6 (ele%taylor, orb_in%vec, ele%vec0, ele%mat6)

else
  call init_coord (orb0, orb_in, ele, upstream_end$, orb_in%species)
  call offset_particle (ele, param, set$, orb0, set_multipoles = .false., set_hvkicks = .false.)
  call taylor_to_mat6 (ele%taylor, orb0%vec, ele%vec0, ele%mat6, orb_out%vec)
  call init_coord (orb_out, orb_out%vec, ele, downstream_end$, orb_in%species)
  call offset_particle (ele, param, unset$, orb_out, set_multipoles = .false., set_hvkicks = .false.)

  if (ele%value(tilt_tot$) /= 0) call tilt_mat6 (ele%mat6, ele%value(tilt_tot$))

  call mat6_add_pitch (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%orientation, ele%mat6)

  ele%vec0 = orb_out%vec - matmul(ele%mat6, orb_in%vec)

endif

end subroutine

