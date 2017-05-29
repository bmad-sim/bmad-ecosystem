!+
! Subroutine make_mat6_bmad (ele, param, orb_in, orb_out, end_in, err)
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
!   end_in -- Logical, optional: If present and True then the end coords orb_out
!               will be taken as input. Not output as normal.
!
! Output:
!   ele       -- Ele_struct: Element with transfer matrix.
!     %vec0     -- 0th order map component
!     %mat6     -- 6x6 transfer matrix.
!   orb_out   -- Coord_struct: Coordinates at the end of element.
!   err       -- Logical, optional: Set True if there is an error. False otherwise.
!-

subroutine make_mat6_bmad (ele, param, orb_in, orb_out, end_in, err)

use track1_mod, dummy1 => make_mat6_bmad

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: orb_in, orb_out, c00
type (lat_param_struct)  param

integer key, tm

logical, optional :: end_in, err
logical err_flag
character(*), parameter :: r_name = 'make_mat6_bmad'

!--------------------------------------------------------
! init

if (present(err)) err = .false.

call mat_make_unit (ele%mat6)

! If element is off.

key = ele%key

if (.not. ele%is_on) then
  select case (key)
  case (taylor$, match$, fiducial$, floor_shift$)
    if (.not. logic_option (.false., end_in)) call set_orb_out (orb_out, orb_in)
    return
  case (ab_multipole$, multipole$, lcavity$, sbend$, patch$)
    ! Nothing to do here
  case default
    key = drift$  
  end select
endif

if (key == sol_quad$ .and. ele%value(k1$) == 0) key = solenoid$

!---------------------------------------------------------

select case (key)
case (ab_multipole$, sad_mult$, beambeam$, sbend$, patch$, quadrupole$, drift$, &
      rcollimator$, ecollimator$, monitor$, instrument$, pipe$, kicker$, hkicker$, vkicker$, &
      elseparator$, rfcavity$, lcavity$, match$, multipole$, octupole$, sextupole$, &
      sol_quad$, solenoid$, taylor$, wiggler$, undulator$)
  tm = ele%tracking_method
  if (key /= wiggler$ .or. ele%sub_key /= map_type$)   ele%tracking_method = bmad_standard$
  call track1 (orb_in, ele, param, c00, mat6 = ele%mat6, make_matrix = .true.)
  ele%tracking_method = tm

  ele%vec0 = c00%vec - matmul(ele%mat6, orb_in%vec)
  if (.not. logic_option (.false., end_in)) call set_orb_out (orb_out, c00)

!--------------------------------------------------------
! Marker, branch, photon_branch, etc.

case (marker$, detector$, fork$, photon_fork$, floor_shift$, fiducial$, mask$) 
  if (.not. logic_option (.false., end_in)) call set_orb_out (orb_out, orb_in)
  

!--------------------------------------------------------
! rbends are not allowed internally

case (rbend$)

  if (present(err)) err = .true.
  call out_io (s_fatal$, r_name,  'RBEND ELEMENTS NOT ALLOWED INTERNALLY!')
  if (global_com%exit_on_error) call err_exit

!--------------------------------------------------------
! Custom

case (custom$)

  if (present(err)) err = .true.
  call out_io (s_fatal$, r_name,  'MAT6_CALC_METHOD = BMAD_STANDARD IS NOT ALLOWED FOR A CUSTOM ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit

!--------------------------------------------------------
! unrecognized element

case default

  if (present(err)) err = .true.
  call out_io (s_fatal$, r_name,  'UNKNOWN ELEMENT KEY: \i0\ ', &
                                  'FOR ELEMENT: ' // ele%name, i_array = [ele%key])
  if (global_com%exit_on_error) call err_exit

end select

!--------------------------------------------------------
contains

subroutine set_orb_out (orb_out, c00)

type (coord_struct) orb_out, c00

orb_out = c00
if (orb_out%direction == 1) then
  orb_out%s = ele%s
else
  orb_out%s = ele%s_start
endif

end subroutine set_orb_out

end subroutine make_mat6_bmad
