!+
! Subroutine make_mat6_bmad (ele, param, start_orb, end_orb, err)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Input:
!   ele       -- Ele_struct: Element to track through.
!   param     -- lat_param_struct: Parameters are needed for some elements.
!   start_orb -- coord_struct: Starting coords.
!
! Output:
!   ele       -- Ele_struct: Element with transfer matrix.
!     %vec0     -- 0th order map component
!     %mat6     -- 6x6 transfer matrix.
!   end_orb   -- Coord_struct: Coordinates at the end of element.
!   err       -- Logical, optional: Set True if there is an error. False otherwise.
!-

subroutine make_mat6_bmad (ele, param, start_orb, end_orb, err)

use bmad_interface, dummy1 => make_mat6_bmad

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: field_ele
type (coord_struct) :: start_orb, end_orb
type (lat_param_struct)  param

integer key, tm

logical, optional :: err
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
    end_orb = start_orb
    call set_end_orb
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
case (ab_multipole$, crab_cavity$, sad_mult$, beambeam$, sbend$, patch$, quadrupole$, drift$, &
      rcollimator$, ecollimator$, monitor$, instrument$, pipe$, kicker$, hkicker$, vkicker$, &
      elseparator$, rfcavity$, lcavity$, match$, multipole$, octupole$, sextupole$, &
      sol_quad$, solenoid$, taylor$, wiggler$, undulator$, ac_kicker$, gkicker$, foil$)
  tm = ele%tracking_method
  field_ele => pointer_to_field_ele(ele, 1)
  if (ele%tracking_method == linear$) ele%tracking_method = bmad_standard$
  !!if (key /= wiggler$ .or. field_ele%field_calc /= fieldmap$)   ele%tracking_method = bmad_standard$
  call track1 (start_orb, ele, param, end_orb, make_map1 = .true.)
  ele%tracking_method = tm

  ele%vec0 = end_orb%vec - matmul(ele%mat6, start_orb%vec)
  call set_end_orb

!--------------------------------------------------------
! Marker, branch, photon_branch, etc.

case (marker$, detector$, fork$, photon_fork$, floor_shift$, fiducial$, mask$, converter$) 
  end_orb = start_orb
  call set_end_orb

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

subroutine set_end_orb ()

type (coord_struct) end_orb, c00

if (end_orb%direction == 1) then
  end_orb%s = ele%s
else
  end_orb%s = ele%s_start
endif

end subroutine set_end_orb

end subroutine make_mat6_bmad
