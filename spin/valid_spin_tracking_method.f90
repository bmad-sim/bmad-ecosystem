!+
! Function valid_spin_tracking_method (ele, spin_tracking_method) result (is_valid)
!
! Routine to return whether a given spin_tracking method is valid for a given element.
!
! Input:
!   ele                   -- ele_struct: Lattice element.
!   spin_tracking_method  -- integer: bmad_standard$, etc.
!
! Output:
!   is_valid              -- logical: True if a valid method. False otherwise.
!-

function valid_spin_tracking_method (ele, spin_tracking_method) result (is_valid)

use bmad_struct

implicit none

type (ele_struct) ele
integer spin_tracking_method
logical is_valid

! 

is_valid = .false.

select case (ele%key)

case (ab_multipole$)
  select case (spin_tracking_method)
  case (tracking$, custom$)
    is_valid = .true.
  end select

case (ac_kicker$)
  select case (spin_tracking_method)
  case (tracking$, custom$)
    is_valid = .true.
  end select

case (capillary$, crystal$, mirror$, multilayer_mirror$, taylor$)
  ! Always False

case (crab_cavity$)
  select case (spin_tracking_method)
  case (tracking$, custom$)
    is_valid = .true.
  end select

case (custom$, rf_bend$)
  select case (spin_tracking_method)
  case (tracking$, custom$)
    is_valid = .true.
  end select

case (hybrid$)
  ! Always False

case (sad_mult$, patch$)
  select case (spin_tracking_method)
  case (custom$, symp_lie_ptc$)
    is_valid = .true.
  end select

case (group$, overlay$, girder$, ramper$, null_ele$)
  ! No valid methods

case (sbend$, quadrupole$, solenoid$, sextupole$, octupole$, drift$, &
      rcollimator$, ecollimator$, monitor$, instrument$, pipe$, kicker$, hkicker$, vkicker$)
  select case (spin_tracking_method)
  case (custom$, symp_lie_ptc$, tracking$, sprint$)
    is_valid = .true.
  end select

case (rfcavity$)
  select case (spin_tracking_method)
  case (custom$, symp_lie_ptc$, tracking$, sprint$)
    is_valid = .true.
  end select

case default
  select case (spin_tracking_method)
  case (custom$, symp_lie_ptc$, tracking$)
    is_valid = .true.
  end select
end select

end function valid_spin_tracking_method 

