!+
! Subroutine attribute_bookkeeper (ele, param)
!
! Subroutine to make sure the attributes of an element are self-consistant.
! Cases:
!   BEAMBEAM:     bbi_const$ = param%n_part * e_mass * charge$ * r_e / 
!                                   (2 * pi * param%energy * (sig_x$ + sig_y$)
!   RFCAVITY:     rf_wavelength$ = param%total_length / harmon$
!   SBEND, RBEND: angle$ = length$ / rho_design$
!   WIGGLER:      k1$ = -0.5 * (0.2997 * b_max$ / param%energy)**2
!                 rho$ = 3.3356 * param%energy / b_max$
!
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ele   -- Ele_struct: Element with attributes 
!   param -- Param_struct: 
!
! Output:
!   ele  -- Ele_struct: Element with self-consistant attributes.
!-

!$Id$
!$Log$
!Revision 1.1  2002/01/08 21:44:36  dcs
!Aligned with VMS version  -- DCS
!

#include "CESR_platform.inc"

subroutine attribute_bookkeeper (ele, param)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct) ele
  type (param_struct) param

  real r

!

  select case (ele%key)

! Bends

  case (sbend$, rbend$)
    ele%value(angle$) = ele%value(l$) / ele%value(rho_design$)

! RFcavity

  case (rfcavity$)
    if (ele%value(harmon$) /= 0) ele%value(rf_wavelength$) =  &
                                   param%total_length / ele%value(harmon$)

! BeamBeam

  case (beambeam$)

    if (ele%value(n_slice$) == 0) ele%value(n_slice$) = 1.0 ! revert to default

    if (ele%value(charge$) == 0 .or. param%n_part == 0) then
      ele%value(bbi_const$) = 0
      return
    endif

    if (ele%value(sig_x$) == 0 .or. ele%value(sig_y$) == 0) then
      type *, 'ERROR IN ATTRIBUTE_BOOKKEEPER: ZERO SIGMA IN BEAMBEAM ELEMENT!'
      call type_ele(ele, .true., 0, .false., .false.)
      call exit
    endif

    ele%value(bbi_const$) = -param%n_part * e_mass * ele%value(charge$) * &
      r_e /  (2 * pi * param%energy * (ele%value(sig_x$) + ele%value(sig_y$)))


! Wiggler

  case (wiggler$) 

    if (param%energy == 0) then
      ele%value(k1$) = 0
    else
      ele%value(k1$) = -0.5 * (0.2997 * ele%value(b_max$) / param%energy)**2
    endif

    if (ele%value(b_max$) == 0) then
      ele%value(rho$) = 0
    else
      ele%value(rho$) = 3.3356 * param%energy / ele%value(b_max$)
    endif
                       
! Loop

  case (loop$)

    r = ele%value(radius$)
    ele%value(diameter$) = 2 * r
    ele%value(r2$) = r**2
    ele%value(ri$) = r * ele%value(current$)
    ele%value(r2i$) = ele%value(r2$) * ele%value(current$)

!

  end select

end subroutine
