!+
! Subroutine attribute_bookkeeper (ele, param)
!
! Subroutine to make sure the attributes of an element are self-consistant.
!
!   BEAMBEAM:     bbi_const$ = param%n_part * e_mass * charge$ * r_e /
!                                   (2 * pi * param%energy * (sig_x$ + sig_y$)
!
!   RFCAVITY:     rf_wavelength$ = param%total_length / harmon$
!
!   SBEND:        angle$   = length$ * G_design$
!                 l_chord$ = 2 * sin(angle$/2) / G_design$
!                 rho = 1 / G_design
!
!   WIGGLER:      k1$       = -0.5 * (0.2997 * b_max$ / param%energy)**2
!                 rho$ = 3.3356 * param%energy / b_max$
!
!
! Modules needed:
!   use bmad
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
!Revision 1.6  2002/07/16 20:44:00  dcs
!*** empty log message ***
!
!Revision 1.5  2002/06/13 14:54:21  dcs
!Interfaced with FPP/PTC
!
!Revision 1.4  2002/02/23 20:32:10  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/23 16:52:52  dcs
!*** empty log message ***
!
!Revision 1.1  2002/01/08 21:44:36  dcs
!Aligned with VMS version  -- DCS
!

#include "CESR_platform.inc"

subroutine attribute_bookkeeper (ele, param)

  use bmad

  implicit none

  type (ele_struct) ele
  type (param_struct) param

  real(rdef) r
  
!

  select case (ele%key)

! Bends

  case (sbend$)
    ele%value(angle$) = ele%value(l$) * ele%value(g_design$)
    if (ele%value(l$) == 0 .or. ele%value(g_design$) == 0) then
      ele%value(l_chord$) = 0
    else
      ele%value(l_chord$) = 2 * sin(ele%value(angle$)/2) / ele%value(g_design$)
    endif
    if (ele%value(g_design$) == 0) then
      ele%value(rho$) = 0
    else
      ele%value(rho$) = 1 / ele%value(g$)
    endif


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
      call type_ele(ele, .true., 0, .false., 0, .false.)
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
                       
  end select

end subroutine
