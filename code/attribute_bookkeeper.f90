!+
! Subroutine attribute_bookkeeper (ele, param)
!
! Subroutine to make sure the attributes of an element are self-consistant.
!
!   BEAMBEAM:     bbi_const$ = param%n_part * m_electron * charge$ * r_e /
!                           (2 * pi * param%beam_energy * (sig_x$ + sig_y$)
!
!   RFCAVITY:     rf_wavelength$ = param%total_length / harmon$
!
!   SBEND:        angle$   = L$ * G$
!                 l_chord$ = 2 * sin(Angle$/2) / G$
!                 rho$     = 1 / G$
!
!   WIGGLER:      k1$  = -0.5 * (c_light * b_max$ / param%beam_energy)**2
!                 rho$ = param%beam_energy / (c_light * b_max$)
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

#include "CESR_platform.inc"

subroutine attribute_bookkeeper (ele, param)

  use bmad_struct
  use bmad_interface
  use ptc_interface_mod

  implicit none

  type (ele_struct) ele
  type (param_struct) param

  real(rdef) r, factor, p
  
! b_field_master

  if (ele%b_field_master) then

    if (ele%value(energy$) == 0) then
      factor = 0
    else
      call energy_to_kinetic (ele%value(energy$), param%particle, p0c = p)
      factor = c_light / p
    endif

    select case (ele%key)
    case (quadrupole$)
      ele%value(k1$) = factor * ele%value(B_gradiant$)
    case (sextupole$)
      ele%value(k2$) = factor * ele%value(B_gradiant$)
    case (octupole$)
      ele%value(k3$) = factor * ele%value(B_gradiant$)
    case (solenoid$)
      ele%value(ks$) = factor * ele%value(B_field$)
    case (sbend$)
      ele%value(g$) = factor * ele%value(B_field$)
    case default
      print *, 'ERROR IN ATTRIBUTE_BOOKKEEPER: ', &
                      '"B_FIELD_MASTER" NOT IMPLEMENTED FOR: ', trim(ele%name)
      call err_exit
    end select

  endif

!

  select case (ele%key)

! Bends

  case (sbend$)
    ele%value(angle$) = ele%value(l$) * ele%value(g$)
    if (ele%value(l$) == 0 .or. ele%value(g$) == 0) then
      ele%value(l_chord$) = 0
    else
      ele%value(l_chord$) = 2 * sin(ele%value(angle$)/2) / ele%value(g$)
    endif
    if (ele%value(g$) == 0) then
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

    ele%value(bbi_const$) = &
        -param%n_part * m_electron * ele%value(charge$) * r_e /  &
        (2 * pi * param%beam_energy * (ele%value(sig_x$) + ele%value(sig_y$)))


! Wiggler

  case (wiggler$) 

    if (param%beam_energy == 0) then
      ele%value(k1$) = 0
    else
      ele%value(k1$) = -0.5 * &
                    (c_light * ele%value(b_max$) / param%beam_energy)**2
    endif

    if (ele%value(b_max$) == 0) then
      ele%value(rho$) = 0
    else
      ele%value(rho$) = param%beam_energy / (c_light * ele%value(b_max$))
    endif
                       
  end select

end subroutine
