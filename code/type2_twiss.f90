!+
! Subroutine type2_twiss (ele, lines, n_lines, frequency_units)
!
! Subroutine to encode Twiss information in an element in an array of strings.
! See also the subroutine: type_twiss.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele          -- Ele_struct: Element containing the Twiss parameters.
!   frequency_units 
!                -- Integer, optional: Units for phi:
!                       = radians$  => Type Twiss, use radians for phi (Default).
!                       = degrees$  => Type Twiss, use degrees for phi.
!                       = cycles$   => Type Twiss, use cycles (1 = 2pi) units.
!
! Output
!   lines(:)     -- Character(*): Character array to hold the output.
!   n_lines      -- Number of lines used
!-

#include "CESR_platform.inc"

subroutine type2_twiss (ele, lines, n_lines, frequency_units)

  use bmad_struct
  use bmad_interface, except => type2_twiss

  implicit none

  type (ele_struct)  ele

  integer, optional :: frequency_units
  integer n_lines

  real(rp) coef

  character(*) lines(:)
  character fmt*80, str*80

! Encode twiss info

  select case (integer_option(radians$, frequency_units))
  case (0)
    n_lines = 0
    return
  case (radians$)
    str = '           (m)       (-)     (1/m)     (rad)'
    fmt = '(a, f11.4, 2f10.3, f10.4, 2f11.4)'
    coef = 1
  case (degrees$)
    str = '           (m)       (-)     (1/m)     (deg)'
    fmt = '(a, f11.4, 2f10.3, f10.2, 2f11.4)'
    coef = 180 / pi
  case (cycles$)
    str = '           (m)       (-)     (1/m)  (cycles)'
    fmt = '(a, f11.4, 2f10.3, f10.4, 2f11.4)'
    coef = 1 / twopi                   
  case default
   lines(1) = 'ERROR IN TYPE2_TWISS: BAD "FREQUENCY_UNITS"'
   n_lines = 1
   return
  end select

  write (lines(1), '(10x, a)')  &
            'Beta     Alpha     Gamma       Phi        Eta       Etap'
  lines(2) = trim(str) // '        (m)        (-)'
  write (lines(3), fmt) ' X:', ele%x%beta,  &
          ele%x%alpha, ele%x%gamma, coef*ele%x%phi, ele%x%eta, ele%x%etap
  write (lines(4), fmt) ' Y:', ele%y%beta,  &
          ele%y%alpha, ele%y%gamma, coef*ele%y%phi, ele%y%eta, ele%y%etap
  n_lines = 4

end subroutine
