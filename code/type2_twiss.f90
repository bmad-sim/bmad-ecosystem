!+
! Subroutine type2_twiss (ele, frequency_units, lines, n_lines)
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
!                -- Integer: Units for phi:
!                       = radians$  => Type Twiss, use radians for phi.
!                       = degrees$  => Type Twiss, use degrees for phi.
!                       = cycles$   => Type Twiss, use cycles (1 = 2pi) units.
!
! Output
!   lines(:)     -- Character*(*): Character array to hold the output.
!   n_lines      -- Number of lines used
!-

#include "CESR_platform.inc"

subroutine type2_twiss (ele, frequency_units, lines, n_lines)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct)  ele

  integer frequency_units, n_lines

  real(rp) coef

  character*(*) lines(:)
  character fmt*80, str*80

! Encode twiss info

  if (frequency_units > 0) then
    if (frequency_units == radians$) then
      str = '           (m)       (-)     (1/m)     (rad)'
      fmt = '(a, f11.4, 2f10.3, f10.4, 2f11.4)'
      coef = 1
    elseif (frequency_units == degrees$) then
      str = '           (m)       (-)     (1/m)     (deg)'
      fmt = '(a, f11.4, 2f10.3, f10.2, 2f11.4)'
      coef = 180 / pi
    elseif (frequency_units == cycles$) then
      str = '           (m)       (-)     (1/m)  (cycles)'
      fmt = '(a, f11.4, 2f10.3, f10.4, 2f11.4)'
      coef = 1 / twopi                   
    else
     lines(1) = 'ERROR IN TYPE2_TWISS: BAD "FREQUENCY_UNITS"'
     n_lines = 1
     return
    endif    

    write (lines(1), *) ' '
    write (lines(2), '(10x, a)')  &
            'Beta     Alpha     Gamma       Phi        Eta       Etap'
    lines(3) = trim(str) // '        (m)        (-)'
    write (lines(4), fmt) ' X:', ele%x%beta,  &
          ele%x%alpha, ele%x%gamma, coef*ele%x%phi, ele%x%eta, ele%x%etap
    write (lines(5), fmt) ' Y:', ele%y%beta,  &
          ele%y%alpha, ele%y%gamma, coef*ele%y%phi, ele%y%eta, ele%y%etap
    n_lines = 5
  endif

end subroutine
