!+
! Subroutine type2_twiss (ele, lines, n_lines, frequency_units, compact_format)
!
! Subroutine to encode Twiss information in an element in an array of strings.
! See also the subroutine: type_twiss.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele             -- Ele_struct: Element containing the Twiss parameters.
!   frequency_units -- Integer, optional: Units for phi:
!                       = radians$  => Type Twiss, use radians for phi (Default).
!                       = degrees$  => Type Twiss, use degrees for phi.
!                       = cycles$   => Type Twiss, use cycles (1 = 2pi) units.
!   compact_format  -- Logical, optional: If present and True then output looks like:
!
!           Beta     Alpha     Gamma       Phi        Eta       Etap
!            (m)       (-)     (1/m)     (rad)        (m)        (-)
!  X:    29.8929    -2.953     0.325   11.9116     1.4442     0.1347
!  Y:     1.3982     0.015     0.715   11.6300    -0.0006     0.0033     
!
! Else the default is for a format like:
!                                X                   Y
! Beta (m)              29.89292748          1.39825638
! Alpha (-)             -2.95314539          0.01539874
! Gamma (1/m)            0.32495843          0.71532874
! Phi (rad)             11.91163456         11.63002398  
! Eta (m)                1.44429482         -0.00066948
! Etap (-)               0.13477010          0.00337943
!
! Output:
!   lines(:)     -- Character(*): Character array to hold the output.
!   n_lines      -- Number of lines used
!-

#include "CESR_platform.inc"

subroutine type2_twiss (ele, lines, n_lines, frequency_units, compact_format)

use bmad_struct
use bmad_interface, except_dummy => type2_twiss

implicit none

type (ele_struct)  ele

integer, optional :: frequency_units
integer n_lines

real(rp) coef

character(*) lines(:)
character(80) fmt, str, freq_str

logical, optional :: compact_format

! Encode twiss info

select case (integer_option(radians$, frequency_units))
case (0)
  n_lines = 0
  return
case (radians$)
  str = '           (m)       (-)     (1/m)     (rad)'
  fmt = '(a, f11.4, 2f10.3, f10.4, 2f11.4)'
  freq_str = 'Phi (rad)'
  coef = 1
case (degrees$)
  str = '           (m)       (-)     (1/m)     (deg)'
  fmt = '(a, f11.4, 2f10.3, f10.2, 2f11.4)'
  freq_str = 'Phi (deg)'
  coef = 180 / pi
case (cycles$)
  str = '           (m)       (-)     (1/m)  (cycles)'
  fmt = '(a, f11.4, 2f10.3, f10.4, 2f11.4)'
  freq_str = 'Phi (cycles)'
  coef = 1 / twopi                   
case default
 lines(1) = 'ERROR IN TYPE2_TWISS: BAD "FREQUENCY_UNITS"'
 n_lines = 1
 return
end select

!

if (logic_option (.false., compact_format)) then
  write (lines(1), '(10x, a)')  &
            'Beta     Alpha     Gamma       Phi        Eta       Etap'
  lines(2) = trim(str) // '        (m)        (-)'
  write (lines(3), fmt) ' X:', ele%a%beta,  &
          ele%a%alpha, ele%a%gamma, coef*ele%a%phi, ele%a%eta, ele%a%etap
  write (lines(4), fmt) ' Y:', ele%b%beta,  &
          ele%b%alpha, ele%b%gamma, coef*ele%b%phi, ele%b%eta, ele%b%etap
  n_lines = 4

else
  write (lines(1), '(9x, 2(19x, a))') 'A', 'B'
  write (lines(2), '(a11, 2f20.8)') 'Beta (m)      ', ele%a%beta, ele%b%beta
  write (lines(3), '(a11, 2f20.8)') 'Alpha (-)     ', ele%a%alpha, ele%b%alpha
  write (lines(4), '(a11, 2f20.8)') 'Gamma (1/m)   ', ele%a%gamma, ele%b%gamma
  write (lines(5), '(a11, 2f20.8)') freq_str(:14),    ele%a%phi*coef, ele%b%phi*coef
  lines(6) = ''
  write (lines(7), '(9x, 2(19x, a))') 'X', 'Y'  
  write (lines(8), '(a11, 2f20.8)') 'Eta (m)       ', ele%x%eta, ele%y%eta
  write (lines(9), '(a11, 2f20.8)') 'Etap (-)      ', ele%x%etap, ele%y%etap
  n_lines = 9
endif

end subroutine
