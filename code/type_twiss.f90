!+
! Subroutine type_twiss (ele, frequency_units, compact_format, lines, n_lines)
!
! Subroutine to print or put in a string array Twiss information from an element.
! If the lines argument is not present, the element information is printed to the terminal.
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
!                                A                   B
! Beta (m)              29.89292748          1.39825638
! Alpha                 -2.95314539          0.01539874
! Gamma (1/m)            0.32495843          0.71532874
! Phi (rad)             11.91163456         11.63002398  
! Eta (m)                1.44429482         -0.00066948
! Etap                   0.13477010          0.00337943
!
! Output:
!   lines(7)     -- Character(200), allocatable, optional :: Character array to hold the output. 
!                     If not present, the information is printed to the terminal.
!   n_lines      -- Integer, optional: Number of lines in lines(:) that hold valid output.
!                     n_lines must be present if lines(:) is. 
!-

subroutine type_twiss (ele, frequency_units, compact_format, lines, n_lines)

use bmad_interface, except_dummy => type_twiss

implicit none

type (ele_struct)  ele

integer, optional :: frequency_units
integer, optional :: n_lines
integer i, nl

real(rp) coef, cbar(2,2)

character(*), optional :: lines(:)
character(200) li(7)
character(80) fmt, str, freq_str

logical, optional :: compact_format

! Encode twiss info

select case (integer_option(radians$, frequency_units))
case (0)
  nl = 0
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
 li(1) = 'ERROR IN TYPE_TWISS: BAD "FREQUENCY_UNITS"'
 nl = 1
 return
end select

!

if (logic_option (.false., compact_format)) then
  write (li(1), '(10x, a)')  &
            'Beta     Alpha     Gamma       Phi        Eta       Etap'
  li(2) = trim(str) // '        (m)        (-)'
  write (li(3), fmt) ' X:', ele%a%beta,  &
          ele%a%alpha, ele%a%gamma, coef*ele%a%phi, ele%a%eta, ele%a%etap
  write (li(4), fmt) ' Y:', ele%b%beta,  &
          ele%b%alpha, ele%b%gamma, coef*ele%b%phi, ele%b%eta, ele%b%etap
  nl = 4

else
  call c_to_cbar (ele, cbar)
  write (li(1), '(12x, 2(14x, a), 12x, a, 24x, a)') 'A', 'B', 'Cbar', 'C_mat'
  write (li(2), '(2x, a12, 2a, 2(3x, 2a))') 'Beta (m)    ', v(ele%a%beta), v(ele%b%beta), v2(cbar(1,1)), v2(cbar(1,2)), v2(ele%c_mat(1,1)), v2(ele%c_mat(1,2))
  write (li(3), '(2x, a12, 2a, 2(3x, 2a))') 'Alpha       ', v(ele%a%alpha), v(ele%b%alpha), v2(cbar(2,1)), v2(cbar(2,2)), v2(ele%c_mat(2,1)), v2(ele%c_mat(2,2))
  write (li(4), '(2x, a12, 2a)') 'Gamma (1/m) ', v(ele%a%gamma), v(ele%b%gamma)
  write (li(5), '(2x, a12, 2a, 12x, a, 3(14x, a))') freq_str, v(ele%a%phi*coef), v(ele%b%phi*coef), 'X', 'Y', 'Z'
  write (li(6), '(2x, a12, 5a)') 'Eta (m)     ', v(ele%a%eta),  v(ele%b%eta),  v(ele%x%eta),  v(ele%y%eta),  v(ele%z%eta)
  write (li(7), '(2x, a12, 5a)') 'Etap        ', v(ele%a%etap), v(ele%b%etap), v(ele%x%etap), v(ele%y%etap), v(ele%z%etap)
  nl = 7
endif

! finish

if (present(lines)) then
  n_lines = nl
  lines(1:nl) = li(1:nl)
else
  do i = 1, nl
    print '(1x, a)', trim(li(i))
  enddo
endif

!--------------------------------------------
contains

function v(val) result (str)
real(rp) val
character(15) str

!

if (abs(val) < 9999) then
  write (str, '(f15.8)') val
else
  write (str, '(es15.5)') val
endif

end function v

!--------------------------------------------
!contains

function v2(val) result (str)
real(rp) val
character(13) str

!

if (abs(val) < 9999) then
  write (str, '(f13.8)') val
else
  write (str, '(es13.5)') val
endif

end function v2

end subroutine
