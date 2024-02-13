!+
! Subroutine type_twiss (ele, frequency_units, compact_format, lines, n_lines)
!
! Subroutine to print or put in a string array Twiss information from an element.
! If the lines argument is not present, the element information is printed to the terminal.
!
! The compact form looks like:
!           Beta     Alpha     Gamma       Phi        Eta       Etap    dEta/ds
!            (m)       (-)     (1/m)     (rad)        (m)        (-)        (-)
!  X:    29.8929    -2.953     0.325   11.9116     1.4442     0.1347     0.1347
!  Y:     1.3982     0.015     0.715   11.6300    -0.0006     0.0033     0.0033
!
! The default verbose form looks like:
!                          A              B            Cbar                        C_mat
!  Beta (m)         0.95312906     0.01777742  |   0.00010435  -0.00377511      0.00103683  -0.00049140
!  Alpha            0.01355730    -0.00986789  |  -0.01199208   0.00067403     -0.09219227   0.00009904
!  Gamma (1/m)      1.04936870    56.25660601  |   Gamma_c =   1.00002260       Mode_Flip = F
!  Phi (rad)        0.00000000     0.00000000            X              Y              Z
!  Eta (m)         -0.00212550     0.00142557    -0.00212346     0.00144302     0.00000000
!  Etap            -0.03500656    -0.00125557    -0.03513890    -0.00102335     1.00000000
!  dEta/ds         -0.03500656    -0.00125557    -0.03513890    -0.00102335     1.00000000
!  Sigma            0.00034547     0.00005437     0.00034547     0.00005437
!
! Input:
!   ele             -- Ele_struct: Element containing the Twiss parameters.
!   frequency_units -- Integer, optional: Units for phi:
!                       = radians$  => Type Twiss, use radians for phi (Default).
!                       = degrees$  => Type Twiss, use degrees for phi.
!                       = cycles$   => Type Twiss, use cycles (1 = 2pi) units.
!   compact_format  -- Logical, optional: If present and True then use a compact output form.
!
! Output:
!   lines(:)     -- Character(*), optional :: Character array to hold the output.
!                     The string length should be at least 120 characters. 8 lines are needed for the verbose form.
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
character(200) li(9)
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
  li(2) = trim(str) // '        (m)        (-)        (-)'
  write (li(3), fmt) ' X:', ele%a%beta,  &
          ele%a%alpha, ele%a%gamma, coef*ele%a%phi, ele%a%eta, ele%a%etap, ele%a%deta_ds
  write (li(4), fmt) ' Y:', ele%b%beta,  &
          ele%b%alpha, ele%b%gamma, coef*ele%b%phi, ele%b%eta, ele%b%etap, ele%b%deta_ds
  nl = 4

else
  call c_to_cbar (ele, cbar)
  write (li(1), '(12x, 2(14x, a), 12x, a, 24x, a)') 'A', 'B', 'Cbar', 'C_mat'
  write (li(2), '(2x, a12, 2a, 3a, 3x, 2a)')  'Beta (m)    ', v(ele%a%beta), v(ele%b%beta), '  |', v2(cbar(1,1)), v2(cbar(1,2)), v2(ele%c_mat(1,1)), v2(ele%c_mat(1,2))
  write (li(3), '(2x, a12, 2a, 3a, 3x, 2a)')   'Alpha       ', v(ele%a%alpha), v(ele%b%alpha), '  |', v2(cbar(2,1)), v2(cbar(2,2)), v2(ele%c_mat(2,1)), v2(ele%c_mat(2,2))
  write (li(4), '(2x, a12, 2a, a, 3x, 2a, 7x, a, l1)') 'Gamma (1/m) ', v(ele%a%gamma), v(ele%b%gamma), '  |', 'Gamma_c =', v2(ele%gamma_c), 'Mode_Flip = ', ele%mode_flip
  write (li(5), '(2x, a12, 2a, 12x, a, 3(14x, a))') freq_str, v(ele%a%phi*coef), v(ele%b%phi*coef), 'X', 'Y', 'Z'
  write (li(6), '(2x, a12, 5a)')               'Eta (m)     ', v(ele%a%eta),  v(ele%b%eta),  v(ele%x%eta),  v(ele%y%eta),  v(ele%z%eta)
  write (li(7), '(2x, a12, 5a)')               'Etap        ', v(ele%a%etap), v(ele%b%etap), v(ele%x%etap), v(ele%y%etap), v(ele%z%etap)
  write (li(8), '(2x, a12, 5a)')               'dEta/ds     ', v(ele%a%deta_ds), v(ele%b%deta_ds), v(ele%x%deta_ds), v(ele%y%deta_ds), v(ele%z%deta_ds)
  nl = 8
  if (ele%a%sigma /= 0 .or. ele%b%sigma /= 0) then
    nl=nl+1; write (li(nl), '(2x, a12, 5a)')   'Sigma       ', v(ele%a%sigma), v(ele%b%sigma), v(ele%x%sigma), v(ele%y%sigma)
  endif
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
