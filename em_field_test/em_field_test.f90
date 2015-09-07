!+
! Program em_field_test
!
! This program is part of the Bmad regression testing suite.
!-

program em_field_test

use bmad
use nr

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct) :: orb, dorb
type (em_field_struct) field0, fp, fm, ff
type (em_potential_struct) p0, pp, pm

real(rp) del
integer i, j
logical err

!

call bmad_parser ('em_field_test.bmad', lat)

!

ele => lat%ele(1)
orb = lat%beam_start
del = orb%vec(2)

call em_field_calc (ele, lat%param, orb%vec(5), 0.0_rp, orb, .true., field0, .true., err, p0)

ff = em_field_struct()

do i = 1, 3
  j = 2*i - 1
  dorb = orb
  dorb%vec(j) = orb%vec(j) + del
  call em_field_calc (ele, lat%param, dorb%vec(5), 0.0_rp, dorb, .true., fp, .false., err, pp)
  dorb%vec(j) = orb%vec(j) - del
  call em_field_calc (ele, lat%param, dorb%vec(5), 0.0_rp, dorb, .true., fm, .false., err, pm)
  ff%dE(:,i) = (fp%e - fm%e) / (2 * del)
  ff%dB(:,i) = (fp%b - fm%b) / (2 * del)
  select case (i)
  case (1)
    ff%B(2) = ff%B(2) - (pp%A(3) - pm%A(3)) / (2 * del)
    ff%B(3) = ff%B(3) + (pp%A(2) - pm%A(2)) / (2 * del)
  case (2)
    ff%B(3) = ff%B(3) - (pp%A(1) - pm%A(1)) / (2 * del)
    ff%B(1) = ff%B(1) + (pp%A(3) - pm%A(3)) / (2 * del)
  case (3)
    ff%B(1) = ff%B(1) - (pp%A(2) - pm%A(2)) / (2 * del)
    ff%B(2) = ff%B(2) + (pp%A(1) - pm%A(1)) / (2 * del)
  end select

enddo

print '(a, 3f10.3)', 'At:', orb%vec(1:5:2)

print *
print '(2x, a, t49, a, t95, a)', 'dB(theory)', 'dB(actual)', 'dB(theory) - dB(actual)'
do i = 1, 3
  print '(3(3es14.6, 4x))', field0%dB(i,:), ff%dB(i,:), field0%dB(i,:) - ff%dB(i,:)
enddo

print *
print *, ' B(actual)     B(theory)     B(actual) - B(theory)'
do i = 1, 3
  print '(3es14.6)', field0%B(i), ff%B(i), field0%B(i) - ff%B(i)
enddo

print *
print '(a, 3es14.6)', 'B:      ', field0%b
print '(a, 3es14.6)', 'B Grad: ', ff%dB(1,1),  ff%dB(2,2),  ff%dB(3,3)
print '(a, 3es14.6)', 'B Curl: ', ff%dB(2,3) - ff%dB(3,2), ff%dB(3,1) - ff%dB(1,3), ff%dB(1,2) - ff%dB(2,1)
print '(a, 3es14.6)', 'B Div:  ', ff%dB(1,1) + ff%dB(2,2) + ff%dB(3,3)

!

close (1)

end program
