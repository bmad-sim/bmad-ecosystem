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
type (em_field_struct) field0, fp, fm

real(rp) del, de(3,3), db(3,3)
integer i, j
logical err

!

call bmad_parser ('em_field_test.bmad', lat)

!

ele => lat%ele(1)
orb = lat%beam_start
del = orb%vec(2)

call em_field_calc (ele, lat%param, orb%vec(5), 0.0_rp, orb, .true., field0, err_flag = err)

do i = 1, 3
  j = 2*i - 1
  dorb = orb
  dorb%vec(j) = orb%vec(j) + del
  call em_field_calc (ele, lat%param, dorb%vec(5), 0.0_rp, dorb, .true., fp, err_flag = err)
  dorb%vec(j) = orb%vec(j) - del
  call em_field_calc (ele, lat%param, dorb%vec(5), 0.0_rp, dorb, .true., fm, err_flag = err)
  de(:,i) = (fp%e - fm%e) / (2 * del)
  db(:,i) = (fp%b - fm%b) / (2 * del)
enddo

print '(a, 3f10.3)', 'At:', orb%vec(1:5:2)
print *
print '(a, 3es14.6)', 'E:      ', field0%e
print '(a, 3es14.6)', 'E Grad: ', de(1,1),  de(2,2),  de(3,3)
print '(a, 3es14.6)', 'E Curl: ', de(2,3) - de(3,2), de(3,1) - de(1,3), de(1,2) - de(2,1)
print '(a, 3es14.6)', 'E Div:  ', de(1,1) + de(2,2) + de(3,3)
print *
do i = 1, 3
  print '(3es14.6)', de(i,:)
enddo


print *
print '(a, 3es14.6)', 'B:      ', field0%b
print '(a, 3es14.6)', 'B Grad: ', db(1,1),  db(2,2),  db(3,3)
print '(a, 3es14.6)', 'B Curl: ', db(2,3) - db(3,2), db(3,1) - db(1,3), db(1,2) - db(2,1)
print '(a, 3es14.6)', 'B Div:  ', db(1,1) + db(2,2) + db(3,3)


!

close (1)

end program
