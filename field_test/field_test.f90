program field_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (em_taylor_struct) emt(3)
type (coord_struct) orb
type (em_field_struct) field

real(rp) em_field(3), err
integer i, j

!

open (1, file = 'output.now')

call bmad_parser ('field_test.bmad', lat)

!

ele => lat%ele(1)
orb%vec = 0
call init_coord (orb, orb, ele, downstream_end$)

call gen_grad1_to_em_taylor(ele, ele%gen_grad_map(1), 0, emt)
err = 0

do i = -1, 2
do j = -1, 2
  orb%vec(1) = 1 * i
  orb%vec(3) = 1 * j
  call em_field_calc (ele, lat%param, 0.01_rp, orb, .true., field)
  call evaluate_em_taylor ([orb%vec(1), orb%vec(3)], emt, em_field)
  err = err + sum(abs(field%b-em_field))
  print '(2i4, 2f8.2, 3(2x, 3f12.6))', i, j, orb%vec(1), orb%vec(3), field%b, em_field, field%b-em_field
enddo
enddo

print *, 'Err:', err

!

write (1, '(a, es20.12)') '"Symp_Err" ABS 1e-10', mat_symp_error(lat%branch(1)%ele(1)%mat6)

close (1)

end program
