!+
! Program superimpose_test
!
! This program is part of the Bmad regression testing suite.
!-

program superimpose_test

use bmad

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, lord
type (coord_struct), allocatable :: orb(:)

real(rp) :: mat6(6,6), vec0(6)
integer i, j, iarr(10)

character(200) line, names(20)
!

call bmad_parser ('superimpose_test.bmad', lat)
open (1, file = 'output.now')

call transfer_matrix_calc (lat, mat6, vec0)
do i = 1, 6
  write (1, '(a, i0, a, 6f14.8)') '"mat_row', i, '" ABS 1E-7', mat6(i,:)
enddo

call reallocate_coord (orb, lat)
orb(0) = lat%particle_start
call track_all (lat, orb)
write (1, '(a, 6f14.8)') '"track" ABS 1E-7', orb(lat%n_ele_track)%vec

!--------------------------

ele => lat%branch(1)%ele(2)
write (1, '(2a)') '"Reversed-Name" STR ', quote(ele%name)
write (1, '(a, es12.4)') '"Reversed-S" ABS 1E-8 ', ele%s

ele => lat%branch(2)%ele(2)
write (1, '(2a)') '"Reflected-Name" STR ', quote(ele%name)
write (1, '(a, es12.4)') '"Reflected-S" ABS 1E-8 ', ele%s

!--------------------------

branch => lat%branch(3)
line = ''
do i = 1, branch%n_ele_track
  line = trim(line) // ' ' // branch%ele(i)%name
enddo

write (1, '(2a)') '"B3-names" STR ', quote(line)

do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  if (ele%slave_status /= super_slave$) cycle
  do j = 1, ele%n_lord
    lord => pointer_to_lord(ele, j)
    iarr(j) = lord%ix_ele
    names(j) = lord%name
  enddo

  write (1, '(a, i0, a, 99i4)') '"B3-ix-', i, '" ABS 1E-8', ele%n_lord, iarr(ele%n_lord)
  write (1, '(a, i0, a, 99es16.8)') '"B3-s-', i, '" ABS 1E-8', ele%s, ele%ref_time
enddo

!--------------------------

branch => pointer_to_branch('LT', lat)
write (1, '(2a)') '"M-markers" STR ', quote(trim(branch%ele(2)%name) // ' ' // &
                                          trim(branch%ele(3)%name) // ' ' // trim(branch%ele(4)%name))
write (1, '(2a)') '"P-markers" STR ', quote(trim(branch%ele(6)%name) // ' ' // &
                                          trim(branch%ele(7)%name) // ' ' // trim(branch%ele(8)%name))

!--------------------------

branch => lat%branch(5)
do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  write (1, '()') '"Lcav-', i, '" REL 1E-8', ele%s, ele%ref_time, ele%value(p0c$)
enddo

end program
