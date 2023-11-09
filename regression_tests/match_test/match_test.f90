program match_test

use bmad

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct) orb0, orb1

integer i, ib, ie
character(100) lat_file

!

open (1, file = 'output.now')

lat_file = 'match_test.bmad'
call bmad_parser(lat_file, lat)

!

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  ele => branch%ele(1)
  call init_coord(orb0, lat%particle_start, ele, upstream_end$)
  call make_mat6(ele, branch%param, lat%particle_start, orb0)

  write (1, '(a, i0, a, 6es16.8)') '"Vec0-', ib, '" ABS 1E-10', ele%vec0
  do i = 1, 6
    write (1, '(2(a, i0), a, 6es16.8)') '"Mat', i, '-', ib, '" ABS 1E-10', ele%mat6(i,:)
  enddo
enddo

close(1)

end program
