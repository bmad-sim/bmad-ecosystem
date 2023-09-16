program fringe_test

use bmad
use fringe_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct) orbit0, orbit

real(rp) mat6(6,6), mat0(6,6)
integer i

!

open (1, file = 'output.now')

call bmad_parser('fringe_test.bmad', lat)

ele => lat%ele(1)
call init_coord (orbit0, lat%particle_start, ele, upstream_end$)
call mat_make_unit(mat0)

orbit = orbit0
mat6 = mat0
call hwang_bend_edge_kick(ele, lat%param, first_track_edge$, orbit, mat6, .true.)

write (1, '(a, 6es14.6)') '"Up-Orb" ABS 1E-12', orbit%vec - orbit0%vec
do i = 1, 6
  write (1, '(a, i0, a, 6es14.6)') '"Up-Mat', i, '" ABS 1E-12', mat6(i,:) - mat0(i,:) 
enddo
write (1, '(a, es12.4)') '"Up-Symp-Err" ABS 1e-10', mat_symp_error(mat6)

write (1, '(a)')
orbit = orbit0
mat6 = mat0
call hwang_bend_edge_kick(ele, lat%param, second_track_edge$, orbit, mat6, .true.)

write (1, '(a, 6es14.6)') '"Dn-Orb" ABS 1E-12', orbit%vec - orbit0%vec
do i = 1, 6
  write (1, '(a, i0, a, 6es14.6)') '"Dn-Mat', i, '" ABS 1E-12', mat6(i,:) - mat0(i,:) 
enddo
write (1, '(a, es12.4)') '"Dn-Symp-Err" ABS 1e-10', mat_symp_error(mat6)

close (1)

end program fringe_test
