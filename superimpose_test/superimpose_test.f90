!+
! Program superimpose_test
!
! This program is part of the Bmad regression testing suite.
!-

program superimpose_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct), allocatable :: orb(:)

real(rp) :: mat6(6,6), vec0(6)

integer i

!

call bmad_parser ('superimpose_test.bmad', lat)
open (1, file = 'output.now')

call transfer_matrix_calc (lat, mat6, vec0)
do i = 1, 6
  write (1, '(a, i0, a, 6f14.8)') '"mat_row', i, '" ABS 1E-7', mat6(i,:)
enddo

call reallocate_coord (orb, lat)
orb(0) = lat%beam_start
call track_all (lat, orb)
write (1, '(a, 6f14.8)') '"track" ABS 1E-7', orb(lat%n_ele_track)%vec

end program
