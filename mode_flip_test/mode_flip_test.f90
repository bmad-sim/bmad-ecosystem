!+
! Program mode_flip_test
!
! This program is part of the Bmad regression testing suite.
!-

program mode_flip_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct), allocatable :: co(:)

real(rp) mat_track(6,6), mat_twiss(6,6)
real(rp) vec0(6), vec1(6)
integer i
character(40) fmt

!

call bmad_parser ('ah04n-coupled-v0.2.bmad', lat)
call set_on_off(rfcavity$, lat, off$)
allocate(co(0:lat%n_ele_track))

call closed_orbit_calc(lat, co, 4)
call lat_make_mat6(lat, -1, co)
call twiss_at_start(lat)
call twiss_propagate_all(lat)

!

open (1, file = 'output.now')

call this_flip_test (689, 710, 'no-flip to flip')
call this_flip_test (689, 730, 'through flip region')
call this_flip_test (710, 730, 'flip to no-flip')
call this_flip_test (600, 650, 'flip free')

close (1)


!---------------------------
contains

subroutine this_flip_test (n1, n2, who)

integer n1, n2

character(*) who

!

call mat_make_unit (mat_track)
do i = n1+1, n2
  mat_track = matmul(lat%ele(i)%mat6, mat_track)
enddo

vec0 = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
vec1 = [0.07, 0.08, 0.09, 0.10, 0.11, 0.12]
call transfer_mat_from_twiss (lat%ele(n1), lat%ele(n2), vec0, vec0, mat_twiss)

fmt = '(3a, 3f20.15, 5x, 3f20.15)'
do i = 1, 6
  write (1, '(3a, i0, a, 6f14.8)') '"', who, ':', i, '" ABS 1e-8 ', mat_track(i,:)-mat_twiss(i,:)
enddo

!call mat_type (mat_track - mat_twiss, 0, who)

end subroutine this_flip_test

end program
