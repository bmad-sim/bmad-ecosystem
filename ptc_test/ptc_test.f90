!+
! Program ptc_test
!
! This program is part of the Bmad regression testing suite.
!-

program ptc_test

use bmad

implicit none

type (lat_struct), target :: lat, lat2
type (ele_struct), pointer :: ele
type (coord_struct) start_orb, end_orb1, end_orb2, end_orb1p, end_orb2p

real(rp) diff_mat(6,6), diff_vec(6)

! ele(1) is a wiggler, ele(2:4) is the same wiggler split by a marker.

call bmad_parser ('ptc_test.bmad', lat)
lat%ele(3)%key = -1
call remove_eles_from_lat (lat)

start_orb%vec = [0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006]
call track1 (start_orb, lat%ele(1), lat%param, end_orb1)
call track1 (start_orb, lat%ele(2), lat%param, end_orb2)
call track1 (end_orb2, lat%ele(3), lat%param, end_orb2)
call lat_make_mat6(lat, -1)
lat%ele(3)%vec0 = lat%ele(3)%vec0 + matmul(lat%ele(3)%mat6, lat%ele(2)%vec0)
lat%ele(3)%mat6 = matmul(lat%ele(3)%mat6, lat%ele(2)%mat6)

diff_mat = lat%ele(3)%mat6 - lat%ele(1)%mat6
diff_vec = lat%ele(3)%vec0 - lat%ele(1)%vec0

open (1, file = 'output.now')
write (1, '(a, es20.10)') '"Bmad:vec(1)" REL  1E-10', end_orb1%vec(1)
write (1, '(a, es20.10)') '"Bmad:vec(2)" REL  1E-10', end_orb1%vec(2)
write (1, '(a, es20.10)') '"Bmad:vec(3)" REL  1E-10', end_orb1%vec(3)
write (1, '(a, es20.10)') '"Bmad:vec(4)" REL  1E-10', end_orb1%vec(4)
write (1, '(a, es20.10)') '"Bmad:vec(5)" REL  1E-10', end_orb1%vec(5)
write (1, '(a, es20.10)') '"Bmad:vec(6)" REL  1E-10', end_orb1%vec(6)

write (1, *)
write (1, '(a, 6es10.2)') '"Bmad:dvec"  ABS  1E-17', end_orb1%vec - end_orb2%vec
write (1, '(a, es10.2)')  '"Bmad:dmat"  ABS  2E-15', maxval(abs(diff_mat))
write (1, '(a, es10.2)')  '"Bmad:dvec0" ABS  1E-17', maxval(abs(diff_vec))

!

lat2 = lat
lat2%ele%tracking_method = symp_lie_ptc$
lat2%ele%mat6_calc_method = symp_lie_ptc$

call track1 (start_orb, lat2%ele(1), lat2%param, end_orb1p)
call track1 (start_orb, lat2%ele(2), lat2%param, end_orb2p)
call track1 (end_orb2p, lat2%ele(3), lat2%param, end_orb2p)
call lat_make_mat6(lat2, -1)
lat2%ele(3)%vec0 = lat2%ele(3)%vec0 + matmul(lat2%ele(3)%mat6, lat2%ele(2)%vec0)
lat2%ele(3)%mat6 = matmul(lat2%ele(3)%mat6, lat2%ele(2)%mat6)

diff_mat = lat2%ele(3)%mat6 - lat%ele(1)%mat6
diff_vec = lat2%ele(3)%vec0 - lat%ele(1)%vec0

write (1, '(a, 6es10.2)') '"PTC1:dvec"  ABS  3E-9', end_orb1%vec - end_orb1p%vec
write (1, '(a, 6es10.2)') '"PTC2:dvec"  ABS  3E-9', end_orb1%vec - end_orb2p%vec
write (1, '(a, es10.2)')  '"PTC2:dmat"  ABS  2E-05', maxval(abs(diff_mat))
write (1, '(a, es10.2)')  '"PTC2:dvec0" ABS  2E-11', maxval(abs(diff_vec))




end program
