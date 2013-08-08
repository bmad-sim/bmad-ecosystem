!+
! Program ptc_test
!
! This program is part of the Bmad regression testing suite.
!-

program ptc_test

use bmad
use polymorphic_taylor
use ptc_layout_mod

implicit none

type (lat_struct), target :: lat, lat2, lat3
type (ele_struct), pointer :: ele
type (coord_struct) start_orb, end_orb1, end_orb2, end_orb1p, end_orb2p
type (coord_struct) end_orb1t, end_orb2t
type (taylor_struct) bmad_taylor(6)
type (real_8) y8(6)
type (branch_struct), pointer :: branch

real(rp) diff_mat(6,6), diff_vec(6)
real(rp) vec_bmad(6), vec_ptc(6), vec_bmad2(6), beta0 
real(rp) m6_to_ptc(6,6), m6_to_bmad(6,6), m6(6,6)

integer i, j

namelist / params / start_orb

!----------------------------------------------------------
! Check information passing between bmad element and associated ptc fibre

branch => lat%branch(1)
!!call branch_to_ptc_layout (branch)



!------------------------------------------------------------------------
! Tracking tests
! ele(1) is a wiggler, ele(2:4) is the same wiggler split by a marker.

call bmad_parser ('ptc_test.bmad', lat)

open (1, file = 'ptc_test.bmad')
read (1, nml = params)
close (1)

lat%ele(3)%key = -1
call remove_eles_from_lat (lat)

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
write (1, '(a, es20.10)') '"Bmad:orb%t " REL  1E-10', end_orb1%t

write (1, *)
write (1, '(a, 6es10.2)') '"Bmad2:dvec"   ABS  1E-11', end_orb1%vec - end_orb2%vec
write (1, '(a, es10.2)')  '"Bmad2:dmat"   ABS  1E-11', maxval(abs(diff_mat))
write (1, '(a, es10.2)')  '"Bmad2:dvec0"  ABS  1E-11', maxval(abs(diff_vec))
write (1, '(a, es10.2)')  '"Bmad2:dorb%t" ABS  1E-22', &
          (end_orb1%t - lat%ele(1)%ref_time) - (end_orb2%t - lat%ele(3)%ref_time)

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

write (1, *)
write (1, '(a, 6es10.2)') '"PTC1:dvec"   ABS  5E-09', end_orb1%vec - end_orb1p%vec
write (1, '(a, es10.2)')  '"PTC1:dorb%t" ABS  1E-20', end_orb1%t - end_orb1p%t
write (1, '(a, 6es10.2)') '"PTC2:dvec"   ABS  5E-09', end_orb1%vec - end_orb2p%vec
write (1, '(a, es10.2)')  '"PTC2:dmat"   ABS  2E-06', maxval(abs(diff_mat))
write (1, '(a, es10.2)')  '"PTC2:dvec0"  ABS  2E-11', maxval(abs(diff_vec))
write (1, '(a, es10.2)')  '"PTC2:dorb%t" ABS  1E-22', &
          (end_orb1%t - lat%ele(1)%ref_time) - (end_orb2%t - lat%ele(3)%ref_time)

!

lat3 = lat2
lat3%ele%tracking_method = taylor$
lat3%ele%mat6_calc_method = taylor$

call track1 (start_orb, lat3%ele(1), lat3%param, end_orb1t)
call track1 (start_orb, lat3%ele(2), lat3%param, end_orb2t)
call track1 (end_orb2t, lat3%ele(3), lat3%param, end_orb2t)
call lat_make_mat6(lat3, -1)
lat3%ele(3)%vec0 = lat3%ele(3)%vec0 + matmul(lat3%ele(3)%mat6, lat3%ele(2)%vec0)
lat3%ele(3)%mat6 = matmul(lat3%ele(3)%mat6, lat3%ele(2)%mat6)

diff_mat = lat3%ele(3)%mat6 - lat2%ele(1)%mat6
diff_vec = lat3%ele(3)%vec0 - lat2%ele(1)%vec0

write (1, *)
write (1, '(a, 6es10.2)') '"TAYLOR1:dvec"   ABS  1E-9', end_orb1p%vec - end_orb1t%vec
write (1, '(a, es10.2)')  '"TAYLOR1:dorb%t" ABS  1E-17', end_orb1p%t - end_orb1t%t
write (1, '(a, 6es10.2)') '"TAYLOR2:dvec"   ABS  1E-09', end_orb1p%vec - end_orb2t%vec
write (1, '(a, es10.2)')  '"TAYLOR2:dmat"   ABS  1E-11', maxval(abs(diff_mat))
write (1, '(a, es10.2)')  '"TAYLOR2:dvec0"  ABS  1E-11', maxval(abs(diff_vec))
write (1, '(a, es10.2)')  '"TAYLOR2:dorb%t" ABS  1E-22', &
          (end_orb1p%t - lat2%ele(1)%ref_time) - (end_orb2p%t - lat2%ele(3)%ref_time)

! Vector translation

vec_bmad = [0.01, 0.02, 0.03, 0.04, 1.0, 2.0]
beta0 = 0.6
 
call vec_bmad_to_ptc (vec_bmad, beta0, vec_ptc, m6_to_ptc)
call vec_ptc_to_bmad (vec_ptc, beta0, vec_bmad2, m6_to_bmad)

m6 = matmul(m6_to_ptc, m6_to_bmad)
forall (j = 1:6) m6(j,j) = m6(j,j) - 1

write (1, *) 
write (1, '(a, 6es10.2)') '"vec_convert" ABS 1E-15', vec_bmad2 - vec_bmad
write (1, '(a, es10.2)')  '"mat_convert" ABS 1E-15', maxval(abs(m6))

! Map translation

call taylor_to_real_8 (lat3%ele(1)%taylor, beta0, y8)
call real_8_to_taylor (y8, beta0, bmad_taylor)
bmad_taylor = bmad_taylor - lat3%ele(1)%taylor

do i = 1, 6
  diff_vec = maxval(abs(bmad_taylor(i)%term(:)%coef))
enddo

write (1, '(a, 6es10.2)') '"map_convert" ABS 1E-15', diff_vec

end program
