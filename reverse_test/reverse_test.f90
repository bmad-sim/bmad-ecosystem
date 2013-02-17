program reverse_test

use bmad
use write_lat_file_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele1, ele2
type (coord_struct) orb0, orb1, orb2

real(rp) mat1(6,6), vec1(6)
logical :: err_flag

! Init

open (1, file = 'output.now')

! Basic idea is that if the particle trajectory is reversed and the particle charge
! is reversed then the particle should come back to its starting point.

call bmad_parser ('reverse.bmad', lat)
call write_bmad_lattice_file ('lat.bmad', lat)
call bmad_parser ('lat.bmad', lat)

lat%branch(1)%ele(0)%value(e_tot$) = lat%branch(0)%ele(1)%value(e_tot$)
lat%branch(1)%ele(0)%value(p0c$) = lat%branch(0)%ele(1)%value(p0c$)

call lat_compute_ref_energy_and_time (lat, err_flag)

orb0 = lat%beam_start
ele1 => lat%branch(0)%ele(1)
call init_coord (orb0, orb0%vec, ele1, .false.)

call track1 (orb0, ele1, ele1%branch%param, orb1)
call make_mat6 (ele1, ele1%branch%param, orb0)

orb1%vec(2) = -orb1%vec(2)
orb1%vec(4) = -orb1%vec(4)
orb1%vec(5) = -orb1%vec(5)
orb1%t      = -orb1%t

ele2 => lat%branch(1)%ele(1)
if (ele2%key == rfcavity$) ele2%value(rf_frequency$) = -ele2%value(rf_frequency$)

call track1 (orb1, ele2, ele2%branch%param, orb2)
call make_mat6 (ele2, ele2%branch%param, orb1)

!

call mat_inverse (ele1%mat6, mat1)
mat1(1:6, 2) = -mat1(1:6, 2)
mat1(1:6, 4) = -mat1(1:6, 4)
mat1(1:6, 5) = -mat1(1:6, 5)
mat1(2, 1:6) = -mat1(2, 1:6) 
mat1(4, 1:6) = -mat1(4, 1:6) 
mat1(5, 1:6) = -mat1(5, 1:6) 
vec1 = ele1%vec0
vec1(2) = -vec1(2)
vec1(4) = -vec1(4)
vec1(5) = -vec1(5)
vec1 = -matmul(mat1, vec1)

!

orb2%vec(2) = -orb2%vec(2)
orb2%vec(4) = -orb2%vec(4)
orb2%vec(5) = -orb2%vec(5)
orb2%t = -orb2%t

print '(6es10.2, 5x, es10.2)', orb1%vec, orb1%t
write (1, '(a, es11.3)') '"dorb(1)" ABS 1d-14 ', orb2%vec(1) - orb0%vec(1)
write (1, '(a, es11.3)') '"dorb(2)" ABS 1d-14 ', orb2%vec(2) - orb0%vec(2)
write (1, '(a, es11.3)') '"dorb(3)" ABS 1d-14 ', orb2%vec(3) - orb0%vec(3)
write (1, '(a, es11.3)') '"dorb(4)" ABS 1d-14 ', orb2%vec(4) - orb0%vec(4)
write (1, '(a, es11.3)') '"dorb(5)" ABS 1d-14 ', orb2%vec(5) - orb0%vec(5)
write (1, '(a, es11.3)') '"dorb(6)" ABS 1d-14 ', orb2%vec(6) - orb0%vec(6)
write (1, '(a, es11.3)') '"c*dt"    ABS 1d-14 ', (orb2%t     - orb0%t) * c_light

write (1, *)
write (1, '(a, 6es11.3)') '"mat1(1,:)" ABS 1d-14 ', ele1%mat6(1,:)
write (1, '(a, 6es11.3)') '"mat1(2,:)" ABS 1d-14 ', ele1%mat6(2,:)
write (1, '(a, 6es11.3)') '"mat1(3,:)" ABS 1d-14 ', ele1%mat6(3,:)
write (1, '(a, 6es11.3)') '"mat1(4,:)" ABS 1d-14 ', ele1%mat6(4,:)
write (1, '(a, 6es11.3)') '"mat1(5,:)" ABS 1d-14 ', ele1%mat6(5,:)
write (1, '(a, 6es11.3)') '"mat1(6,:)" ABS 1d-14 ', ele1%mat6(6,:)

write (1, *)
write (1, '(a, 6es11.3)') '"mat2(1,:)" ABS 1d-14 ', ele2%mat6(1,:)
write (1, '(a, 6es11.3)') '"mat2(2,:)" ABS 1d-14 ', ele2%mat6(2,:)
write (1, '(a, 6es11.3)') '"mat2(3,:)" ABS 1d-14 ', ele2%mat6(3,:)
write (1, '(a, 6es11.3)') '"mat2(4,:)" ABS 1d-14 ', ele2%mat6(4,:)
write (1, '(a, 6es11.3)') '"mat2(5,:)" ABS 1d-14 ', ele2%mat6(5,:)
write (1, '(a, 6es11.3)') '"mat2(6,:)" ABS 1d-14 ', ele2%mat6(6,:)

write (1, *)
write (1, '(a, 6es11.3)') '"dmat(1,:)" ABS 1d-14 ', ele2%mat6(1,:) - mat1(1,:)
write (1, '(a, 6es11.3)') '"dmat(2,:)" ABS 1d-14 ', ele2%mat6(2,:) - mat1(2,:)
write (1, '(a, 6es11.3)') '"dmat(3,:)" ABS 1d-14 ', ele2%mat6(3,:) - mat1(3,:)
write (1, '(a, 6es11.3)') '"dmat(4,:)" ABS 1d-14 ', ele2%mat6(4,:) - mat1(4,:)
write (1, '(a, 6es11.3)') '"dmat(5,:)" ABS 1d-14 ', ele2%mat6(5,:) - mat1(5,:)
write (1, '(a, 6es11.3)') '"dmat(6,:)" ABS 1d-14 ', ele2%mat6(6,:) - mat1(6,:)

write (1, *)
write (1, '(a, 6es11.3)') '"vec1"  ABS 1d-14 ', vec1(:)
write (1, '(a, 6es11.3)') '"vec2"  ABS 1d-14 ', ele2%vec0(:)
write (1, '(a, 7es11.3)') '"dvec0" ABS 1d-14 ', ele2%vec0(:) - vec1(:)

write (1, *)
write (1, '(a,  3es11.3)') '"max(dmat, vec0)" ', maxval(abs([orb2%vec-orb0%vec, (orb2%t - orb0%t) * c_light])), &
                                  maxval(abs(ele2%mat6 - mat1)), maxval(abs(ele2%vec0 - vec1))

! And close

close (1)

end program
