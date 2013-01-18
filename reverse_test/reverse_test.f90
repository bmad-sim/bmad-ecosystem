program reverse_test

use bmad
use write_lat_file_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct) orb0, orb1, orb2
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
ele => lat%branch(0)%ele(1)
call init_coord (orb0, orb0%vec, ele, .false.)

call track1 (orb0, ele, ele%branch%param, orb1)

orb1%vec(2) = -orb1%vec(2)
orb1%vec(4) = -orb1%vec(4)
orb1%vec(5) = -orb1%vec(5)
orb1%t      = -orb1%t

ele => lat%branch(1)%ele(1)
if (ele%key == rfcavity$) ele%value(rf_frequency$) = -ele%value(rf_frequency$)
call track1 (orb1, ele, ele%branch%param, orb2)

print '(6es10.2, 5x, es10.2)', orb1%vec, orb1%t
write (1, '(a, es11.3)') '"dorb(1)" ABS 1d-14 ', orb2%vec(1) - orb0%vec(1)
write (1, '(a, es11.3)') '"dorb(2)" ABS 1d-14 ', orb2%vec(2) + orb0%vec(2)
write (1, '(a, es11.3)') '"dorb(3)" ABS 1d-14 ', orb2%vec(3) - orb0%vec(3)
write (1, '(a, es11.3)') '"dorb(4)" ABS 1d-14 ', orb2%vec(4) + orb0%vec(4)
write (1, '(a, es11.3)') '"dorb(5)" ABS 1d-14 ', orb2%vec(5) + orb0%vec(5)
write (1, '(a, es11.3)') '"dorb(6)" ABS 1d-14 ', orb2%vec(6) - orb0%vec(6)
write (1, '(a, es11.3)') '"dt"      ABS 1d-14 ', orb2%t      + orb0%t

! And close

close (1)

end program
