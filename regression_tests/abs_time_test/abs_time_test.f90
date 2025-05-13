program abs_time_test

use bmad

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (coord_struct), allocatable :: orb1(:), orb2(:)

integer n
logical err_flag, ok

!

bmad_com%auto_bookkeeper = .false.
call bmad_parser ('abs_time_test.bmad', lat)

call reallocate_coord (orb1, lat%n_ele_max)
call init_coord (orb1(0), lat%particle_start, lat%ele(0), downstream_end$)
call track_all (lat, orb1)

bmad_com%absolute_time_tracking = .true.
!!call autoscale_phase_and_amp (lat%ele(2), lat%param, err_flag)
call lattice_bookkeeper (lat)

call reallocate_coord (orb2, lat%n_ele_max)
call init_coord (orb2(0), lat%particle_start, lat%ele(0), downstream_end$)
call track_all (lat, orb2)


open (1, file = 'output.now')

n = lat%n_ele_track
write (1, '(a, es22.12)') '"vec(1)" REL  1E-10', orb1(n)%vec(1)
write (1, '(a, es22.12)') '"vec(2)" REL  1E-10', orb1(n)%vec(2)
write (1, '(a, es22.12)') '"vec(3)" REL  1E-10', orb1(n)%vec(3)
write (1, '(a, es22.12)') '"vec(4)" REL  1E-10', orb1(n)%vec(4)
write (1, '(a, es22.12)') '"vec(5)" REL  1E-10', orb1(n)%vec(5)
write (1, '(a, es22.12)') '"vec(6)" REL  1E-10', orb1(n)%vec(6)
write (1, '(a, es22.12)') '"t"      REL  1E-10', orb1(n)%t

write (1, *)
write (1, '(a, es22.12)') '"dvec(1)" ABS  2E-19', orb2(n)%vec(1) - orb1(n)%vec(1)
write (1, '(a, es22.12)') '"dvec(2)" ABS  1E-19', orb2(n)%vec(2) - orb1(n)%vec(2)
write (1, '(a, es22.12)') '"dvec(3)" ABS  2E-19', orb2(n)%vec(3) - orb1(n)%vec(3)
write (1, '(a, es22.12)') '"dvec(4)" ABS  1E-19', orb2(n)%vec(4) - orb1(n)%vec(4)
write (1, '(a, es22.12)') '"dvec(5)" ABS  2E-15', orb2(n)%vec(5) - orb1(n)%vec(5)
write (1, '(a, es22.12)') '"dvec(6)" ABS  5E-15', orb2(n)%vec(6) - orb1(n)%vec(6)
write (1, '(a, es22.12)') '"c*dt"    ABS  1E-15', c_light * (orb2(n)%t - orb1(n)%t)

!----------------------------

branch => lat%branch(1)
call reallocate_coord (orb1, branch%n_ele_max)

call init_coord (orb1(0), lat%particle_start, branch%ele(0), downstream_end$)
call track_all (lat, orb1, branch%ix_branch)

n = lat%n_ele_track
write (1, *)
write (1, '(a, es22.12)') '"RF-vec(1)" REL  1E-10', orb1(n)%vec(1)
write (1, '(a, es22.12)') '"RF-vec(2)" REL  1E-10', orb1(n)%vec(2)
write (1, '(a, es22.12)') '"RF-vec(3)" REL  1E-10', orb1(n)%vec(3)
write (1, '(a, es22.12)') '"RF-vec(4)" REL  1E-10', orb1(n)%vec(4)
write (1, '(a, es22.12)') '"RF-vec(5)" REL  1E-10', orb1(n)%vec(5)
write (1, '(a, es22.12)') '"RF-vec(6)" REL  1E-10', orb1(n)%vec(6)
write (1, '(a, es22.12)') '"RF-t"      REL  1E-10', orb1(n)%t

close (1)


end program 
