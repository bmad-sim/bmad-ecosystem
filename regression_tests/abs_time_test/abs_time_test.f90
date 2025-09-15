program abs_time_test

use bmad

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (coord_struct), allocatable :: orb1(:), orb2(:)

integer n
logical err_flag, ok

!

open (1, file = 'output.now')

bmad_com%auto_bookkeeper = .false.
call bmad_parser ('abs_time_test.bmad', lat)
n = lat%n_ele_track

call reallocate_coord (orb1, lat%n_ele_max)
call init_coord (orb1(0), lat%particle_start, lat%ele(0), downstream_end$)

call reallocate_coord (orb2, lat%n_ele_max)
call init_coord (orb2(0), lat%particle_start, lat%ele(0), downstream_end$)

! Multipass tracking

bmad_com%absolute_time_tracking = .false.
call track_all (lat, orb1)

write (1, '(a, es22.12)') '"mrel-vec(1)" REL  1E-10', orb1(n)%vec(1)
write (1, '(a, es22.12)') '"mrel-vec(2)" REL  1E-10', orb1(n)%vec(2)
write (1, '(a, es22.12)') '"mrel-vec(3)" REL  1E-10', orb1(n)%vec(3)
write (1, '(a, es22.12)') '"mrel-vec(4)" REL  1E-10', orb1(n)%vec(4)
write (1, '(a, es22.12)') '"mrel-vec(5)" REL  1E-10', orb1(n)%vec(5)
write (1, '(a, es22.12)') '"mrel-vec(6)" REL  1E-10', orb1(n)%vec(6)
write (1, '(a, es22.12)') '"mrel-t"      REL  1E-10', c_light*orb1(n)%t*orb1(n)%beta - lat%ele(n)%s

bmad_com%absolute_time_tracking = .true.
bmad_com%absolute_time_ref_shift = .true.
call lattice_bookkeeper (lat)
call track_all (lat, orb2)

write (1, '(a, es22.12)') '"mshift-vec(1)" REL  1E-10', orb2(n)%vec(1)
write (1, '(a, es22.12)') '"mshift-vec(2)" REL  1E-10', orb2(n)%vec(2)
write (1, '(a, es22.12)') '"mshift-vec(3)" REL  1E-10', orb2(n)%vec(3)
write (1, '(a, es22.12)') '"mshift-vec(4)" REL  1E-10', orb2(n)%vec(4)
write (1, '(a, es22.12)') '"mshift-vec(5)" REL  1E-10', orb2(n)%vec(5)
write (1, '(a, es22.12)') '"mshift-vec(6)" REL  1E-10', orb2(n)%vec(6)
write (1, '(a, es22.12)') '"mshift-t"      REL  1E-10', c_light*orb2(n)%t*orb2(n)%beta - lat%ele(n)%s

!

bmad_com%absolute_time_ref_shift = .false.
call lattice_bookkeeper (lat)
call track_all (lat, orb1)

write (1, '(a, es22.12)') '"abs-vec(1)" REL  1E-10', orb1(n)%vec(1)
write (1, '(a, es22.12)') '"abs-vec(2)" REL  1E-10', orb1(n)%vec(2)
write (1, '(a, es22.12)') '"abs-vec(3)" REL  1E-10', orb1(n)%vec(3)
write (1, '(a, es22.12)') '"abs-vec(4)" REL  1E-10', orb1(n)%vec(4)
write (1, '(a, es22.12)') '"abs-vec(5)" REL  1E-10', orb1(n)%vec(5)
write (1, '(a, es22.12)') '"abs-vec(6)" REL  1E-10', orb1(n)%vec(6)
write (1, '(a, es22.12)') '"abs-t"      REL  1E-10', c_light*orb1(n)%t*orb1(n)%beta - lat%ele(n)%s

! Non-multipass tracking

bmad_com%absolute_time_tracking = .false.
call lattice_bookkeeper (lat)
call track_all (lat, orb1, ix_branch = 1)

write (1, '(a, es22.12)') '"non-rel-vec(1)" REL  1E-10', orb1(n)%vec(1)
write (1, '(a, es22.12)') '"non-rel-vec(2)" REL  1E-10', orb1(n)%vec(2)
write (1, '(a, es22.12)') '"non-rel-vec(3)" REL  1E-10', orb1(n)%vec(3)
write (1, '(a, es22.12)') '"non-rel-vec(4)" REL  1E-10', orb1(n)%vec(4)
write (1, '(a, es22.12)') '"non-rel-vec(5)" REL  1E-10', orb1(n)%vec(5)
write (1, '(a, es22.12)') '"non-rel-vec(6)" REL  1E-10', orb1(n)%vec(6)
write (1, '(a, es22.12)') '"non-rel-t"      REL  1E-10', c_light*orb1(n)%t*orb1(n)%beta - lat%ele(n)%s

bmad_com%absolute_time_tracking = .true.
bmad_com%absolute_time_ref_shift = .true.
call lattice_bookkeeper (lat)
call track_all (lat, orb2, ix_branch = 1)

write (1, '(a, es22.12)') '"dvec(1)" ABS  1E-9', orb2(n)%vec(1) - orb1(n)%vec(1)
write (1, '(a, es22.12)') '"dvec(2)" ABS  1E-9', orb2(n)%vec(2) - orb1(n)%vec(2)
write (1, '(a, es22.12)') '"dvec(3)" ABS  1E-9', orb2(n)%vec(3) - orb1(n)%vec(3)
write (1, '(a, es22.12)') '"dvec(4)" ABS  1E-9', orb2(n)%vec(4) - orb1(n)%vec(4)
write (1, '(a, es22.12)') '"dvec(5)" ABS  1E-9', orb2(n)%vec(5) - orb1(n)%vec(5)
write (1, '(a, es22.12)') '"dvec(6)" ABS  1E-9', orb2(n)%vec(6) - orb1(n)%vec(6)
write (1, '(a, es22.12)') '"dt"      ABS  1E-9', c_light * (orb1(n)%t*orb1(n)%beta - orb1(n)%t*orb2(n)%beta)

close (1)


end program 
