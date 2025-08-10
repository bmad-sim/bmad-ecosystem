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

bmad_com%absolute_time_tracking = .false.
call reallocate_coord (orb1, lat%n_ele_max)
call init_coord (orb1(0), lat%particle_start, lat%ele(0), downstream_end$)
call track_all (lat, orb1)

write (1, '(a, es22.12)') '"vec(1)" REL  1E-10', orb1(n)%vec(1)
write (1, '(a, es22.12)') '"vec(2)" REL  1E-10', orb1(n)%vec(2)
write (1, '(a, es22.12)') '"vec(3)" REL  1E-10', orb1(n)%vec(3)
write (1, '(a, es22.12)') '"vec(4)" REL  1E-10', orb1(n)%vec(4)
write (1, '(a, es22.12)') '"vec(5)" REL  1E-10', orb1(n)%vec(5)
write (1, '(a, es22.12)') '"vec(6)" REL  1E-10', orb1(n)%vec(6)
write (1, '(a, es22.12)') '"t"      REL  1E-10', c_light*orb1(n)%t*orb1(n)%beta - lat%ele(n)%s

bmad_com%absolute_time_tracking = .true.
bmad_com%absolute_time_ref_shift = .true.
call lattice_bookkeeper (lat)

call track_all (lat, orb1)

write (1, '(a, es22.12)') '"vec(1)" REL  1E-10', orb1(n)%vec(1)
write (1, '(a, es22.12)') '"vec(2)" REL  1E-10', orb1(n)%vec(2)
write (1, '(a, es22.12)') '"vec(3)" REL  1E-10', orb1(n)%vec(3)
write (1, '(a, es22.12)') '"vec(4)" REL  1E-10', orb1(n)%vec(4)
write (1, '(a, es22.12)') '"vec(5)" REL  1E-10', orb1(n)%vec(5)
write (1, '(a, es22.12)') '"vec(6)" REL  1E-10', orb1(n)%vec(6)
write (1, '(a, es22.12)') '"t"      REL  1E-10', c_light*orb1(n)%t*orb1(n)%beta - lat%ele(n)%s

bmad_com%absolute_time_ref_shift = .false.
call lattice_bookkeeper (lat)

call track_all (lat, orb1)

write (1, '(a, es22.12)') '"vec(1)" REL  1E-10', orb1(n)%vec(1)
write (1, '(a, es22.12)') '"vec(2)" REL  1E-10', orb1(n)%vec(2)
write (1, '(a, es22.12)') '"vec(3)" REL  1E-10', orb1(n)%vec(3)
write (1, '(a, es22.12)') '"vec(4)" REL  1E-10', orb1(n)%vec(4)
write (1, '(a, es22.12)') '"vec(5)" REL  1E-10', orb1(n)%vec(5)
write (1, '(a, es22.12)') '"vec(6)" REL  1E-10', orb1(n)%vec(6)
write (1, '(a, es22.12)') '"t"      REL  1E-10', c_light*orb1(n)%t*orb1(n)%beta - lat%ele(n)%s

close (1)


end program 
