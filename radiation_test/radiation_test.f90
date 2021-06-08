program radiation_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct) orb_start, orb_end

!

open (1, file = 'output.now', recl = 200)

call bmad_parser('radiation_test.bmad', lat)
ele => lat%ele(1)

call init_coord (orb_start, lat%particle_start, ele, upstream_end$)

call track1 (orb_start, ele, lat%param, orb_end)
write (1, '(a, 7es16.8)') '"No-Rad"   ABS 1e-16', orb_end%vec

bmad_com%radiation_damping_on = .true.
call track1 (orb_start, ele, lat%param, orb_end)
write (1, '(a, 7es16.8)') '"Rad-Damp" ABS 1e-16', orb_end%vec

bmad_com%radiation_zero_average = .true.
call track1 (orb_start, ele, lat%param, orb_end)
write (1, '(a, 7es16.8)') '"Rad-Zero" ABS 1e-16', orb_end%vec

close(1)

end program
