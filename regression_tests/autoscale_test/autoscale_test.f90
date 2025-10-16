program test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct) orb_start, orb_end

!

call bmad_parser ('autoscale.bmad', lat)
open (1, file = 'output.now')

ele => lat%ele(1)
call init_coord(orb_start, ele, upstream_end$)
call track1(orb_start, ele, lat%param, orb_end)

write (1, '(a, f16.12)') '"Auto_Phase" ABS 2e-4', ele%value(phi0_autoscale$)
write (1, '(a, es18.10)') '"Auto_Amp"   REL 1e-6', ele%value(field_autoscale$)
write (1, '(a, es18.10)') '"z" ABS 4e-10', orb_end%vec(5)
write (1, '(a, es18.10)') '"pz" ABS 1e-8', orb_end%vec(6)

end program
