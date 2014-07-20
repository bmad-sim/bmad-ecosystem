program geometry_test

use bmad
use test_mod

implicit none

type (lat_struct), target :: lat
type (floor_position_struct) local, floor
type (ele_struct), pointer :: ele1

integer status

!

call bmad_parser ('geometry_test.bmad', lat)

open (1, file = 'output.now')

floor%r = [1.0_rp, 0.1_rp, 2.1_rp]
local = coords_floor_to_curvilinear (floor, lat%ele(2), ele1, status)
write (1, '(a, 2i6, 3f12.6)') '"G1" ABS 0 ', ele1%ix_ele, status, local%r

close (1)

end program
