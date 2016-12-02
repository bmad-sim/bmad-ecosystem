program geometry_test

use bmad
use geometry_mod

implicit none

type (lat_struct), target :: lat
type (coord_struct) orbit
type (floor_position_struct) local, floor, pos_ele
type (ele_struct), pointer :: ele
type (em_field_struct) f1, f2

real(rp) w_mat(3,3), s

integer status

!

call bmad_parser ('geometry_test.bmad', lat)

open (1, file = 'output.now')

floor%r = [-0.11_rp, 0.2_rp, 2.0_rp]
local = coords_floor_to_local_curvilinear(floor, lat%branch(1)%ele(1), status)
write (1, '(a, 3f12.8, i6)') '"Floor-to-loc" ABS 0 ', local%r, status

!

floor%r = [1.0_rp, 0.1_rp, 2.1_rp]
local = coords_floor_to_curvilinear (floor, lat%ele(2), ele, status)
write (1, '(a, 2i6, 3f12.6)') '"G1" ABS 0 ', ele%ix_ele, status, local%r

ele => lat%ele(4)

call init_coord (orbit, lat%beam_start, ele, inside$)
s = orbit%vec(5)
call offset_particle (ele, lat%param, unset$, orbit, ds_pos = s)

call em_field_calc (ele, lat%param, s, orbit, .false., f1)

call em_field_calc (ele, lat%param, s, orbit, .true., f2)

pos_ele%r = [orbit%vec(1), orbit%vec(3), s]
local = coords_element_frame_to_local(pos_ele, ele, w_mat)
f2%e = matmul(f2%e, w_mat)
f2%b = matmul(f2%b, w_mat)

write (1, '(a, 3f12.8)') '"Sbend-B1"  ABS 1e-8', f1%b
!write (1, '(a, 3es12.4)') '"Sbend-B2"  ABS 1e-8', f2%b
!write (1, '(a, 3es12.4)') '"Sbend-Orb"  ABS 1e-8', [orbit%vec(1), orbit%vec(3), s]
!write (1, '(a, 3es12.4)') '"Sbend-Loc"  ABS 1e-8', local%r
!write (1, '(a, 3es12.4)') '"Sbend-dB" ABS 1e-8', f1%b - f2%b

close (1)

end program
