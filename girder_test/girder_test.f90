!+
! Program girder_test
!
! This program is part of the Bmad regression testing suite.
!-

program girder_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: girder, ele
type (ele_pointer_struct), allocatable :: eles(:)

real(rp) w_mat(3,3), w_mat_inv(3,3), mat3(3,3)
real(rp), pointer :: v(:)
integer n_loc
logical err

!

call bmad_parser ('girder_test.bmad', lat)

!

open (1, file = 'output.now')


call lat_ele_locator ('G1', lat, eles, n_loc, err)
girder => eles(1)%ele

ele => pointer_to_slave(girder, 1)
v => ele%value
write (1, '(a, 3es17.9)') '"Offset1"    ABS 1E-10', v(x_offset_tot$), v(y_offset_tot$), v(z_offset_tot$)
write (1, '(a, 3es17.9)') '"Angle1"     ABS 1E-10', v(x_pitch_tot$), v(y_pitch_tot$), v(tilt_tot$)

ele => pointer_to_slave(girder, 2)
v => ele%value
write (1, '(a, 3es17.9)') '"Offset2"    ABS 1E-10', v(x_offset_tot$), v(y_offset_tot$), v(z_offset_tot$)
write (1, '(a, 3es17.9)') '"Angle2"     ABS 1E-10', v(x_pitch_tot$), v(y_pitch_tot$), v(tilt_tot$)

!

call floor_angles_to_w_mat (0.1_rp, 0.2_rp, 0.3_rp, w_mat, w_mat_inv)
call mat_make_unit (mat3)
mat3 = matmul(w_mat, w_mat_inv) - mat3
write (1, '(a, es10.2)') '"W_mat_inv"  ABS 1E-14', maxval(abs(mat3))

!

close (1)

end program
