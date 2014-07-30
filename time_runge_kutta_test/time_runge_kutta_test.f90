!+
! Program time_runge_kutta_test
!
! This program is part of the Bmad regression testing suite.
!-

program time_runge_kutta_test

use bmad
use time_tracker_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct) start_orb, end_orb1, end_orb2

namelist / params / start_orb 

!

call bmad_parser ('lat.bmad', lat)

open (1, file = 'lat.bmad')
read (1, nml = params)

ele => lat%ele(1)
call init_coord (start_orb, start_orb%vec, ele, upstream_end$, lat%param%particle)

ele%tracking_method = time_runge_kutta$
call track1 (start_orb, ele, lat%param, end_orb1)

lat%ele%tracking_method = runge_kutta$
call track1 (start_orb, ele, lat%param, end_orb2)

open (1, file = 'output.now')
write (1, '(a, es20.10)') '"time_runge_kutta:diff:vec(1)" ABS  1E-10', end_orb1%vec(1) - end_orb2%vec(1)
write (1, '(a, es20.10)') '"time_runge_kutta:diff:vec(2)" ABS  1E-10', end_orb1%vec(2) - end_orb2%vec(2)
write (1, '(a, es20.10)') '"time_runge_kutta:diff:vec(3)" ABS  1E-10', end_orb1%vec(3) - end_orb2%vec(3)
write (1, '(a, es20.10)') '"time_runge_kutta:diff:vec(4)" ABS  1E-10', end_orb1%vec(4) - end_orb2%vec(4)
write (1, '(a, es20.10)') '"time_runge_kutta:diff:vec(5)" ABS  6E-09', end_orb1%vec(5) - end_orb2%vec(5)
write (1, '(a, es20.10)') '"time_runge_kutta:diff:vec(6)" ABS  4E-09', end_orb1%vec(6) - end_orb2%vec(6)
write (1, '(a, es20.10)') '"time_runge_kutta:diff:s"      ABS  2E-11', end_orb1%s - end_orb2%s
write (1, '(a, es20.10)') '"time_runge_kutta:diff:t"      ABS  1E-18', end_orb1%t - end_orb2%t

write(1, *) 
write (1, '(a, es20.10)') '"runge_kutta:vec(1)" ABS  1E-10', end_orb2%vec(1)
write (1, '(a, es20.10)') '"runge_kutta:vec(2)" ABS  1E-10', end_orb2%vec(2)
write (1, '(a, es20.10)') '"runge_kutta:vec(3)" ABS  1E-10', end_orb2%vec(3)
write (1, '(a, es20.10)') '"runge_kutta:vec(4)" ABS  1E-10', end_orb2%vec(4)
write (1, '(a, es20.10)') '"runge_kutta:vec(5)" ABS  1E-10', end_orb2%vec(5)
write (1, '(a, es20.10)') '"runge_kutta:vec(6)" ABS  1E-10', end_orb2%vec(6)
write (1, '(a, es20.10)') '"runge_kutta:s"      ABS  1E-10', end_orb2%s
write (1, '(a, es20.10)') '"runge_kutta:t"      ABS  1E-10', end_orb2%t



end program
