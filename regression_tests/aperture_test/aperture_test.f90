program aperture_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct) start_orb, end_orb

integer :: i, j, ib, nargs, isn

logical print_extra
 
character(100) lat_file

!

global_com%exit_on_error = .false.
lat_file = 'aperture_test.bmad'

print_extra = .false.
nargs = command_argument_count()
if (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit

elseif (nargs > 0)then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
  print_extra = .true.
endif

open (1, file = 'output.now')

!

call bmad_parser (lat_file, lat)

ele => lat%ele(2)
call init_coord (start_orb, lat%particle_start, ele, upstream_end$)
ele%tracking_method = runge_kutta$
call track1(start_orb, ele, lat%param, end_orb)
write (1, '(a, 6f10.6, i4)') '"RK"   ABS 0', end_orb%vec, end_orb%state

ele%tracking_method = time_runge_kutta$
call track1(start_orb, ele, lat%param, end_orb)
write (1, '(a, 6f10.6, i4)') '"TRK"  ABS 0', end_orb%vec, end_orb%state

do i = 3, lat%n_ele_track-1  ! Do not include END marker element.
  call check_this_aperture (lat%ele(i))
enddo

!

close (1)

!----------------------------------------------------------------------
contains

subroutine check_this_aperture (ele)

type (ele_struct) ele
type (coord_struct) orbit
integer i, j, k, state(9)
real(rp) xw2, x0, yw2, y0, eps, unstable(9)

!

call init_coord(orbit, ele = ele, element_end = downstream_end$)

!----------

xw2 = 0.01;   x0 = 0.01
yw2 = 0.015;  y0 = 0
eps = 1e-5

k = 0
do i = -1, 1;  do j = -1, 1
  orbit%vec = [0.0_rp, x0+i*(xw2-eps), 0.0_rp, y0+j*(yw2-eps), -0.6_rp, 0.02_rp]
  orbit%state = alive$
  call check_aperture_limit(orbit, ele, second_track_edge$, lat%param)
  k = k + 1
  state(k) = orbit%state
  unstable(k) = lat%param%unstable_factor
enddo;  enddo
write (1, '(a, 4x, a, 9i3)') quote(trim(ele%name) // '_S_Pxy-'), 'ABS  1E-8', (state(k), k = 1, 9) 
write (1, '(a, 4x, a, 9es10.2)') quote(trim(ele%name) // '_U_Pxy-'), 'ABS  1E-8', (unstable(k), k = 1, 9) 

k = 0
do i = -1, 1;  do j = -1, 1
  orbit%vec = [0.0_rp, x0+i*(xw2+eps), 0.0_rp, y0+j*(yw2+eps), -0.6_rp, 0.02_rp]
  orbit%state = alive$
  call check_aperture_limit(orbit, ele, second_track_edge$, lat%param)
  k = k + 1
  state(k) = orbit%state
  unstable(k) = lat%param%unstable_factor
enddo;  enddo
write (1, '(a, 4x, a, 9i3)') quote(trim(ele%name) // '_S_Pxy+'), 'ABS  1E-8', (state(k), k = 1, 9) 
write (1, '(a, 4x, a, 9es10.2)') quote(trim(ele%name) // '_U_Pxy+'), 'ABS  1E-8', (unstable(k), k = 1, 9) 

!------------

xw2 = 0.01;   x0 = 0.01
yw2 = 0.015;  y0 = 0
eps = 1e-5

k = 0
do i = -1, 1;  do j = -1, 1
  orbit%vec = [x0+i*(xw2-eps), 0.1_rp, y0+j*(yw2-eps), 0.0_rp, -0.6_rp, 0.02_rp]
  orbit%state = alive$
  call check_aperture_limit(orbit, ele, second_track_edge$, lat%param)
  k = k + 1
  state(k) = orbit%state
  unstable(k) = lat%param%unstable_factor
enddo;  enddo
write (1, '(a, 4x, a, 9i3)') quote(trim(ele%name) // '_S_xy-'), 'ABS  1E-8', (state(k), k = 1, 9) 
write (1, '(a, 4x, a, 9es10.2)') quote(trim(ele%name) // '_U_xy-'), 'ABS  1E-8', (unstable(k), k = 1, 9) 

k = 0
do i = -1, 1;  do j = -1, 1
  orbit%vec = [x0+i*(xw2+eps), 0.0_rp, y0+j*(yw2+eps), 0.0_rp, -0.6_rp, 0.02_rp]
  orbit%state = alive$
  call check_aperture_limit(orbit, ele, second_track_edge$, lat%param)
  k = k + 1
  state(k) = orbit%state
  unstable(k) = lat%param%unstable_factor
enddo;  enddo
write (1, '(a, 4x, a, 9i3)') quote(trim(ele%name) // '_S_xy+'), 'ABS  1E-8', (state(k), k = 1, 9) 
write (1, '(a, 4x, a, 9es10.2)') quote(trim(ele%name) // '_U_xy+'), 'ABS  1E-8', (unstable(k), k = 1, 9) 

xw2 = 0.2;   x0 = -0.6
yw2 = 0.03;   y0 = 0.02
eps = 1e-5

!-------------

k = 0
do i = -1, 1;  do j = -1, 1
  orbit%vec = [0.0_rp, 0.01_rp, 0.0_rp, 0.0_rp, x0+i*(xw2-eps), y0+j*(yw2-eps)]
  orbit%state = alive$
  call check_aperture_limit(orbit, ele, second_track_edge$, lat%param)
  k = k + 1
  state(k) = orbit%state
  unstable(k) = lat%param%unstable_factor
enddo;  enddo
write (1, '(a, 4x, a, 9i3)') quote(trim(ele%name) // '_S_zPz-'), 'ABS  1E-8', (state(k), k = 1, 9) 
write (1, '(a, 4x, a, 9es10.2)') quote(trim(ele%name) // '_U_zPz-'), 'ABS  1E-8', (unstable(k), k = 1, 9) 

k = 0
do i = -1, 1;  do j = -1, 1
  orbit%vec = [0.0_rp, 0.01_rp, 0.0_rp, 0.0_rp, x0+i*(xw2+eps), y0+j*(yw2+eps)]
  orbit%state = alive$
  call check_aperture_limit(orbit, ele, second_track_edge$, lat%param)
  k = k + 1
  state(k) = orbit%state
  unstable(k) = lat%param%unstable_factor
enddo;  enddo
write (1, '(a, 4x, a, 9i3)') quote(trim(ele%name) // '_S_zPz+'), 'ABS  1E-8', (state(k), k = 1, 9) 
write (1, '(a, 4x, a, 9es10.2)') quote(trim(ele%name) // '_U_zPz+'), 'ABS  1E-8', (unstable(k), k = 1, 9) 

end subroutine check_this_aperture

end program
