program aperture_test

use bmad

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct) start_orb, end_orb, orbit

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
branch => lat%branch(1)
do i = 1, branch%n_ele_track
  call check_this_aperture2(branch%ele(i))
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

!----------------------------------------------------------------------
! contains

subroutine check_this_aperture2 (ele)

type (ele_struct) ele
type (coord_struct) orbit

real(rp) :: x_set(6) =                   [0.000_rp, 0.009_rp, 0.011_rp, 0.09_rp, 0.11_rp, 1.0_rp]
character(8) :: ap_at(6) = [character(8):: '???',   '<0.01?', '0.01',   '>0.1?', '0.1',   '>0.1?']   

integer i, j
integer :: particle_at(3) = [first_track_edge$, in_between$, second_track_edge$]
integer :: where(3) = [upstream_end$, inside$, downstream_end$]

! The possible aperture limits are one of: none, 0.01, 0.1

where_loop: do i = 1, 3
  call init_coord(orbit, ele, where(i))
  orbit%state = alive$

  do j = 1, size(x_set)
    orbit%vec(1) = x_set(j)
    call check_aperture_limit(orbit, ele, particle_at(i), ele%branch%param)
    if (orbit%state == alive$) cycle
    write (1, '(a, i0, 1x, a, 1x, i0, a, 2x, a)') '"', ele%ix_ele, trim(ele%name), i, '" STR', quote(ap_at(j))
    cycle where_loop
  enddo

  write (1, '(a, i0, 1x, a, 1x, i0, a)') '"', ele%ix_ele, trim(ele%name), i, '" STR "No Aperture"'
enddo where_loop

end subroutine check_this_aperture2

end program
