module dynamic_aperture_mod

use bmad_struct
use bmad_interface

type aperture_struct
  type (coord_struct)  closed_orbit
  real(rp) x, y
  integer plane
  integer ix_lat
  integer i_turn
end type

type track_input_struct
  integer n_turn
  real(rp) x_init, y_init
  real(rp) accuracy
end type

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine dynamic_aperture (lat, orb0, theta_xy, track_input, aperture)
!
! Subroutine to determine the dynamic aperture of a lattice by tracking.
! The subroutine works by determining where on a radial line y = const * x
! the aperture is. Here x and y are deviations from the closed orbit.
!
! Modules Needed:
!   use dynamic_aperture_mod
!
! Input:
!   lat        -- lat_struct: Lat containing the lattice.
!   orb0        -- Coord_struct: Closed orbit at the start.
!     %vec(6)      -- Energy offset of closed orbit.
!   theta_xy    -- Real(rp): Angle of radial line (in radians) in x-y space.
!                    Angle is "normalized" by %x_init, %y_init.
!   track_input -- Track_input_struct: Structure holding the input data:
!     %n_turn     -- Number of turns tracked.
!     %x_init     -- Initial x coordinate to start with for theta_xy = 0.
!     %y_init     -- Initial y coordinate to start with for theta_xy = pi/2.
!     %accuracy   -- Accuracy needed of aperture results.
!
! Output:
!     aperture  -- Aperture_struct:
!       %closed_orbit -- Closed orbit coordinates
!       %a            -- X at aperture limit
!       %b            -- Y at aperture limit
!       %plane        -- Plane in which lost (X_PLANE$ or Y_PLANE$)
!       %ix_lat      -- Index where lost
!       %i_turn       -- Turn where lost
!
! Note: The radial lines are spaced equally in angle using coordinates
!       normalized by %X_INIT and %Y_INIT
!-

subroutine dynamic_aperture (lat, orb0, theta_xy, track_input, aperture)

implicit none

type (lat_struct)  lat
type (lat_param_struct)  param_save
type (coord_struct)  orb0
type (coord_struct), allocatable :: orbit(:)
type (aperture_struct)  aperture
type (track_input_struct)  track_input

integer it, turn_lost, track_state

real(rp) theta_xy, x0, x1, x2, y0, y1, y2

logical aperture_bracketed

! init setup

if (track_input%x_init == 0 .or. track_input%y_init == 0) then
  print *, 'ERROR IN DYNAMIC_APERTURE: TRACK_INPUT.X_INIT OR',  &
                                           ' TRACK_INPUT.Y_INIT = 0'
  if (global_com%exit_on_error) call err_exit
endif

param_save = lat%param
lat%param%aperture_limit_on = .true.

call reallocate_coord (orbit, lat%n_ele_max)

! Find starting point

x0 = 0
y0 = 0
x1 = cos(theta_xy) * track_input%x_init
y1 = sin(theta_xy) * track_input%y_init

aperture_bracketed = .false.

! use a binary search to find where the aparture is along the line

test_loop: do

  orbit(0) = orb0
  orbit(0)%vec(1) = orbit(0)%vec(1) + x1 
  orbit(0)%vec(3) = orbit(0)%vec(3) + y1 


  ! Make sure eta values have been calculated. 
  call twiss_at_start(lat)

  ! track n_turns

  do it = 1, track_input%n_turn
    call track_all (lat, orbit, 0, track_state)
    if (track_state /= moving_forward$) then
      aperture%plane = orbit(track_state)%state
      exit
    endif
    orbit(0) = orbit(lat%n_ele_track)
  enddo

  ! change search interval end

  if (track_state /= moving_forward$) then
    x2 = x1
    y2 = y1
    turn_lost = it
    aperture_bracketed = .true.
    aperture%ix_lat = track_state 
  else
    x0 = x1
    y0 = y1
  endif

  ! calculate new starting point

  if (aperture_bracketed) then
    x1 = (x0 + x2) / 2
    y1 = (y0 + y2) / 2
  else
    x1 = 2 * x1
    y1 = 2 * y1
  endif

  ! check to see if there is an aperture

  if (x1 > 1000*track_input%x_init .or. y1 > 1000*track_input%y_init) then
    print *, 'ERROR IN DYNAMIC_APERTURE: CANNOT FIND APERTURE LIMIT'
    if (global_com%exit_on_error) call err_exit
  endif

  ! see if we are accurate enough

  if (aperture_bracketed) then
    if (sqrt((x2-x0)**2 + (y2-y0)**2) <= 2*track_input%accuracy) exit test_loop
  endif

enddo test_loop

! fill in the info

aperture%x = x1
aperture%y = y1
aperture%i_turn = turn_lost
aperture%closed_orbit = orb0

lat%param = param_save

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine dynamic_aperture_parallel (lat, orb0, theta_xy_list, track_input, aperture_list)
!
! Parallel version of subroutine dynamic_aperture using OpenMP  
! to process a list of angles theta_xy_list(:) 
! and return a list of apertures aperture_list(:) 
!
!
! See subroutine dynamic_aperture
!
!----------------------------------------------------------------------


subroutine dynamic_aperture_parallel(lat, orb0, angle_list, track_input, aperture_list)

!$ use omp_lib

implicit none

type (lat_struct) :: lat
type (coord_struct) :: orb0
type (lat_struct), allocatable :: omp_lat(:)
type (track_input_struct) :: track_input
type (aperture_struct) :: aperture
type (aperture_struct) :: aperture_list(:)
type (aperture_struct), allocatable :: omp_aperture(:)

real(rp) :: angle_list(:)

integer :: i
integer :: omp_n, omp_i

!

! Parallel memory setup
omp_n = 1
!$ omp_n = omp_get_max_threads()
allocate(omp_lat(omp_n))
allocate(omp_aperture(omp_n))
do i=1, omp_n
  omp_lat(i) = lat
end do

!$OMP parallel &
!$OMP default(private), &
!$OMP shared(omp_lat, orb0, omp_aperture, angle_list, aperture_list, track_input)
!$OMP do schedule(dynamic)
do i=lbound(angle_list, 1), ubound(angle_list, 1)
  omp_i = 1
  !$ omp_i = omp_get_thread_num()+1
  call dynamic_aperture (omp_lat(omp_i), orb0, angle_list(i), track_input, omp_aperture(omp_i))
  aperture_list(i) = omp_aperture(omp_i)
end do  
!$OMP end do
!$OMP end parallel
  
! Cleanup
do i=1, omp_n
  call deallocate_lat_pointers(omp_lat(i))
end do
deallocate(omp_lat, omp_aperture)
  
end subroutine


end module
