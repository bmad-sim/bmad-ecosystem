module dynamic_aperture_mod

use bmad_struct
use bmad_interface

implicit none


type aperture_data_struct
  real(rp) x, y     ! aperture
  integer plane     ! plane determining loss
  integer ix_lat    ! ele index lost at
  integer i_turn    ! turn lost at
end type

type aperture_param_struct
  integer :: n_turn = 100                   ! Number of turns a particle must survive
  real(rp) :: x_init = 1e-3_rp              ! initial estimate for horizontal aperture
  real(rp) :: y_init = 1e-3_rp              ! initial estimate for vertical aperture
  real(rp) :: accuracy = 1e-5_rp            ! resolution of bracketed aperture
end type

type aperture_scan_struct
  type(aperture_data_struct), allocatable :: aperture(:) ! set of apertures at different angles
  type(aperture_param_struct) :: param                   ! parameters used for the scan            
  real(rp)           :: min_angle = 0
  real(rp)           :: max_angle = pi
  integer            :: n_angle   = 9
end type

contains


!----------------------------------------------------------------------
!+
! Subroutine dynamic_aperture_scan(lat, orb0, aperture_scan, parallel)
!
! Driver routine for dynamic_aperture. 
! 
! If aperture_scan%param%x_init or %y_init == 0, 
! then a separate scan will be done to set them starting with 0.001 m 
!
! Modules Needed:
!   use dynamic_aperture_mod
!
! Input:
!   lat                 -- lat_struct: Lat containing the lattice.
!   orb0                -- coord_struct: Closed orbit at IP.
!   aperture_scan       -- aperture_scan_struct: 
!     %param            -- aperture_param_struct: input parameters
!     %min_angle, max_angle, n_angle -- integer: angle scan parameters
!   parallel            -- logical, optional :: Use OpenMP parallel routine. 
!                                               Default: False
!
! Output:
!   aperture_scan       -- aperture_scan_struct 
!     %aperture(:)      -- aperture_data_struct: apertures for each angle
!-
subroutine dynamic_aperture_scan(lat, orb0, aperture_scan, parallel)

!$ use omp_lib

implicit none
type (lat_struct) :: lat
type (aperture_scan_struct) :: aperture_scan
type (aperture_data_struct) :: aperture
type (coord_struct) orb0
real(rp), allocatable  :: angle_list(:)
real(rp) :: delta_angle, x_init_temp, y_init_temp
integer :: i, omp_n
logical, optional :: parallel
character(40) :: r_name = 'dynamic_aperture_scan'

!

! Angle preparation
allocate(angle_list(aperture_scan%n_angle))
delta_angle = (aperture_scan%max_angle - aperture_scan%min_angle)/(aperture_scan%n_angle -1)
do i=1, aperture_scan%n_angle
  angle_list(i) = (i-1)*delta_angle + aperture_scan%min_angle
enddo

! Array control
if ( allocated(aperture_scan%aperture)) then
  if (size(aperture_scan%aperture) /= aperture_scan%n_angle) deallocate(aperture_scan%aperture)
endif
if (.not. allocated(aperture_scan%aperture)) allocate(aperture_scan%aperture(aperture_scan%n_angle))


! Auto-set x_init and y_init if they are zero
if (aperture_scan%param%x_init == 0 ) then
  aperture_scan%param%x_init = 0.001_rp
  y_init_temp = aperture_scan%param%y_init
  aperture_scan%param%y_init = 0.001_rp
  call dynamic_aperture (lat, orb0, 0.0_rp, aperture_scan%param, aperture)
  aperture_scan%param%x_init = aperture%x
  aperture_scan%param%y_init = y_init_temp
endif
if (aperture_scan%param%y_init == 0 ) then
  aperture_scan%param%y_init = 0.001_rp
  x_init_temp = aperture_scan%param%x_init
  aperture_scan%param%x_init = 0.001_rp
  call dynamic_aperture (lat, orb0, pi/2, aperture_scan%param, aperture)
  aperture_scan%param%y_init = aperture%y
  aperture_scan%param%x_init = x_init_temp
endif

! Only call the parallel routine if there is more than one thread available 
omp_n = 1
!$ omp_n = omp_get_max_threads()
if (logic_option(.false., parallel) .and. omp_n > 1) then
  call dynamic_aperture_parallel(lat, orb0, angle_list, aperture_scan%param, aperture_scan%aperture)
else
  do i=1, aperture_scan%n_angle
    call dynamic_aperture (lat, orb0, angle_list(i), aperture_scan%param, aperture_scan%aperture(i))
  enddo
endif

end subroutine dynamic_aperture_scan

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine dynamic_aperture (lat, orb0, theta_xy, aperture_param, aperture)
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
!   aperture_param -- aperture_param_struct: Structure holding the input data:
!     %n_turn     -- Number of turns tracked.
!     %x_init     -- Initial x coordinate to start with for theta_xy = 0.
!     %y_init     -- Initial y coordinate to start with for theta_xy = pi/2.
!     %accuracy   -- Accuracy needed of aperture results.
!
! Output:
!     aperture  -- aperture_data_struct:
!       %x            -- X at aperture limit
!       %y            -- Y at aperture limit
!       %plane        -- Plane in which lost (X_PLANE$ or Y_PLANE$)
!       %ix_lat       -- Index where lost
!       %i_turn       -- Turn where lost
!
! Note: The radial lines are spaced equally in angle using coordinates
!       normalized by %X_INIT and %Y_INIT
!-

subroutine dynamic_aperture (lat, orb0, theta_xy, aperture_param, aperture)

implicit none

type (lat_struct)  lat
type (lat_param_struct)  param_save
type (coord_struct)  orb0
type (coord_struct), allocatable :: orbit(:)
type (aperture_data_struct)  aperture
type (aperture_param_struct)  aperture_param

real(rp) theta_xy, x0, x1, x2, y0, y1, y2

integer it, turn_lost, track_state

character(40) :: r_name = 'dynamic_aperture'

logical aperture_bracketed

! init setup

if (aperture_param%x_init == 0) then
  call out_io(s_error$, r_name, 'aperture_param.x_init == 0') 
  if (global_com%exit_on_error) call err_exit
endif


if (aperture_param%y_init == 0) then
  call out_io(s_error$, r_name, 'aperture_param.y_init == 0') 
  if (global_com%exit_on_error) call err_exit
endif

param_save = lat%param
lat%param%aperture_limit_on = .true.

call reallocate_coord (orbit, lat%n_ele_max)

! Find starting point

x0 = 0
y0 = 0
x1 = cos(theta_xy) * aperture_param%x_init
y1 = sin(theta_xy) * aperture_param%y_init

aperture_bracketed = .false.

! use a binary search to find where the aparture is along the line

test_loop: do

  orbit(0) = orb0
  orbit(0)%vec(1) = orbit(0)%vec(1) + x1 
  orbit(0)%vec(3) = orbit(0)%vec(3) + y1 


  ! Make sure eta values have been calculated. 
  call twiss_at_start(lat)

  ! track n_turns

  do it = 1, aperture_param%n_turn
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

  if (x1 > 1000*aperture_param%x_init .or. y1 > 1000*aperture_param%y_init) then
    call out_io(s_error$, r_name, 'CANNOT FIND APERTURE LIMIT')
    if (global_com%exit_on_error) call err_exit
  endif

  ! see if we are accurate enough

  if (aperture_bracketed) then
    if (sqrt((x2-x0)**2 + (y2-y0)**2) <= 2*aperture_param%accuracy) exit test_loop
  endif

enddo test_loop

! fill in the info

aperture%x = x1
aperture%y = y1
aperture%i_turn = turn_lost

lat%param = param_save

end subroutine dynamic_aperture

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine dynamic_aperture_parallel (lat, orb0, theta_xy_list, aperture_param, aperture_list)
!
! Parallel version of subroutine dynamic_aperture using OpenMP  
! to process a list of angles theta_xy_list(:) 
! and return a list of apertures aperture_list(:) 
!
!
! See subroutine dynamic_aperture
!
!----------------------------------------------------------------------


subroutine dynamic_aperture_parallel(lat, orb0, angle_list, aperture_param, aperture_list)

!$ use omp_lib

implicit none

type (lat_struct) :: lat
type (coord_struct) :: orb0
type (lat_struct), allocatable :: omp_lat(:)
type (aperture_param_struct) :: aperture_param
type (aperture_data_struct) :: aperture
type (aperture_data_struct) :: aperture_list(:)
type (aperture_data_struct), allocatable :: omp_aperture(:)

real(rp) :: angle_list(:)

integer :: i
integer :: omp_n, omp_i

character(40) :: r_name = 'dynamic_aperture_parallel'

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
!$OMP shared(omp_lat, orb0, omp_aperture, angle_list, aperture_list, aperture_param)
!$OMP do schedule(dynamic)
do i=lbound(angle_list, 1), ubound(angle_list, 1)
  omp_i = 1
  !$ omp_i = omp_get_thread_num()+1
  call dynamic_aperture (omp_lat(omp_i), orb0, angle_list(i), aperture_param, omp_aperture(omp_i))
  aperture_list(i) = omp_aperture(omp_i)
end do  
!$OMP end do
!$OMP end parallel
  
! Cleanup
do i=1, omp_n
  call deallocate_lat_pointers(omp_lat(i))
end do
deallocate(omp_lat, omp_aperture)
  
end subroutine dynamic_aperture_parallel

end module
