module custom_dynamic_aperture_mod

use bmad

implicit none

type custom_aperture_data_struct
  real(rp) x, y     ! aperture
  integer plane     ! plane determining loss
  integer ix_lat    ! ele index lost at
  integer i_turn    ! turn lost at
end type

type custom_aperture_param_struct
  type (coord_struct)  :: closed_orbit      ! about which the aperture is scanned
  integer :: n_turn = 100                   ! Number of turns a particle must survive
  integer :: n_adts = -1                    ! If positive, will accumulate and average the phase advance per turn over n_adts turns at the 
                                            ! aperture along min_angle, max_angle, and +y.  This is a calculation of the amplitude-dependent tune
                                            ! shift at the dynamic aperture.
  real(rp) :: init_len = 1e-3_rp              ! initial estimate for horizontal aperture
  real(rp) :: step_len = 5e-4_rp
  real(rp) :: accuracy = 1e-5_rp            ! resolution of bracketed aperture
  real(rp) :: adts_x_max
  real(rp) :: adts_x_min
  real(rp) :: adts_y_max
  real(rp) :: adts_y_min
end type

type custom_aperture_scan_struct
  type(custom_aperture_data_struct), allocatable :: aperture(:) ! set of apertures at different angles
  type(custom_aperture_param_struct) :: param                   ! parameters used for the scan            
  real(rp)           :: min_angle = 0
  real(rp)           :: max_angle = pi
  integer            :: n_angle   = 9
  real(rp)           :: Sx = 1.0
  real(rp)           :: Sy = 1.0
end type

contains

!----------------------------------------------------------------------
!+
! Subroutine dynamic_aperture_scan(lat, aperture_scan, parallel)
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
subroutine custom_dynamic_aperture_scan(lat, aperture_scan, parallel)

implicit none
type (lat_struct) :: lat
type (custom_aperture_scan_struct) :: aperture_scan
type (custom_aperture_data_struct) :: aperture
real(rp), allocatable  :: angle_list(:)
real(rp) :: delta_angle, x_init_temp, y_init_temp
integer :: i
logical, optional :: parallel
character(40) :: r_name = 'custom_dynamic_aperture_scan'
logical do_adts
logical mask_x(1:lat%n_ele_track)
logical mask_y(1:lat%n_ele_track)
integer n_angle

!
n_angle = aperture_scan%n_angle

! Angle preparation
allocate(angle_list(n_angle))
delta_angle = (aperture_scan%max_angle - aperture_scan%min_angle)/(aperture_scan%n_angle -1)
mask_x = abs(lat%ele(1:lat%n_ele_track)%value(x1_limit$)) > 0.0001
mask_y = abs(lat%ele(1:lat%n_ele_track)%value(y1_limit$)) > 0.0001
aperture_scan%Sx = minval( lat%ele(1:lat%n_ele_track)%value(x1_limit$) / sqrt(lat%ele(1:lat%n_ele_track)%a%beta), mask_x ) * sqrt(lat%ele(1)%a%beta)
aperture_scan%Sy = minval( lat%ele(1:lat%n_ele_track)%value(y1_limit$) / sqrt(lat%ele(1:lat%n_ele_track)%b%beta), mask_y ) * sqrt(lat%ele(1)%b%beta)
do i=1, n_angle
  angle_list(i) = (i-1)*delta_angle + aperture_scan%min_angle
enddo

! Array control
if ( allocated(aperture_scan%aperture)) then
  if (size(aperture_scan%aperture) /= n_angle) deallocate(aperture_scan%aperture)
endif
if (.not. allocated(aperture_scan%aperture)) allocate(aperture_scan%aperture(n_angle))

if( aperture_scan%param%n_adts .gt. 0 ) then
  call custom_dynamic_aperture (lat, angle_list(1), aperture_scan, aperture_scan%aperture(1), .true., .false.)
  call custom_dynamic_aperture (lat, angle_list(n_angle), aperture_scan, aperture_scan%aperture(n_angle), .true., .false.)
  call custom_dynamic_aperture (lat, angle_list((n_angle+1)/2), aperture_scan, aperture_scan%aperture((n_angle+1)/2), .false., .true.)
  do i=1, n_angle
    if( i==1 .or. i==n_angle .or. i==(n_angle+1)/2 ) then
      cycle
    endif
    call custom_dynamic_aperture (lat, angle_list(i), aperture_scan, aperture_scan%aperture(i), .false., .false.)
    if( i .lt. (n_angle+1)/2 ) then
      if( aperture_scan%aperture(i)%x .gt. aperture_scan%aperture(1)%x ) aperture_scan%aperture(i)%x = aperture_scan%aperture(1)%x
    else
      if( aperture_scan%aperture(i)%x .lt. aperture_scan%aperture(n_angle)%x ) aperture_scan%aperture(i)%x = aperture_scan%aperture(n_angle)%x
    endif
    if( aperture_scan%aperture(i)%y .gt. aperture_scan%aperture((n_angle+1)/2)%y ) aperture_scan%aperture(i)%y = aperture_scan%aperture((n_angle+1)/2)%y
  enddo
else
  do i=1, n_angle
    call custom_dynamic_aperture (lat, angle_list(i), aperture_scan, aperture_scan%aperture(i), .false., .false.)
  enddo
endif

end subroutine custom_dynamic_aperture_scan

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine dynamic_aperture (lat, theta_xy, aperture_scan, aperture)
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
!   theta_xy    -- Real(rp): Angle of radial line (in radians) in x-y space.
!                    Angle is "normalized" by %x_init, %y_init.
!   aperture_scan -- aperture_param_struct: Structure holding the input data:
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
!     aperture_param     -- aperture_param_struct:
!       %closed_orbit -- Closed orbit coordinates
!
!-

subroutine custom_dynamic_aperture (lat, theta_xy, aperture_scan, aperture, do_adts_x, do_adts_y)
use dynap_mod
use adts_mod

implicit none

type (lat_struct)  lat
type (coord_struct), allocatable :: orbit(:)
type (coord_struct), allocatable :: orbit_temp(:,:)
type (custom_aperture_data_struct)  aperture
type (custom_aperture_scan_struct)  aperture_scan
type (custom_aperture_param_struct)  aperture_param
logical do_adts_x, do_adts_y

real(rp) theta_xy, x0, x1, x2, y0, y1, y2

integer it, turn_lost, track_state
integer ix

character(40) :: r_name = 'dynamic_aperture'

logical aperture_bracketed

real(rp) adts_x, adts_y
real(rp) d_nu_x, d_nu_y

! init setup

aperture_param = aperture_scan%param

bmad_com%aperture_limit_on = .true.

!call reallocate_coord (orbit, lat%n_ele_max)
allocate(orbit(0:lat%n_ele_track))
if(do_adts_x .or. do_adts_y) then
  allocate(orbit_temp(aperture_param%n_adts,0:lat%n_ele_track))
endif

! Find starting point

x0 = 0
y0 = 0
x1 = local_cosphi(theta_xy,aperture_scan%Sx,aperture_scan%Sy) * aperture_param%init_len
y1 = local_sinphi(theta_xy,aperture_scan%Sx,aperture_scan%Sy) * aperture_param%init_len

aperture_bracketed = .false.

! use a binary search to find where the aparture is along the line

test_loop: do
  orbit(0) = aperture_param%closed_orbit
  orbit(0)%vec(1) = orbit(0)%vec(1) + x1 
  orbit(0)%vec(3) = orbit(0)%vec(3) + y1 

  do it = 1, aperture_param%n_turn
    call track_all (lat, orbit, track_state=track_state)

    if (track_state /= moving_forward$) then
      aperture%plane = orbit(track_state)%state
      exit
    endif
    if(do_adts_x .or. do_adts_y) then
      if( it .ge. aperture_param%n_turn-aperture_param%n_adts+1 ) then
        ix = it-(aperture_param%n_turn-aperture_param%n_adts)
        orbit_temp(ix,:) = orbit(:)
      endif
    endif
    orbit(0) = orbit(lat%n_ele_track)
  enddo

  if(do_adts_x .or. do_adts_y) then
    if( track_state == moving_forward$ ) then
      adts_x = 0.0d0
      adts_y = 0.0d0
      do it=1,aperture_param%n_adts
        call accumulate_phase_advance(lat, orbit_temp(it,:), d_nu_x, d_nu_y)
        adts_x = adts_x + d_nu_x
        adts_y = adts_y + d_nu_y
      enddo
      adts_x = adts_x/twopi/aperture_param%n_adts
      adts_y = adts_y/twopi/aperture_param%n_adts
      if(do_adts_x) then
        if( adts_x .gt. aperture_param%adts_x_max) track_state = 1
        if( adts_x .lt. aperture_param%adts_x_min) track_state = 1
        if( adts_y .gt. aperture_param%adts_y_max) track_state = 1
        if( adts_y .lt. aperture_param%adts_y_min) track_state = 1
      endif
      if(do_adts_y) then
        if( adts_y .gt. aperture_param%adts_y_max) track_state = 1
        if( adts_y .lt. aperture_param%adts_y_min) track_state = 1
      endif
    endif
  endif

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
    x1 = x1 + local_cosphi(theta_xy,aperture_scan%Sx,aperture_scan%Sy) * aperture_param%step_len
    y1 = y1 + local_sinphi(theta_xy,aperture_scan%Sx,aperture_scan%Sy) * aperture_param%step_len
  endif

  ! check to see if there is an aperture

  if (x1 > 1.0d0 .or. y1 > 1.0d0) then
    call out_io(s_error$, r_name, 'CANNOT FIND APERTURE LIMIT')
    if (global_com%exit_on_error) call err_exit
  endif

  ! see if we are accurate enough

  if (aperture_bracketed) then
    if (sqrt((x2-x0)**2 + (y2-y0)**2) <= 2*aperture_param%accuracy) exit test_loop
  endif

enddo test_loop

! fill in the info

aperture%x = aperture_param%closed_orbit%vec(1) + x0
aperture%y = aperture_param%closed_orbit%vec(3) + y0
aperture%i_turn = turn_lost

deallocate(orbit)
if(allocated(orbit_temp)) deallocate(orbit_temp)

end subroutine

function local_cosphi(th,Sx,Sy) result(x)
  real(rp) th, Sx, Sy, x
  x = Sx * cos(th) / sqrt(Sx**2 * cos(th)**2 + Sy**2 * sin(th)**2 )
end function

function local_sinphi(th,Sx,Sy) result(x)
  real(rp) th, Sx, Sy, x
  x = Sy * sin(th) / sqrt(Sx**2 * cos(th)**2 + Sy**2 * sin(th)**2 )
end function

end module
