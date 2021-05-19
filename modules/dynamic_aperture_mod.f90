module dynamic_aperture_mod

use coord_mod
use equal_mod

implicit none

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine dynamic_aperture_scan(branch, aperture_scan, parallel)
!
! Driver routine for dynamic_aperture1. 
! 
! If aperture_scan%param%x_init or %y_init == 0, 
! then a separate scan will be done to set them starting with 0.001 m 
!
! Input:
!   branch              -- branch_struct: Lattice branch to track through.
!   aperture_scan       -- aperture_scan_struct: 
!     %ref_orb            -- Reference orbit.
!   parallel            -- logical, optional :: Use OpenMP parallel routine. 
!                                               Default: False
!
! Output:
!   aperture_scan       -- aperture_scan_struct 
!     %aperture(:)      -- aperture_data_struct: apertures for each angle
!-

subroutine dynamic_aperture_scan(branch, aperture_scan, parallel)

!$ use omp_lib

type (branch_struct) :: branch
type (aperture_scan_struct), target :: aperture_scan
type (aperture_data_struct) :: aperture
type (aperture_data_struct), allocatable :: apertures(:)
type (aperture_param_struct), pointer :: ap_param

real(rp), allocatable  :: angle_list(:)
real(rp) :: delta_angle, time0, time1, time2
real(rp) x_lims(1:branch%n_ele_track)
real(rp) y_lims(1:branch%n_ele_track)
real(rp) Sx, Sy
integer :: i, omp_n

logical, optional :: parallel
character(*), parameter :: r_name = 'dynamic_aperture_scan'

! Angle preparation. Since the vertical aperture can be quite different from the horizontal one, the
! scan angles are scaled to give roughly equally spaced points.

ap_param => aperture_scan%param

if (ap_param%n_angle < 1) then
  call out_io (s_error$, r_name, 'N_ANGLE MUST BE AT LEAST 1')
  return
endif

allocate(angle_list(ap_param%n_angle))

if (ap_param%n_angle == 1) then
  aperture_scan%S_xy = 1.0
  angle_list(1) = (ap_param%min_angle + ap_param%max_angle) / 2
else
  delta_angle = (ap_param%max_angle - ap_param%min_angle)/(ap_param%n_angle-1)
  where (branch%ele(1:branch%n_ele_track)%value(x1_limit$) > 0)
    x_lims(:) = branch%ele(1:)%value(x1_limit$)
  elsewhere
    x_lims(:) = 1.0
  endwhere
  where (branch%ele(1:branch%n_ele_track)%value(y1_limit$) > 0)
    y_lims(:) = branch%ele(1:)%value(y1_limit$)
  elsewhere
    y_lims(:) = 1.0
  endwhere
  Sx = minval(x_lims / sqrt(branch%ele(1:branch%n_ele_track)%a%beta)) * sqrt(branch%ele(1)%a%beta)
  Sy = minval(y_lims / sqrt(branch%ele(1:branch%n_ele_track)%b%beta)) * sqrt(branch%ele(1)%b%beta)
  aperture_scan%S_xy = Sx/Sy
  do i=1, ap_param%n_angle
    angle_list(i) = (i-1)*delta_angle + ap_param%min_angle
  enddo
endif

! Array control

if ( allocated(aperture_scan%aperture)) then
  if (size(aperture_scan%aperture) /= ap_param%n_angle) deallocate(aperture_scan%aperture)
endif
if (.not. allocated(aperture_scan%aperture)) allocate(aperture_scan%aperture(ap_param%n_angle))

! Auto-set x_init and y_init if they are zero

call run_timer('ABS', time0)
call out_io (s_info$, r_name, 'Angle scan setup...')

!$OMP parallel sections
!$OMP section

if (ap_param%x_init == 0) then
  ap_param%x_init = 0.001_rp
  call dynamic_aperture1 (branch, aperture_scan%ref_orb, 0.0_rp, aperture_scan%S_xy, ap_param, aperture, .false.)
  ap_param%x_init = aperture%x
endif

!$OMP section

if (ap_param%y_init == 0) then
  ap_param%y_init = 0.001_rp
  call dynamic_aperture1 (branch, aperture_scan%ref_orb, pi/2, aperture_scan%S_xy, ap_param, aperture, .false.)
  ap_param%y_init = aperture%y
endif

!$OMP end parallel sections

call run_timer('ABS', time1)
call out_io (s_info$, r_name, 'Angle scale factors calculated. dTime(min): \f8.2\ ', &
                              '  Starting angle scan...', r_array = [(time1-time0)/60])

! Only call the parallel routine if there is more than one thread available 

omp_n = 1
!$ omp_n = omp_get_max_threads()

if (logic_option(.false., parallel) .and. omp_n > 1) then
  allocate (apertures(ap_param%n_angle))
  !$OMP parallel do
  do i = 1, ap_param%n_angle
    call dynamic_aperture1 (branch, aperture_scan%ref_orb, angle_list(i), aperture_scan%S_xy, ap_param, apertures(i))
  end do  
  !$OMP end parallel do

  aperture_scan%aperture = apertures
  deallocate (apertures)

else
  do i = 1, ap_param%n_angle
    call dynamic_aperture1 (branch, aperture_scan%ref_orb, angle_list(i), aperture_scan%S_xy, ap_param, aperture_scan%aperture(i))
    call run_timer('ABS', time2)
    if (time2 - time1 < 60) cycle
    call out_io (s_info$, r_name, 'Finished angle scan #' // int_str(i) // ' of ' // int_str(ap_param%n_angle) // ' in (min): \f8.2\ ', &
    r_array = [(time2 - time0) / 60])
    time1 = time2
  enddo
endif

end subroutine dynamic_aperture_scan

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine dynamic_aperture1 (branch, orb0, theta_xy, S_xy aperture_param, aperture, check_xy_init)
!
! Subroutine to determine one dynamic aperture point by tracking.
! This routine works by determining where on a radial line y = const * x
! the aperture is. Here x and y are deviations from the reference orbit.
!
! Input:
!   branch         -- branch_struct: Lattice branch to track through.
!   orb0           -- Coord_struct: reference orbit at the start.
!   theta_xy       -- Real(rp): Angle of radial line (in radians) in x-y space.
!                         Angle is "normalized" by %x_init, %y_init.
!   S_xy           -- real(rp): used to scale angles for linear aperture aspect ratio.
!   aperture_param -- aperture_param_struct: Structure holding the input data:
!   check_xy_init  -- logical, optional: If True, do not check that aperture_param%x_init 
!                         and %y_init are non-zero. Default is True.
!
! Output:
!     aperture  -- aperture_data_struct:
!
! Note: The radial lines are spaced equally in angle using coordinates
!       normalized by %X_INIT and %Y_INIT
!-

subroutine dynamic_aperture1 (branch, orb0, theta_xy, S_xy, aperture_param, aperture, check_xy_init)

type (branch_struct)  branch
type (bmad_common_struct)  com_save
type (coord_struct)  orb0
type (coord_struct), allocatable :: orbit(:)
type (aperture_data_struct)  aperture
type (aperture_param_struct)  aperture_param

real(rp) theta_xy, x0, x1, x2, y0, y1, y2
real(rp) S_xy
real(rp) init_len

integer it, turn_lost, track_state

character(*), parameter :: r_name = 'dynamic_aperture1'

logical, optional :: check_xy_init
logical aperture_bracketed

! init setup

if (logic_option(.true., check_xy_init)) then
  if (aperture_param%x_init == 0) then
    call out_io(s_fatal$, r_name, 'aperture_param.x_init == 0') 
    if (global_com%exit_on_error) call err_exit
  endif


  if (aperture_param%y_init == 0) then
    call out_io(s_fatal$, r_name, 'aperture_param.y_init == 0') 
    if (global_com%exit_on_error) call err_exit
  endif
endif

com_save = bmad_com
bmad_com%aperture_limit_on = .true.

call reallocate_coord (orbit, branch%n_ele_max)

! Find starting point

x0 = 0
y0 = 0
init_len = sqrt(aperture_param%x_init**2 + aperture_param%y_init**2)
x1 = cosphi(theta_xy,S_xy) * init_len
y1 = sinphi(theta_xy,S_xy) * init_len

aperture_bracketed = .false.

! use a binary search to find where the aparture is along the line

test_loop: do

  orbit(0) = orb0
  orbit(0)%vec(1) = orbit(0)%vec(1) + x1 
  orbit(0)%vec(3) = orbit(0)%vec(3) + y1 

  ! track n_turns

  do it = 1, aperture_param%n_turn
    call track_all (branch%lat, orbit, branch%ix_branch, track_state)
    if (track_state /= moving_forward$) then
      aperture%plane = orbit(track_state)%state
      exit
    endif
    orbit(0) = orbit(branch%n_ele_track)
  enddo

  ! change search interval end

  if (track_state /= moving_forward$) then
    x2 = x1
    y2 = y1
    turn_lost = it
    aperture_bracketed = .true.
    aperture%ix_ele = track_state 
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

  if (x1 > 1.0 .or. y1 > 1.0) then
    call out_io(s_error$, r_name, 'DA exceeds 1 meter ... giving up.')
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

bmad_com = com_save

!-----------------------------------------------------------------
contains

function cosphi(th,S_xy) result(x)
  real(rp) th, S_xy, x
  !x = Sx * cos(th) / sqrt(Sx**2 * cos(th)**2 + Sy**2 * sin(th)**2 )
  x = S_xy * cos(th) / sqrt(S_xy**2 * cos(th)**2 + sin(th)**2 )
end function cosphi

function sinphi(th,S_xy) result(x)
  real(rp) th, S_xy, x
  !x = Sy * sin(th) / sqrt(Sx**2 * cos(th)**2 + Sy**2 * sin(th)**2 )
  x = sin(th) / sqrt(S_xy**2 * cos(th)**2 + sin(th)**2 )
end function sinphi

end subroutine dynamic_aperture1

end module
