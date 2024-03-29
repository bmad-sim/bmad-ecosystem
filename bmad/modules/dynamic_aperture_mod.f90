module dynamic_aperture_mod

use coord_mod
use equal_mod

implicit none

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine dynamic_aperture_scan(aperture_scan, aperture_param, pz_start, lat, print_timing)
!
! Routine to do a set of dynamic aperture scans.
! 
! Input:
!   aperture_param      -- aperture_param_struct: Scan parameters
!   lat                 -- lat_struct: Lattice.
!   pz_start(:)         -- real(rp): Starting phase space pz values.
!   print_timing        -- logical, optional: If True print info on calculation time. Default is True.
!
! Output:
!   aperture_scan(:)    -- aperture_scan_struct, allocatable: Set of scans. One for each pz_start(:).
!-

subroutine dynamic_aperture_scan(aperture_scan, aperture_param, pz_start, lat, print_timing)

!$ use omp_lib

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (coord_struct), allocatable :: orbit(:)
type (coord_struct) orb0
type (ele_struct), pointer :: ele0
type (ele_pointer_struct), allocatable :: eles(:)
type (lat_pointer_struct), allocatable :: lats(:)
type (lat_ele_loc_struct) ele0_loc
type (aperture_scan_struct), allocatable, target :: aperture_scan(:)
type (aperture_scan_struct), pointer :: ap_scan
type (aperture_point_struct) :: ap_point, x_point, y_point
type (aperture_point_struct), allocatable :: ap_points(:)
type (aperture_param_struct), target :: aperture_param
type (aperture_param_struct), pointer :: ap_param

real(rp) pz_start(:)
real(rp), allocatable  :: angle_list(:)
real(rp) :: delta_angle, time0, time1, time2
real(rp) Sx, Sy

integer :: i, omp_n, n_loc, ns, ix_lat
logical err, thread_safe
logical, optional :: print_timing

character(*), parameter :: r_name = 'dynamic_aperture_scan'

! Angle preparation. Since the vertical aperture can be quite different from the horizontal one, the
! scan angles are scaled to give roughly equally spaced points.

ap_param => aperture_param

if (ap_param%n_angle < 1) then
  call out_io (s_error$, r_name, 'N_ANGLE MUST BE AT LEAST 1')
  return
endif

if (ap_param%start_ele == '') then
  branch => lat%branch(0)
  ele0 => branch%ele(0)
else
  call lat_ele_locator(ap_param%start_ele, lat, eles, n_loc, err)
  if (n_loc == 0) then
    call out_io (s_error$, r_name, 'NO LATTICE ELEMENT FOUND MATCHING: ' // ap_param%start_ele)
    return
  endif
  if (n_loc > 1) then
    call out_io (s_error$, r_name, 'MULTIPLE LATTICE ELEMENTS FOUND MATCHING: ' // ap_param%start_ele)
    return
  endif
  ele0 => eles(1)%ele
  branch => lat%branch(ele0%ix_branch)
endif

ele0_loc = lat_ele_loc_struct(ele0%ix_ele, ele0%ix_branch)

allocate(angle_list(ap_param%n_angle))

if (ap_param%n_angle == 1) then
  angle_list(1) = (ap_param%min_angle + ap_param%max_angle) / 2
else
  delta_angle = (ap_param%max_angle - ap_param%min_angle)/(ap_param%n_angle-1)
  do i=1, ap_param%n_angle
    angle_list(i) = (i-1)*delta_angle + ap_param%min_angle
  enddo
endif

!

thread_safe = .true.
do i = 1, branch%n_ele_track
  if (branch%ele(i)%tracking_method /= symp_lie_ptc$) cycle
  thread_safe = .false.
  !$ call out_io (s_warn$, r_name, 'Element with symp_lie_ptc tracking detected: ' // branch%ele(i)%name, &
  !$                               'PTC is not thread safe so will NOT run multithreaded!')
  exit
enddo

! Need separate lat instances when ramping since ramping will modify a lattice.

allocate (lats(max(2,ap_param%n_angle)))
lats(1)%lat => lat
lats(2)%lat => lat
!$ allocate(lats(2)%lat)
!$ lats(2)%lat = lat

!---------------

if (allocated(aperture_scan)) then
  if (size(aperture_scan) /= size(pz_start)) deallocate (aperture_scan)
endif
if (.not. allocated(aperture_scan)) allocate(aperture_scan(size(pz_start)))

allocate (ap_points(ap_param%n_angle))

!

call reallocate_coord(orbit, branch%n_ele_max)

do ns = 1, size(pz_start)
  ap_scan => aperture_scan(ns)

  if (rf_is_on(branch)) then
    if (ns == 1) call closed_orbit_calc (lat, orbit, 6, ix_branch = branch%ix_branch, err_flag = err)
    ap_scan%ref_orb = orbit(ele0%ix_ele)
    ap_scan%ref_orb%vec = ap_scan%ref_orb%vec + pz_start(ns) * [ele0%x%eta, ele0%x%etap, ele0%y%eta, ele0%y%etap, 0.0_rp, 1.0_rp]
  else
    orbit(0)%vec(6) = pz_start(ns)
    call closed_orbit_calc (lat, orbit, 4, ix_branch = branch%ix_branch, err_flag = err)
    ap_scan%ref_orb = orbit(ele0%ix_ele)
  endif

  if (err) then
    call out_io (s_error$, r_name, 'CANNOT FIND CLOSED ORBIT FOR STARTING ENERGY: ' // real_str(pz_start(ns), 4))
    if (allocated(ap_scan%point)) deallocate (ap_scan%point)
    cycle
  endif

  ! Allocate point array.

  if (allocated(ap_scan%point)) then
    if (size(ap_scan%point) /= ap_param%n_angle) deallocate(ap_scan%point)
  endif
  if (.not. allocated(ap_scan%point)) allocate(ap_scan%point(ap_param%n_angle))

  ! Auto-set x_init and y_init if they are zero

  call run_timer('ABS', time0)
  if (logic_option(.true., print_timing)) call out_io (s_info$, r_name, 'Dynamic aperture scan setup for pz: \f6.4\ ', r_array = [pz_start(ns)])

  if (ap_param%x_init == 0) ap_param%x_init = 0.001_rp
  if (ap_param%y_init == 0) ap_param%y_init = 0.001_rp

  !$OMP parallel sections if (thread_safe) private(branch, ele0)
  !$OMP section

  call set_branch_and_ele_for_omp(1, ele0_loc, lats, lat, branch, ele0)
  call dynamic_aperture_point (branch, ele0, ap_scan%ref_orb, 0.0_rp, ap_param, x_point, .false.)
  ap_param%x_init = x_point%x

  !$OMP section

  call set_branch_and_ele_for_omp(2, ele0_loc, lats, lat, branch, ele0)
  call dynamic_aperture_point (branch, ele0, ap_scan%ref_orb, pi/2, ap_param, y_point, .false.)
  ap_param%y_init = y_point%y

  !$OMP end parallel sections

  call run_timer('ABS', time1)
  if (logic_option(.true., print_timing)) call out_io (s_blank$, r_name, '  Scale factors calculated. dTime(min): \f8.2\ ', &
                                                                                r_array = [(time1-time0)/60])

  !$OMP parallel do if (thread_safe) private(branch, ele0, ix_lat)
  do i = 1, ap_param%n_angle
    ix_lat = 1
    !$ ix_lat = omp_get_thread_num() + 1
    call set_branch_and_ele_for_omp(ix_lat, ele0_loc, lats, lat, branch, ele0)

    if (abs(angle_list(i)) < 1d-6) then
      ap_points(i) = x_point
    elseif (abs(angle_list(i) - pi/2) < 1d-6) then
      ap_points(i) = y_point
    else
      call dynamic_aperture_point (branch, ele0, ap_scan%ref_orb, angle_list(i), ap_param, ap_points(i))
    endif
  end do
  !$OMP end parallel do

  ap_scan%point = ap_points

  call run_timer('ABS', time2)
  if (logic_option(.true., print_timing)) call out_io (s_blank$, r_name, '  Finished scan at pz in (min): \f8.2\ ', &
                                                                               r_array = [(time2 - time0) / 60])
enddo

deallocate (ap_points)

end subroutine dynamic_aperture_scan

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine set_branch_and_ele_for_omp(ix_lat, ele0_loc, lats, lat, branch, ele0)
!
! Routine to set branch and ele for the lattice to be used to track through.
!
! Different lattices are used to avoid problems when there is ramping since ramping
! will modify lattice element parameters.
!
! Input:
!   ix_lat      -- integer: Lattice index.
!   ele0_loc    -- lat_ele_loc_struct: ele0 location in lattice.
!   lats(:)     -- lat_pointer_struct: Pointers to lattices to track through.
!   lat         -- lat_struct: Original lattice.
!
! Output:
!   lats(:)     -- lat_pointer_struct: Properly initialized if needed.
!   branch      -- branch_struct, pointer: Branch to track through
!   ele0        -- starting element for tracking.
!-

subroutine set_branch_and_ele_for_omp(ix_lat, ele0_loc, lats, lat, branch, ele0)

type (lat_ele_loc_struct) ele0_loc
type (lat_pointer_struct) lats(:)
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele0

integer ix_lat
character(*), parameter :: r_name = 'set_branch_and_ele_for_omp'

!

if (ix_lat > size(lats)) then
  call out_io (s_fatal$, r_name, 'IX_LAT INDEX OUT OF RANGE! PLEASE REPORT THIS! ' // int_str(ix_lat))
  stop
endif

if (.not. associated(lats(ix_lat)%lat)) then
  allocate (lats(ix_lat)%lat)
  lats(ix_lat)%lat = lat
endif

branch => lats(ix_lat)%lat%branch(ele0_loc%ix_branch)
ele0 => branch%ele(ele0_loc%ix_ele)

end subroutine set_branch_and_ele_for_omp

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine dynamic_aperture_point (branch, ele0, orb0, theta_xy, ap_param, ap_point, check_xy_init)
!
! Subroutine to determine one dynamic aperture point by tracking.
! This routine works by determining where on a radial line y = const * x the aperture is.
! Here x and y are deviations from the reference orbit.
!
! Input:
!   branch          -- branch_struct: Lattice branch to track through.
!   ele0            -- ele_struct: Lattice element at start of tracking
!   orb0            -- Coord_struct: reference orbit at the start of tracking.
!   theta_xy        -- Real(rp): Angle of radial line (in radians) in x-y space.
!                         Angle is "normalized" by %x_init, %y_init.
!   ap_param        -- aperture_param_struct: Structure holding the input data:
!   check_xy_init   -- logical, optional: If True, do not check that aperture_param%x_init 
!                         and %y_init are non-zero. Default is True.
!
! Output:
!     ap_point      -- aperture_point_struct:
!-

subroutine dynamic_aperture_point (branch, ele0, orb0, theta_xy, ap_param, ap_point, check_xy_init)

type (branch_struct)  branch
type (ele_struct) ele0
type (coord_struct)  orb0
type (coord_struct), allocatable :: orbit(:)
type (aperture_point_struct)  ap_point
type (aperture_param_struct)  ap_param
type (bmad_common_struct) bmad_com_save

real(rp) theta_xy, x0, x1, x2, y0, y1, y2, r

integer n, it, turn_lost, track_state, n_search

character(*), parameter :: r_name = 'dynamic_aperture_point'

logical, optional :: check_xy_init
logical aperture_bracketed

! init setup

if (ap_param%x_init == 0 .and. abs(cos(theta_xy)) > 1e-10) then
  call out_io(s_fatal$, r_name, 'ap_param.x_init == 0') 
  if (global_com%exit_on_error) call err_exit
endif


if (ap_param%y_init == 0 .and. abs(sin(theta_xy)) > 1e-10) then
  call out_io(s_fatal$, r_name, 'ap_param.y_init == 0') 
  if (global_com%exit_on_error) call err_exit
endif

! Turn spin tracking off for calc.

bmad_com_save = bmad_com
bmad_com%spin_tracking_on  = .false.
bmad_com%aperture_limit_on = .true.

call reallocate_coord (orbit, branch%n_ele_max)

! Find starting point

x0 = 0
y0 = 0
x1 = cos(theta_xy) * ap_param%x_init
y1 = sin(theta_xy) * ap_param%y_init

aperture_bracketed = .false.

! use a binary search to find where the aparture is along the line

do n_search = 1, 1000000
  n = ele0%ix_ele
  orbit(n) = orb0
  orbit(n)%vec(1) = orbit(n)%vec(1) + x1 
  orbit(n)%vec(3) = orbit(n)%vec(3) + y1 

  ! track n_turns

  do it = 1, ap_param%n_turn
    call track_many (branch%lat, orbit, ele0%ix_ele, ele0%ix_ele, +1, branch%ix_branch, track_state)
    if (track_state /= moving_forward$) then
      ap_point%plane = orbit(track_state)%state
      exit
    endif
  enddo

  ! change search interval end

  if (track_state /= moving_forward$) then
    x2 = x1
    y2 = y1
    turn_lost = it
    aperture_bracketed = .true.
    ap_point%ix_ele = track_state 
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
    call out_io(s_error$, r_name, 'DA exceeds 1 meter ... giving up (set an aperture in the lattice file to avoid this).')
    if (global_com%exit_on_error) call err_exit
  endif

  ! See if we are accurate enough

  if (aperture_bracketed) then
    r = 0.5 * sqrt((x0+x2)**2 + (y0+y2)**2)
    if (sqrt((x2-x0)**2 + (y2-y0)**2) <= 2*(ap_param%abs_accuracy + r * ap_param%rel_accuracy)) exit
  endif

enddo

! fill in the info

ap_point%x = x1
ap_point%y = y1
ap_point%i_turn = turn_lost

bmad_com = bmad_com_save

end subroutine dynamic_aperture_point

end module
