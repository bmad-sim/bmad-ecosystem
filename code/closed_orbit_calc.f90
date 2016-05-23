!+                           
! Subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, ix_branch, err_flag)
!
! Subroutine to calculate the closed orbit for a circular machine.
! Closed_orbit_calc uses the 1-turn transfer matrix to converge upon a  
! solution. 
!
! For i_dim = 4 this routine tracks through the lattice with the RF turned off.
! and the particle energy will be determined by closed_orb(0)%vec(6) for direction = 1
! and closed_orb(n0)%vec(6) with n0 = lat%n_ele_track for direction = -1.
! The z component (closed_orb(i)%vec(5)) will be set to zero at the beginning
! of tracking and will generally not be zero at the end.
!
! i_dim = 5 simulates the affect of the RF that makes the beam change 
! its energy until the change of path length in the closed orbit over 
! one turn is zero. Tracking is done with RF off. This can be useful in
! determining the affect of kicks on the beam energy (this can, of course,
! also be done with i_dim = 6).
!
! i_dim = 6 finds the closed orbit with energy variation. The RF needs to be
! turned on in this case. Additionally, to simulate cases where the RF frequency
! is not a multiple of the revolution harmonic (EG in a dispersion measurement), 
! lat%absolute_time_tracking needs to be set to True and the phi0_fieldmap attributes 
! of the RF cavities should be adjusted using autoscale_phase_and_amp.
!
! Note: This routine uses the 1-turn matrix lat%param%t1_no_RF or 
! lat%param%t1_with_RF in the computations. If you have changed conditions 
! significantly enough you might want to force a remake of the 1-turn matrices
! by calling clear_lat_1turn_mats.
!
! The closed orbit calculation stops when the following condition is satisfied:
!   amp_del < amp_co * bmad_com%rel_tol_tracking + bmad_com%abs_tol_tracking
! Where:
!   amp_co = abs(closed_orb[at beginning of lattice])
!   amp_del = abs(closed_orb[at beginning] - closed_orb[at end])
!   closed_orb = vector: (x, px, y, py, z, pz)
! closed_orb[at end] is calculated by tracking through the lattice starting 
! from closed_orb[at start].
!
! Note: See also closed_orbit_from_tracking as an alternative method
! of finding the closed orbit.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat            -- lat_struct: Lat to track through.
!   closed_orb(0:) -- Coord_struct, allocatable: closed_orb(n0) 
!                      is the initial guess where n0 = 0 for direction = 1 and 
!                      n0 = lat%n_ele_track for direction = -1. Additionally, 
!                      if i_dim = 4, then closed_orb(n0)%vec(6) is used as the energy 
!                      around which the closed orbit is calculated.
!   i_dim          -- Integer, optional: Phase space dimensions to use:
!                     = 4  Transverse closed orbit at constant energy (RF off).
!                          (dE/E = closed_orb(0)%vec(6))
!                     = 5 Transverse closed orbit at constant energy (RF off) with 
!                          the energy adjusted so that vec(5) is the same 
!                          at the beginning and at the end.
!                     = 6 True closed orbit.
!                     Default: Use 4 or 6 depending upon if RF is on or off.
!   direction      -- Integer, optional: Direction of tracking. 
!                       +1 --> forwad (default), -1 --> backward.
!                       The closed orbit will be dependent on direction only
!                       in the case that radiation damping is turned on.
!   ix_branch      -- Integer, optional: Lattice branch to find the closed orbit of. 
!                       Default is 0 (main branch).
!
!   bmad_com       -- Bmad_common_struct: Bmad common block.
!     %rel_tol_tracking -- Relative error. See above. Default = 1d-8
!     %abs_tol_tracking -- Absolute error. See above. Default = 1d-10
!
! Output:
!   closed_orb(0:) -- Coord_struct, allocatable: Closed orbit. closed_orb(i)
!                      is the orbit at the exit end of the ith element.
!   err_flag       -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, ix_branch, err_flag)

use bmad_interface, except_dummy => closed_orbit_calc
use bookkeeper_mod, only: set_on_off, restore_state$, off_and_save$
use reverse_mod, only: lat_reverse
use super_recipes_mod
use eigen_mod

implicit none

type (lat_struct), target ::  lat
type (lat_struct), target, save :: rev_lat
type (lat_struct), pointer :: this_lat
type (ele_struct), pointer :: ele, ele_start
type (branch_struct), pointer :: branch
type (coord_struct), allocatable, target ::  closed_orb(:), co_saved(:)
type (coord_struct), pointer :: start, end
type (bmad_common_struct) bmad_com_saved

type matrix_save
  real(rp) mat6(6,6)
  real(rp) vec0(6)
end type
type (matrix_save), allocatable :: m(:)

real(rp) t1(6,6), del_orb(6)
real(rp) :: amp_co(6), amp_del(6), dt, amp, dorb(6), old_start(6), old_end(6)
real(rp) z0, dz, z_here, this_amp, dz_norm, max_del, this_del, max_eigen
real(rp) a_lambda, chisq, old_chisq, rf_freq
real(rp), allocatable :: on_off_state(:), vec0(:), weight(:), a(:), covar(:,:), alpha(:,:)

complex(rp) eigen_val(6), eigen_vec(6,6)

integer, optional :: direction, ix_branch, i_dim
integer j, ie, i_loop, n_dim, n_ele, i_max, dir, track_state, n, status, j_max

logical, optional, intent(out) :: err_flag
logical err, error, allocate_m_done, stable_orbit_found, t1_needs_checking
logical, allocatable :: maska(:)

character(20) :: r_name = 'closed_orbit_calc'

!----------------------------------------------------------------------
! init
! Random fluctuations must be turned off to find the closed orbit.

allocate_m_done = .false.

dir = integer_option(+1, direction)
if (dir /= 1 .and. dir /= -1) then
  call out_io (s_error$, r_name, 'BAD DIRECTION ARGUMENT.')
  return
endif

if (dir == 1) then
  this_lat => lat
else
  call lat_reverse(lat, rev_lat)
  this_lat => rev_lat
endif  

if (present(err_flag)) err_flag = .true.
branch => this_lat%branch(integer_option(0, ix_branch))

call reallocate_coord (closed_orb, branch%n_ele_max)  ! allocate if needed

bmad_com_saved = bmad_com
bmad_com%radiation_fluctuations_on = .false.  
bmad_com%aperture_limit_on = .false.

n_ele = branch%n_ele_track

start => closed_orb(0)
end   => closed_orb(n_ele)
ele_start => branch%ele(0)

!----------------------------------------------------------------------
! Further init

if (present(i_dim)) then
  n_dim = i_dim ! dimension of transfer matrix
elseif (rf_is_on(branch)) then
  n_dim = 6
else
  n_dim = 4
endif

select case (n_dim)

! Constant energy case
! Turn off RF voltage if i_dim == 4 (for constant delta_E)

case (4, 5)

  ! Check if rf is on and if so issue a warning message

  if (rf_is_on(branch)) then
    call out_io (s_warn$, r_name, 'Inconsistant calculation: RF ON with i_dim = \i4\ ', i_dim)
  endif

  !

  bmad_com%radiation_damping_on = .false.  ! Want constant energy

  if (all(branch%param%t1_no_RF == 0)) &
              call transfer_matrix_calc (this_lat, branch%param%t1_no_RF, ix_branch = branch%ix_branch)
  t1 = branch%param%t1_no_RF
  start%vec(5) = 0

  !

  call set_on_off (rfcavity$, this_lat, off_and_save$, ix_branch = branch%ix_branch, saved_values = on_off_state)

! Variable energy case: i_dim = 6

case (6)
  if (all(branch%param%t1_with_RF == 0)) &
              call transfer_matrix_calc (this_lat, branch%param%t1_with_RF, ix_branch = branch%ix_branch)
  t1 = branch%param%t1_with_RF

  if (t1(6,5) == 0) then
    call out_io (s_error$, r_name, 'CANNOT DO FULL 6-DIMENSIONAL', &
                                   'CALCULATION WITH NO RF VOLTAGE!')
    return
  endif

  ! Assume that frequencies are comensurate otherwise a closed orbit does not exist.

  rf_freq = 1d30
  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (ele%key /= rfcavity$) cycle
    if (.not. ele%is_on) cycle
    rf_freq = min(rf_freq, abs(ele%value(rf_frequency$)))
  enddo

! Error

case default
  call out_io (s_error$, r_name, 'BAD I_DIM ARGUMENT: \i4\ ', i_dim)
  return
end select
        
! Orbit correction = (T-1)^-1 * (orbit_end - orbit_start)
!                  = t11_inv  * (orbit_end - orbit_start)


!--------------------------------------------------------------------------
! Because of nonlinearities we may need to iterate to find the solution

allocate(co_saved(0:ubound(closed_orb, 1)))
call init_coord (start, start, ele_start, start_end$, start%species)

allocate(vec0(n_dim), weight(n_dim), a(n_dim), maska(n_dim), covar(n_dim, n_dim), alpha(n_dim, n_dim))
vec0 = 0
maska = .false.
stable_orbit_found = .false.
a_lambda = -1
old_chisq = 1d30   ! Something large
old_start = 1d30
old_end = 1d30
a = start%vec(1:n_dim)
if (n_dim == 5) a(5) = start%vec(6)

!

t1_needs_checking = .true.  
i_max = 100

do i_loop = 1, i_max

  weight(1:n_dim) = 1 / (start%vec(1:n_dim) * bmad_com%rel_tol_tracking + bmad_com%abs_tol_tracking)**2
  if (n_dim == 5) weight(5) = 1 / (start%vec(6) * bmad_com%rel_tol_tracking + bmad_com%abs_tol_tracking)**2

  call super_mrqmin (vec0, weight, a, covar, alpha, chisq, co_func, a_lambda, status)

  if (a_lambda < 1d-10) a_lambda = 1d-10

  if (status < 0) then  
    call out_io (s_error$, r_name, 'Singular matrix Encountered!')
    call end_cleanup
    return
  endif

  if (i_loop == 1 .and. .not. stable_orbit_found) then
    a(1:4) = 0  ! Desperation move.
  elseif (i_loop == 2 .and. .not. stable_orbit_found) then
    call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING!!', 'ABORTING CLOSED ORBIT SEARCH.')
    call end_cleanup
    return
  else
    amp_co = abs(start%vec)
    dorb = end%vec - start%vec

    if (n_dim == 6 .and. branch%lat%absolute_time_tracking) then
      dt = (end%t - start%t) - nint((end%t - start%t) * rf_freq) / rf_freq
      dorb(5) = -end%beta * c_light * dt
    endif

    amp_del = abs(dorb)
    if (all(amp_del(1:n_dim) < amp_co(1:n_dim) * bmad_com%rel_tol_tracking + bmad_com%abs_tol_tracking)) exit
  endif

  if (track_state == moving_forward$ .and. chisq < old_chisq .and. status /= 1) co_saved = closed_orb

  ! If not converging fast enough, remake the transfer matrix.
  ! This is computationally intensive so only do this if the orbit has shifted significantly.
  ! Status == 1 means that the change in "a" was too large but t1 does not need to be remade.

  if (chisq > old_chisq/2 .and. status /= 1 .and. &
                maxval(abs(old_start(1:n_dim)-co_saved(0)%vec(1:n_dim))) > 1d-6) then ! If not converging
    if (.not. allocate_m_done) then
      allocate (m(branch%n_ele_max))
      allocate_m_done = .true.
    endif

    do n = 1, branch%n_ele_max
      m(n)%mat6 = branch%ele(n)%mat6
      m(n)%vec0 = branch%ele(n)%vec0
    enddo

    call lat_make_mat6 (lat, -1, co_saved, branch%ix_branch)
    call transfer_matrix_calc (lat, t1, ix_branch = branch%ix_branch)
    t1_needs_checking = .true.  ! New t1 matrix needs to be checked for stability.

    do n = 1, branch%n_ele_max
      branch%ele(n)%mat6 = m(n)%mat6
      branch%ele(n)%vec0 = m(n)%vec0
    enddo

    old_start = co_saved(0)%vec
  endif

  ! The super_mrqmin routine  will converge to the nearest fixed point even if that fixed point is an unstable one.
  ! For n_dim = 6, there are longitudinally unstable fixed points which we want to avoid.
  ! If we are near an unstable fixed point look for a better spot by shifting the particle in z in steps of pi/4.
  ! Note: due to inaccuracies, the maximum eigen value may be slightly over 1 at the stable fixed point.

  if (n_dim == 6 .and. t1_needs_checking .and. stable_orbit_found) then
    call mat_eigen (t1, eigen_val, eigen_vec, error)
    if (maxval(abs(eigen_val)) - 1 > 1d-5) then    ! Is unstable
      max_del = 1d10  ! Something large
      z0 = a(5)
      dz = start%beta * c_light / (8 * rf_freq)
      do j = -3, 4
        z_here = z0 + j * dz
        call track_this_lat(z_here, this_del, max_eigen)
        if (max_eigen - 1 > 1d-5) cycle
        if (this_del > max_del) cycle
        j_max = j
        max_del = this_del
      enddo

      ! Reset
      a(5) = z0 + j_max * dz
      call track_this_lat(a(5), this_del, max_eigen)

      ! If z needs to be shifted, reset super_mrqmin
      if (j_max /= 0) then
        a_lambda = -1               ! Signal super_mrqmin reset
      endif
    endif
    t1_needs_checking = .false.
  endif

  !

  old_chisq = chisq

  if (i_loop == i_max) then
    call out_io (s_error$, r_name, &
              'Closed orbit not converging! error in closed orbit: \es10.2\ ', &
              'If this error is acceptable, change bmad_com%rel_tol_tracking (\es10.2\) and/or', &
              'bmad_com%abs_tol_tracking (\es10.2\)', &
              r_array = [maxval(amp_del(1:n_dim)), bmad_com%rel_tol_tracking, bmad_com%abs_tol_tracking])
    call end_cleanup
    return
  endif

enddo

! Cleanup

call end_cleanup
if (present(err_flag)) err_flag = .false.

!------------------------------------------------------------------------------
! return rf cavities to original state

contains

subroutine end_cleanup

bmad_com = bmad_com_saved  ! Restore

if (n_dim == 4 .or. n_dim == 5) then
  call set_on_off (rfcavity$, this_lat, restore_state$, ix_branch = branch%ix_branch, saved_values = on_off_state)
endif

if (dir == -1) then
  closed_orb(1:n_dim) = closed_orb(n_dim:1:-1)
endif

end subroutine

!------------------------------------------------------------------------------
! contains

subroutine track_this_lat(z_set, del, max_eigen)

real(rp) z_set, dt, del, dorb(6), max_eigen

!

start%vec(5) = z_set
call track_all (this_lat, closed_orb, branch%ix_branch, track_state)

if (track_state /= moving_forward$) then
  max_eigen = 10
  return
endif

dorb = end%vec - start%vec
if (branch%lat%absolute_time_tracking) then
  dt = (end%t - start%t) - nint((end%t - start%t) * rf_freq) / rf_freq
  dorb(5) = -end%beta * c_light * dt
endif
del = maxval(abs(dorb))

call lat_make_mat6 (this_lat, -1, closed_orb, branch%ix_branch)
call transfer_matrix_calc (this_lat, t1, ix_branch = branch%ix_branch)
call mat_eigen (t1, eigen_val, eigen_vec, error)

max_eigen = maxval(abs(eigen_val))

end subroutine track_this_lat

!------------------------------------------------------------------------------
! contains

subroutine co_func (a_try, y_fit, dy_da, status)

real(rp), intent(in) :: a_try(:)
real(rp), intent(out) :: y_fit(:)
real(rp), intent(out) :: dy_da(:, :)
real(rp) del_orb(6), dz_norm

integer status, i

! For i_dim = 6, if at peak of RF then delta z may be singularly large. 
! To avoid this, veto any step where z changes by more than lambda_rf/10.

if (n_dim == 6) then
  dz_norm = abs(a_try(5)-a(5)) / (start%beta * c_light / (10 * rf_freq))
  if (dz_norm > 1) then
    status = 1  ! Veto step
    return
  endif
endif

!

if (n_dim == 5) then
  start%vec(1:4) = a_try(1:4)
  start%vec(6)   = a_try(5)
else
  start%vec(1:n_dim) = a_try
endif

call init_coord (start, start, ele_start, start_end$, start%species)
call track_all (this_lat, closed_orb, branch%ix_branch, track_state)

status = 0

!

del_orb = end%vec - start%vec

if (n_dim == 6 .and. branch%lat%absolute_time_tracking) then
  dt = (end%t - start%t) - nint((end%t - start%t) * rf_freq) / rf_freq
  del_orb(5) = -end%beta * c_light * dt
endif

if (track_state == moving_forward$) then
  stable_orbit_found = .true.
  y_fit = del_orb(1:n_dim)
else
  y_fit = track_state   ! Some large number
endif

dy_da = t1(1:n_dim,1:n_dim)
forall (i = 1:n_dim) dy_da(i,i) = dy_da(i,i) - 1

if (n_dim == 5) then
  dy_da(:,5) = t1(1:5,6)
endif

end subroutine co_func

end subroutine
