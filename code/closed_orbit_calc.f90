!+                           
! Subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, ix_branch, err_flag)
!
! Subroutine to calculate the closed orbit for a circular machine.
! Closed_orbit_calc uses the 1-turn transfer matrix to converge upon a  
! solution. 
!
! If bmad_com%spin_tracking_on = True, the invariant spin direction will be calculated and
! put in closed_orb(:)%spin.
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
use super_recipes_mod
use eigen_mod
use rotation_3d_mod

implicit none

type (lat_struct), target ::  lat
type (ele_struct), pointer :: ele, ele_start
type (branch_struct), pointer :: branch
type (coord_struct), allocatable, target ::  closed_orb(:), co_saved(:)
type (coord_struct), pointer :: orb_start, orb_end
type (bmad_common_struct) bmad_com_saved
type (super_mrqmin_storage_struct) storage

type matrix_save
  real(rp) mat6(6,6)
  real(rp) vec0(6)
  type (coord_struct) map_ref_orb_in
  type (coord_struct) map_ref_orb_out
end type
type (matrix_save), allocatable :: m(:)

real(rp) t1(6,6), del_co(6), t11_inv(6,6), i1_int
real(rp) :: amp_co(6), amp_del(6), dt, amp, dorb(6), old_start(6), old_end(6)
real(rp) z0, dz, z_here, this_amp, dz_norm, max_del, this_del, max_eigen
real(rp) a_lambda, chisq, old_chisq, rf_freq, svec(3), mat3(3,3)
real(rp), allocatable :: on_off_state(:), vec0(:), weight(:), a(:)

complex(rp) eigen_val(6), eigen_vec(6,6)

integer, optional :: direction, ix_branch, i_dim
integer i, j, ie, i_loop, n_dim, n_ele, i_max, dir, track_state, n, status, j_max
integer ix_ele_start, ix_ele_end

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
  bmad_com = bmad_com_saved  ! Restore
  return
endif

if (present(err_flag)) err_flag = .true.
branch => lat%branch(integer_option(0, ix_branch))

call reallocate_coord (closed_orb, branch%n_ele_max)  ! allocate if needed

bmad_com_saved = bmad_com
bmad_com%radiation_fluctuations_on = .false.  
bmad_com%aperture_limit_on = .false.
bmad_com%spin_tracking_on = .false.

n_ele = branch%n_ele_track

if (dir == 1) then
  ix_ele_start = 0
  ix_ele_end   = n_ele
else
  ix_ele_start = n_ele
  ix_ele_end   = 0
endif

orb_start => closed_orb(ix_ele_start)
orb_end   => closed_orb(ix_ele_end)
ele_start => branch%ele(ix_ele_start)

call init_coord (orb_start, orb_start, ele_start, start_end$, orb_start%species, dir)

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
    call out_io (s_warn$, r_name, 'Inconsistant calculation: RF ON with i_dim = \i4\ ', &
                  'Using branch: ' // branch_name(branch), i_array = [i_dim])
  endif

  !

  bmad_com%radiation_damping_on = .false.  ! Want constant energy
  call set_on_off (rfcavity$, lat, off_and_save$, ix_branch = branch%ix_branch, saved_values = on_off_state)

  if (any(branch%param%t1_no_RF /= 0) .and. dir == 1) then
    t1 = branch%param%t1_no_RF
  else
    call this_t1_calc(branch, dir, .false., t1, err); if (err) return
    if (dir == 1) branch%param%t1_no_RF = t1
  endif

  orb_start%vec(5) = 0

  if (n_dim == 5) then  ! crude I1 integral calculation
    i1_int = 0
    do ie = 1, branch%n_ele_track
      ele => branch%ele(ie)
      if (ele%key == sbend$) then
        i1_int = i1_int + ele%value(l$) * ele%value(g$) * (branch%ele(ie-1)%x%eta + ele%x%eta) / 2
      endif
    enddo
  endif

! Variable energy case: i_dim = 6

case (6)
  if (any(branch%param%t1_with_RF /= 0) .and. dir == 1) then
    t1 = branch%param%t1_with_RF
  else
    call this_t1_calc(branch, dir, .false., t1, err); if (err) return
    if (dir == 1)  branch%param%t1_with_RF = t1
  endif

  if (t1(6,5) == 0) then
    call out_io (s_error$, r_name, 'CANNOT DO FULL 6-DIMENSIONAL', &
                                   'CALCULATION WITH NO RF VOLTAGE!', &
                                    'USING BRANCH: ' // branch_name(branch))
    bmad_com = bmad_com_saved  ! Restore
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
  bmad_com = bmad_com_saved  ! Restore
  return
end select
        
! Orbit correction = (T-1)^-1 * (orbit_end - orbit_start)
!                  = t11_inv  * (orbit_end - orbit_start)

!--------------------------------------------------------------------------
! Because of nonlinearities we may need to iterate to find the solution

allocate(co_saved(0:ubound(closed_orb, 1)))
allocate(vec0(n_dim), weight(n_dim), a(n_dim), maska(n_dim))
vec0 = 0
maska = .false.
stable_orbit_found = .false.
a_lambda = -1
old_chisq = 1d30   ! Something large
old_start = 1d30
old_end = 1d30
a = orb_start%vec(1:n_dim)
if (n_dim == 5) a(5) = orb_start%vec(6)

!

t1_needs_checking = .true.  
i_max = 100

do i_loop = 1, i_max

  weight(1:n_dim) = 1 / (abs(a) * bmad_com%rel_tol_tracking + bmad_com%abs_tol_tracking)**2

  call super_mrqmin (vec0, weight, a, chisq, co_func, storage, a_lambda, status)

  if (a_lambda < 1d-10) a_lambda = 1d-10

  if (status < 0) then  
    call out_io (s_error$, r_name, 'Singular matrix Encountered!', 'Using branch: ' // branch_name(branch))
    call end_cleanup
    return
  endif

  if (i_loop == 1 .and. .not. stable_orbit_found) then
    a(1:4) = 0  ! Desperation move.
  elseif (i_loop == 2 .and. .not. stable_orbit_found) then
    call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING!!', 'ABORTING CLOSED ORBIT SEARCH.', &
                                   'USING BRANCH: ' // branch_name(branch))
    call end_cleanup
    return
  else
    amp_co = abs(orb_start%vec)
    dorb = orb_end%vec - orb_start%vec

    if (n_dim == 6 .and. branch%lat%absolute_time_tracking) then
      dt = (orb_end%t - orb_start%t) - nint((orb_end%t - orb_start%t) * rf_freq) / rf_freq
      dorb(5) = -orb_end%beta * c_light * dt
    endif

    amp_del = abs(dorb)
    if (all(amp_del(1:n_dim) < amp_co(1:n_dim) * bmad_com%rel_tol_tracking + bmad_com%abs_tol_tracking)) exit
  endif

  if (track_state == moving_forward$ .and. chisq < old_chisq .and. status /= 1) co_saved = closed_orb

  ! If not converging fast enough, remake the transfer matrix.
  ! This is computationally intensive so only do this if the orbit has shifted significantly.
  ! Status == 1 means that the change in "a" was too large but t1 does not need to be remade.

  if (chisq > old_chisq/2 .and. status /= 1) then
    if (maxval(abs(old_start(1:n_dim)-co_saved(0)%vec(1:n_dim))) > 1d-6) then ! If not converging
      if (.not. allocate_m_done) then
        allocate (m(branch%n_ele_max))
        allocate_m_done = .true.
      endif

      do n = 1, branch%n_ele_max
        m(n)%mat6 = branch%ele(n)%mat6
        m(n)%vec0 = branch%ele(n)%vec0
        m(n)%map_ref_orb_in  = branch%ele(n)%map_ref_orb_in
        m(n)%map_ref_orb_out = branch%ele(n)%map_ref_orb_out
      enddo

      call this_t1_calc (branch, dir, .true., t1, err); if (err) return
      t1_needs_checking = .true.  ! New t1 matrix needs to be checked for stability.

      do n = 1, branch%n_ele_max
        branch%ele(n)%mat6 = m(n)%mat6
        branch%ele(n)%vec0 = m(n)%vec0
        branch%ele(n)%map_ref_orb_in  = m(n)%map_ref_orb_in
        branch%ele(n)%map_ref_orb_out = m(n)%map_ref_orb_out
      enddo

      old_start = co_saved(0)%vec

    else
      call make_t11_inv (err)
      if (err) then
        call end_cleanup
        return
      endif

      del_co(1:n_dim) = matmul(t11_inv(1:n_dim,1:n_dim), dorb(1:n_dim)) 

      ! For n_dim = 6, if at peak of RF then del_co(5) may be singularly large. 
      ! To avoid this, limit z step to be no more than lambda_rf/10

      if (n_dim == 5) then
        del_co(5) = dorb(5) / i1_int      

      elseif (n_dim == 6) then
        dz_norm = abs(del_co(5)) / (orb_start%beta * c_light / (10 * rf_freq))
        if (dz_norm > 1) then
          del_co = del_co / dz_norm
        endif
      endif

      a = a + del_co(1:n_dim)

    endif
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
      dz = orb_start%beta * c_light / (8 * rf_freq)
      do j = -3, 4
        z_here = z0 + j * dz
        call track_this_lat(z_here, this_del, max_eigen)
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
              'Using branch: ' // branch_name(branch), &
              r_array = [maxval(amp_del(1:n_dim)), bmad_com%rel_tol_tracking, bmad_com%abs_tol_tracking])
    call end_cleanup
    return
  endif

enddo

! Calc invarient spin axis.

if (bmad_com_saved%spin_tracking_on) then
  bmad_com%spin_tracking_on = .true.
  do i = 1, 3
    orb_start%spin = 0
    orb_start%spin(i) = 1
    call track_many (lat, closed_orb, ix_ele_start, ix_ele_end, dir, branch%ix_branch, track_state)
    mat3(:,i) = orb_end%spin
  enddo
  call w_mat_to_axis_angle(mat3, orb_start%spin, branch%param%spin_tune)
  if (branch%param%spin_tune < 0) branch%param%spin_tune = branch%param%spin_tune + twopi
  call track_many (lat, closed_orb, ix_ele_start, ix_ele_end, dir, branch%ix_branch, track_state)
endif

! Cleanup

call end_cleanup
if (present(err_flag)) err_flag = .false.

!------------------------------------------------------------------------------
! return rf cavities to original state

contains

subroutine end_cleanup()

bmad_com = bmad_com_saved  ! Restore

if (n_dim == 4 .or. n_dim == 5) then
  call set_on_off (rfcavity$, lat, restore_state$, ix_branch = branch%ix_branch, saved_values = on_off_state)
endif

end subroutine

!------------------------------------------------------------------------------
! contains

subroutine track_this_lat(z_set, del, max_eigen)

real(rp) z_set, dt, del, dorb(6), max_eigen
logical err

!

orb_start%vec(5) = z_set
call track_many (lat, closed_orb, ix_ele_start, ix_ele_end, dir, branch%ix_branch, track_state)

if (track_state /= moving_forward$) then
  max_eigen = 10
  return
endif

dorb = orb_end%vec - orb_start%vec
if (branch%lat%absolute_time_tracking) then
  dt = (orb_end%t - orb_start%t) - nint((orb_end%t - orb_start%t) * rf_freq) / rf_freq
  dorb(5) = -orb_end%beta * c_light * dt
endif
del = maxval(abs(dorb))

call this_t1_calc (branch, dir, .true., t1, err); if (err) return
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
  dz_norm = abs(a_try(5)-a(5)) / (orb_start%beta * c_light / (10 * rf_freq))
  if (dz_norm > 1) then
    status = 1  ! Veto step
    return
  endif
endif

!

if (n_dim == 5) then
  orb_start%vec(1:4) = a_try(1:4)
  orb_start%vec(6)   = a_try(5)
else
  orb_start%vec(1:n_dim) = a_try
endif

call init_coord (orb_start, orb_start, ele_start, start_end$, orb_start%species, dir)
call track_many (lat, closed_orb, ix_ele_start, ix_ele_end, dir, branch%ix_branch, track_state)

status = 0

!

del_orb = orb_end%vec - orb_start%vec

if (n_dim == 6 .and. branch%lat%absolute_time_tracking) then
  dt = (orb_end%t - orb_start%t) - nint((orb_end%t - orb_start%t) * rf_freq) / rf_freq
  del_orb(5) = -orb_end%beta * c_light * dt
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

!------------------------------------------------------------------------------
! contains

subroutine this_t1_calc (branch, dir, make_mat6, t1, err)

type (branch_struct) branch
type (coord_struct) start_saved, end_saved

real(rp) t1(6,6), delta
integer dir, i, track_state
logical make_mat6, err

!

err = .false.

if (dir == 1) then
  if (make_mat6) call lat_make_mat6 (branch%lat, -1, closed_orb, branch%ix_branch)
  call transfer_matrix_calc (branch%lat, t1, ix_branch = branch%ix_branch)

else
  call track_many (branch%lat, closed_orb, branch%n_ele_track, 0, dir, branch%ix_branch, track_state)
  start_saved = orb_start
  end_saved   = orb_end

  delta = 1d-6
  do i = 1, 6
    orb_start%vec(i) = start_saved%vec(i) + delta
    call track_many (branch%lat, closed_orb, branch%n_ele_track, 0, dir, branch%ix_branch, track_state)
    if (track_state /= moving_forward$) then 
      call out_io (s_error$, r_name, 'PARTICLE LOST TRACKING BACKWARDS. [POSSIBLE CAUSE: WRONG PARTICLE SPECIES.]', &
                                     'USING BRANCH: ' // branch_name(branch))
      call end_cleanup
      err = .true.
      return
    endif 
    t1(:,i) = (orb_end%vec - end_saved%vec) / delta
    orb_start%vec = start_saved%vec
  enddo
endif

end subroutine this_t1_calc

!------------------------------------------------------------------------------
! contains

subroutine make_t11_inv (err)

real(rp) mat(6,6)
logical ok1, ok2, err

!

err = .true.

ok1 = .true.
call mat_make_unit (mat(1:n_dim,1:n_dim))
mat(1:n_dim,1:n_dim) = mat(1:n_dim,1:n_dim) - t1(1:n_dim,1:n_dim)
call mat_inverse(mat(1:n_dim,1:n_dim), t11_inv(1:n_dim,1:n_dim), ok2)

if (.not. ok1 .or. .not. ok2) then 
  call out_io (s_error$, r_name, 'MATRIX INVERSION FAILED!')
  return
endif

err = .false.

end subroutine

end subroutine closed_orbit_calc
