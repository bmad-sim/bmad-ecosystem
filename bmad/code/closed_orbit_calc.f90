!+                           
! Subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, ix_branch, err_flag, print_err)
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
! and closed_orb(nt)%vec(6) with nt = lat%n_ele_track for direction = -1.
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
! bmad_com%absolute_time_tracking needs to be set to True and the phi0_fieldmap attributes 
! of the RF cavities should be adjusted using autoscale_phase_and_amp.
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
! Input:
!   lat            -- lat_struct: Lat to track through.
!   closed_orb(0:) -- Coord_struct, allocatable: closed_orb(nt) 
!                      is the initial guess where nt = 0 for direction = 1 and 
!                      nt = lat%n_ele_track for direction = -1. Additionally, 
!                      if i_dim = 4, then closed_orb(nt)%vec(6) is used as the energy 
!                      around which the closed orbit is calculated.
!   i_dim          -- Integer, optional: Phase space dimensions to use:
!                     = 4  Transverse closed orbit at constant energy (RF off).
!                          (dE/E = closed_orb(0)%vec(6))
!                     = 5 Transverse closed orbit at constant energy (RF off) with 
!                          the energy adjusted so that vec(5) is the same 
!                          at the beginning and at the end.
!                     = 6 True closed orbit.
!                     Default: 4 if RF is off, 6 if RF is on.
!   direction      -- Integer, optional: Direction of tracking. 
!                       +1 --> forwad (default), -1 --> backward.
!                       The closed orbit will be dependent on direction only
!                       in the case that radiation damping is turned on.
!   ix_branch      -- Integer, optional: Lattice branch to find the closed orbit of. 
!                       Default is 0 (main branch).
!   print_err      -- logical, optional: Print error message if calc does not converge?
!                       Default is True. Note: Condition messages like no RF voltage with 
!                       i_dim = 6 will always be printed. 
!
!   bmad_com       -- Bmad_common_struct: Bmad common block.
!     %rel_tol_tracking -- Relative error. See above. Default = 1d-8
!     %abs_tol_tracking -- Absolute error. See above. Default = 1d-10
!     %spin_tracking_on -- If True, the closed orbit invariant spin will be calculated and put in closed_orb%spin
!
! Output:
!   closed_orb(0:) -- Coord_struct, allocatable: Closed orbit. closed_orb(i)
!                      is the orbit at the exit end of the ith element.
!     %vec(6)             -- Closed orbit phase space
!     %spin(3)            -- Closed orbit invariant spin if bmad_com%spin_tracking_on = T.
!   lat            -- lat_struct:
!     %branch(ix_branch)%param%spin_tune -- Spin tune if bmad_com%spin_tracking_on = T.
!   err_flag       -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, ix_branch, err_flag, print_err)

use bmad_interface, except_dummy => closed_orbit_calc
use super_recipes_mod
use lmdif_mod

implicit none

type (lat_struct), target ::  lat
type (ele_struct), pointer :: ele, ele_start
type (branch_struct), pointer :: branch
type (coord_struct), allocatable, target ::  closed_orb(:), co_best(:)
type (coord_struct), pointer :: start_orb, end_orb
type (coord_struct), pointer :: c_orb(:), c_best(:), s_orb, e_orb
type (coord_struct) original_start_orb
type (bmad_common_struct) bmad_com_saved
type (super_mrqmin_storage_struct) storage

type matrix_save
  real(rp) mat6(6,6)
  real(rp) vec0(6)
  type (coord_struct) map_ref_orb_in
  type (coord_struct) map_ref_orb_out
end type
type (matrix_save), allocatable :: m(:)

type closed_orb_com_struct   ! Common block for
  real(rp) :: t1(6,6)
  real(rp), allocatable :: a_vec(:)
  real(rp) rf_wavelen, dz_step
end type

type (closed_orb_com_struct), target :: coc
type (closed_orb_com_struct), pointer :: cocp

real(rp) del_co(6), t11_inv(6,6), i1_int
real(rp) :: amp_co(6), amp_del(6), dt, amp, dorb(6), start_orb_t1(6)
real(rp) z0, dz, z_here, this_amp, dz_norm, max_eigen, min_max_eigen, z_original, betas(2)
real(rp) a_lambda, chisq, old_chisq, rf_freq, svec(3), mat3(3,3), this(5)
real(rp), allocatable, target :: on_off_state(:), scatter_state(:), vec0(:), dvec(:), weight(:), merit_vec(:), dy_da(:,:)

integer, optional :: direction, ix_branch, i_dim
integer i, j, ie, i2, i_loop, n_dim, n_ele, i_max, dir, track_state, n, status, j_max
integer ix_ele_start, ix_ele_end

logical, optional, intent(out) :: err_flag
logical, optional, intent(in) :: print_err
logical err, error, allocate_m_done, has_been_alive_at_end, check_for_z_stability, printit, at_end, try_lmdif
logical invert_step_tried
logical, allocatable :: maska(:)

character(20) :: r_name = 'closed_orbit_calc'

!----------------------------------------------------------------------
! init
! Random fluctuations must be turned off to find the closed orbit.

allocate_m_done = .false.
printit = logic_option(.true., print_err)
cocp => coc   ! To get around Ifort bug.

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
bmad_com%spin_sokolov_ternov_flipping_on = .false.
bmad_private%random_on = .false.

n_ele = branch%n_ele_track
betas = 1

if (dir == 1) then
  ix_ele_start = 0
  ix_ele_end   = n_ele
else
  ix_ele_start = n_ele
  ix_ele_end   = 0
endif

start_orb => closed_orb(ix_ele_start)
end_orb   => closed_orb(ix_ele_end)
ele_start => branch%ele(ix_ele_start)

c_orb => closed_orb   ! Used to get around ifort bug
s_orb => start_orb    ! Used to get around ifort bug
e_orb => end_orb      ! Used to get around ifort bug

call init_coord (start_orb, start_orb, ele_start, start_end$, start_orb%species, dir)
original_start_orb = start_orb

call set_on_off (foil$, lat, off_and_save$, ix_branch = branch%ix_branch, saved_values = scatter_state, &
                                                                      attribute = 'SCATTER_METHOD', set_val = off$)

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
    coc%t1 = branch%param%t1_no_RF
  else
    call this_t1_calc(branch, dir, .false., coc%t1, betas, start_orb_t1, err)
    if (err) then
      call end_cleanup(branch, .false.)
      return
    endif
    if (dir == 1) branch%param%t1_no_RF = coc%t1
  endif

  start_orb%vec(5) = 0

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
  ! Assume that frequencies are comensurate otherwise a closed orbit does not exist.

  rf_freq = 1d30
  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (ele%key /= rfcavity$) cycle
    if (.not. ele%is_on) cycle
    rf_freq = min(rf_freq, abs(ele%value(rf_frequency$)))
  enddo

  z_original = start_orb%vec(5) 
  coc%rf_wavelen = c_light * (branch%ele(ix_ele_start)%value(p0c$) / branch%ele(ix_ele_start)%value(e_tot$)) / rf_freq

  !

  if (.not. rf_is_on(branch)) then
    call out_io (s_error$, r_name, 'CANNOT DO FULL 6-DIMENSIONAL', &
                                   'CALCULATION WITH NO RF VOLTAGE!', &
                                   'USING BRANCH: ' // branch_name(branch))
    call end_cleanup(branch, .false.)
    return
  endif

  call transfer_matrix_calc (branch%lat, coc%t1, ix_branch = branch%ix_branch)

  if (all(coc%t1 == 0)) then  ! Might be that matrices have not yet been made
    call this_t1_calc(branch, dir, .false., coc%t1, betas, start_orb_t1, err)
    if (err) then
      call end_cleanup(branch, .false.)
      return
    endif
  endif

! Error

case default
  call out_io (s_error$, r_name, 'BAD I_DIM ARGUMENT: \i4\ ', i_dim)
  call end_cleanup(branch, .false.)
  return
end select

! Orbit correction = (T-1)^-1 * (orbit_end - orbit_start)
!                  = t11_inv  * (orbit_end - orbit_start)

!--------------------------------------------------------------------------
! Because of nonlinearities we may need to iterate to find the solution

allocate(co_best(0:ubound(closed_orb, 1)))
allocate(vec0(n_dim), dvec(n_dim), weight(n_dim), coc%a_vec(n_dim), maska(n_dim), merit_vec(n_dim), dy_da(n_dim,n_dim))
vec0 = 0
maska = .false.
has_been_alive_at_end = .false.
check_for_z_stability = (n_dim == 6)
invert_step_tried = .false.
a_lambda = -1
old_chisq = 1d30   ! Something large
coc%a_vec = start_orb%vec(1:n_dim)
if (n_dim == 5) coc%a_vec(5) = start_orb%vec(6)
try_lmdif = .true.
co_best = closed_orb
c_best => co_best

!

i_max = 100

do i_loop = 1, i_max

  weight(1:n_dim) = 1 / (abs(coc%a_vec) * bmad_com%rel_tol_tracking + bmad_com%abs_tol_tracking)**2
  weight(1:4) = weight(1:4) * [1/sqrt(betas(1)), sqrt(betas(1)), 1/sqrt(betas(2)), sqrt(betas(2))]

  call super_mrqmin (vec0, weight, coc%a_vec, chisq, co_func, storage, a_lambda, status, print_err = .false.)
  ! If super_mrqmin rejects a trial step then will need to reset closed_orb.
  if (chisq < old_chisq .and. status /= 1 .and. end_orb%state /= not_set$) then
    co_best = closed_orb
  endif

  if (any(coc%a_vec /= start_orb%vec(1:n_dim)) .and. all(coc%a_vec == co_best(ix_ele_start)%vec(1:n_dim))) then
    closed_orb = co_best
    branch%param%unstable_factor = 0
  endif

  !

  if (.not. has_been_alive_at_end) then
    select case (i_loop)
    case (1)
      n = min(5, n_dim)
      coc%a_vec(1:n) = 0  ! Desperation move.
      a_lambda = -1       ! Reset super_mrqmin
      cycle
    case (2)
      n = min(5, n_dim)
      this = start_orb%vec(6) * [ele_start%x%eta, ele_start%x%etap, ele_start%y%eta, ele_start%y%etap, 0.0_rp]
      coc%a_vec(1:n) = this(1:n)  ! Another desperation move.
      a_lambda = -1  ! Reset super_mrqmin
      cycle
    case default
      ! Thought: Put in code using lmdif and minimizing branch%param%unstable_factor to see if a stable orbit
      ! can be found.
      if (printit) call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING!!', 'ABORTING CLOSED ORBIT SEARCH.', &
                                                  'TRACKING BRANCH: ' // branch_name(branch))
      call end_cleanup(branch, .true.)
      return
    end select
  endif

  !

  if ((branch%param%unstable_factor /= 0 .and. try_lmdif) .or. a_lambda > 1d8) then  ! Particle lost. Try lmdif since lmdif does not need derivatives
    call initial_lmdif()
    do i2 = 1, i_max/4
      call co_func(coc%a_vec, dvec, dy_da, status)
      do j = 1, n_dim
        merit_vec(j) = sqrt(weight(j)) * dvec(j)
      enddo
      call suggest_lmdif (coc%a_vec, merit_vec, bmad_com%rel_tol_tracking, i_max, at_end)
    enddo

    call co_func(coc%a_vec, dvec, dy_da, status)
    if (status /= 0) then
      if (printit) call out_io (s_error$, r_name, 'CANNOT FIND STABLE ORBIT.', 'TRACKING BRANCH: ' // branch_name(branch))
      call end_cleanup(branch, .true.)
      return
    endif

    call this_t1_calc (branch, dir, .true., coc%t1, betas, start_orb_t1, err)
    if (err) then
      call end_cleanup(branch, .true.)
      return
    endif

    check_for_z_stability = (n_dim == 6)  ! New t1 matrix needs to be checked for stability.
    a_lambda = -1    ! Reset super_mrqmin
    try_lmdif = .false.
    cycle
  endif

  if (status == 1) then ! too large step in z
    a_lambda = a_lambda * coc%dz_step  ! Scale next step 
  endif

  if (a_lambda < 1d-10) a_lambda = 1d-10

  if (status < 0) then  
    if (printit) call out_io (s_error$, r_name, 'SINGULAR MATRIX ENCOUNTERED!', 'ABORTING CLOSED ORBIT SEARCH.', &
                                                'TRACKING BRANCH: ' // branch_name(branch))
    call end_cleanup(branch, .true.)
    return
  endif

  !

  amp_co = abs(start_orb%vec)
  dorb = end_orb%vec - start_orb%vec

  if (n_dim == 6) then
    coc%a_vec(5) = modulo2 (coc%a_vec(5), coc%rf_wavelen / 2)   ! Keep z within 1/2 wavelength of original value

    if (bmad_com%absolute_time_tracking) then
      dt = (end_orb%t - start_orb%t) - nint((end_orb%t - start_orb%t) * rf_freq) / rf_freq
      dorb(5) = -end_orb%beta * c_light * dt
    endif
  endif

  amp_del = abs(dorb)
  if (all(amp_del(1:n_dim) < amp_co(1:n_dim) * bmad_com%rel_tol_tracking + bmad_com%abs_tol_tracking)) then
    if (n_dim == 6) then
      call this_t1_calc (branch, dir, .true., coc%t1, betas, start_orb_t1, err)
      if (err) then
        call end_cleanup(branch, .true.)
        return
      endif

      if (max_z_eigen(coc%t1) < 1.000001_rp) exit
      check_for_z_stability = .true.   ! Is unstable
    else
      exit
    endif
  endif

  ! If not converging fast enough, remake the transfer matrix.
  ! This is computationally intensive so only do this if the orbit has shifted significantly or the
  ! 1-turn matrix inverting step was a bust.
  ! Status == 1 means that the change in "a" was too large but coc%t1 does not need to be remade.

  if (chisq > old_chisq/2 .and. status == 0 .and. branch%param%unstable_factor == 0) then ! If not converging
    if (maxval(abs(start_orb_t1(1:n_dim)-co_best(0)%vec(1:n_dim))) > 1d-6 .or. invert_step_tried) then 
      invert_step_tried = .false.
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

      call this_t1_calc (branch, dir, .true., coc%t1, betas, start_orb_t1, err)
      if (err) then
        call end_cleanup(branch, .true.)
        return
      endif
      check_for_z_stability = (n_dim == 6)  ! New t1 matrix needs to be checked for stability.

      do n = 1, branch%n_ele_max
        branch%ele(n)%mat6 = m(n)%mat6
        branch%ele(n)%vec0 = m(n)%vec0
        branch%ele(n)%map_ref_orb_in  = m(n)%map_ref_orb_in
        branch%ele(n)%map_ref_orb_out = m(n)%map_ref_orb_out
      enddo

    ! Else try a step by inverting the 1-turn matrix.
    else
      invert_step_tried = .true.
      call make_t11_inv (err)
      if (err) then
        call end_cleanup(branch, .true.)
        return
      endif

      del_co(1:n_dim) = matmul(t11_inv(1:n_dim,1:n_dim), dorb(1:n_dim)) 

      ! For n_dim = 6, if at peak of RF then del_co(5) may be singularly large. 
      ! To avoid this, limit z step to be no more than lambda_rf/10

      if (n_dim == 5) then
        del_co(5) = dorb(5) / i1_int      

      elseif (n_dim == 6) then
        dz_norm = abs(del_co(5)) / (coc%rf_wavelen / 10)
        if (dz_norm > 1) then
          del_co = del_co / dz_norm
        endif
      endif

      if (all(abs(del_co(1:n_dim)) < 1e-2)) then
        call co_func (coc%a_vec + del_co(1:n_dim), dvec, dy_da, status)
        if (track_state == moving_forward$) coc%a_vec = coc%a_vec + del_co(1:n_dim)
      endif
    endif
  endif

  ! The super_mrqmin routine  will converge to the nearest fixed point even if that fixed point is an unstable one.
  ! For n_dim = 6, there are longitudinally unstable fixed points which we want to avoid.
  ! If we are near an unstable fixed point look for a better spot by shifting the particle in z in steps of pi/4.
  ! Note: due to inaccuracies, the maximum eigen value may be slightly over 1 at the stable fixed point.

  if (check_for_z_stability .and. has_been_alive_at_end) then
    min_max_eigen = max_z_eigen(coc%t1)
    if (min_max_eigen > 1.000001_rp) then    ! Is unstable
      j_max = 0
      z0 = coc%a_vec(5)
      dz = coc%rf_wavelen / 8
      do j = -3, 4
        z_here = z0 + j * dz
        call max_eigen_calc(branch, z_here, max_eigen, start_orb_t1, betas)
        if (max_eigen > min_max_eigen) cycle
        j_max = j
        min_max_eigen = max_eigen
      enddo

      ! Reset
      coc%a_vec(5) = z0 + j_max * dz
      call max_eigen_calc(branch, coc%a_vec(5), max_eigen, start_orb_t1, betas)

      ! If z needs to be shifted, reset super_mrqmin
      if (j_max /= 0) then
        a_lambda = -1               ! Signal super_mrqmin reset
        try_lmdif = .true.
      endif
    endif
    check_for_z_stability = .false.
  endif

  !

  old_chisq = chisq

  if (i_loop == i_max) then
    if (maxval(amp_del(1:n_dim)) < 1d-4) then
      if (printit) call out_io (s_error$, r_name, &
              'Closed orbit not converging! error in closed orbit: \es10.2\ ', &
              'Using branch: ' // branch_name(branch), &
              'If this error is acceptable, change bmad_com%rel_tol_tracking (\es10.2\) and/or', &
              'bmad_com%abs_tol_tracking (\es10.2\)', &
              r_array = [maxval(amp_del(1:n_dim)), bmad_com%rel_tol_tracking, bmad_com%abs_tol_tracking])
    else
      if (printit) call out_io (s_error$, r_name, &
              'Closed orbit not converging! error in closed orbit: \es10.2\ ', &
              'Using branch: ' // branch_name(branch), &
              r_array = [maxval(amp_del(1:n_dim))])
    endif
    call end_cleanup(branch)
    return
  endif

enddo

if (i_loop == i_max + 1) then
  if (printit) call out_io (s_error$, r_name, 'CLOSED ORBIT CANNOT BE FOUND FOR BRANCH: ' // branch_name(branch))
  call end_cleanup(branch)
  return
endif

! Calc invarient spin axis.

if (bmad_com_saved%spin_tracking_on) then
  bmad_com%spin_tracking_on = .true.
  do i = 1, 3
    start_orb%spin = 0
    start_orb%spin(i) = 1
    call track_many (lat, closed_orb, ix_ele_start, ix_ele_end, dir, branch%ix_branch, track_state)
    mat3(:,i) = end_orb%spin
  enddo
  call w_mat_to_axis_angle(mat3, start_orb%spin, branch%param%spin_tune)
  if (bmad_com%spin_n0_direction_user_set) then
    if (dot_product(start_orb%spin, lat%particle_start%spin) < 0) then
      start_orb%spin = -start_orb%spin
      branch%param%spin_tune = 1.0_rp - branch%param%spin_tune
    endif
  endif
  call track_many (lat, closed_orb, ix_ele_start, ix_ele_end, dir, branch%ix_branch, track_state)
endif

! Cleanup

call end_cleanup(branch)
if (present(err_flag)) err_flag = .false.

!------------------------------------------------------------------------------
! return rf cavities to original state

contains

subroutine end_cleanup(branch, reset_orb)

type (branch_struct) branch
logical, optional :: reset_orb

!

bmad_com = bmad_com_saved  ! Restore
bmad_private%random_on = .true.

if (n_dim == 4 .or. n_dim == 5) then
  call set_on_off (rfcavity$, branch%lat, restore_state$, ix_branch = branch%ix_branch, saved_values = on_off_state)
endif

call set_on_off (foil$, lat, restore_state$, ix_branch = branch%ix_branch, saved_values = scatter_state, attribute = 'SCATTER_METHOD')

if (logic_option(.false., reset_orb)) closed_orb(0) = original_start_orb

end subroutine

!------------------------------------------------------------------------------
! contains

subroutine max_eigen_calc(branch, z_set, max_eigen, start_orb_t1, betas)

type (branch_struct) branch
real(rp) z_set, max_eigen, start_orb_t1(6), betas(2)
logical err

!

start_orb%vec(5) = z_set
call init_coord (start_orb, start_orb, ele_start, start_end$, start_orb%species, dir)

call track_many (lat, closed_orb, ix_ele_start, ix_ele_end, dir, branch%ix_branch, track_state)

if (track_state /= moving_forward$) then
  max_eigen = 10
  return
endif

call this_t1_calc (branch, dir, .true., coc%t1, betas, start_orb_t1, err, .false.); if (err) return
max_eigen = max_z_eigen(coc%t1)

end subroutine max_eigen_calc

!------------------------------------------------------------------------------
! contains

function max_z_eigen(t1) result (z_eigen)

real(rp) t1(6,6), z_eigen
complex(rp) eigen_val(6), eigen_vec(6,6)
logical error

!

call mat_eigen (t1, eigen_val, eigen_vec, error)
z_eigen = max(abs(eigen_val(5)), abs(eigen_val(6)))

end function max_z_eigen

!------------------------------------------------------------------------------
! contains

subroutine co_func (a_try, y_fit, dy_da, status)

type (branch_struct), pointer :: branch

real(rp), intent(in) :: a_try(:)
real(rp), intent(out) :: y_fit(:)
real(rp), intent(out) :: dy_da(:, :)
real(rp) del_orb(6)

integer status, i

! For i_dim = 6, if at peak of RF then delta z may be singularly large. 
! To avoid this, veto any step where z changes by more than lambda_rf/10.

branch => lat%branch(integer_option(0, ix_branch))  ! To get around ifort compiler bug

if (n_dim == 6) then
  coc%dz_step = abs(a_try(5)-coc%a_vec(5)) / (coc%rf_wavelen / 10)
  if (coc%dz_step > 1) then
    status = 1  ! Veto step
    return
  endif
endif

!

if (n_dim == 5) then
  start_orb%vec(1:4) = a_try(1:4)
  start_orb%vec(6)   = a_try(5)
else
  start_orb%vec(1:n_dim) = a_try
endif

call init_coord (start_orb, start_orb, ele_start, start_end$, start_orb%species, dir)
call track_many (lat, closed_orb, ix_ele_start, ix_ele_end, dir, branch%ix_branch, track_state)

!

status = 0
del_orb = end_orb%vec - start_orb%vec

if (n_dim == 6 .and. bmad_com%absolute_time_tracking) then
  dt = (end_orb%t - start_orb%t) - nint((end_orb%t - start_orb%t) * rf_freq) / rf_freq
  del_orb(5) = -end_orb%beta * c_light * dt
endif

if (track_state == moving_forward$) then
  if (.not. has_been_alive_at_end) co_best = closed_orb
  has_been_alive_at_end = .true.
  y_fit = del_orb(1:n_dim)
else
  y_fit = 1d10 * branch%param%unstable_factor   ! Some large number
endif

dy_da = coc%t1(1:n_dim,1:n_dim)
forall (i = 1:n_dim) dy_da(i,i) = dy_da(i,i) - 1

if (n_dim == 5) then
  dy_da(:,5) = coc%t1(1:5,6)
endif

end subroutine co_func

!------------------------------------------------------------------------------
! contains

subroutine this_t1_calc (branch, dir, make_mat6, t1, betas, start_orb_t1, err, cleanup_on_err)

type (branch_struct) branch
type (coord_struct) start_saved, end_saved
type (ele_struct) ele0

real(rp) t1(6,6), delta, betas(2), growth_rate, start_orb_t1(6)
integer dir, i, track_state, stat
logical make_mat6, err, stable
logical, optional :: cleanup_on_err

!

err = .false.

if (dir == 1) then
  if (make_mat6) call lat_make_mat6 (branch%lat, -1, closed_orb, branch%ix_branch)
  call transfer_matrix_calc (branch%lat, t1, ix_branch = branch%ix_branch)
  if (all(t1 == 0)) then ! Something is really wrong so try a hard reset of t1 around the zero orbit.
    call lat_make_mat6 (branch%lat, -1, ix_branch = branch%ix_branch, err_flag = err)
    if (err) then
      call out_io (s_error$, r_name, 'PARTICLE LOST. NO CLOSED ORBIT CALCULATED.', &
                                     'USING BRANCH: ' // branch_name(branch))
      if (logic_option(.true., cleanup_on_err)) call end_cleanup(branch, .true.)
      return
    endif
    call transfer_matrix_calc (branch%lat, t1, ix_branch = branch%ix_branch)
  endif

else
  call track_many (branch%lat, closed_orb, branch%n_ele_track, 0, dir, branch%ix_branch, track_state)
  start_saved = start_orb
  end_saved   = end_orb

  delta = 1d-6
  do i = 1, 6
    start_orb%vec(i) = start_saved%vec(i) + delta
    call track_many (branch%lat, closed_orb, branch%n_ele_track, 0, dir, branch%ix_branch, track_state)
    if (track_state /= moving_forward$) then 
      call out_io (s_error$, r_name, 'PARTICLE LOST TRACKING BACKWARDS. [POSSIBLE CAUSE: WRONG PARTICLE SPECIES.]', &
                                     'NO CLOSED ORBIT CALCULATED.', &
                                     'USING BRANCH: ' // branch_name(branch))
      if (logic_option(.true., cleanup_on_err)) call end_cleanup(branch, .true.)
      err = .true.
      return
    endif 
    t1(:,i) = (end_orb%vec - end_saved%vec) / delta
    start_orb%vec = start_saved%vec
  enddo
endif

start_orb_t1 = start_orb%vec

!

call twiss_from_mat6 (t1, start_orb%vec, ele0, stable, growth_rate, stat, .false.)
if (stat == ok$) then
  betas = [ele0%a%beta, ele0%b%beta]
else
  betas = 1
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
mat(1:n_dim,1:n_dim) = mat(1:n_dim,1:n_dim) - coc%t1(1:n_dim,1:n_dim)
call mat_inverse(mat(1:n_dim,1:n_dim), t11_inv(1:n_dim,1:n_dim), ok = ok2)

if (.not. ok1 .or. .not. ok2) then 
  if (printit) call out_io (s_error$, r_name, 'MATRIX INVERSION FAILED!')
  call end_cleanup(branch, .true.)
  return
endif

err = .false.

end subroutine

end subroutine closed_orbit_calc
