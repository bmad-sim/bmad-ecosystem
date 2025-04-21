module space_charge_mod

use csr_and_space_charge_mod

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine sc_field_calc (bunch, branch, include_image, t_end, sc_field)
!
! Routine to calculate the space charge field of a bunch
!
! Input:
!   bunch           -- bunch_struct: Starting bunch position in time-based coordinates
!   branch          -- branch_struct: Lattice branch being tracked through.
!   include_image   -- logical: True if cathode image charge fields are to be included.
!   t_end           -- real(rp): Calculate space charge field for preborn particles emitted before t_end.
!   bunch_params    -- bunch_params_struct, optional: If present, particles too far from the bunch will be ignored.
!                        The threashold is set by space_charge_com%particle_sigma_cutoff. 
! Output:
!   include_image   -- logical: Set False if image charge calc no longer needed (Note: never set True).
!   sc_field        -- em_field_struct: Space charge field at particle positions
!-

subroutine sc_field_calc (bunch, branch, include_image, t_end, sc_field, bunch_params)

implicit none

type this_w_struct
  real(rp) w_mat(3,3)
  integer ixp
end type

type (bunch_struct), target :: bunch
type (branch_struct), target :: branch
type (coord_struct), pointer :: p
type (csr_particle_position_struct) :: position(size(bunch%particle))
type (this_w_struct) :: w(size(bunch%particle))
type (em_field_struct) :: sc_field(size(bunch%particle))
type (mesh3d_struct) :: mesh3d, mesh3d_image
type (bunch_params_struct), optional :: bunch_params
type (floor_position_struct) pos, pos0

integer :: n, n_alive, i, imin(1), k
real(rp) :: beta, ratio, t_end, s_ave, w_mat_inv(3,3)
real(rp) :: Evec(3), Bvec(3), Evec_image(3), sigma(3)
logical :: include_image, err, bend_here

! Initialize variables
mesh3d%nhi = space_charge_com%space_charge_mesh_size
if (present(bunch_params)) then
  sigma(1) = sqrt(bunch_params%sigma(1,1))
  sigma(2) = sqrt(bunch_params%sigma(3,3))
  sigma(3) = sqrt(bunch_params%sigma(5,5))
endif

! Gather particles
n = 0
n_alive = 0
beta = 0
bend_here = .false.

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state == pre_born$ .and. p%t <= t_end) then ! Particles to be emitted
    n = n + 1
    position(n)%r = [p%vec(1), p%vec(3), p%s]
    position(n)%charge = 0
  else if (p%state /= alive$) then  ! Lost particles
    cycle
  else  ! Living particles
    if (present(bunch_params)) then
      if (out_of_sigma_cutoff(p)) cycle  ! Ignore particles too far off the bunch
    endif
    n = n + 1
    position(n)%r = [p%vec(1), p%vec(3), p%s]
    position(n)%charge = p%charge * charge_of(p%species)
    beta = beta + p%beta
    n_alive = n_alive + 1
    if (branch%ele(p%ix_ele)%key == sbend$) bend_here = .true.
    w(i)%ixp = n
  endif
enddo

! Return if not enough living particles
if (n_alive<2) return
beta = beta/n_alive

! If there is a bend then need to transform from curvilinear to Carteasion coords

if (bend_here .and. .false.) then
  s_ave = sum(position%r(3)*position%charge) / sum(position%charge)
  pos0 = coords_curvilinear_to_floor([0.0_rp, 0.0_rp, s_ave], branch, err)
  w_mat_inv = transpose(pos0%w)

  do i = 1, n
    pos = coords_curvilinear_to_floor(position(i)%r, branch, err)
    position(i)%r = matmul(w_mat_inv, pos%r - pos0%r)
    w(i)%w_mat = matmul(transpose(pos%w), pos0%w)
  enddo
endif

! Calculate space charge field
mesh3d%gamma = 1/sqrt(1- beta**2)
call deposit_particles (position(1:n)%r(1), position(1:n)%r(2), position(1:n)%r(3), mesh3d, qa=position(1:n)%charge)
call space_charge_3d(mesh3d, at_cathode=include_image, calc_bfield=.true., image_efield=mesh3d_image%efield)

! Determine if cathode image field should be turned off
if (include_image) then
  ! Copy mesh3d dimensions and allocate image_efield
  mesh3d_image%nlo = mesh3d%nlo
  mesh3d_image%nhi = mesh3d%nhi
  mesh3d_image%min = mesh3d%min
  mesh3d_image%max = mesh3d%max
  mesh3d_image%delta = mesh3d%delta

  ! Compare efield and image_efield on the last particle
  imin = minloc(bunch%particle%s,mask=(bunch%particle%state == alive$))
  p => bunch%particle(imin(1))
  call interpolate_field(p%vec(1), p%vec(3), p%s,  mesh3d, E=Evec)
  call interpolate_field(p%vec(1), p%vec(3), p%s,  mesh3d_image, E=Evec_image)
  ratio = maxval(abs(Evec_image/Evec))
  ! If image field is small compared to the bunch field, turn it off from here on
  if (ratio <= space_charge_com%cathode_strength_cutoff) include_image = .false.
endif

! Calculate field at particle locations

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state == pre_born$ .and. p%t <= t_end) then  ! Particles to be emitted
    call interpolate_field(p%vec(1), p%vec(3), p%s,  mesh3d, E=sc_field(i)%E, B=sc_field(i)%B)
  else if (p%state /= alive$) then ! Lost particles
    sc_field(i)%E = 0
    sc_field(i)%B = 0
  else  ! Living particles
    if (present(bunch_params)) then
      if (out_of_sigma_cutoff(p)) then  ! Ignore particles too far off the bunch
        sc_field(i)%E = 0
        sc_field(i)%B = 0
      endif
    endif
    call interpolate_field(p%vec(1), p%vec(3), p%s,  mesh3d, E=sc_field(i)%E, B=sc_field(i)%B)
    if (bend_here .and. .false.) then
      k = w(i)%ixp
      sc_field(i)%E = matmul(w(k)%w_mat, sc_field(i)%E)
      sc_field(i)%B = matmul(w(k)%w_mat, sc_field(i)%B)
    endif
  end if
enddo

!-------------------------------------------------------------------
contains
! Check if a particle is too far from the bunch
!
function out_of_sigma_cutoff(p) result (out_of_cutoff)

type (coord_struct) :: p
logical :: out_of_cutoff
real(rp) :: sig_cut

out_of_cutoff = .false.
sig_cut = space_charge_com%particle_sigma_cutoff
if (sig_cut .le. 0) return
out_of_cutoff = (abs(p%vec(1) - bunch_params%centroid%vec(1)) > sigma(1) * sig_cut .or. &
                 abs(p%vec(3) - bunch_params%centroid%vec(3)) > sigma(2) * sig_cut .or. &
                 abs(p%s - bunch_params%centroid%vec(5)) > sigma(3) * sig_cut)

end function out_of_sigma_cutoff
end subroutine sc_field_calc

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine sc_step(bunch, ele, include_image, t_end, n_emit)
!
! Subroutine to track a bunch through a given time step with space charge
!
! Input:
!   bunch         -- bunch_struct: Starting bunch position in t-based coordinates
!   ele           -- ele_struct: Nominal element being tracked through.
!   include_image -- logical: Include image charge forces?
!   t_end         -- real(rp): Time at which the tracking ends.
!   sc_field      -- em_field_struct(:): Array to hold space charge fields. 
!                       Its length should be the number of particles.
!
! Output:
!   bunch         -- bunch_struct: Ending bunch position in t-based coordinates after space charge kick.
!   include_image -- logical: Set False if image charge calc no longer needed (Note: never set True).
!   n_emit        -- integer, optional: The number of particles emitted in this step.
!-

subroutine sc_step(bunch, ele, include_image, t_end, sc_field, n_emit)

implicit none

type (bunch_struct), target :: bunch
type (ele_struct) :: ele
type (em_field_struct) :: extra_field(size(bunch%particle))
type (coord_struct), pointer :: p
type (em_field_struct) :: sc_field(:)
type (bunch_params_struct) :: bunch_params

real(rp) t_end, sum_z
logical include_image, error
integer i, n
integer, optional :: n_emit

! Calculate space charge field
if (space_charge_com%particle_sigma_cutoff > 0) then
  call calc_bunch_params(bunch, bunch_params, error, .false., is_time_coords=.true., ele=ele)
  if (error) then
    call sc_field_calc(bunch, ele%branch, include_image, t_end, sc_field)
  else
    call sc_field_calc(bunch, ele%branch, include_image, t_end, sc_field, bunch_params)
  endif
else
  call sc_field_calc(bunch, ele%branch, include_image, t_end, sc_field)
endif

! Generate particles at the cathode
do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state == pre_born$ .and. p%t <= t_end) then
    p%state = alive$
    if (present(n_emit)) n_emit = n_emit + 1
  endif
enddo 

! Update bunch info
bunch%charge_live = sum(bunch%particle%charge, (bunch%particle%state == alive$))
bunch%n_live = count((bunch%particle%state == alive$))

! And track
call track_bunch_time(bunch, ele%branch, t_end, 1e30_rp, extra_field=sc_field)
end subroutine sc_step

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine sc_adaptive_step(bunch, ele, include_image, t_now, dt_step, dt_next)
!
! Routine to track a bunch of particles with space charge for one step using
! adaptive step size control and determine appropriate step size for the next step
!
! Input:
!   bunch         -- bunch_struct: Starting bunch position in t-based coordinates
!   ele           -- ele_struct: Nominal lattice element being tracked through.
!   include_image -- logical: Include image charge forces?
!   t_now         -- real(rp): Current time at the beginning of tracking
!   dt_step       -- real(rp): Initial SC time step to take
!   sc_field      -- em_field_struct(:): Array to hold space charge fields. 
!                       Its length should be the number of particles.
!
! Output:
!   bunch         -- bunch_struct: Ending bunch position in t-based coordinates.
!   include_image -- logical: Set False if image charge calc no longer needed (Note: never set True).
!   dt_next       -- real(rp): Next SC time step the tracker would take based on the error tolerance
!   dt_step       -- real(rp): Step done.
!-

subroutine sc_adaptive_step(bunch, ele, include_image, t_now, dt_step, dt_next, sc_field)

implicit none

type (bunch_struct) :: bunch, bunch_full, bunch_half
type (ele_struct) ele
type (coord_struct), pointer :: p
type (em_field_struct) :: sc_field(:)

real(rp) :: t_now, dt_step, dt_next, sqrt_N
real(rp) :: r_err(6), r_scale(6), rel_tol, abs_tol, err_max
real(rp), parameter :: tiny = 1.0e-30_rp

integer i, N, n_emit(2)
logical include_image

!

sqrt_N = sqrt(abs(1/(c_light*dt_step)))  ! number of steps we would take to cover 1 meter
rel_tol = space_charge_com%rel_tol_tracking / sqrt_N
abs_tol = space_charge_com%abs_tol_tracking / sqrt_N

bunch_full = bunch
bunch_half = bunch
dt_next = dt_step

do
  n_emit = [0,0]
  ! Full step
  call sc_step(bunch_full, ele, include_image, t_now+dt_step, sc_field)
  ! Two half steps
  call sc_step(bunch_half, ele, include_image, t_now+dt_step/2, sc_field, n_emit=n_emit(1))
  call sc_step(bunch_half, ele, include_image, t_now+dt_step, sc_field, n_emit = n_emit(2))

  r_scale = abs(bunch_rms_vec(bunch))+abs(bunch_rms_vec(bunch_half)) + tiny
  r_err = [0,0,0,0,0,0]
  N = 0
  ! Calculate error from the difference
  do i = 1, size(bunch%particle)
    ! If going backwards then particle is considered dead.
    if (bunch_half%particle(i)%direction == -1 .or. &
                      bunch_full%particle(i)%direction == -1) bunch_half%particle(i)%state = lost_z$
    if (bunch_half%particle(i)%state /= alive$) cycle ! Only count living particles
    r_err(:) = r_err(:) + abs(bunch_full%particle(i)%vec(:)-bunch_half%particle(i)%vec(:))
    N = N +1
  enddo
  ! If no living particle, finish step
  if (N==0) then
    bunch = bunch_half
    return
  end if
  
  ! Compare error to the tolerance
  r_err = r_err/N
  
  err_max = maxval(r_err(:)/(r_scale*rel_tol + abs_tol))
  if (space_charge_com%debug) print *, dt_step, err_max, n_emit
  
  ! If error is larger than tolerance, try again with a smaller step
  if (err_max <= 1.0) exit
  dt_step = next_step_size(dt_step, err_max)
  bunch%n_bad = bunch%n_bad + 1
  bunch_full = bunch
  bunch_half = bunch
enddo

bunch%n_good = bunch%n_good + 1
! Adjust next step size
dt_next = next_step_size(dt_step, err_max)

bunch = bunch_half

!-------------------------------------------------------------------
contains

! Use RMS of coordinates to estimate the scale of motion
! It combine the centroid position and the width of the bunch

function bunch_rms_vec(bunch) result (rms_vec)

type (bunch_struct) bunch
real(rp) rms_vec(6)
integer i, N

!

N = 0
rms_vec = [0,0,0,0,0,0]
do i = 1,size(bunch%particle)
  if (bunch%particle(i)%state /= alive$) cycle
  rms_vec(:) = rms_vec(:) + bunch%particle(i)%vec(:)**2
  N = N +1
enddo
if (N==0) return
rms_vec = sqrt(rms_vec/N)

end function bunch_rms_vec

! Determines next time step size
function next_step_size(dt_step, err_max) result (dt_next)

real(rp) :: dt_step, err_max, dt_next, t_end
real(rp), parameter :: safety = 0.9_rp, p_grow = -0.5_rp
real(rp), parameter :: p_shrink = -0.5_rp, err_con = 1.89d-4
integer ix_next

! If early during emission, emit single particle


! shrink or grow step size based on err_max
if (err_max > 1) then
  dt_next = safety * dt_step * (err_max**p_shrink)
else if (err_max > err_con) then
  dt_next = safety * dt_step * (err_max**p_grow)
else
  dt_next = 5.0_rp * dt_step
endif

end function next_step_size

end subroutine sc_adaptive_step

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine track_bunch_to_s (bunch, s, branch)
!
! Drift a bunch of particles to the same s coordinate
!
! Input:
!   bunch     -- bunch_struct: Input bunch position in s-based coordinate.
!   s         -- real(rp): Target coordinate.
!   branch    -- branch_struct: Branch being tracked through.
!
! Output:
!   bunch     -- bunch_struct: Output bunch position in s-based coordinate. Particles will be at the same s coordinate
!-

subroutine track_bunch_to_s (bunch, s, branch)

implicit none

type (bunch_struct), target :: bunch
type (coord_struct), pointer :: p
type (branch_struct) :: branch
type (coord_struct) :: position

integer i, track_state
real(rp) s, ds, s0, s_end, s_begin

!

if (bunch%drift_between_t_and_s) then
  do i = 1, size(bunch%particle)
    p => bunch%particle(i)
    if (p%state /= alive$) cycle
    call drift_particle_to_s(p, s, branch)
  enddo
  return
endif

!

s_end = branch%ele(branch%n_ele_track)%s
s_begin = branch%ele(0)%s

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state /= alive$) cycle

  ! Can happen that particle has been "drifted" past the end of the branch.
  ! Track_from_s_to_s cannot handle such a case so use drift_particle_to_s instead.
  if (p%s > s_end) then
    call drift_particle_to_s (p, s_end, branch)
  elseif (p%s < s_begin) then
    call drift_particle_to_s (p, s_begin, branch)
  endif

  ds = s - p%s
  if (ds == 0) cycle
  p%time_dir = sign_of(ds) * p%direction
  s0 = p%s
  call track_from_s_to_s (branch%lat, s0, s, p, p, ix_branch = branch%ix_branch, track_state = track_state)
  if (track_state /= moving_forward$) p%state = lost$
  p%time_dir = 1
enddo

end subroutine track_bunch_to_s

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine track_bunch_to_t (bunch, t_target, branch)
!
! Drift a bunch of particles to the same t coordinate
!
! Input:
!   bunch     -- bunch_struct: Input bunch position in s-based coordinate.
!   t_target  -- real(rp): Target t coordinate.
!   branch    -- branch_struct: Lattice branch being tracked through.
!
! Output:
!   bunch     -- bunch_struct: Output bunch position in s-based coordinate. Particles will be at the same t coordinate
!-

subroutine track_bunch_to_t (bunch, t_target, branch)

use super_recipes_mod, only: super_zbrent

implicit none

type (bunch_struct), target :: bunch
type (coord_struct), pointer :: p0, p_ptr
type (coord_struct), target :: p
type (branch_struct) :: branch
type (coord_struct) :: position

integer i, status
real(rp) ds, t_target, dt, dt2, s1, s_end, s_target, s_begin

!

if (bunch%drift_between_t_and_s) then
  do i = 1, size(bunch%particle)
    p0 => bunch%particle(i)
    if (p0%state /= alive$) cycle
    call drift_particle_to_t(p0, t_target, branch)
  enddo
  return
endif

! Convert bunch to s-based coordinates

s_begin = branch%ele(0)%s
s_end = branch%ele(branch%n_ele_track)%s
p_ptr => p    ! To get around ifort debug info bug

particle_loop: do i = 1, size(bunch%particle)
  p0 => bunch%particle(i)
  if (p0%state /= alive$) cycle
  p = p0

  ! bracket solution
  do
    dt = t_target - p0%t
    ds = 1.01 * dt * c_light * p0%beta * p0%direction
    s_target = ds + p0%s
    if (abs(ds) < 0.1_rp * bmad_com%significant_length) then
      bunch%particle(i)%time_dir = 1
      cycle particle_loop
    endif
    dt2 = track_func(s_target, status)
    if (status /= 0) then
      p0%state = lost$
      cycle particle_loop
    endif
    if (dt*dt2 <= 0) exit
    p0 = p
  enddo

  ! Use zbrent to find bracketed solution
  s1 = p%s
  s1 = super_zbrent(track_func, p0%s, s1, 1e-10_rp, bmad_com%significant_length, status)
  if (status /= 0) then
    p0%state = lost$
    cycle particle_loop
  endif
  dt = track_func(s1, status)
  p%time_dir = 1
  bunch%particle(i) = p
enddo particle_loop

!--------------------------------------------------------------------
contains

function track_func (s_target, status) result (dt)

real(rp), intent(in) :: s_target
integer status, track_state
real(rp) :: dt

!

p0%time_dir = sign_of(s_target - p0%s)
p%time_dir  = sign_of(s_target - p0%s)
status = 0

! Can happen that particle needs to "drift" past the end of the branch.
! Track_from_s_to_s cannot handle such a case so use drift_particle_to_t instead.

if (s_target < s_begin .and. p0%s >= s_begin) then
  call track_from_s_to_s (branch%lat, p0%s, s_begin, p0, p, ix_branch = p%ix_branch, track_state = track_state)
  call drift_particle_to_s(p, s_target, branch)
elseif (s_target > s_end .and. p0%s <= s_end) then
  call track_from_s_to_s (branch%lat, p0%s, s_end, p0, p, ix_branch = p%ix_branch, track_state = track_state)
  call drift_particle_to_s(p, s_target, branch)
else
  call track_from_s_to_s (branch%lat, p0%s, s_target, p0, p, ix_branch = p%ix_branch, track_state = track_state)
endif

if (track_state /= moving_forward$) then
  status = 1
  dt = real_garbage$
  return
endif

dt = t_target - p%t 

end function track_func

end subroutine track_bunch_to_t
  
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine drift_particle_to_s (p, s, branch)
!
! Drift a particle to a given s-coordinate
!
! Input:
!   p         -- coord_struct: Init particle position.
!   s         -- real(rp): Target s coordinate.
!   branch    -- branch_struct: Branch being tracked through.
!
! Output:
!   p         -- coord_struct: Final particle position.
!-

subroutine drift_particle_to_s (p, s, branch)

implicit none

type (coord_struct) :: p
type (branch_struct) :: branch
type (coord_struct) :: position

real(rp) s, ds

!

if (p%state /= alive$) return
ds = s - p%s

call track_a_drift(p, ds)

if (p%s > branch%ele(branch%n_ele_track)%s) then
  p%ix_ele = branch%n_ele_track
  p%location = downstream_end$
elseif (p%s < branch%ele(0)%s) then
  p%ix_ele = 1
  p%location = upstream_end$
else
  p%ix_ele = element_at_s(branch, p%s, (ds < 0), position=position)
  p%location = position%location
endif

end subroutine drift_particle_to_s

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine drift_particle_to_t (p, t, branch)
!
! Drift a particle to a given t-coordinate
!
! Input:
!   p         -- coord_struct: Init particle position.
!   t         -- real(rp): Target t coordinate.
!   branch    -- branch_struct: Lattice branch being tracked through.
!
! Output:
!   p         -- coord_struct: Final particle position.
!-

subroutine drift_particle_to_t (p, t, branch)

implicit none

type (coord_struct) :: p
type (branch_struct) :: branch
type (coord_struct) :: position

real(rp) t, pz0, E_tot, dt, ds

!

if (p%state /= alive$) return

pz0 = sqrt((1.0_rp + p%vec(6))**2 - p%vec(2)**2 - p%vec(4)**2 ) ! * p0 
E_tot = sqrt((1.0_rp + p%vec(6))**2 + (mass_of(p%species)/p%p0c)**2) ! * p0
dt = t - p%t
ds = dt*(c_light*pz0/E_tot)
call track_a_drift(p,ds)

if (p%s > branch%ele(branch%n_ele_track)%s) then
  p%ix_ele = branch%n_ele_track
  p%location = downstream_end$
else
  p%ix_ele = element_at_s(branch, p%s, (ds < 0), position=position)
  p%location = position%location
endif

end subroutine drift_particle_to_t

end module
