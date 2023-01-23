module space_charge_mod

use bmad

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine sc_field_calc (bunch, include_image, sc_field)
!
! Routine to calculate the space charge field of a bunch
!
! Input:
!   bunch           -- bunch_struct: Starting bunch position in time-based coordinates
!   include_image   -- logical: True if cathode image charge fields are to be included.
!
! Output:
!   include_image   -- logical: Set False if image charge calc no longer needed (Note: never set True).
!   sc_field        -- em_field_struct: space charge field at particle positions
!-

subroutine sc_field_calc (bunch, include_image, sc_field)

use csr_and_space_charge_mod

implicit none

type (bunch_struct), target :: bunch
type (coord_struct), pointer :: p
type (csr_particle_position_struct) :: position(size(bunch%particle))
type (em_field_struct) :: sc_field(size(bunch%particle))
type (mesh3d_struct) :: mesh3d, mesh3d_image

integer :: n, i, imin(1)
real(rp) :: beta, ratio
real(rp) :: Evec(3), Bvec(3), Evec_image(3)
logical :: include_image, err

! Initialize variables
mesh3d%nhi = space_charge_com%space_charge_mesh_size

! Gather alive particles

n = 0
beta = 0
do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state /= alive$) cycle
  n = n + 1
  position(n)%r = [p%vec(1), p%vec(3), p%s]
  position(n)%charge = p%charge * charge_of(p%species)
  beta = beta + p%beta
enddo
if (n<2) return
beta = beta/n

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
  if (p%state /= alive$) then
    sc_field(i)%E = [0,0,0]
    sc_field(i)%B = [0,0,0]
  else
    call interpolate_field(p%vec(1), p%vec(3), p%s,  mesh3d, E=sc_field(i)%E, B=sc_field(i)%B)
  end if
enddo

end subroutine sc_field_calc

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine sc_step(bunch, ele, include_image, t_end)
!
! Subroutine to track a bunch through a given time step with space charge
!
! Input:
!   bunch         -- bunch_struct: Starting bunch position in t-based coordinates
!   ele           -- ele_struct: Element being tracked through.
!   include_image -- logical: Include image charge forces?
!   t_end         -- real(rp): Time at which the tracking ends.
!
! Output:
!   bunch         -- bunch_struct: Ending bunch position in t-based coordinates after space charge kick.
!   include_image -- logical: Set False if image charge calc no longer needed (Note: never set True).
!-

subroutine sc_step(bunch, ele, include_image, t_end)

implicit none

type (bunch_struct), target :: bunch
type (ele_struct) :: ele
type (em_field_struct) :: extra_field(size(bunch%particle))
type (coord_struct), pointer :: p

real(rp) t_end, sum_z
logical include_image
integer i

!

call sc_field_calc (bunch, include_image, extra_field)

! Generate particles at the cathode

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state == pre_born$ .and. p%t <= t_end) p%state = alive$
enddo 

! And track

call track_bunch_time(bunch, ele, t_end, 1e30_rp, extra_field=extra_field)

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
!
! Output:
!   bunch         -- bunch_struct: Ending bunch position in t-based coordinates.
!   include_image -- logical: Set False if image charge calc no longer needed (Note: never set True).
!   dt_next       -- real(rp): Next SC time step the tracker would take based on the error tolerance
!   dt_step       -- real(rp): Step done.
!-

subroutine sc_adaptive_step(bunch, ele, include_image, t_now, dt_step, dt_next)

implicit none

type (bunch_struct) :: bunch, bunch_full, bunch_half
type (ele_struct) ele
type (coord_struct), pointer :: p

real(rp) :: t_now, dt_step, dt_next, sqrt_N
real(rp) :: r_err(6), r_scale(6), rel_tol, abs_tol, err_max
real(rp), parameter :: safety = 0.9_rp, p_grow = -0.25_rp
real(rp), parameter :: p_shrink = -0.25_rp, err_con = 1.89d-4
real(rp), parameter :: tiny = 1.0e-30_rp

integer i, N
logical include_image

!

sqrt_N = sqrt(abs(1/(c_light*dt_step)))  ! number of steps we would take to cover 1 meter
rel_tol = space_charge_com%rel_tol_tracking / sqrt_N
abs_tol = space_charge_com%abs_tol_tracking / sqrt_N

bunch_full = bunch
bunch_half = bunch
dt_next = dt_step

do
  ! Full step
  call sc_step(bunch_full, ele, include_image, t_now+dt_step)
  ! Two half steps
  call sc_step(bunch_half, ele, include_image, t_now+dt_step/2)
  call sc_step(bunch_half, ele, include_image, t_now+dt_step)

  r_err = [0,0,0,0,0,0]
  N = 0
  ! Calculate error from the difference
  do i = 1, size(bunch%particle)
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
  r_scale = abs(bunch_rms_vec(bunch))+abs(bunch_rms_vec(bunch_half)) + tiny
  err_max = maxval(r_err(:)/(r_scale*rel_tol + abs_tol))
  ! If error is larger than tolerance, try again with a smaller step
  if (err_max <= 1.0) exit
  dt_step = safety * dt_step * (err_max**p_shrink)
  bunch%n_bad = bunch%n_bad + 1
  bunch_full = bunch
  bunch_half = bunch
enddo

bunch%n_good = bunch%n_good + 1
! Adjust next step size
! Copied from Runge-Kutta
if (err_max > err_con) then
  dt_next = safety * dt_step * (err_max**p_grow)
else
  dt_next = 5.0_rp * dt_step
endif

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

end subroutine sc_adaptive_step

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine track_to_s (bunch, s, branch)
!
! Drift a bunch of particles to the same s coordinate
!
! Input:
!   bunch     -- bunch_struct: Input bunch position in s-based coordinate.
!   s         -- real(rp): Target s coordinate.
!   branch    -- branch_struct: Branch being tracked through.
!
! Output:
!   bunch     -- bunch_struct: Output bunch position in s-based coordinate. Particles will be at the same s coordinate
!-

subroutine track_to_s (bunch, s, branch)

implicit none

type (bunch_struct), target :: bunch
type (coord_struct), pointer :: p
type (branch_struct) :: branch
type (coord_struct) :: position

integer i, dir_save
real(rp) s, s_end

! Convert bunch to s-based coordinates

s_end = min(s, branch%ele(branch%n_ele_track)%s)

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state /= alive$) cycle

  dir_save = p%direction
  if (s_end-s < 0) then
    p%direction = -1
  else
    p%direction = 1
  endif

  call track_from_s_to_s (branch%lat, p%s, s_end, p, p, ix_branch = branch%ix_branch)

  p%ix_ele = element_at_s(branch, p%s, (p%direction < 0), position=position)
  p%location = position%location
  p%direction = dir_save
enddo

end subroutine track_to_s

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine track_to_t (bunch, t, branch)
!
! Drift a bunch of particles to the same t coordinate
!
! Input:
!   bunch     -- bunch_struct: Input bunch position in s-based coordinate.
!   t         -- real(rp): Target t coordinate.
!   branch    -- branch_struct: Lattice branch being tracked through.
!
! Output:
!   bunch     -- bunch_struct: Output bunch position in s-based coordinate. Particles will be at the same t coordinate
!-

subroutine track_to_t (bunch, t, branch)

implicit none

type (bunch_struct), target :: bunch
type (coord_struct), pointer :: p
type (branch_struct) :: branch
type (coord_struct) :: position

integer i, dir_save
real(rp) t, pz0, E_tot, dt, s_end

! Convert bunch to s-based coordinates

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  if (p%state /= alive$) cycle

  dir_save = p%direction
  pz0 = sqrt((1.0_rp + p%vec(6))**2 - p%vec(2)**2 - p%vec(4)**2 ) ! * p0 
  E_tot = sqrt((1.0_rp + p%vec(6))**2 + (mass_of(p%species)/p%p0c)**2) ! * p0
  dt = t - p%t
  s_end = min(p%s + dt*(c_light*pz0/E_tot), branch%ele(branch%n_ele_track)%s)
  call track_from_s_to_s (branch%lat, p%s, s_end, p, p, ix_branch = branch%ix_branch)

  p%ix_ele = element_at_s(branch, p%s, (p%direction < 0), position=position)
  p%location = position%location
  p%direction = dir_save
enddo

end subroutine track_to_t
  
end module
