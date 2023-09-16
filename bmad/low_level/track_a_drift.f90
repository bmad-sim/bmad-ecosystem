!+
! Subroutine track_a_drift (orb, length, mat6, make_matrix, ele_orientation, include_ref_motion, time)
!
! Bmad_standard tracking of particle as through a drift. Not to be used with photons.
!
! Input:
!   orb                 -- coord_struct: Orbit at start of the drift.
!   length              -- real(rp): Length to drift through in body coordinates.
!                       --    If orb%direction = 1, positive length is in +z direction and vice versa.
!   mat6(6,6)           -- real(rp), optional: Transfer matrix up to the drift.
!   make_matrix         -- logical, optional: Propagate the transfer matrix? Default is false.
!   ele_orientation     -- integer, optional: Element orientation. Default is orb%direction.
!   include_ref_motion  -- logical, optional: Include effect of the motion of the reference particle?
!                           Default is True. False is basically only used by offset_particle.
!                           Additionally, if False, orb%s is not changed.
!   time                -- real(rp), optional: Particle time before drifting. Typically this is an RF
!                           clock time which may not be equal to orb%t
!
! Output:
!   orb         -- coord_struct: Orbit at end of the drift.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix including the drift.
!   time        -- real(rp), optional: Updated time.
!-

subroutine track_a_drift (orb, length, mat6, make_matrix, ele_orientation, include_ref_motion, time)

use bmad_routine_interface, dummy => track_a_drift

implicit none

type (coord_struct) orb
type (ele_struct) ele
type (lat_param_struct) param

real(rp), optional :: mat6(6,6), time
real(rp) matd(6,6), e_tot_ref, e_particle, rel_len, dt
real(rp) length, rel_pc, dz, px, py, ps, delta, pxy2, mc2, beta_ref

integer, optional :: ele_orientation
integer rel_z_vel
logical, optional :: make_matrix, include_ref_motion

! If the element orientation is opposite the particle direction, px and py are reversed.

if (length == 0) return
delta = orb%vec(6)
rel_pc = 1 + delta
rel_z_vel = integer_option(orb%direction, ele_orientation) * orb%direction

px = rel_z_vel * orb%vec(2) / rel_pc
py = rel_z_vel * orb%vec(4) / rel_pc
pxy2 = px**2 + py**2
if (pxy2 >= 1) then
  orb%state = lost_pz_aperture$
  return
endif
ps = sqrt(1 - pxy2)

orb%vec(1) = orb%vec(1) + length * px / ps
orb%vec(3) = orb%vec(3) + length * py / ps

! Length is the length in body coordinates

if (orb%beta > 0) then
  if (logic_option(.true., include_ref_motion)) then
    mc2 = mass_of(orb%species)
    ! dz = length * ([beta/beta_ref - 1] - [1/ps - 1])
    if (orb%direction == 1) then
      dz = length * (sqrt_one((mc2**2 * (2*delta+delta**2))/((orb%p0c*rel_pc)**2 + mc2**2)) + sqrt_one(-pxy2)/ps)
    else
      beta_ref = orb%p0c / sqrt(orb%p0c**2 + mc2**2)
      dz = length * ((-orb%beta/beta_ref - 1.0_rp) - (1.0_rp / ps - 1.0_rp))
    endif
    orb%s = orb%s + orb%direction * length
  else
    dz = -length /ps
  endif

  dt = rel_z_vel * length / (orb%beta * ps * c_light)
  orb%t = orb%t + dt
  if (present(time)) time = time + dt

else
  if (logic_option(.true., include_ref_motion)) then
    dz = length * (1 - 1/ps)
    orb%s = orb%s + orb%direction * length
  else
    dz = -length /ps
  endif
endif

orb%vec(5) = orb%vec(5) + rel_z_vel * dz

if (logic_option(.false., make_matrix)) then
  call mat_make_unit(matd)
  rel_len = length / (rel_pc * ps)
  matd(1,2) =  rel_len * (px**2 / ps**2 + 1)
  matd(3,4) =  rel_len * (py**2 / ps**2 + 1)
  matd(1,4) =  rel_len * px*py / ps**2
  matd(3,2) =  rel_len * px*py / ps**2
  matd(1,6) = -rel_len * px / ps**2
  matd(3,6) = -rel_len * py / ps**2
  matd(5,2) = -rel_len * px / ps**2 
  matd(5,4) = -rel_len * py / ps**2
  if (logic_option(.true., include_ref_motion)) then
    e_tot_ref = sqrt(orb%p0c**2 + mass_of(orb%species)**2)
    e_particle = orb%p0c * (1 + orb%vec(6)) / orb%beta
    matd(5,6) =  rel_len * (px**2 + py**2) / ps**2 + length * mass_of(orb%species)**2 * e_tot_ref / e_particle**3
  else
    matd(5,6) =  rel_len * (px**2 + py**2) / ps**2
  endif

  mat6 = matmul(matd, mat6)
endif

end subroutine track_a_drift

