!+
! Function low_energy_z_correction (orbit, ele, ds, mat6, make_matrix) result (dz)
! 
! Routine to add a correction to the change in the z-coordinate due to speed < c_light.
! This routine assumes a constant particle velocity. The change in z is
!   dz = (ds - d_particle) + ds * (beta - beta_ref) / beta_ref
! where d_particle is the distance the particle travels. 
! The second term in the above equation is the correction calculated by this routine.
!
! Input:
!   orbit   -- coord_struct: Position before correction
!   ele     -- ele_sturct: Element being tracked through.
!   ds      -- real(rp): Longitudinal distance traveled by reference particle.
!   mat6(6,6)    -- Real(rp), optional: Transfer matrix before the multipole.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   dz      -- real(rp): Change in z.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix transfer matrix including multipole.
!-

function low_energy_z_correction (orbit, ele, ds, mat6, make_matrix) result (dz)

use bmad_struct

implicit none

type (coord_struct) orbit
type (ele_struct) ele

real(rp) ds, dz, f
real(rp) p0c, beta0, mass, e_tot, pz
real(rp), optional :: mat6(6,6)
logical, optional :: make_matrix

!

mass = mass_of(orbit%species)
p0c = orbit%p0c
e_tot = ele%value(e_tot$)
beta0 = p0c / e_tot
pz = orbit%vec(6)

if (mass * (beta0*pz)**2 < 3d-7 * e_tot) then
  f = beta0**2 * (2 * beta0**2 - (mass / e_tot)**2 / 2)
  dz = ds * pz * (1 - 3 * (pz * beta0**2) / 2 + pz**2 * f) * (mass / e_tot)**2
  if (logic_option(.false., make_matrix)) then
    mat6(5,:) = mat6(5,:) + ds * (1 - 3 * pz * beta0**2 + 3 * pz**2 * f) * (mass / e_tot)**2 * mat6(6,:)
  endif

else
  dz = ds * (orbit%beta - beta0) / beta0
  if (logic_option(.false., make_matrix)) then
    mat6(5,:) = mat6(5,:) + ds * (mass**2 * e_tot / ((1 + pz) * p0c / orbit%beta)**3) * mat6(6,:)
  endif
endif

end function low_energy_z_correction

