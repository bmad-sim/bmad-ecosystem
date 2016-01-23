!+
! Subroutine track_a_drift (orb, length)
!
! Subroutine to track a particle as through a drift.
!
! Input:
!   orb      -- coord_struct: Orbit at start of the drift.
!   length   -- Real(rp): Length to drift through.
!
! Output:
!   orb      -- coord_struct: Orbit at end of the drift
!-

subroutine track_a_drift (orb, length)

use bmad_struct, only: coord_struct, ele_struct, lat_param_struct, rp, c_light, mass_of, lost_z_aperture$

implicit none

type (coord_struct) orb
type (ele_struct) ele
type (lat_param_struct) param
real(rp) length, rel_pc, dz, px, py, pz, pxy2, beta_ref

! Everything but photons

rel_pc = 1 + orb%vec(6)
px = orb%vec(2) / rel_pc
py = orb%vec(4) / rel_pc
pxy2 = px**2 + py**2
if (pxy2 >= 1) then
  orb%state = lost_z_aperture$
  return
endif
pz = sqrt(1 - pxy2)
beta_ref = orb%p0c / sqrt(mass_of(orb%species)**2 + orb%p0c**2)

orb%vec(1) = orb%vec(1) + length * px / pz
orb%vec(3) = orb%vec(3) + length * py / pz

if (orb%beta > 0) then
  dz = length * (orb%beta /beta_ref - 1/pz)
  orb%t = orb%t + length / (orb%beta * pz * c_light)
else
  dz = length * (1 - 1/pz)
endif

orb%vec(5) = orb%vec(5) + dz
orb%s = orb%s + length

end subroutine track_a_drift

