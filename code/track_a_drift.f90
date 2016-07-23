!+
! Subroutine track_a_drift (orb, length, mat6, make_matrix)
!
! Subroutine to track a particle as through a drift.
!
! Input:
!   orb         -- coord_struct: Orbit at start of the drift.
!   length      -- Real(rp): Length to drift through.
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix before fringe.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orb        -- coord_struct: Orbit at end of the drift.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix transfer matrix including fringe.
!-

subroutine track_a_drift (orb, length, mat6, make_matrix)

use bmad_struct

implicit none

type (coord_struct) orb
type (ele_struct) ele
type (lat_param_struct) param

real(rp), optional :: mat6(6,6)
real(rp) matd(6,6), e_tot, e_particle, rel_len
real(rp) length, rel_pc, dz, px, py, pz, pxy2, beta_ref

logical, optional :: make_matrix

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

if (logic_option(.false., make_matrix)) then
  call mat_make_unit(matd)
  e_tot = sqrt(orb%p0c**2 + mass_of(orb%species)**2)
  e_particle = orb%p0c * (1 + orb%vec(6)) / orb%beta
  rel_len = length / (rel_pc * pz)
  matd(1,2) =  rel_len * (px**2 / pz**2 + 1)
  matd(3,4) =  rel_len * (py**2 / pz**2 + 1)
  matd(1,4) =  rel_len * px*py / pz**2
  matd(3,2) =  rel_len * px*py / pz**2
  matd(1,6) = -rel_len * px / pz**2
  matd(3,6) = -rel_len * py / pz**2
  matd(5,2) = -rel_len * px / pz**2 
  matd(5,4) = -rel_len * py / pz**2
  matd(5,6) =  rel_len * (px**2 + py**2) / pz**2 + &
                    length * mass_of(orb%species)**2 * e_tot / e_particle**3
  mat6 = matmul(matd, mat6)
endif

end subroutine track_a_drift

