!+
! Subroutine apply_energy_kick (dE, orbit, mat6, make_matrix)
! 
! Routine to change the energy of a particle by an amount dE.
! Appropriate changes to z and beta will be made.
!
! Note: the matrix propagation assumes a thin lens energy kick which is independent of px, py and pz.
!
! Input:
!   dE          -- real(rp): Energy change
!   orbit       -- coord_struct: Beginning coordinates
!   ddE_dr(2)   -- real(rp), Derivatives of dE [ddE_dx, ddE_dy].
!   mat6(6,6)   -- real(rp), optional: Transfer matrix before fringe.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: coordinates with added dE energy kick.
!     %vec(6)    -- Set to -1 if particle energy becomes negative. 
!     %state     -- Set to lost_pz$ is particle energy becomes negative.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix transfer matrix including energy kick.
!-

subroutine apply_energy_kick (dE, orbit, ddE_dr, mat6, make_matrix)

use equal_mod, except_dummy => apply_energy_kick

implicit none

type (coord_struct) orbit
real(rp) dE, ddE_dr(2), mc2, pc_new, beta_new, p0c, E_new, pc, beta_old, t3, f, E_old, kmat(6,6)
real(rp), optional :: mat6(6,6)
logical, optional :: make_matrix

!

mc2 = mass_of(orbit%species)
p0c = orbit%p0c
pc = (1 + orbit%vec(6)) * p0c
beta_old = orbit%beta

E_old = pc / orbit%beta
E_new = E_old + dE

if (E_new < mc2) then
  orbit%vec(6) = -1
  orbit%beta = 0
  orbit%state = lost_pz$
  return
endif

orbit%vec(6) = orbit%vec(6) + (1 + orbit%vec(6)) * sqrt_one((2*E_old*dE + dE**2)/pc**2)
pc_new = p0c * (1 + orbit%vec(6))
beta_new = pc_new / E_new

!! t3 = sqrt(E_new**2 - mc2**2) / p0c - (1 + orbit%vec(6))  !! Should be zero

if (logic_option(.false., make_matrix)) then
  call mat_make_unit(kmat)
  f = 1 / (p0c * beta_new)
  kmat(6,1) = f * ddE_dr(1) 
  kmat(6,3) = f * ddE_dr(2)
  kmat(6,6) = beta_old / beta_new
  f = orbit%vec(5) * mc2**2 / (beta_old * pc_new * E_new**2)
  kmat(5,1) = f * ddE_dr(1)
  kmat(5,3) = f * ddE_dr(2) 
  kmat(5,5) = beta_new / beta_old
  kmat(5,6) = orbit%vec(5) * p0c * (1 - beta_new**2 * (1 + (mc2/pc)**2 * E_new / E_old)) / pc_new
  mat6 = matmul (kmat, mat6)
endif

orbit%vec(5) = orbit%vec(5) * beta_new / beta_old
orbit%beta = beta_new

end subroutine apply_energy_kick

